//! A deterministic, reproducible wind model of a real rifle range.
//!
//! Hibbing Rifle & Pistol Club's 600 m lane, firing line to target, over the
//! checked-in MnTOPO DEM (`tests/data/hibbing_lane.tif`). Everything is fixed
//! — the DEM bytes, monitor positions and readings, the background wind, and
//! the sample points — so the printed digest is stable run to run and can be
//! pinned in CI.
//!
//! Run with: cargo run --example range --features tiff
use nalgebra::Vector3;
use winds_engine::{geotiff_grid, ConstantField, Monitor, TerrainFlowField, WindField};

/// Lane endpoints: firing line and target (degrees).
const FIRING_LAT: f64 = 47.403947;
const FIRING_LON: f64 = -92.995347;
const TARGET_LAT: f64 = 47.407748;
const TARGET_LON: f64 = -92.990340;

/// Firing → target over the flat Earth: the lane's downrange unit vector in
/// (east, north), and its length in metres.
fn lane_geometry() -> ((f64, f64), (f64, f64), f64) {
    const R: f64 = 6_371_000.0;
    let latm = ((FIRING_LAT + TARGET_LAT) / 2.0).to_radians();
    let east = (TARGET_LON - FIRING_LON).to_radians() * R * latm.cos();
    let north = (TARGET_LAT - FIRING_LAT).to_radians() * R;
    let length = (east * east + north * north).sqrt();
    let (ex, nx) = (east / length, north / length);
    ((FIRING_LON, FIRING_LAT), (ex, nx), length)
}

fn web_mercator(lat: f64, lon: f64) -> (f64, f64) {
    let x = lon * 20_037_508.34 / 180.0;
    let y = (lat.to_radians() / 2.0 + std::f64::consts::FRAC_PI_4)
        .tan()
        .ln()
        * 20_037_508.34
        / std::f64::consts::PI;
    (x, y)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dem = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/hibbing_lane.tif"
    ))?;
    let grid = geotiff_grid(&dem)?.downsample(30.0)?;
    println!(
        "down-sampled grid {}x{} @ {:.1} m",
        grid.ncols, grid.nrows, grid.dx
    );

    let ((f_lon, f_lat), (_ex, _nx), length) = lane_geometry();
    let (x0, y0) = web_mercator(f_lat, f_lon);
    let (x1, y1) = web_mercator(TARGET_LAT, TARGET_LON);
    println!("lane length {:.1} m", length);

    // Four fixed wind monitors at deterministic distances down the lane.
    // Each reading is (fraction along lane, speed m/s, FROM-direction deg
    // clockwise from north). We rotate the compass reading into the
    // easting/northing frame the grid uses.
    let readings: [(f64, f64, f64); 4] = [
        (0.25, 3.1, 240.0),
        (0.40, 3.8, 225.0),
        (0.60, 2.9, 250.0),
        (0.75, 3.5, 220.0),
    ];
    let mut monitors = Vec::new();
    for &(t, speed, from_north_deg) in &readings {
        let x = x0 + (x1 - x0) * t;
        let z = y0 + (y1 - y0) * t;
        // A wind FROM `from_north_deg` blows TOWARD `toward = from_north_deg + 180`.
        let toward = (from_north_deg + 180.0_f64).to_radians();
        let east_mps = speed * toward.sin();
        let north_mps = speed * toward.cos();
        // Monitor uses "direction_from_deg" with the same convention as
        // horizontal(s,d): FROM-direction measured from +x (east) toward +z.
        // from-x is such that the toward vector is at (d+180): east=s cos(d+180),
        // north=s sin(d+180) => d = atan2(north, east) - 180.
        let d = (north_mps.atan2(east_mps).to_degrees() - 180.0).rem_euclid(360.0);
        monitors.push(Monitor::at(x, z, speed, d));
    }

    // Background wind: 6 m/s from the west (270° compass = blowing toward east).
    let base = ConstantField(Vector3::new(6.0, 0.0, 0.0));

    let field = TerrainFlowField::adjust(
        grid.clone(),
        &base,
        &monitors,
        &[10.0, 50.0, 150.0, 400.0],
        800.0,
        200.0,
        300.0,
    )?;

    // Sample the adjusted field at fixed stations, 10 m above ground.
    let mut digest = String::new();
    for (label, t) in [
        ("firing", 0.00),
        ("25%", 0.25),
        ("mid", 0.50),
        ("75%", 0.75),
        ("target", 1.00),
    ] {
        let x = x0 + (x1 - x0) * t;
        let z = y0 + (y1 - y0) * t;
        let h = grid.height_at(x, z) + 10.0;
        let w = field.wind(0.0, Vector3::new(x, h, z));
        let (wx, wz) = (w.x, w.z);
        let speed = (wx * wx + wz * wz).sqrt();
        println!("{label:>6}  {wx:+.3} {wz:+.3}  {speed:.2} m/s");
        digest.push_str(&format!("{t:.2}:{wx:.6}:{wz:.6};"));
    }
    println!("\nDIGEST {digest}");

    Ok(())
}
