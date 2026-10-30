//! A reproducible, CI-pinned wind model of the Hibbing Rifle & Pistol Club
//! 600 m lane over the checked-in MnTOPO DEM.
//!
//! Everything is fixed — DEM bytes, monitor positions/readings, the
//! background wind and the sample points — so the digest is a golden value
//! that must not change. If it does, either the terrain model changed or the
//! fixture did; investigate before re-pinning.
#![cfg(feature = "tiff")]

use nalgebra::Vector3;
use winds_engine::{geotiff_grid, ConstantField, Monitor, TerrainFlowField, WindField};

const FIRING_LAT: f64 = 47.403947;
const FIRING_LON: f64 = -92.995347;
const TARGET_LAT: f64 = 47.407748;
const TARGET_LON: f64 = -92.990340;

fn web_mercator(lat: f64, lon: f64) -> (f64, f64) {
    let x = lon * 20_037_508.34 / 180.0;
    let y = (lat.to_radians() / 2.0 + std::f64::consts::FRAC_PI_4)
        .tan()
        .ln()
        * 20_037_508.34
        / std::f64::consts::PI;
    (x, y)
}

fn digest() -> String {
    let dem = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/hibbing_lane.tif"
    ))
    .unwrap();
    let grid = geotiff_grid(&dem).unwrap().downsample(30.0).unwrap();
    let (x0, y0) = web_mercator(FIRING_LAT, FIRING_LON);
    let (x1, y1) = web_mercator(TARGET_LAT, TARGET_LON);

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
        let toward = (from_north_deg + 180.0_f64).to_radians();
        let north_mps = speed * toward.cos();
        let east_mps = speed * toward.sin();
        let d = (north_mps.atan2(east_mps).to_degrees() - 180.0).rem_euclid(360.0);
        monitors.push(Monitor::at(x, z, speed, d));
    }

    let base = ConstantField(Vector3::new(6.0, 0.0, 0.0));
    let field = TerrainFlowField::adjust(
        grid.clone(),
        &base,
        &monitors,
        &[10.0, 50.0, 150.0, 400.0],
        800.0,
        200.0,
        300.0,
    )
    .unwrap();

    let mut out = String::new();
    for t in [0.00_f64, 0.25, 0.50, 0.75, 1.00] {
        let x = x0 + (x1 - x0) * t;
        let z = y0 + (y1 - y0) * t;
        let h = grid.height_at(x, z) + 10.0;
        let w = field.wind(0.0, Vector3::new(x, h, z));
        out.push_str(&format!("{t:.2}:{:.6}:{:.6};", w.x, w.z));
    }
    out
}

#[test]
fn lane_wind_digest_is_stable() {
    let got = digest();
    // Golden value pinned 2026-10-30 from the example (see examples/range.rs).
    let expected = "0.00:2.735458:1.881235;0.25:2.714291:1.887971;0.50:2.694908:1.899095;0.75:2.683987:1.916732;1.00:2.688016:1.941466;";
    assert_eq!(got, expected, "lane wind digest changed");
}

#[test]
fn lane_wind_is_deterministic() {
    assert_eq!(digest(), digest());
}
