use winds_engine::terrain::{HeightGrid, TerrainFlowField};
use winds_engine::{ConstantField, Wind, WindField};

/// A 41×21 grid (100 m cells) with a smooth ridge running north–south.
fn ridge_grid() -> HeightGrid {
    let (ncols, nrows) = (41usize, 21usize);
    let mut heights = Vec::with_capacity(ncols * nrows);
    for i in 0..nrows {
        for j in 0..ncols {
            let x = (j as f64 - 20.0) * 100.0;
            let _z = (i as f64 - 10.0) * 100.0;
            heights.push(400.0 * (-(x * x) / (500.0_f64 * 500.0)).exp());
        }
    }
    HeightGrid::new(-2000.0, -1000.0, 100.0, 100.0, ncols, nrows, heights).unwrap()
}

fn east_wind(speed: f64) -> ConstantField {
    // FROM 0° = from +X = blows toward −X (a headwind). For an eastward
    // frame we want air moving +X, i.e. FROM 180°.
    ConstantField {
        wind: Wind::from_speed_direction(speed, 180.0),
    }
}

#[test]
fn ridge_speeds_up_over_the_crest_and_returns_to_free_stream() {
    let grid = ridge_grid();
    let base = east_wind(10.0);
    let field =
        TerrainFlowField::adjust(grid.clone(), &base, &[], &[50.0], 2000.0, 300.0, 400.0).unwrap();
    let speed_at = |x: f64| {
        let w = field.wind(0.0, [x, grid.height_at(x, 0.0) + 50.0, 0.0]);
        w.speed_mps()
    };
    // A symmetric hill under uniform wind speeds up monotonically toward
    // the crest and relaxes to the free stream far upwind — no shelter on a
    // symmetric ridge (Jackson-Hunt behaviour).
    let upwind = speed_at(-1500.0);
    let windward = speed_at(-500.0);
    let crest = speed_at(0.0);
    assert!(
        crest > windward && windward > upwind,
        "crest {crest} > windward {windward} > upwind {upwind}"
    );
    assert!((upwind - 10.0).abs() < 1.5, "far upwind {upwind}");
    // Below the linearized speed-up bound 2H/L.
    assert!(crest < 10.0 * (1.0 + 2.0 * 400.0 / 500.0), "crest {crest}");
    // Symmetric so the lee mirrors the windward slope.
    assert!((speed_at(500.0) - windward).abs() < 0.05);
}

#[test]
fn flat_ground_with_uniform_wind_is_unchanged() {
    let grid = HeightGrid::new(0.0, 0.0, 100.0, 100.0, 6, 6, vec![100.0; 36]).unwrap();
    let base = east_wind(8.0);
    let field =
        TerrainFlowField::adjust(grid.clone(), &base, &[], &[30.0], 1000.0, 200.0, 400.0).unwrap();
    for x in [50.0, 250.0, 450.0] {
        for z in [50.0, 250.0, 450.0] {
            let w = field.wind(0.0, [x, grid.height_at(x, z) + 30.0, z]);
            assert!(
                (w.speed_mps() - 8.0).abs() < 1e-3,
                "speed {} at ({x}, {z})",
                w.speed_mps()
            );
        }
    }
}

#[test]
fn monitors_bend_the_wind_near_them_and_relax_away() {
    let grid = HeightGrid::new(0.0, 0.0, 100.0, 100.0, 21, 21, vec![0.0; 21 * 21]).unwrap();
    // Background east at 5 m/s; a monitor at (1000, 1000) says FROM 90°
    // (blowing toward −z, a crosswind).
    let monitors = [(1000.0, 1000.0, 5.0, 90.0)];
    let field = TerrainFlowField::adjust(
        grid,
        &east_wind(5.0),
        &monitors,
        &[20.0],
        400.0,
        500.0,
        400.0,
    )
    .unwrap();
    let near = field.wind(0.0, [1000.0, 20.0, 1000.0]);
    // Near the monitor the wind approaches the observed crosswind direction.
    assert!(
        near.z_mps.abs() > near.x_mps.abs(),
        "near the monitor should blow crosswind: {near:?}"
    );
    // Away from the monitor it relaxes toward the east background (the
    // domain edge holds the background wind).
    let far = field.wind(0.0, [0.0, 20.0, 0.0]);
    assert!(
        far.x_mps.abs() > far.z_mps.abs(),
        "far away should return toward eastward: {far:?}"
    );
}

#[test]
fn adjustment_is_deterministic() {
    let grid = ridge_grid();
    let a = TerrainFlowField::adjust(
        grid.clone(),
        &east_wind(10.0),
        &[(0.0, 0.0, 12.0, 270.0)],
        &[40.0],
        2000.0,
        300.0,
        400.0,
    )
    .unwrap();
    let b = TerrainFlowField::adjust(
        grid,
        &east_wind(10.0),
        &[(0.0, 0.0, 12.0, 270.0)],
        &[40.0],
        2000.0,
        300.0,
        400.0,
    )
    .unwrap();
    for x in [-1500.0, -500.0, 0.0, 500.0, 1500.0] {
        assert_eq!(a.wind(0.0, [x, 140.0, 0.0]), b.wind(0.0, [x, 140.0, 0.0]));
    }
}

#[test]
fn ascii_dem_parses_with_row_reversal() {
    let text = "ncols 3\nnrows 2\nxllcorner 100.0\nyllcorner 200.0\ncellsize 50.0\n100.0 101.0 102.0\n200.0 201.0 202.0\n";
    let grid = HeightGrid::from_ascii(text).unwrap();
    assert_eq!((grid.ncols, grid.nrows), (3, 2));
    assert_eq!(grid.origin_x, 100.0);
    assert_eq!(grid.dx, 50.0);
    // Rows reversed: north row is the file's second line (200,201,202).
    assert_eq!(&grid.heights[0..3], &[200.0, 201.0, 202.0]);
    assert_eq!(&grid.heights[3..6], &[100.0, 101.0, 102.0]);
}

#[test]
fn validates_its_inputs() {
    let grid = HeightGrid::new(0.0, 0.0, 100.0, 100.0, 6, 6, vec![0.0; 36]).unwrap();
    assert!(TerrainFlowField::adjust(
        grid.clone(),
        &east_wind(5.0),
        &[],
        &[],
        1000.0,
        300.0,
        400.0
    )
    .is_err());
    assert!(
        TerrainFlowField::adjust(grid, &east_wind(5.0), &[], &[30.0], 0.0, 300.0, 400.0).is_err()
    );
    assert!(HeightGrid::new(0.0, 0.0, 100.0, 100.0, 2, 2, vec![0.0; 3]).is_err());
}
