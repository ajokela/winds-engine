#![cfg(feature = "grib")]

use winds_engine::grib::profile_from_grib2;

fn golden() -> Vec<u8> {
    std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/gfs_20261026_00z_mpx.grib2"
    ))
    .unwrap()
}

#[test]
fn golden_gfs_subset_decodes_to_a_profile() {
    // Golden file: NOMADS filter_gfs_0p25.pl, cycle 2026-10-26 00Z, f000,
    // levels 850/700/500 mb, vars HGT/UGRD/VGRD, subregion 44-46 N, 95-92 W.
    let profile = profile_from_grib2(&golden(), 44.883, -93.233, 0.0).unwrap();
    let levels = profile.levels();
    assert_eq!(levels.len(), 3);
    // Isobaric order: 500 mb highest, 850 mb lowest.
    assert!(levels[0].altitude_m < levels[2].altitude_m);
    for level in levels {
        println!(
            "PIN alt={:.4} speed={:.4} dir={:.4}",
            level.altitude_m, level.speed_mps, level.direction_from_deg
        );
    }
    // Pinned decode of the golden file at the nearest grid point (0.1%
    // tolerance for float formatting).
    let expected = [
        (1583.0547, 2.5750, 63.5293),
        (3167.6299, 15.0881, 249.2842),
        (5809.4062, 21.6409, 259.7409),
    ];
    for (level, &(alt, speed, dir)) in levels.iter().zip(&expected) {
        assert!(
            (level.altitude_m - alt).abs() / alt < 0.001,
            "alt {} vs {alt}",
            level.altitude_m
        );
        assert!(
            (level.speed_mps - speed).abs() / speed.max(1.0) < 0.001,
            "speed {} vs {speed}",
            level.speed_mps
        );
        assert!(
            (level.direction_from_deg - dir).abs() < 0.001,
            "dir {} vs {dir}",
            level.direction_from_deg
        );
    }
    // Wind speeds are meteorologically sane.
    for level in levels {
        assert!(level.speed_mps >= 0.0 && level.speed_mps < 100.0);
        assert!(level.direction_from_deg >= 0.0 && level.direction_from_deg < 360.0);
    }
}

#[test]
fn direction_rotates_with_launch_azimuth() {
    let north = profile_from_grib2(&golden(), 44.883, -93.233, 0.0).unwrap();
    let east = profile_from_grib2(&golden(), 44.883, -93.233, 90.0).unwrap();
    for (n, e) in north.levels().iter().zip(east.levels()) {
        let rotated = (n.direction_from_deg - 90.0).rem_euclid(360.0);
        assert!(
            (rotated - e.direction_from_deg).abs() < 1e-9,
            "{rotated} vs {}",
            e.direction_from_deg
        );
        assert_eq!(n.speed_mps, e.speed_mps);
    }
}

#[test]
fn validates_its_inputs() {
    let data = golden();
    assert!(profile_from_grib2(&data, 91.0, 0.0, 0.0).is_err());
    assert!(profile_from_grib2(&data, 45.0, f64::NAN, 0.0).is_err());
    assert!(profile_from_grib2(b"not grib", 45.0, 0.0, 0.0).is_err());
}
