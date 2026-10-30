use nalgebra::Vector3;
use winds_engine::validate;
use winds_engine::{
    ConstantField, EkmanProfile, LogProfile, ProfileField, ProfileLevel, Sounding, WindField,
};

fn level(altitude_m: f64, speed_mps: f64, direction_from_deg: f64) -> ProfileLevel {
    ProfileLevel {
        altitude_m,
        speed_mps,
        direction_from_deg,
    }
}

/// Horizontal wind from a speed and a wind-FROM direction in the frame
/// convention (0° = headwind, 90° = from the right).
fn horizontal(speed_mps: f64, direction_from_deg: f64) -> Vector3<f64> {
    let toward = (direction_from_deg + 180.0).to_radians();
    Vector3::new(speed_mps * toward.cos(), 0.0, speed_mps * toward.sin())
}

/// Horizontal speed of a wind vector.
fn speed(w: Vector3<f64>) -> f64 {
    w.x.hypot(w.z)
}

/// Wind-FROM direction of a wind vector, degrees in [0, 360).
fn direction_from(w: Vector3<f64>) -> f64 {
    (w.z.atan2(w.x).to_degrees() + 180.0).rem_euclid(360.0)
}

/// Signed difference a - b, wrapped to (-180, 180].
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 {
        d - 360.0
    } else {
        d
    }
}

#[test]
fn wind_direction_round_trips() {
    // FROM 0° blows toward −X (a headwind); FROM 90° toward −Z.
    let headwind = horizontal(10.0, 0.0);
    assert!((headwind.x + 10.0).abs() < 1e-12);
    assert!(headwind.z.abs() < 1e-12);
    assert!(direction_from(headwind).abs() < 1e-9);
    let cross = horizontal(10.0, 90.0);
    assert!(cross.x.abs() < 1e-12);
    assert!((cross.z + 10.0).abs() < 1e-12);
    for deg in [1.0, 90.0, 179.9, 181.0, 270.0, 359.9] {
        let w = horizontal(7.0, deg);
        assert!((speed(w) - 7.0).abs() < 1e-9);
        assert!(angle_diff(direction_from(w), deg).abs() < 1e-9, "deg {deg}");
    }
}

#[test]
fn constant_field_is_constant() {
    let field = ConstantField(horizontal(5.0, 45.0));
    assert_eq!(
        field.wind(0.0, Vector3::zeros()),
        field.wind(99.0, Vector3::new(1.0, 2.0, 3.0))
    );
}

#[test]
fn profile_interpolates_speed_linearly() {
    let profile = ProfileField::new(vec![level(0.0, 2.0, 0.0), level(1000.0, 12.0, 0.0)]).unwrap();
    assert!((speed(profile.wind_at(500.0)) - 7.0).abs() < 1e-12);
    // Below and above the levels, the nearest holds.
    assert!((speed(profile.wind_at(-50.0)) - 2.0).abs() < 1e-12);
    assert!((speed(profile.wind_at(5000.0)) - 12.0).abs() < 1e-12);
}

#[test]
fn profile_interpolates_direction_along_the_short_arc() {
    // 350° to 10°: through 0°, never through 180°.
    let profile =
        ProfileField::new(vec![level(0.0, 5.0, 350.0), level(1000.0, 5.0, 10.0)]).unwrap();
    let mid = direction_from(profile.wind_at(500.0));
    assert!(!(10.0..=350.0).contains(&mid), "mid {mid}");
    let quarter = direction_from(profile.wind_at(250.0));
    assert!(quarter > 350.0, "quarter {quarter}");
}

#[test]
fn profile_validates_its_levels() {
    assert!(ProfileField::new(vec![]).is_err());
    assert!(ProfileField::new(vec![level(0.0, -1.0, 0.0)]).is_err());
    assert!(ProfileField::new(vec![level(0.0, 1.0, 0.0), level(0.0, 2.0, 0.0)]).is_err());
    // Unsorted input is accepted and sorted.
    let profile = ProfileField::new(vec![level(1000.0, 12.0, 0.0), level(0.0, 2.0, 0.0)]).unwrap();
    assert!((speed(profile.wind_at(500.0)) - 7.0).abs() < 1e-12);
}

#[test]
fn log_profile_hits_its_reference_and_roughness() {
    let log = LogProfile::new(8.0, 10.0, 0.03, 270.0).unwrap();
    assert!(
        (log.speed_at(10.0) - 8.0).abs() < 1e-12,
        "{}",
        log.speed_at(10.0)
    );
    assert_eq!(log.speed_at(0.03), 0.0);
    assert_eq!(log.speed_at(0.01), 0.0);
    // Doubling height adds a constant increment in a log law.
    let increment = log.speed_at(20.0) - log.speed_at(10.0);
    assert!((increment - (log.speed_at(40.0) - log.speed_at(20.0))).abs() < 1e-12);
    assert!(LogProfile::new(8.0, 10.0, 20.0, 0.0).is_err());
    assert!(LogProfile::new(-1.0, 10.0, 0.03, 0.0).is_err());
}

#[test]
fn ekman_profile_has_the_right_endpoints() {
    let ekman = EkmanProfile::new(12.0, 45.0, 5.0, 0.0).unwrap();
    // Calm at the surface, geostrophic aloft.
    assert!(speed(ekman.wind_at(0.0)) < 1e-12);
    let aloft = ekman.wind_at(5_000.0);
    assert!(
        (speed(aloft) - 12.0).abs() / 12.0 < 1e-3,
        "{}",
        speed(aloft)
    );
    // Aloft it comes FROM the requested geostrophic direction.
    assert!(angle_diff(direction_from(aloft), 0.0).abs() < 1.0);
    // The boundary layer turns with height: speed rises monotonically-ish
    // from zero, and the flow is rotated away from geostrophic near the ground.
    let low = ekman.wind_at(100.0);
    assert!(speed(low) > 0.0 && speed(low) < 12.0);
    assert!(EkmanProfile::new(12.0, 90.0, 5.0, 0.0).is_err());
    assert!(EkmanProfile::new(12.0, 45.0, 0.0, 0.0).is_err());
}

#[test]
fn sounding_parses_the_golden_file() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/mpx_20261023_00z.txt"
    ))
    .unwrap();
    let sounding = Sounding::parse(&text).unwrap();
    assert_eq!(sounding.station.as_deref(), Some("KMPX (Minneapolis, MN)"));
    assert_eq!(sounding.valid.as_deref(), Some("2026-10-23T00:00:00Z"));
    assert_eq!(sounding.levels.len(), 115);
    let first = sounding.levels[0];
    assert_eq!(first.altitude_m, 287.0);
    assert_eq!(first.speed_mps, 0.0);
    let profile = sounding.to_profile_field().unwrap();
    assert_eq!(profile.levels().len(), 115);
}

#[test]
fn sounding_rejects_bad_lines() {
    assert!(Sounding::parse("").is_err());
    assert!(Sounding::parse("# only a comment").is_err());
    let err = Sounding::parse("100 5").unwrap_err();
    assert!(err.reason.contains("line 1"), "{err}");
    let err = Sounding::parse("100 fast 90").unwrap_err();
    assert!(err.reason.contains("not numeric"), "{err}");
}

#[test]
fn validation_is_zero_for_a_perfect_profile() {
    let sounding = Sounding::parse("0 5 90\n1000 10 100\n2000 15 110\n").unwrap();
    let profile = sounding.to_profile_field().unwrap();
    let report = validate::validate(&profile, &sounding);
    assert_eq!(report.levels, 3);
    assert!(report.speed_rmse_mps < 1e-9, "{}", report.speed_rmse_mps);
    assert!(report.direction_rmse_deg < 1e-9);
    assert_eq!(report.layers.len(), 2); // 0–1 km and 1–3 km
}

#[test]
fn validation_measures_a_uniform_offset() {
    let sounding = Sounding::parse("0 5 90\n500 5 90\n1500 5 90\n").unwrap();
    let levels: Vec<ProfileLevel> = sounding
        .levels
        .iter()
        .map(|l| level(l.altitude_m, l.speed_mps + 2.0, l.direction_from_deg))
        .collect();
    let profile = ProfileField::new(levels).unwrap();
    let report = validate::validate(&profile, &sounding);
    assert!(
        (report.speed_rmse_mps - 2.0).abs() < 1e-9,
        "{}",
        report.speed_rmse_mps
    );
    assert!(report.direction_rmse_deg < 1e-9);
    for layer in &report.layers {
        assert!((layer.mean_abs_speed_error_mps - 2.0).abs() < 1e-9);
    }
}

#[test]
fn golden_sounding_validates_against_a_smoothed_profile() {
    // A coarse 6-level profile fit to the golden sounding should land within
    // a few m/s RMS — a smoke test that the metrics are sane on real data.
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/mpx_20261023_00z.txt"
    ))
    .unwrap();
    let sounding = Sounding::parse(&text).unwrap();
    let coarse: Vec<ProfileLevel> = sounding
        .levels
        .iter()
        .step_by(20)
        .map(|l| level(l.altitude_m, l.speed_mps, l.direction_from_deg))
        .collect();
    let profile = ProfileField::new(coarse).unwrap();
    let report = validate::validate(&profile, &sounding);
    // Actual: 4.14 m/s — sane for 6 levels against 115 real measurements.
    assert!(
        report.speed_rmse_mps < 5.0,
        "rmse {}",
        report.speed_rmse_mps
    );
}

#[test]
fn sounding_heights_stay_at_sea_level() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/mpx_20261023_00z.txt"
    ))
    .unwrap();
    let sounding = Sounding::parse(&text).unwrap();
    // The station's surface report is the lowest level: KMPX sits 287 m above sea level.
    assert_eq!(sounding.surface_altitude_m(), Some(287.0));
    // The profile keeps the measured sea-level heights, level for level.
    let msl = sounding.to_profile_field().unwrap();
    assert_eq!(msl.levels().len(), sounding.levels.len());
    assert_eq!(msl.levels()[0].altitude_m, 287.0);
    let first = sounding.levels[0];
    let wind = msl.wind_at(287.0);
    let expected = horizontal(first.speed_mps, first.direction_from_deg);
    assert!((wind - expected).norm() < 1e-9, "{wind:?} vs {expected:?}");
}

#[test]
fn an_empty_sounding_has_no_surface() {
    let sounding = Sounding {
        station: None,
        valid: None,
        levels: vec![],
    };
    assert_eq!(sounding.surface_altitude_m(), None);
}

#[test]
fn ekman_backs_near_the_surface_in_the_north_and_veers_in_the_south() {
    // Holton: in the northern hemisphere the surface wind is turned toward low pressure, 45°
    // to the LEFT of (backed from) the geostrophic wind, and the wind veers with height. With
    // +Y up and +Z to the right of +X (FROM 90° is FROM +Z), "left" is counterclockwise seen
    // from above, which lowers the FROM angle. The southern hemisphere mirrors it.
    let geo = 200.0;
    let north = EkmanProfile::new(10.0, 45.0, 5.0, geo).unwrap();
    let south = EkmanProfile::new(10.0, -45.0, 5.0, geo).unwrap();
    assert!(
        (angle_diff(direction_from(north.wind_at(0.01)), geo) + 45.0).abs() < 0.01,
        "{}",
        angle_diff(direction_from(north.wind_at(0.01)), geo)
    );
    assert!(
        (angle_diff(direction_from(south.wind_at(0.01)), geo) - 45.0).abs() < 0.01,
        "{}",
        angle_diff(direction_from(south.wind_at(0.01)), geo)
    );
    // Veering with height in the north: the FROM angle rises from about geo - 45°.
    let mut previous = angle_diff(direction_from(north.wind_at(1.0)), geo);
    for z in [50.0, 100.0, 200.0, 400.0] {
        let now = angle_diff(direction_from(north.wind_at(z)), geo);
        assert!(now > previous, "z {z}: {now} <= {previous}");
        previous = now;
    }
}
