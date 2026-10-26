#![cfg(feature = "download")]

use winds_engine::download::{Model, NomadsRequest};

fn request() -> NomadsRequest {
    NomadsRequest {
        model: Model::Gfs,
        date_utc: 20261026,
        cycle_hour_utc: 0,
        forecast_hour: 0,
        latitude_deg: 44.883,
        longitude_deg: -93.233,
        levels_mb: vec![850, 700, 500],
        launch_azimuth_deg: 0.0,
    }
}

#[test]
fn url_is_the_nomads_subset_form() {
    let url = request().url();
    let expected = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.20261026%2F00%2Fatmos&file=gfs.t00z.pgrb2.0p25.f000&lev_850_mb=on&lev_700_mb=on&lev_500_mb=on&var_HGT=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=-94.733&rightlon=-91.733&toplat=46.383&bottomlat=43.383";
    assert_eq!(url, expected);
}

#[test]
fn hrrr_index_parses_message_ranges() {
    // Golden index: noaa-hrrr-bdp-pds hrrr.20261025/conus
    // hrrr.t23z.wrfprsf00.grib2.idx.
    let idx = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/hrrr_20261025_23z_wrfprsf00.grib2.idx"
    ))
    .unwrap();
    let ranges = winds_engine::download::hrrr_parse_ranges(&idx, &[850, 700, 500]).unwrap();
    // 3 variables x 3 levels = 9 messages, each ended by the next message's
    // offset except the file's last.
    assert_eq!(ranges.len(), 9);
    assert!(ranges
        .iter()
        .all(|(start, end)| end.is_none_or(|e| e > *start)));
    // Pinned offsets from the golden index.
    assert!(ranges.contains(&(203029096, Some(203639003)))); // UGRD 850 mb
    assert!(ranges.contains(&(203639004, Some(204255049)))); // VGRD 850 mb
    assert!(winds_engine::download::hrrr_parse_ranges(&idx, &[999]).is_err());
}

#[test]
fn hrrr_urls_have_the_mirror_form() {
    let mut req = request();
    req.model = Model::Hrrr;
    req.date_utc = 20261025;
    req.cycle_hour_utc = 23;
    let (idx, data) = req.hrrr_urls();
    assert_eq!(
        idx,
        "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.20261025/conus/hrrr.t23z.wrfprsf00.grib2.idx"
    );
    assert_eq!(data, idx.strip_suffix(".idx").unwrap());
    assert!(req.url().starts_with("https://noaa-hrrr-bdp-pds"));
}

/// Live NOAA round-trip; run explicitly with
/// `cargo test --features download -- --include-ignored`.
#[test]
#[ignore]
fn live_fetch_decodes_a_profile() {
    let profile = winds_engine::download::fetch_profile(&request()).unwrap();
    assert!(!profile.levels().is_empty());
    println!("live profile: {:?}", profile.levels());
}

/// Live HRRR round-trip over Minneapolis.
#[test]
#[ignore]
fn live_hrrr_fetch_decodes_a_profile() {
    let mut req = request();
    req.model = Model::Hrrr;
    req.date_utc = 20261025;
    req.cycle_hour_utc = 23;
    let profile = winds_engine::download::fetch_profile(&req).unwrap();
    assert_eq!(profile.levels().len(), 3);
    println!("live hrrr profile: {:?}", profile.levels());
}
