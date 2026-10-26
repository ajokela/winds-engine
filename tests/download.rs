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

/// Live NOAA round-trip; run explicitly with
/// `cargo test --features download -- --include-ignored`.
#[test]
#[ignore]
fn live_fetch_decodes_a_profile() {
    let profile = winds_engine::download::fetch_profile(&request()).unwrap();
    assert!(!profile.levels().is_empty());
    println!("live profile: {:?}", profile.levels());
}
