//! NOAA NOMADS subset downloads, behind the `download` feature.
//!
//! This is the only networked module in the crate. The caller supplies the
//! cycle date and forecast hour explicitly — no wall clock here either, so
//! requests are reproducible.

use crate::field::ProfileField;
use crate::grib::profile_from_grib2;
use crate::{Error, Result};
use std::io::Read;

/// The forecast model to subset. GFS 0.25° pressure-level files carry
/// isobaric levels in one file; HRRR pressure files are not served by the
/// NOMADS subset scripts (documented limit: only GFS is offered here).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Model {
    /// NOAA GFS at 0.25° (~25 km), global coverage, cycles every 6 h.
    Gfs,
}

/// One subset request against the NOMADS filter API.
#[derive(Debug, Clone)]
pub struct NomadsRequest {
    /// Forecast model.
    pub model: Model,
    /// Cycle date, yyyymmdd (UTC).
    pub date_utc: u32,
    /// Cycle hour (0, 6, 12 or 18 for GFS).
    pub cycle_hour_utc: u32,
    /// Forecast hour offset (0 = analysis).
    pub forecast_hour: u32,
    /// Point latitude, degrees north.
    pub latitude_deg: f64,
    /// Point longitude, degrees east (negative for west).
    pub longitude_deg: f64,
    /// Isobaric levels to request, in millibars.
    pub levels_mb: Vec<u32>,
    /// Compass direction of the trajectory's downrange axis, rotating the
    /// profile into the crate's frame convention.
    pub launch_azimuth_deg: f64,
}

impl NomadsRequest {
    /// The subset URL for the NOMADS filter script.
    pub fn url(&self) -> String {
        let Model::Gfs = self.model;
        let mut url = format!(
            "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.{}%2F{:02}%2Fatmos&file=gfs.t{:02}z.pgrb2.0p25.f{:03}",
            self.date_utc, self.cycle_hour_utc, self.cycle_hour_utc, self.forecast_hour
        );
        for level in &self.levels_mb {
            url.push_str(&format!("&lev_{level}_mb=on"));
        }
        url.push_str("&var_HGT=on&var_UGRD=on&var_VGRD=on&subregion=");
        url.push_str(&format!(
            "&leftlon={}&rightlon={}&toplat={}&bottomlat={}",
            self.longitude_deg - 1.5,
            self.longitude_deg + 1.5,
            self.latitude_deg + 1.5,
            self.latitude_deg - 1.5
        ));
        url
    }
}

/// Download the subset file and decode it into a profile at the requested
/// point.
pub fn fetch_profile(request: &NomadsRequest) -> Result<ProfileField> {
    let url = request.url();
    let response = ureq::get(&url)
        .call()
        .map_err(|e| Error::invalid("download", format!("request failed: {e}")))?;
    let mut bytes = Vec::new();
    response
        .into_body()
        .into_reader()
        .read_to_end(&mut bytes)
        .map_err(|e| Error::invalid("download", format!("read failed: {e}")))?;
    if !bytes.starts_with(b"GRIB") {
        return Err(Error::invalid(
            "download",
            "response is not GRIB2 (bad cycle or unavailable data)",
        ));
    }
    profile_from_grib2(
        &bytes,
        request.latitude_deg,
        request.longitude_deg,
        request.launch_azimuth_deg,
    )
}
