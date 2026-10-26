//! NOAA NOMADS subset downloads, behind the `download` feature.
//!
//! This is the only networked module in the crate. The caller supplies the
//! cycle date and forecast hour explicitly — no wall clock here either, so
//! requests are reproducible.

use crate::field::ProfileField;
use crate::grib::profile_from_grib2;
use crate::{Error, Result};
use std::io::Read;

/// The forecast model to subset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Model {
    /// NOAA GFS at 0.25° (~25 km), global coverage, cycles every 6 h,
    /// via the NOMADS filter script (server-side subsetting).
    Gfs,
    /// NOAA HRRR at 3 km, CONUS coverage, hourly cycles, via byte-range
    /// requests against the AWS mirror's GRIB2 index files.
    Hrrr,
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
    /// The primary URL for the request's model: the NOMADS filter script
    /// for GFS, the AWS mirror's data file for HRRR (its index lives at
    /// [`Self::hrrr_urls`]).
    pub fn url(&self) -> String {
        match self.model {
            Model::Gfs => {
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
            Model::Hrrr => self.hrrr_urls().1,
        }
    }

    /// The (index, data) URLs for the HRRR AWS mirror.
    pub fn hrrr_urls(&self) -> (String, String) {
        let data = format!(
            "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.{}/conus/hrrr.t{:02}z.wrfprsf{:02}.grib2",
            self.date_utc, self.cycle_hour_utc, self.forecast_hour
        );
        (format!("{data}.idx"), data)
    }
}

/// Download the subset file and decode it into a profile at the requested
/// point.
pub fn fetch_profile(request: &NomadsRequest) -> Result<ProfileField> {
    let bytes = match request.model {
        Model::Gfs => http_get(&request.url(), None)?,
        Model::Hrrr => fetch_hrrr_messages(request)?,
    };
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

/// Byte ranges `(start, end inclusive)` of the messages a profile needs.
/// Public for testing and for callers doing their own fetching.
pub fn hrrr_parse_ranges(idx: &str, levels_mb: &[u32]) -> Result<Vec<(u64, Option<u64>)>> {
    let lines: Vec<&str> = idx.lines().collect();
    let mut ranges = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        let parts: Vec<&str> = line.split(':').collect();
        if parts.len() < 5 {
            continue;
        }
        if !matches!(parts[3], "HGT" | "UGRD" | "VGRD") {
            continue;
        }
        let Some(mb) = parts[4]
            .strip_suffix(" mb")
            .and_then(|s| s.trim().parse::<u32>().ok())
        else {
            continue;
        };
        if !levels_mb.contains(&mb) {
            continue;
        }
        let start: u64 = parts[1]
            .parse()
            .map_err(|_| Error::invalid("hrrr index", "non-numeric offset"))?;
        let end = lines
            .get(i + 1)
            .and_then(|next| next.split(':').nth(1))
            .and_then(|s| s.parse::<u64>().ok())
            .map(|next| next - 1);
        ranges.push((start, end));
    }
    if ranges.is_empty() {
        return Err(Error::invalid(
            "hrrr index",
            "no HGT/UGRD/VGRD messages at the requested levels",
        ));
    }
    Ok(ranges)
}

fn fetch_hrrr_messages(request: &NomadsRequest) -> Result<Vec<u8>> {
    let (idx_url, data_url) = request.hrrr_urls();
    let idx = http_get(&idx_url, None)?;
    let idx =
        String::from_utf8(idx).map_err(|_| Error::invalid("hrrr index", "index is not text"))?;
    let ranges = hrrr_parse_ranges(&idx, &request.levels_mb)?;
    let mut bytes = Vec::new();
    for (start, end) in ranges {
        let range = match end {
            Some(end) => format!("bytes={start}-{end}"),
            None => format!("bytes={start}-"),
        };
        bytes.extend(http_get(&data_url, Some(&range))?);
    }
    Ok(bytes)
}

fn http_get(url: &str, range: Option<&str>) -> Result<Vec<u8>> {
    let config = ureq::Agent::config_builder()
        .timeout_global(Some(std::time::Duration::from_secs(60)))
        .build();
    let agent = ureq::Agent::new_with_config(config);
    let mut call = agent.get(url);
    if let Some(range) = range {
        call = call.header("Range", range);
    }
    let response = call
        .call()
        .map_err(|e| Error::invalid("download", format!("request failed: {e}")))?;
    let mut bytes = Vec::new();
    response
        .into_body()
        .into_reader()
        .read_to_end(&mut bytes)
        .map_err(|e| Error::invalid("download", format!("read failed: {e}")))?;
    Ok(bytes)
}
