//! GRIB2 forecast grids to wind profiles, via the `grib` crate.
//!
//! Decoded offline and deterministic: a byte slice in, a [`ProfileField`]
//! out. Acquisition (downloading those bytes) lives behind the `download`
//! feature in [`crate::download`].

use crate::{Error, ProfileField, ProfileLevel, Result};
use grib::{Grib2SubmessageDecoder, LatLons};

/// GRIB2 parameter (discipline, category, number) for the variables a wind
/// profile needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Variable {
    /// Geopotential height, gpm (3/5).
    Height,
    /// U-component of wind, m/s (2/2).
    Ugrd,
    /// V-component of wind, m/s (2/3).
    Vgrd,
}

impl Variable {
    fn of(category: Option<u8>, number: Option<u8>) -> Option<Self> {
        match (category, number) {
            (Some(3), Some(5)) => Some(Self::Height),
            (Some(2), Some(2)) => Some(Self::Ugrd),
            (Some(2), Some(3)) => Some(Self::Vgrd),
            _ => None,
        }
    }
}

/// Quantize a coordinate for map keys (1e-5° buckets).
fn quantize(deg: f64) -> i64 {
    (deg * 100_000.0).round() as i64
}

/// Bilinear interpolation over the enclosing grid cell, falling back to the
/// nearest point when the target sits on the grid edge or outside it.
fn bilinear_sample(
    latlons: &[(f64, f64)],
    values: &[f32],
    latitude_deg: f64,
    longitude_deg: f64,
) -> Result<f32> {
    use std::collections::HashMap;
    // Only points near the target matter; prefilter to a small window so
    // million-point grids stay cheap.
    const WINDOW_DEG: f64 = 0.6;
    let mut by_point: HashMap<(i64, i64), f32> = HashMap::new();
    let mut lats: Vec<f64> = Vec::new();
    let mut lons: Vec<f64> = Vec::new();
    for (i, &(lat, lon)) in latlons.iter().enumerate() {
        let mut dlon = (lon - longitude_deg).abs();
        if dlon > 180.0 {
            dlon = 360.0 - dlon;
        }
        if (lat - latitude_deg).abs() > WINDOW_DEG || dlon > WINDOW_DEG {
            continue;
        }
        by_point.insert((quantize(lat), quantize(lon)), values[i]);
        lats.push(lat);
        lons.push(lon);
    }
    lats.sort_by(f64::total_cmp);
    lats.dedup();
    lons.sort_by(f64::total_cmp);
    lons.dedup();
    let cell = |axis: &[f64], x: f64| -> Option<(usize, usize, f64)> {
        if x < axis[0] || x > axis[axis.len() - 1] {
            return None;
        }
        let hi = axis
            .iter()
            .position(|&v| v >= x)
            .unwrap_or(axis.len() - 1)
            .max(1);
        let (lo_v, hi_v) = (axis[hi - 1], axis[hi]);
        let f = if hi_v > lo_v {
            (x - lo_v) / (hi_v - lo_v)
        } else {
            0.0
        };
        Some((hi - 1, hi, f))
    };
    let nearest = |by_point: &HashMap<(i64, i64), f32>| -> Result<f32> {
        if by_point.is_empty() {
            // Target outside the window: full nearest scan (rare path).
            let mut best = f32::NAN;
            let mut best_dist = f64::INFINITY;
            for (i, &(lat, lon)) in latlons.iter().enumerate() {
                let dlat = lat - latitude_deg;
                let mut dlon = (lon - longitude_deg).abs();
                if dlon > 180.0 {
                    dlon = 360.0 - dlon;
                }
                let dist = dlat * dlat + dlon * dlon;
                if dist < best_dist {
                    best_dist = dist;
                    best = values[i];
                }
            }
            if best_dist.is_infinite() {
                return Err(Error::invalid("grib2", "grid has no points"));
            }
            return Ok(best);
        }
        // Nearest within the window (grid edge: no enclosing cell).
        let mut best = f32::NAN;
        let mut best_dist = f64::INFINITY;
        for &(lat, lon) in latlons {
            let key = (quantize(lat), quantize(lon));
            let Some(&value) = by_point.get(&key) else {
                continue;
            };
            let dlat = lat - latitude_deg;
            let mut dlon = (lon - longitude_deg).abs();
            if dlon > 180.0 {
                dlon = 360.0 - dlon;
            }
            let dist = dlat * dlat + dlon * dlon;
            if dist < best_dist {
                best_dist = dist;
                best = value;
            }
        }
        if best_dist.is_infinite() {
            return Err(Error::invalid("grib2", "grid has no points"));
        }
        Ok(best)
    };
    let (Some((j0, j1, fj)), Some((i0, i1, fi))) =
        (cell(&lats, latitude_deg), cell(&lons, longitude_deg))
    else {
        return nearest(&by_point);
    };
    let corner = |j: usize, i: usize| {
        by_point
            .get(&(quantize(lats[j]), quantize(lons[i])))
            .copied()
    };
    let (Some(v00), Some(v01), Some(v10), Some(v11)) = (
        corner(j0, i0),
        corner(j0, i1),
        corner(j1, i0),
        corner(j1, i1),
    ) else {
        return nearest(&by_point);
    };
    Ok(v00 * (1.0 - fi) as f32 * (1.0 - fj) as f32
        + v01 * fi as f32 * (1.0 - fj) as f32
        + v10 * (1.0 - fi) as f32 * fj as f32
        + v11 * fi as f32 * fj as f32)
}

/// Decode GRIB2 messages with geopotential height, U and V components on
/// isobaric surfaces into a [`ProfileField`] at one point.
///
/// Sampling is bilinear over the enclosing grid cell, falling back to the
/// nearest point at the grid edges. Heights come from the geopotential-height
/// field, so the profile is on the same metres-above-sea-level datum as the
/// rest of the crate.
///
/// `launch_azimuth_deg` is the compass direction of the downrange +X axis,
/// converting meteorological FROM-directions (clockwise from north) into the
/// crate's frame convention (0° = headwind, 90° = from the right).
pub fn profile_from_grib2(
    data: &[u8],
    latitude_deg: f64,
    longitude_deg: f64,
    launch_azimuth_deg: f64,
) -> Result<ProfileField> {
    if !latitude_deg.is_finite() || latitude_deg.abs() > 90.0 {
        return Err(Error::invalid(
            "latitude_deg",
            "must be finite within [-90, 90]",
        ));
    }
    if !longitude_deg.is_finite() || !launch_azimuth_deg.is_finite() {
        return Err(Error::invalid("longitude_deg", "must be finite"));
    }
    let grib2 = grib::from_bytes(data)
        .map_err(|e| Error::invalid("grib2", format!("could not decode: {e}")))?;

    let mut levels: Vec<(i32, Variable, f64, f64, f64)> = Vec::new();
    for (_index, submessage) in grib2.iter() {
        let prod = submessage.prod_def();
        let Some(variable) = Variable::of(prod.parameter_category(), prod.parameter_number())
        else {
            continue;
        };
        let Some((surface, _)) = prod.fixed_surfaces() else {
            continue;
        };
        // 100 = isobaric surface, value in Pa (Code table 4.5).
        if surface.surface_type != 100 {
            continue;
        }
        let latlons: Vec<(f64, f64)> = submessage
            .latlons()
            .map_err(|e| Error::invalid("grib2", format!("grid has no lat/lons: {e:?}")))?
            .collect();
        let values: Vec<f32> = Grib2SubmessageDecoder::from(submessage)
            .map_err(|e| Error::invalid("grib2", format!("cannot decode values: {e}")))?
            .dispatch()
            .map_err(|e| Error::invalid("grib2", format!("cannot decode values: {e}")))?
            .collect();
        let value = bilinear_sample(&latlons, &values, latitude_deg, longitude_deg)?;
        levels.push((surface.scaled_value, variable, f64::from(value), 0.0, 0.0));
    }

    // Assemble per-level samples: (pressure Pa, height, u, v).
    type LevelParts = (i32, Option<f64>, Option<f64>, Option<f64>);
    let mut assembled: Vec<LevelParts> = Vec::new();
    for (level_pa, variable, value, _lat, _lon) in levels {
        let entry = match assembled.iter_mut().find(|(l, _, _, _)| *l == level_pa) {
            Some(entry) => entry,
            None => {
                assembled.push((level_pa, None, None, None));
                assembled.last_mut().expect("just pushed")
            }
        };
        match variable {
            Variable::Height => entry.1 = Some(value),
            Variable::Ugrd => entry.2 = Some(value),
            Variable::Vgrd => entry.3 = Some(value),
        }
    }

    let mut profile_levels = Vec::new();
    for (_pa, height, u, v) in assembled {
        let (Some(altitude_m), Some(u_mps), Some(v_mps)) = (height, u, v) else {
            return Err(Error::invalid(
                "grib2",
                "levels lack geopotential height, UGRD or VGRD",
            ));
        };
        let speed_mps = u_mps.hypot(v_mps);
        // Meteorological FROM (clockwise from north) rotated into the frame:
        // FROM = atan2(-u, -v), then subtract the launch azimuth.
        let from_met = (-u_mps).atan2(-v_mps).to_degrees().rem_euclid(360.0);
        let direction_from_deg = (from_met - launch_azimuth_deg).rem_euclid(360.0);
        profile_levels.push(ProfileLevel {
            altitude_m,
            speed_mps,
            direction_from_deg,
        });
    }

    ProfileField::new(profile_levels).map_err(|e| Error::invalid("grib2", e.to_string()))
}
