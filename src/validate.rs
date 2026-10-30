//! Compare a wind profile against a measured sounding.
//!
//! Metrics are sampled at the sounding's own levels: speed RMSE, direction
//! RMSE (circular), and mean absolute speed error per altitude layer
//! (1, 3, 6 and 10 km boundaries, then the rest). Layers bin on the
//! sounding's altitudes, which are above sea level, not above ground: for a
//! station at 287 m the first layer spans only 713 m above ground.

use crate::sounding::Sounding;
use crate::ProfileField;
use nalgebra::Vector3;

/// Horizontal speed of a wind vector, m/s.
fn horizontal_speed(w: Vector3<f64>) -> f64 {
    w.x.hypot(w.z)
}

/// Wind-FROM direction of a wind vector, degrees in [0, 360) (0° = headwind,
/// 90° = from the right).
fn direction_from_deg(w: Vector3<f64>) -> f64 {
    (w.z.atan2(w.x).to_degrees() + 180.0).rem_euclid(360.0)
}

/// Altitude-layer boundaries, metres above sea level (the sounding's datum,
/// not above ground): 0–1 km, 1–3 km, 3–6 km, 6–10 km, rest.
const LAYER_TOPS_M: [f64; 4] = [1_000.0, 3_000.0, 6_000.0, 10_000.0];

/// Mean absolute speed error within one altitude layer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LayerError {
    /// Top of the layer, metres above sea level (`f64::INFINITY` for the
    /// open layer).
    pub top_m: f64,
    /// Mean absolute speed error, m/s.
    pub mean_abs_speed_error_mps: f64,
    /// Sounding levels in the layer.
    pub levels: usize,
}

/// The comparison result.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidationReport {
    /// Sounding levels compared.
    pub levels: usize,
    /// RMS of (profile speed − measured speed), m/s.
    pub speed_rmse_mps: f64,
    /// RMS of the circular direction difference, degrees.
    pub direction_rmse_deg: f64,
    /// Per-layer mean absolute speed errors.
    pub layers: Vec<LayerError>,
}

/// Sample `profile` at every sounding level and compute the metrics. The
/// profile is read at the sounding's altitudes, above sea level, so its
/// levels must be in that datum too.
pub fn validate(profile: &ProfileField, sounding: &Sounding) -> ValidationReport {
    let mut speed_sq = 0.0;
    let mut direction_sq = 0.0;
    let mut layer_abs = vec![0.0; LAYER_TOPS_M.len() + 1];
    let mut layer_n = vec![0usize; LAYER_TOPS_M.len() + 1];
    for level in &sounding.levels {
        let wind = profile.wind_at(level.altitude_m);
        let speed_err = horizontal_speed(wind) - level.speed_mps;
        speed_sq += speed_err * speed_err;
        let mut dir_err = (direction_from_deg(wind) - level.direction_from_deg) % 360.0;
        if dir_err > 180.0 {
            dir_err -= 360.0;
        } else if dir_err < -180.0 {
            dir_err += 360.0;
        }
        direction_sq += dir_err * dir_err;
        let layer = LAYER_TOPS_M
            .iter()
            .position(|top| level.altitude_m < *top)
            .unwrap_or(LAYER_TOPS_M.len());
        layer_abs[layer] += speed_err.abs();
        layer_n[layer] += 1;
    }
    let n = sounding.levels.len().max(1) as f64;
    ValidationReport {
        levels: sounding.levels.len(),
        speed_rmse_mps: (speed_sq / n).sqrt(),
        direction_rmse_deg: (direction_sq / n).sqrt(),
        layers: layer_abs
            .into_iter()
            .zip(layer_n)
            .enumerate()
            .filter(|(_, (_, count))| *count > 0)
            .map(|(i, (abs, count))| LayerError {
                top_m: LAYER_TOPS_M.get(i).copied().unwrap_or(f64::INFINITY),
                mean_abs_speed_error_mps: abs / count as f64,
                levels: count,
            })
            .collect(),
    }
}
