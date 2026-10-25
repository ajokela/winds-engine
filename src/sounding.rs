//! Radiosonde soundings in the crate's line-based interchange format.
//!
//! The format is deliberately trivial so any archive (IEM RAOB JSON, UWyo
//! text lists, SHARPpy exports) can be converted offline with a one-liner:
//!
//! ```text
//! # station: KMPX
//! # valid: 2026-10-23T00:00:00Z
//! # columns: altitude_m speed_mps direction_from_deg
//! 287.0 0.0000 0.0
//! 305.0 3.0867 125.0
//! ```
//!
//! Blank lines and `#` comments are ignored; `# station:` and `# valid:`
//! header keys are carried along. Directions are meteorological FROM (see
//! [`Wind::from_speed_direction`](crate::Wind::from_speed_direction)).

use crate::field::{ProfileField, ProfileLevel};
use crate::{Error, Result};

/// One measured level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoundingLevel {
    /// Height above mean sea level, metres.
    pub altitude_m: f64,
    /// Measured wind speed, m/s.
    pub speed_mps: f64,
    /// Measured direction the air comes FROM, degrees.
    pub direction_from_deg: f64,
}

/// A parsed sounding.
#[derive(Debug, Clone, PartialEq)]
pub struct Sounding {
    /// Station identifier, if the header carried one.
    pub station: Option<String>,
    /// Observation time string, if the header carried one (kept opaque:
    /// this crate never parses clocks).
    pub valid: Option<String>,
    /// Measured levels in file order.
    pub levels: Vec<SoundingLevel>,
}

impl Sounding {
    /// Parse the interchange format.
    pub fn parse(text: &str) -> Result<Self> {
        let mut station = None;
        let mut valid = None;
        let mut levels = Vec::new();
        for (number, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(comment) = line.strip_prefix('#') {
                let comment = comment.trim();
                if let Some(value) = comment.strip_prefix("station:") {
                    station = Some(value.trim().to_string());
                } else if let Some(value) = comment.strip_prefix("valid:") {
                    valid = Some(value.trim().to_string());
                }
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() != 3 {
                return Err(Error::invalid(
                    "sounding",
                    format!(
                        "line {}: expected 'altitude_m speed_mps direction_from_deg'",
                        number + 1
                    ),
                ));
            }
            let parse = |token: &str| -> Result<f64> {
                token.parse().map_err(|_| {
                    Error::invalid(
                        "sounding",
                        format!("line {}: '{token}' is not numeric", number + 1),
                    )
                })
            };
            levels.push(SoundingLevel {
                altitude_m: parse(fields[0])?,
                speed_mps: parse(fields[1])?,
                direction_from_deg: parse(fields[2])?,
            });
        }
        if levels.is_empty() {
            return Err(Error::invalid("sounding", "no levels"));
        }
        Ok(Self {
            station,
            valid,
            levels,
        })
    }

    /// The lowest measured altitude, metres above mean sea level. In a standard sounding this
    /// is the station's surface report, so it is the station's elevation. `None` when there
    /// are no levels.
    pub fn surface_altitude_m(&self) -> Option<f64> {
        self.levels.iter().map(|l| l.altitude_m).reduce(f64::min)
    }

    /// The levels as a validated [`ProfileField`] (sorted and checked). Heights stay in metres
    /// above mean sea level, as measured, so the field matches a frame whose y is altitude above
    /// sea level (rockets-engine's) with no conversion, wherever that frame's ground sits.
    ///
    /// Convert datums in one place only, never in these levels as well. This crate's staged
    /// types have no offset wrapper: ballistics-engine's `HeightOffset` (from its next release)
    /// wraps only ballistics-engine's own `WindField`, not this crate's (the two crates' wind
    /// types differ; see the [crate docs](crate)). So a caller in a sea-level frame reads a
    /// ground-relative field such as [`LogProfile`](crate::analytic::LogProfile) or
    /// [`EkmanProfile`](crate::analytic::EkmanProfile) at `y - ground elevation` at its single
    /// call site, or builds ballistics-engine's `ProfileField`/`LogProfile`/`EkmanProfile` from
    /// the same values and wraps those in `HeightOffset`. This field, being above sea level,
    /// needs no offset in a sea-level frame.
    pub fn to_profile_field(&self) -> Result<ProfileField> {
        ProfileField::new(
            self.levels
                .iter()
                .map(|l| ProfileLevel {
                    altitude_m: l.altitude_m,
                    speed_mps: l.speed_mps,
                    direction_from_deg: l.direction_from_deg,
                })
                .collect(),
        )
    }
}
