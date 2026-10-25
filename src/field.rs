//! The [`WindField`] trait and the concrete field types.

use crate::{Error, Result};

/// World-frame air velocity, m/s: +X is downrange (tailwind), +Y is up, +Z
/// is the horizontal cross axis to the right of +X (a right-handed frame,
/// ballistics-engine's McCoy frame). Vertical wind is zero for the horizontal
/// field types in this crate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wind {
    /// Downrange component, m/s.
    pub x_mps: f64,
    /// Vertical component, m/s.
    pub y_mps: f64,
    /// Cross-range component, m/s.
    pub z_mps: f64,
}

impl Wind {
    /// Calm air.
    pub const ZERO: Self = Self {
        x_mps: 0.0,
        y_mps: 0.0,
        z_mps: 0.0,
    };

    /// Horizontal wind from a speed and a meteorological direction: the
    /// direction the air comes FROM, degrees clockwise, where 0° means from
    /// downrange (+X, so the wind blows toward −X, a headwind) and 90° means
    /// from +Z. Vertical component is zero.
    pub fn from_speed_direction(speed_mps: f64, direction_from_deg: f64) -> Self {
        let toward = (direction_from_deg + 180.0).to_radians();
        Self {
            x_mps: speed_mps * toward.cos(),
            y_mps: 0.0,
            z_mps: speed_mps * toward.sin(),
        }
    }

    /// Horizontal speed.
    pub fn speed_mps(&self) -> f64 {
        self.x_mps.hypot(self.z_mps)
    }

    /// Meteorological direction the air comes FROM, degrees in [0, 360).
    pub fn direction_from_deg(&self) -> f64 {
        let toward = self.z_mps.atan2(self.x_mps).to_degrees() + 180.0;
        toward.rem_euclid(360.0)
    }
}

/// A deterministic, stateless wind field: same time and position always
/// give the same wind, so integrator substeps, snapshots and seeded
/// reproduction are exact.
pub trait WindField {
    /// World-frame air velocity at a time (seconds) and position
    /// (`[x, y, z]` metres, Y up).
    fn wind(&self, time_s: f64, position_m: [f64; 3]) -> Wind;
}

/// Uniform wind everywhere.
#[derive(Debug, Clone, Copy)]
pub struct ConstantField {
    /// The constant air velocity.
    pub wind: Wind,
}

impl WindField for ConstantField {
    fn wind(&self, _time_s: f64, _position_m: [f64; 3]) -> Wind {
        self.wind
    }
}

/// One level of a vertical profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileLevel {
    /// Height above the profile datum, metres.
    pub altitude_m: f64,
    /// Horizontal wind speed, m/s (nonnegative).
    pub speed_mps: f64,
    /// Meteorological direction the air comes FROM, degrees (see
    /// [`Wind::from_speed_direction`]).
    pub direction_from_deg: f64,
}

/// Piecewise-linear vertical profile: speed interpolates linearly in
/// altitude, direction along the shortest circular arc. Below the first
/// level and above the last, the nearest level holds.
#[derive(Debug, Clone)]
pub struct ProfileField {
    levels: Vec<ProfileLevel>,
}

impl ProfileField {
    /// Build from levels, which must be finite, nonnegative in speed and
    /// strictly ascending in altitude.
    pub fn new(mut levels: Vec<ProfileLevel>) -> Result<Self> {
        if levels.is_empty() {
            return Err(Error::invalid("profile", "needs at least one level"));
        }
        levels.sort_by(|a, b| a.altitude_m.total_cmp(&b.altitude_m));
        for pair in levels.windows(2) {
            if pair[0].altitude_m == pair[1].altitude_m {
                return Err(Error::invalid(
                    "profile",
                    "altitudes must be strictly ascending",
                ));
            }
        }
        for level in &levels {
            if !level.altitude_m.is_finite()
                || !level.speed_mps.is_finite()
                || !level.direction_from_deg.is_finite()
            {
                return Err(Error::invalid("profile", "levels must be finite"));
            }
            if level.speed_mps < 0.0 {
                return Err(Error::invalid("profile", "speeds must be nonnegative"));
            }
        }
        Ok(Self { levels })
    }

    /// The sorted levels.
    pub fn levels(&self) -> &[ProfileLevel] {
        &self.levels
    }

    /// Wind at an altitude, with circular direction interpolation.
    pub fn wind_at(&self, altitude_m: f64) -> Wind {
        let levels = &self.levels;
        if altitude_m <= levels[0].altitude_m {
            let l = levels[0];
            return Wind::from_speed_direction(l.speed_mps, l.direction_from_deg);
        }
        for pair in levels.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            if altitude_m <= hi.altitude_m {
                let f = (altitude_m - lo.altitude_m) / (hi.altitude_m - lo.altitude_m);
                let speed = lo.speed_mps + f * (hi.speed_mps - lo.speed_mps);
                let mut delta = (hi.direction_from_deg - lo.direction_from_deg) % 360.0;
                if delta > 180.0 {
                    delta -= 360.0;
                } else if delta < -180.0 {
                    delta += 360.0;
                }
                let direction = lo.direction_from_deg + f * delta;
                return Wind::from_speed_direction(speed, direction);
            }
        }
        let l = levels[levels.len() - 1];
        Wind::from_speed_direction(l.speed_mps, l.direction_from_deg)
    }
}

impl WindField for ProfileField {
    fn wind(&self, _time_s: f64, position_m: [f64; 3]) -> Wind {
        self.wind_at(position_m[1])
    }
}
