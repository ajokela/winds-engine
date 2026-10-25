//! Analytic boundary-layer profiles: logarithmic and Ekman spiral.

use crate::field::{Wind, WindField};
use crate::{Error, Result};

/// Log-law boundary layer: `u(z) = u_ref · ln(z/z0) / ln(z_ref/z0)` for
/// `z > z0`, zero below the roughness length `z0`. Direction is constant
/// with height.
#[derive(Debug, Clone, Copy)]
pub struct LogProfile {
    /// Wind speed at the reference height, m/s.
    pub reference_speed_mps: f64,
    /// Reference height, metres (above the roughness length).
    pub reference_height_m: f64,
    /// Aerodynamic roughness length z0, metres.
    pub roughness_m: f64,
    /// Meteorological direction the air comes FROM, degrees (constant with
    /// height).
    pub direction_from_deg: f64,
}

impl LogProfile {
    /// Validate the parameters.
    pub fn new(
        reference_speed_mps: f64,
        reference_height_m: f64,
        roughness_m: f64,
        direction_from_deg: f64,
    ) -> Result<Self> {
        let profile = Self {
            reference_speed_mps,
            reference_height_m,
            roughness_m,
            direction_from_deg,
        };
        profile.validate()?;
        Ok(profile)
    }

    fn validate(&self) -> Result<()> {
        if !self.reference_speed_mps.is_finite() || self.reference_speed_mps < 0.0 {
            return Err(Error::invalid(
                "reference_speed_mps",
                "must be a nonnegative number",
            ));
        }
        if !self.reference_height_m.is_finite() || self.reference_height_m <= 0.0 {
            return Err(Error::invalid("reference_height_m", "must be positive"));
        }
        if !self.roughness_m.is_finite()
            || self.roughness_m <= 0.0
            || self.roughness_m >= self.reference_height_m
        {
            return Err(Error::invalid(
                "roughness_m",
                "must be positive and below the reference height",
            ));
        }
        Ok(())
    }

    /// Wind speed at an altitude, m/s.
    pub fn speed_at(&self, altitude_m: f64) -> f64 {
        if altitude_m <= self.roughness_m {
            return 0.0;
        }
        let z0 = self.roughness_m;
        self.reference_speed_mps * (altitude_m / z0).ln() / (self.reference_height_m / z0).ln()
    }
}

impl WindField for LogProfile {
    fn wind(&self, _time_s: f64, position_m: [f64; 3]) -> Wind {
        Wind::from_speed_direction(self.speed_at(position_m[1]), self.direction_from_deg)
    }
}

/// Ekman spiral: the classic solution for a neutrally stratified boundary
/// layer with constant eddy viscosity `K`. At the ground the wind is zero;
/// aloft it approaches the geostrophic speed `G`. See Holton, *An
/// Introduction to Dynamic Meteorology*.
///
/// `u = G(1 − e^(−γz)·cos(γz))`, `v = sign(f)·G·e^(−γz)·sin(γz)` with
/// `γ = sqrt(|f| / 2K)` and Coriolis `f = 2Ω·sin(latitude)`. `u` runs along
/// the geostrophic flow and `v` to its LEFT, so in the northern hemisphere
/// the surface wind is backed about 45° from the geostrophic wind (turned
/// toward low pressure) and veers with height; the southern hemisphere
/// mirrors this. The result is rotated so the geostrophic direction matches
/// `direction_from_deg` aloft.
#[derive(Debug, Clone, Copy)]
pub struct EkmanProfile {
    /// Geostrophic wind speed aloft, m/s.
    pub geostrophic_mps: f64,
    /// Latitude, degrees (sign sets the turning sense).
    pub latitude_deg: f64,
    /// Constant eddy viscosity K, m²/s.
    pub eddy_viscosity_m2s: f64,
    /// Meteorological direction the geostrophic wind comes FROM aloft,
    /// degrees.
    pub direction_from_deg: f64,
}

impl EkmanProfile {
    /// Earth's rotation rate, rad/s.
    const OMEGA: f64 = 7.2921159e-5;

    /// Validate the parameters.
    pub fn new(
        geostrophic_mps: f64,
        latitude_deg: f64,
        eddy_viscosity_m2s: f64,
        direction_from_deg: f64,
    ) -> Result<Self> {
        let profile = Self {
            geostrophic_mps,
            latitude_deg,
            eddy_viscosity_m2s,
            direction_from_deg,
        };
        profile.validate()?;
        Ok(profile)
    }

    fn validate(&self) -> Result<()> {
        if !self.geostrophic_mps.is_finite() || self.geostrophic_mps < 0.0 {
            return Err(Error::invalid(
                "geostrophic_mps",
                "must be a nonnegative number",
            ));
        }
        if !self.latitude_deg.is_finite() || self.latitude_deg.abs() >= 90.0 {
            return Err(Error::invalid(
                "latitude_deg",
                "must be finite and away from the poles",
            ));
        }
        if !self.eddy_viscosity_m2s.is_finite() || self.eddy_viscosity_m2s <= 0.0 {
            return Err(Error::invalid("eddy_viscosity_m2s", "must be positive"));
        }
        Ok(())
    }

    fn coriolis(&self) -> f64 {
        2.0 * Self::OMEGA * self.latitude_deg.to_radians().sin()
    }

    /// Ekman components `(u, v)` in the geostrophic frame at an altitude.
    fn components_at(&self, altitude_m: f64) -> (f64, f64) {
        let f = self.coriolis();
        if f == 0.0 || altitude_m <= 0.0 {
            return (0.0, 0.0);
        }
        let gamma = (f.abs() / (2.0 * self.eddy_viscosity_m2s)).sqrt();
        let decay = (-gamma * altitude_m).exp();
        let g = self.geostrophic_mps;
        (
            g * (1.0 - decay * (gamma * altitude_m).cos()),
            f.signum() * g * decay * (gamma * altitude_m).sin(),
        )
    }

    /// Wind at an altitude, rotated so the geostrophic direction aloft
    /// matches `direction_from_deg`.
    pub fn wind_at(&self, altitude_m: f64) -> Wind {
        let (u, v) = self.components_at(altitude_m);
        // The solution's u points along the geostrophic flow, (cos rot, sin rot)
        // in (x, z), and v to its left. With +Y up and +Z to the right of +X,
        // left of the flow is (sin rot, -cos rot).
        let flow_dir = self.direction_from_deg + 180.0;
        let rot = flow_dir.to_radians();
        let x = u * rot.cos() + v * rot.sin();
        let z = u * rot.sin() - v * rot.cos();
        Wind {
            x_mps: x,
            y_mps: 0.0,
            z_mps: z,
        }
    }
}

impl WindField for EkmanProfile {
    fn wind(&self, _time_s: f64, position_m: [f64; 3]) -> Wind {
        self.wind_at(position_m[1])
    }
}
