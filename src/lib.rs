//! Deterministic wind data for trajectory simulation, consumed by
//! rockets-engine.
//!
//! The wind-field types — the [`WindField`] trait, [`ConstantField`],
//! [`ProfileField`], [`LogProfile`], [`EkmanProfile`] and the [`HeightOffset`]
//! datum wrapper — are ballistics-engine's (its `wind` module, since
//! ballistics-engine 0.49.0), re-exported here so this crate's pipelines and
//! its consumers share one definition. Everything is stateless and
//! reproducible from inputs: no wall clock or unseeded randomness anywhere,
//! and no network except the `download` feature.
//!
//! - [`sounding`]: radiosonde interchange-format parsing.
//! - [`validate`]: profile-vs-sounding comparison metrics.
//! - [`terrain`]: a diagnostic mass-consistent terrain wind.
//! - `grib` (feature): GRIB2 forecast files into a [`ProfileField`].
//! - `download` (feature): live GFS/HRRR acquisition, the only networked
//!   code in the crate.

#[cfg(feature = "download")]
pub mod download;
#[cfg(feature = "grib")]
pub mod grib;
pub mod sounding;
pub mod terrain;
pub mod validate;

pub use ballistics_engine::wind::{
    ConstantField, EkmanProfile, HeightOffset, LogProfile, ProfileField, ProfileLevel, WindField,
    WindFieldError, WindFieldParam, WindFieldRule,
};
pub use sounding::{Sounding, SoundingLevel};
pub use terrain::{local_en_m, monitors_from_geo_csv, HeightGrid, Monitor, TerrainFlowField};
pub use validate::{LayerError, ValidationReport};

/// Crate-wide error: invalid input with a field name and reason.
#[derive(Debug, Clone, PartialEq)]
pub struct Error {
    /// The offending input or field.
    pub field: &'static str,
    /// Why it was rejected.
    pub reason: String,
}

impl Error {
    pub(crate) fn invalid(field: &'static str, reason: impl Into<String>) -> Self {
        Self {
            field,
            reason: reason.into(),
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.field, self.reason)
    }
}

impl std::error::Error for Error {}

/// Crate-wide result.
pub type Result<T> = std::result::Result<T, Error>;
