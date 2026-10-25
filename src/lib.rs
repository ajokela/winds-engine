//! Deterministic wind fields and radiosonde validation for trajectory
//! simulation, consumed by rockets-engine.
//!
//! Everything here is stateless and reproducible from inputs: no network, no
//! wall clock, no unseeded randomness. The [`WindField`] trait, the profile
//! types and the analytic profiles are the staging ground for the modeling
//! core planned for ballistics-engine's `wind` module (see
//! `docs/winds-engine-plan.md`); when that kernel release lands they move
//! there unchanged and this crate re-exports them.
//!
//! - [`field`]: the [`WindField`] trait, [`ConstantField`], [`ProfileField`].
//! - [`analytic`]: [`LogProfile`](analytic::LogProfile) and
//!   [`EkmanProfile`](analytic::EkmanProfile) boundary-layer models.
//! - [`sounding`]: radiosonde interchange-format parsing.
//! - [`validate`]: profile-vs-sounding comparison metrics.

pub mod analytic;
pub mod field;
pub mod sounding;
pub mod validate;

pub use field::{ConstantField, ProfileField, ProfileLevel, Wind, WindField};
pub use sounding::{Sounding, SoundingLevel};
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
