# winds-engine

Deterministic wind profiles and radiosonde validation for trajectory
simulation — the wind *data* pipeline for rocket modeling, consumed by
rockets-engine. Kernel-side wind modeling (the `WindField` trait, profile
types, analytic profiles) is staged here pending promotion to
ballistics-engine's `wind` module; see
[docs/winds-engine-plan.md](docs/winds-engine-plan.md).

```rust
use winds_engine::{Sounding, WindField, validate};

// A measured sounding in the line-based interchange format. Heights are
// metres above sea level; this station's surface report is at 287 m.
let sounding = Sounding::parse("287 0.0 0\n305 3.09 125\n610 3.60 145\n")?;
// The profile keeps those heights, so it is read with y = height above sea
// level, the frame rockets-engine uses.
let profile = sounding.to_profile_field()?;
let wind = profile.wind(0.0, [0.0, 500.0, 0.0]); // 500 m above sea level
// Or score a candidate profile against the measurement.
let report = validate::validate(&profile, &sounding);
```

- `WindField` trait: stateless, deterministic `wind(time, position)` — the
  same sampling pattern trajectory integrators already rely on.
- `ProfileField` (piecewise-linear, circular direction interpolation),
  `ConstantField`, `LogProfile`, `EkmanProfile`.
- `Sounding`: radiosonde interchange format (tested against a real
  115-level KMPX sounding from the IEM RAOB archive). Heights are above sea
  level and stay that way in `to_profile_field()`, so the profile reads
  directly in a frame whose y is altitude above sea level, as in
  rockets-engine (with `BallisticForces` and in its flat-frame CLI, whose
  ground sits at `--ground-altitude-m`, sea level by default).
- Datum conversion: `LogProfile` and `EkmanProfile` take y as height above
  the ground, where they are calm. This crate's staged types have no offset
  wrapper; ballistics-engine's `HeightOffset` (from its next release, which
  rockets-engine's `=0.48.0` pin does not have) wraps only
  ballistics-engine's own `WindField`, not this crate's. So a caller in a
  sea-level frame either reads them at `y - ground elevation` at its single
  call site, or builds
  ballistics-engine's `ProfileField`/`LogProfile`/`EkmanProfile` from the
  same values and wraps those in `HeightOffset` (offset = the height of the
  caller's y = 0 in the field's datum: minus the ground elevation for a log
  or Ekman profile, none for a sea-level sounding). Convert in one place
  only: sounding heights stay above sea level, and shifting them as well
  would move every wind by the site elevation.
- `validate`: speed/direction RMSE and per-layer errors vs a sounding.
- Relationship to ballistics-engine: these types are this crate's own.
  ballistics-engine's wind fields (unreleased) were adapted from them and
  differ: `WindField::wind` takes and returns `nalgebra::Vector3<f64>` and
  there is no `Wind` type, `ConstantField` is a tuple struct over that
  vector, `LogProfile`/`EkmanProfile` keep their fields private behind
  getters, and constructors return `WindFieldError`. Moving winds-engine
  onto them changes this crate's public API (the example above,
  `validate()`'s signature), so it is a deliberate breaking release, not a
  drop-in re-export.

Zero dependencies, no network, no wall clock, no unseeded randomness.
GRIB2/HRRR/GFS acquisition is the next milestone (risk item in the plan).

MIT licensed (see LICENSE).
