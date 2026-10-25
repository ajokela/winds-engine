# winds-engine

Deterministic wind profiles and radiosonde validation for trajectory
simulation — the wind *data* pipeline for rocket modeling, consumed by
rockets-engine. Kernel-side wind modeling (the `WindField` trait, profile
types, analytic profiles) is staged here pending promotion to
ballistics-engine's `wind` module; see
[docs/winds-engine-plan.md](docs/winds-engine-plan.md).

```rust
use winds_engine::{ProfileField, ProfileLevel, Sounding, WindField, validate};

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
  level and stay that way in `to_profile_field()`. The log and Ekman
  profiles measure height from the ground, so in a sea-level frame they are
  read through ballistics-engine's `HeightOffset` with minus the ground
  elevation. Convert in one place only: shifting a sounding as well would
  move every wind by the site elevation.
- `validate`: speed/direction RMSE and per-layer errors vs a sounding.

Zero dependencies, no network, no wall clock, no unseeded randomness.
GRIB2/HRRR/GFS acquisition is the next milestone (risk item in the plan).

MIT licensed (see LICENSE).
