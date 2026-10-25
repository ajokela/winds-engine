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
// A frame whose y = 0 is the ground at a launch site 287 m above sea level:
// convert the heights once, using the site's elevation (not the station's,
// when the balloon went up elsewhere).
let profile = sounding.to_profile_field_above(287.0)?;
let wind = profile.wind(0.0, [0.0, 500.0, 0.0]); // 500 m above the pad
// Score a candidate profile against the measurement on its own sea-level
// heights.
let report = validate::validate(&sounding.to_profile_field()?, &sounding);
```

- `WindField` trait: stateless, deterministic `wind(time, position)` — the
  same sampling pattern trajectory integrators already rely on.
- `ProfileField` (piecewise-linear, circular direction interpolation),
  `ConstantField`, `LogProfile`, `EkmanProfile`.
- `Sounding`: radiosonde interchange format (tested against a real
  115-level KMPX sounding from the IEM RAOB archive). Heights are above sea
  level; `to_profile_field_above(ground_msl_m)` gives heights above a site's
  ground, which is what a ground-relative frame (a launch pad, the log and
  Ekman profiles) needs.
- `validate`: speed/direction RMSE and per-layer errors vs a sounding.

Zero dependencies, no network, no wall clock, no unseeded randomness.
GRIB2/HRRR/GFS acquisition is the next milestone (risk item in the plan).

MIT licensed (see LICENSE).
