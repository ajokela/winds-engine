# winds-engine

Deterministic wind profiles and radiosonde validation for trajectory
simulation — the wind *data* pipeline for rocket modeling, consumed by
rockets-engine. The wind-field types themselves (the `WindField` trait, the
profile types, the analytic profiles, `HeightOffset`) are ballistics-engine's
`wind` module, re-exported here so the pipelines and their consumers share
one definition; see [docs/winds-engine-plan.md](docs/winds-engine-plan.md).

```rust
use nalgebra::Vector3;
use winds_engine::{Sounding, WindField, validate};

// A measured sounding in the line-based interchange format. Heights are
// metres above sea level; this station's surface report is at 287 m.
let sounding = Sounding::parse("287 0.0 0\n305 3.09 125\n610 3.60 145\n")?;
// The profile keeps those heights, so it is read with y = height above sea
// level, the frame rockets-engine uses.
let profile = sounding.to_profile_field()?;
let wind = profile.wind(0.0, Vector3::new(0.0, 500.0, 0.0)); // 500 m above sea level
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
  the ground, where they are calm. Wrap them in `HeightOffset` (re-exported)
  to read them in another frame: the offset is the height of the caller's
  y = 0 in the field's datum — minus the ground elevation in a sea-level
  frame, nothing for a sea-level sounding. Convert in one place only:
  sounding heights stay above sea level, and shifting them as well would
  move every wind by the site elevation.
- `validate`: speed/direction RMSE and per-layer errors vs a sounding.
- Relationship to ballistics-engine: winds-engine depends on it and
  re-exports its wind-field core (since 0.6.0, against ballistics-engine
  0.49.0). The kernel owns gravity, drag, the ICAO atmosphere and the
  integrator; winds-engine owns the wind *data*: sounding/profile parsing,
  forecast ingestion and terrain flow. They meet in rockets-engine, which
  depends on both and samples a `WindField` through ballistics-engine's
  force kernel.

- `grib` (feature): decode GRIB2 forecast files (GFS isobaric HGT/UGRD/VGRD)
  into a `ProfileField` at a point — bilinear-interpolated, offline and
  deterministic.
- `download` (feature): fetch live forecasts — GFS via the NOAA NOMADS
  filter script (server-side subsetting) or HRRR's 3 km CONUS grids via
  byte-range requests against the AWS mirror's index files — the only
  networked code in the crate.

- `terrain`: a diagnostic mass-consistent terrain wind — a `HeightGrid`
  (ASCII DEM + raw arrays) plus `TerrainFlowField`, which adjusts a base
  field to satisfy 2-D horizontal mass conservation over terrain,
  constrained by surface wind `Monitor`s (one per reading; `from_geo` and
  `monitors_from_geo_csv` project GPS lat/lon and compass directions into
  the field's frame).

The default build depends on ballistics-engine (no default features) and
nalgebra; `grib` adds the `grib` crate, `download` adds `ureq`. No wall
clock or unseeded randomness anywhere — cycle date and forecast hour are
explicit parameters.

MIT licensed (see LICENSE).
