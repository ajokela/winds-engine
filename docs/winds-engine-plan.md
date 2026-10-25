# winds-engine: go/no-go analysis and plan

Status: revised 2026-10-24 (architecture decision recorded below)
Consumers: rockets-engine (trajectories, weathercock, dispersion), moving-batteries (artillery, via rockets-engine)

## Decision (2026-10-24)

Wind **modeling** is upgraded *inside* ballistics-engine's existing `wind`
module — **additively**: `WindSock`/`WindSegment` and every existing API stay
exactly as they are. New trait and profile types are pure additions; nothing
is replaced, no API expectations change.

**winds-engine is explicitly scoped to rocket modeling**: data acquisition
(GRIB2/HRRR/GFS ingestion), radiosonde validation, and feeding the kernel's
profile types to rockets-engine. It is not a general ballistics-ecosystem
wind crate and will not grow kernel-side modeling.

Rationale: the kernel's modeling pieces are small, pure physics, and belong
next to the atmosphere tables — but the ballistics-engine release train
(platform binaries, npm, PyPI, RubyGems, WASM terminal, docs) is too heavy
to couple to data-pipeline iteration. One deliberate kernel minor release
for the stable modeling core; everything that iterates lives here.

## What goes into ballistics-engine (one additive kernel release)

All in the existing `wind` module, alongside the untouched `WindSock`:

- `WindField` trait: `wind(time_s, position_m) -> Vec3` — stateless and
  deterministic, the same sampling pattern the rockets-engine turbulence
  work already relies on (RK4 substeps, snapshots, seed reproduction stay
  bit-exact).
- `ConstantField(Vec3)` — today's uniform vector as a field.
- `ProfileField` — piecewise-linear vertical profile
  `Vec<(altitude_m, speed_mps, direction_deg)>`.
- Analytic profiles: `LogProfile` (roughness boundary layer) and
  `EkmanProfile` (geostrophic speed, latitude, a constant eddy viscosity
  and the geostrophic wind's direction), consistent with the kernel's
  atmosphere tables.
- An adapter so `WindSock` satisfies `WindField` — old and new interop;
  no consumer changes required.

rockets-engine then gains an additive `BallisticForces::with_wind_field`
(keeping `with_wind_segments`), and the weathercock samples the same field
(the turbulence wiring already proved that path).

## What stays in winds-engine (rocket-modeling only)

1. GRIB2 ingestion (HRRR/GFS u/v on isobaric levels) → kernel
   `ProfileField` at a launch coordinate: bilinear horizontal
   interpolation, altitude conversion to heights above sea level
   (geopotential height of each isobaric level to metres above sea level,
   the same datum as radiosonde soundings; any ground-relative conversion
   is left to ballistics-engine's `HeightOffset` at the consumer, in one
   place only), offline cache, NOAA NOMADS subset
   HTTP as the acquisition path (feature-gated `download`, off by default
   so the crate stays embeddable/WASM-safe).
2. Radiosonde validation datasets (University of Wyoming sounding
   archives): layer-mean speed/direction tolerances on a handful of
   station-days.
3. rockets-engine-facing glue examples: `--wind-profile`, dispersion
   profile perturbation with per-run seeding.

### GRIB2 risk item

1. Evaluate existing Rust GRIB crates for HRRR/GFS u/v-on-isobaric records.
2. Fall back to the NOAA NOMADS subset API (server-side extraction, small
   downloads) with a checked-in golden file for tests.
3. Last resort: a tiny GRIB2 template decoder for exactly the two products
   needed — explicitly not a general GRIB library.

## Explicitly out

- Tier C (WindNinja clone: FEM solver + weather init + GUI) — permanently.
- Diagnostic mass-consistent terrain flow (formerly Tier B): if ever
  needed, its home is the kernel, not this crate; revisit only with a
  concrete moving-batteries requirement.
- Replacing or wrapping existing ballistics-engine wind APIs.

## Determinism and dependency policy

- Kernel additions: no network, no wall-clock, no unseeded randomness;
  everything reproducible from inputs.
- winds-engine: acquisition behind the non-default `download` feature;
  the core decode/validate path is offline and deterministic.
- No AI/assistant attribution in commits; commit as the human developer.

## Effort and order

| Item | Where | Effort |
| --- | --- | --- |
| `WindField` trait + constant/profile/analytic fields + `WindSock` adapter | ballistics-engine `wind` module | 3–4 days |
| Kernel minor release (run the full train, once) | ballistics-engine | 1 day |
| rockets-engine `with_wind_field` + `--wind-profile` + weathercock wiring | rockets-engine | 2–3 days |
| Radiosonde validation set | winds-engine | 2–3 days |
| GRIB acquisition + decode (one product) | winds-engine | 1–2 weeks (risk item) |

First milestone: kernel `WindField` release + rockets-engine
`--wind-profile` end to end, validated against two sounding days. Then
`winds-engine 0.1.0` with one GRIB product.

## Status: winds-engine 0.1.0 scaffolded (2026-10-24)

The crate exists with the modeling core **staged in it** rather than
waiting on the kernel release, because a useful initial crate needed the
types to exist somewhere. `src/field.rs` (`Wind`/`WindField`/
`ConstantField`/`ProfileField`) and `src/analytic.rs` (`LogProfile`,
`EkmanProfile`) are this crate's own types. ballistics-engine's `wind`
module (unreleased) adapted them rather than copying them, and its
versions differ: `WindField::wind(&self, time_s: f64, position_m:
Vector3<f64>) -> Vector3<f64>` with no `Wind` type,
`ConstantField(pub Vector3<f64>)`, private `LogProfile`/`EkmanProfile`
fields behind getters, and constructors returning `WindFieldError` (its
`EkmanProfile` also rejects latitude 0, which this crate accepts).
Moving winds-engine onto those types changes winds-engine's public API
(the README example's `[f64; 3]` positions and `Wind` results,
`validate()`'s `ProfileField` argument), so it is a deliberate breaking
release, not a drop-in re-export. Zero dependencies, 15 tests
including a real 115-level KMPX sounding (IEM RAOB archive, checked in as
`tests/data/mpx_20261023_00z.txt`), CI on GitHub.

Revised order: (1) promote the staged modeling core into ballistics-engine
(one kernel release), (2) rockets-engine `with_wind_field` /
`--wind-profile` against it, (3) GRIB acquisition here, (4) more sounding
validation days (one is checked in already).

## Status: integrated (2026-10-25)

winds-engine is on crates.io and rockets-engine consumes it:
`BallisticForces::with_wind_field`, `Weathercock::with_wind_field` and
`RecoveryForces::with_wind_field` share an `Arc<dyn WindField>`, and the
CLI flies `--wind-profile alt:speed:dir,...` and `--wind-sounding file`
(the radiosonde interchange format) end to end — forces, weathercock and
recovery all sample the same field, turbulence composing on top. The
kernel's own `WindField` work is on ballistics-engine's
`agent/wind-field-20261025` branch; when that release lands, this crate
moves onto the kernel types in a breaking release (above) and
rockets-engine switches its import. 0.1.0 turned the Ekman spiral the
wrong way; 0.1.1 turns it the right way. Next up here: GRIB acquisition (the remaining risk item).
