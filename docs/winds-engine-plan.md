# winds-engine: go/no-go analysis and Tier A plan

Status: proposal, 2026-10-24
Consumers: rockets-engine (trajectories, weathercock, dispersion), moving-batteries (artillery)

## Verdict

- **Tier A — build.** Real-wind profiles and gridded-model ingestion as a small
  deterministic crate. ~2–3 weeks.
- **Tier B — defer.** Diagnostic mass-consistent terrain flow, gated on a
  concrete moving-batteries requirement. ~3–4 weeks if approved.
- **Tier C — never.** A full WindNinja clone (FEM solver + weather init + GUI)
  is 6–12 months competing with a mature free tool. Out of scope permanently.

## What WindNinja is, so we scope against reality

WindNinja is three products under one name:

1. A mass-conserving terrain flow solver (variational FEM over a DEM).
2. A weather-data initialization pipeline (HRRR/NAM/GFS grids, METAR point
   initialization, thermal effects).
3. An operational wildfire tool (GUI, GeoTIFF/KMZ export, FARSITE/FlamMap
   coupling).

Roughly 100k+ lines of C++ on GDAL/VTK, hardened over 15+ years for fire
behavior prediction. Our domain is not fire; our trajectories are seconds to
minutes, not hours. Only the middle layer — real atmospheric data, made
deterministic — maps to our needs, plus optionally a much simpler terrain
model than theirs.

## The gap this closes

Current wind handling in rockets-engine:

- constant world-frame vector (`--wind`),
- downrange-segmented zones (`--wind-zone`, upstream `WindSock`),
- deterministic frozen-field turbulence (`--turbulence`, 8-component
  multisine, pink spectrum, seed-fixed phases).

What is missing is **truth**: a real vertical profile for a real site on a
real day, and (maybe) terrain-channeled flow. Dispersion currently varies
wind as white noise around a guessed mean; with Tier A the mean and its
altitude structure come from an actual forecast/analysis grid.

## Architecture

New crate `winds-engine`, sibling of ballistics-engine, consumed by
rockets-engine the same way (pinned, read-only across repos). One trait,
stateless and deterministic — the same pattern the turbulence work proved
threads through forces, weathercock, and dispersion without integrator
changes:

```rust
/// A deterministic, stateless wind field.
pub trait WindField {
    /// World-frame air velocity at a time and position. +X tailwind, Y up.
    fn wind(&self, time_s: f64, position_m: Vec3) -> Vec3;
}
```

Provided implementations:

- `ConstantField(Vec3)` — today's behaviour.
- `ProfileField { levels: Vec<(altitude_m, speed_mps, direction_deg)> }`
  — piecewise-linear vertical profile; the target of all providers below.
- `LogProfile { reference_speed, reference_height_m, roughness_m }` and
  `EkmanProfile { geostrophic, latitude_deg, roughness_m }` — analytic
  boundary-layer profiles.
- `CompositeField<A, B>` — additive composition (profile + turbulence).

Integration with rockets-engine:

- `BallisticForces::with_wind_field(Box<dyn WindField>)` (or a generic
  parameter); `effective_wind` delegates, turbulence composes on top.
- `Weathercock` samples the same field (it already takes the turbulence
  gust field; same wiring).
- CLI: `--wind-profile <alt:speed:dir,...>` immediately; `--wind-grib
  <path-or-url> --at <lat,lon>` once ingestion lands.
- Dispersion: profile uncertainty via explicit perturbation parameters
  (profile speed scale / direction jitter per run, seeded), reusing the
  existing per-run reseeding pattern.

## Tier A scope

1. `WindField` trait + `ConstantField` + `ProfileField` + analytic profiles
   (log, Ekman).
2. GRIB2 ingestion of HRRR/GFS grids → `ProfileField` at a coordinate:
   - decode one file/one point (u/v components on isobaric levels),
   - bilinear horizontal interpolation, altitude conversion,
   - offline cache directory; NOAA filter/subset HTTP endpoint as the
     acquisition path so CI and users avoid multi-GB downloads.
3. Determinism and reproducibility tests throughout; no chrono/network in
   the core (acquisition behind a feature flag).

### GRIB2 risk item

Rust GRIB2 decoding is the main technical risk. Mitigation ladder:

1. Evaluate existing Rust GRIB crates for HRRR/GFS u/v-on-isobaric records.
2. Fall back to the NOAA NOMADS subset API (server-side extraction, small
   downloads) with a checked-in golden file for tests.
3. Last resort: ship a tiny GRIB2 template-5.0/5.2 decoder for the two
   products we need; explicitly not a general GRIB library.

## Tier B scope (deferred; decision gate)

Diagnostic mass-consistent wind over terrain (NOABL/CONMESS family, not
WindNinja's FEM solver):

- GeoTIFF DEM ingestion.
- Terrain-following coordinate transform, Poisson adjustment for mass
  consistency (nalgebra sparse; already in the dependency tree).
- Validation against published WindNinja runs on canonical terrain
  (idealized ridge, one real site).
- **Gate:** build only if moving-batteries states a concrete terrain-wind
  requirement. Rockets do not need it (low downrange, short exposure).

## Validation plan

- Profile providers: compare against radiosonde archives (University of
  Wyoming sounding data) for a handful of station-days; assert layer-mean
  speed/direction within stated tolerances.
- GRIB path: checked-in golden subset file; decoded profile pinned to
  stored expected values (bit-stable decoder test).
- Integration: rockets-engine flight with a profile vs equivalent zones —
  sanity ordering, determinism per `--seed`, and dispersion still
  reproducible.
- Tier B (if built): WindNinja comparison runs documented like
  `scripts/openrocket-parity/` — driver, constants, tolerance, CI job.

## Dependency and determinism policy

- Core crate: no network, no wall-clock, no unseeded randomness; everything
  reproducible from inputs (matches the ballistics-engine/rockets-engine
  rule).
- Acquisition (HTTP subset downloads) behind a non-default `download`
  feature so the library stays embeddable (WASM included).
- No AI/assistant attribution in commits; commit as the human developer.

## Effort and order

| Item | Effort | Depends on |
| --- | --- | --- |
| Trait + constant/profile/analytic fields | 3–4 days | — |
| rockets-engine integration (`--wind-profile`) | 2–3 days | trait |
| GRIB acquisition + decode (one product) | 1–2 weeks | trait; risk item above |
| Radiosonde validation set | 2–3 days | profile fields |
| Tier B terrain flow (if gated yes) | 3–4 weeks | DEM source decision |

Suggested first milestone: `winds-engine 0.1.0` = trait + profiles +
`--wind-profile` in rockets-engine, validated against two sounding days.
