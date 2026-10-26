# Terrain-aware wind: scope and plan

Status: plan, 2026-10-26. Gated: build only when a concrete terrain-wind
requirement exists (moving-batteries artillery over a ridge/valley, or
similar). Rockets-engine flights are seconds over profile-scale terrain and
do not need this.

## Where this sits vs WindNinja

WindNinja's three products: a weather-initialization pipeline (we have the
trajectory-relevant slice — GFS/HRRR acquisition, profiles, soundings), an
operational application (wildfire; ours is trajectory consumption), and a
**terrain flow solver**. This plan is about the cheapest defensible version
of that third piece: the 1980s *diagnostic mass-consistent adjustment*
(NOABL / CONMESS / Sherman 1978), not WindNinja's variational FEM and not
CFD.

## The physics, in one paragraph

A mass-conserving wind field satisfies continuity. Take an initial wind
(a `WindField`: profile or GRIB), compute a first-guess 3D velocity field,
then enforce `div(ρu) = 0` by solving one scalar Poisson equation for a
potential `φ`: `u = u0 − ∇φ`. Terrain enters through terrain-following
coordinates (the vertical metric keeps flow along slopes rather than
through them) and a no-flow boundary at the ground. One reference wind at
one point constrains the field. That is the entire algorithm — a hundred
lines of physics on top of a sparse solve.

## Scope: Tier 1 only (diagnostic adjustment)

1. **DEM ingestion.** Read SRTM 1-arcsec (~30 m) GeoTIFF and `.hgt` into a
   height grid. No GDAL: `tiff` + raw `.hgt` parsing. ASCII gridded XYZ as
   a test-friendly fallback.
2. **Initialization field.** Evaluate an existing `WindField` over the
   domain to seed the first guess; terrain is a correction, not a
   replacement. A single `ProfileField` (or the GRIB path) is the natural
   source.
3. **Adjustment.** Terrain-following coordinates, Gauss–Seidel (or SOR)
   relaxation of the scalar potential on the volume grid, Neumann BCs at
   ground and top, one reference-wind constraint. nalgebra's `sparse`
   feature; no new scientific stack.
4. **Output.** A `TerrainFlowField` implementing `WindField`, so it drops
   into rockets-engine's `with_wind_field` unchanged — forces, weathercock,
   recovery and turbulence all see it.

## API sketch

```rust
/// A terrain-adjusted field wrapping a base wind.
pub struct TerrainFlowField { /* height grid, adjusted u/v/w grids */ }

impl TerrainFlowField {
    /// Build from a DEM height grid and a base field, adjusting at the
    /// reference point (site location and measured speed/direction).
    pub fn adjust(dem: &HeightGrid, base: impl WindField, site: Site) -> Result<Self>;
}

impl WindField for TerrainFlowField {
    fn wind(&self, time_s: f64, position_m: [f64; 3]) -> Wind;
}
```

## Determinism and dependency policy (unchanged)

No network in the core, no wall clock, no unseeded randomness — the same
promise as the rest of the crate. DEM decoding is offline and deterministic.
No C dependencies, no GDAL/VTK, stays WASM-safe.

## Validation plan (the real cost driver)

Golden-pattern parity, like everything else in this repo:

1. **Idealized ridge**: a symmetric 2-D ridge DEM; pin the speed-up over the
   crest, the sheltered lee, and `div = 0` to a tolerance.
2. **Mass-conservation invariant**: flux through nested control volumes is
   conserved by construction.
3. **Cross-check vs a reference tool**: run WindNinja (or OpenWind) on one
   canonical terrain, record the flow, pin ours within a stated tolerance —
   same pattern as `scripts/openrocket-parity` in rockets-engine.
4. **Determinism**: byte-identical field for identical inputs; snapshot/replay
   unaffected.

## Effort and order

| Step | Effort |
| --- | --- |
| DEM ingestion (GeoTIFF + .hgt) | 1–2 days |
| Terrain-following adjustment + sparse solve | 3–5 days |
| `TerrainFlowField` + rockets-engine integration seam | 1–2 days |
| Idealized-ridge + conservation tests | 2 days |
| WindNinja/OpenWind cross-check (one terrain) | 2–3 days |

~2–3 weeks to a validated Tier 1. Revisit only with a concrete consumer.

## Status: Tier 1 built (2026-10-26)

The gate opened (moving-batteries: distributed wind monitors down a
shooting range). The core is in, zero new dependencies:

- `HeightGrid`: metres-grid terrain, `new` + an ASCII `.asc` reader
  (GeoTIFF/`.hgt` loaders still to come; `.hgt` needs a degree→metre
  projection step).
- `TerrainFlowField::adjust`: per-level 2-D horizontal mass-consistent
  adjustment with the terrain-induced vertical divergence as the forcing
  (air rising up a slope converges, lee descents diverge), surface wind
  monitors blended near the ground with influence decaying with height and
  distance, and a background-clamped domain boundary.
- Tests pin the honest physics: symmetric-hill speed-up toward the crest,
  return to free stream far upwind, flat-ground invariances, monitor
  localisation, determinism, ASCII parsing.

Known v1 limits (all Tier-2: full terrain-following mesh and vertical
coupling): the adjustment is 2-D per level, monitors are blended into a
first guess rather than hard observation constraints, and the far field
relaxes to the background within ~10%.

- **Tier 2** (full terrain-following mesh, Richardson iteration, steep-canyon
  fidelity — WindNinja's default): months; not planned.
- **Tier 3** (CFD, OpenFOAM-style RANS/LES): never.
- Any UI, raster I/O beyond the two formats above, or wildfire coupling.
