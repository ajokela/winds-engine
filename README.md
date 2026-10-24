# winds-engine

Wind *data* for rocket modeling: GRIB2/HRRR/GFS ingestion and radiosonde
validation feeding ballistics-engine's wind profile types. Kernel-side wind
modeling (the `WindField` trait, profiles, analytics) lives in
ballistics-engine itself — this crate is the acquisition and validation
pipeline, consumed by rockets-engine.

See [docs/winds-engine-plan.md](docs/winds-engine-plan.md) for the
architecture decision and plan.

MIT licensed (see LICENSE).
