# Changelog

## [0.9.3]

- Depend on ballistics-engine 0.50 (was 0.49). No source change was needed:
  the wind-field core winds-engine re-exports (`WindField`, `ProfileField`,
  `ConstantField`, `LogProfile`, `EkmanProfile`, `HeightOffset`,
  `WindFieldError`) is untouched between ballistics-engine 0.49.0 and 0.50.0
  (`git diff v0.49.0 v0.50.0 -- src/wind*` there shows only
  `wind_scenarios.rs`, which this crate does not use), and winds-engine uses no
  humidity, station-air or kernel-tuple API. Patch release: no public
  signature or re-exported type changed.
- No numeric shift for winds-engine users: the full test suite, run unchanged
  against 0.50.0 with every feature, produces the same results. The 0.50.0
  humidity and station-air changes affect trajectory solving, not wind
  profiles.
