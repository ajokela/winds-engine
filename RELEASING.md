# Releasing winds-engine

Library-only crate (no binaries, no browser terminal). The public website
(winds.rs) deploys from a separate checkout — see below.

## Authorship

- Commit as **Alex Jokela <alex@camulus.com>**, author and committer:
  `git -c user.name="Alex Jokela" -c user.email="alex@camulus.com" commit ...`
- No AI/assistant attribution anywhere.

## 1. Version bump

`Cargo.toml`: `version = "0.N.0"`. Then `cargo generate-lockfile`.

## 2. Verify

```sh
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test                      # core: zero-dependency build
cargo test --features grib      # GRIB2 decoding
cargo test --features download  # HTTP acquisition (urllib + range requests)
```

The default build depends only on ballistics-engine (`default-features =
false`: the wind-field core, re-exported) and nalgebra — check with
`cargo tree --depth 1`. `grib` and `download` are the only optional
deps/features.

Live NOAA tests are `#[ignore]`d — run manually if you touch fetch/decoding:

```sh
cargo test --features download -- --include-ignored   # can take ~30 s
```

## 3. Commit, publish, push

```sh
git add -A
git -c user.name="Alex Jokela" -c user.email="alex@camulus.com" commit -m "…"
cargo publish
git push origin main
```

## 4. Deploy the winds.rs site

The site lives in this repo's `site/` and deploys from here (the same
lifecycle as the crate, guarded by CI's version-string check):

```sh
aws s3 sync site/ s3://winds.rs/ --delete
INV=$(aws cloudfront create-invalidation --distribution-id E2RFG1P8MLK47T --paths "/*" --query 'Invalidation.Id' --output text)
aws cloudfront wait invalidation-completed --distribution-id E2RFG1P8MLK47T --id "$INV"
curl -s https://winds.rs/ | grep -o 'winds-engine = "0\.N"'   # confirm
```

Infrastructure: S3 bucket `winds.rs`, CloudFront `E2RFG1P8MLK47T`
(`d3utvj3ubazgyc.cloudfront.net`), aliases winds.rs + www.

Note: an older working copy at `~/projects/winds-engine-agent` used to
carry the site on a side branch; it is stale and no longer the source of
truth — edit `site/` here.

## Determinism / dependency policy (per AGENTS.md)

Core: no network, no wall clock, no unseeded randomness. Acquisition
(HTTP) lives behind the `download` feature only; the decode/validate path
is offline and deterministic. The wind-field core is ballistics-engine's
`wind` module (no default features), re-exported; terrain (`terrain.rs`)
adds nothing beyond nalgebra.