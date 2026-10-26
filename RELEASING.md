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

The default build must stay zero-dependency (`cargo tree --depth 0` shows
only winds-engine). `grib` and `download` are the only optional deps/features.

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

## 4. Deploy the winds.rs site (separate checkout!)

The site source lives in `~/projects/winds-engine-agent` on branch
`agent/sounding-datum-20261025` (it is not on this repo's main). To update
the site copy:

```sh
cd ~/projects/winds-engine-agent
# edit site/index.html
git add site/index.html
git -c user.name="Alex Jokela" -c user.email="alex@camulus.com" commit -m "…"
git push origin agent/sounding-datum-20261025
aws s3 sync site/ s3://winds.rs/ --delete
INV=$(aws cloudfront create-invalidation --distribution-id E2RFG1P8MLK47T --paths "/*" --query 'Invalidation.Id' --output text)
aws cloudfront wait invalidation-completed --distribution-id E2RFG1P8MLK47T --id "$INV"
curl -s https://winds.rs/ | grep -o 'winds-engine = "0\.N"'   # confirm
```

Infrastructure: S3 bucket `winds.rs`, CloudFront `E2RFG1P8MLK47T`
(`d3utvj3ubazgyc.cloudfront.net`), aliases winds.rs + www.

## Determinism / dependency policy (per AGENTS.md)

Core: no network, no wall clock, no unseeded randomness. Acquisition
(HTTP) lives behind the `download` feature only; the decode/validate path
is offline and deterministic. Terrain (`terrain.rs`), profiles and
analytics are zero-dependency.