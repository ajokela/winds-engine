//! Diagnostic mass-consistent terrain wind (NOABL/CONMESS-style).
//!
//! A terrain height grid (in metres) plus a background [`WindField`] and
//! optional surface wind monitors become a terrain-adjusted field: each
//! fixed above-ground level is adjusted independently to satisfy 2-D
//! horizontal mass conservation (`div(u) = 0` with `u = u0 − ∇φ`), so flow
//! is deflected over and around terrain. This is the cheap diagnostic tier
//! of terrain flow, not a full CFD solve. Entirely offline and
//! deterministic.

use crate::field::{Wind, WindField};
use crate::{Error, Result};

/// Earth's mean radius for the local-tangent metre projection.
const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// One surface wind observation, in the terrain grid's frame (metres).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Monitor {
    /// Easting (frame x), metres.
    pub x_m: f64,
    /// Northing (frame z), metres.
    pub z_m: f64,
    /// Wind speed, metres per second.
    pub speed_mps: f64,
    /// Meteorological FROM-direction in the frame convention: 0° = from
    /// +x (a headwind), 90° = from the right (see [`Wind::from_speed_direction`]).
    pub direction_from_deg: f64,
}

impl Monitor {
    /// An observation already in the grid's metre frame.
    pub fn at(x_m: f64, z_m: f64, speed_mps: f64, direction_from_deg: f64) -> Self {
        Self {
            x_m,
            z_m,
            speed_mps,
            direction_from_deg,
        }
    }

    /// A geographic reading, projected into local metres and rotated from
    /// compass degrees into the frame. `ref_lat`/`ref_lon` is the point
    /// that maps to the frame's origin, and `north_azimuth_deg` is the
    /// compass direction of the frame's +x axis.
    pub fn from_geo(
        lat_deg: f64,
        lon_deg: f64,
        speed_mps: f64,
        dir_compass_from_deg: f64,
        ref_lat_deg: f64,
        ref_lon_deg: f64,
        north_azimuth_deg: f64,
    ) -> Self {
        let (x_m, z_m) = local_en_m(lat_deg, lon_deg, ref_lat_deg, ref_lon_deg);
        Self {
            x_m,
            z_m,
            speed_mps,
            direction_from_deg: (dir_compass_from_deg - north_azimuth_deg).rem_euclid(360.0),
        }
    }
}

/// Local equirectangular projection: `(east, north)` metres from the
/// reference point, for the small domains a wind network covers. Not for
/// continental scales.
pub fn local_en_m(lat_deg: f64, lon_deg: f64, ref_lat_deg: f64, ref_lon_deg: f64) -> (f64, f64) {
    let dlat = (lat_deg - ref_lat_deg).to_radians();
    let dlon = (lon_deg - ref_lon_deg).to_radians();
    let east = EARTH_RADIUS_M * dlon * ref_lat_deg.to_radians().cos();
    let north = EARTH_RADIUS_M * dlat;
    (east, north)
}

/// Parse a monitor file of one reading per line:
/// `latitude,longitude,speed_ms,direction_compass_from_deg`, `#` comments
/// and a header line are skipped. Coordinates are projected via
/// [`local_en_m`] and directions rotated into the frame.
pub fn monitors_from_geo_csv(
    text: &str,
    ref_lat_deg: f64,
    ref_lon_deg: f64,
    north_azimuth_deg: f64,
) -> Result<Vec<Monitor>> {
    let mut monitors = Vec::new();
    for (number, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 4 {
            return Err(Error::invalid(
                "monitor csv",
                format!(
                    "line {}: expected lat,lon,speed_ms,dir_compass_deg",
                    number + 1
                ),
            ));
        }
        let num = |token: &str, what: &str| -> Result<f64> {
            token.trim().parse().map_err(|_| {
                Error::invalid(
                    "monitor csv",
                    format!("line {}: {what} '{token}' is not numeric", number + 1),
                )
            })
        };
        let (lat, lon) = (num(fields[0], "lat")?, num(fields[1], "lon")?);
        let speed = num(fields[2], "speed")?;
        let dir = num(fields[3], "direction")?;
        if !lat.is_finite()
            || lat.abs() > 90.0
            || !lon.is_finite()
            || !speed.is_finite()
            || speed < 0.0
            || !dir.is_finite()
        {
            return Err(Error::invalid(
                "monitor csv",
                format!(
                    "line {}: readings must be finite with nonnegative speed",
                    number + 1
                ),
            ));
        }
        monitors.push(Monitor::from_geo(
            lat,
            lon,
            speed,
            dir,
            ref_lat_deg,
            ref_lon_deg,
            north_azimuth_deg,
        ));
    }
    Ok(monitors)
}

/// Terrain heights on an axis-aligned, uniform grid in metres. Row-major:
/// `heights[i * ncols + j]` is column `j` (east), row `i` (north). Position
/// `y` is up; the grid lies in the `x`–`z` plane.
#[derive(Debug, Clone)]
pub struct HeightGrid {
    /// Frame-x of the first column centre.
    pub origin_x: f64,
    /// Frame-z of the first row centre.
    pub origin_y: f64,
    /// Column spacing, metres.
    pub dx: f64,
    /// Row spacing, metres.
    pub dy: f64,
    /// Number of columns (east).
    pub ncols: usize,
    /// Number of rows (north).
    pub nrows: usize,
    /// Row-major height samples, `ncols * nrows`.
    pub heights: Vec<f64>,
}

impl HeightGrid {
    /// Build and validate a grid.
    pub fn new(
        origin_x: f64,
        origin_y: f64,
        dx: f64,
        dy: f64,
        ncols: usize,
        nrows: usize,
        heights: Vec<f64>,
    ) -> Result<Self> {
        if dx <= 0.0 || dy <= 0.0 || !dx.is_finite() || !dy.is_finite() {
            return Err(Error::invalid("height grid", "spacing must be positive"));
        }
        if ncols < 2 || nrows < 2 {
            return Err(Error::invalid("height grid", "needs at least 2x2 cells"));
        }
        if heights.len() != ncols * nrows {
            return Err(Error::invalid(
                "height grid",
                "height count does not match the grid shape",
            ));
        }
        if !heights.iter().all(|h| h.is_finite()) {
            return Err(Error::invalid("height grid", "heights must be finite"));
        }
        Ok(Self {
            origin_x,
            origin_y,
            dx,
            dy,
            ncols,
            nrows,
            heights,
        })
    }

    /// Terrain height at a horizontal position, bilinear between samples;
    /// clamps at the grid edges.
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let fx = ((x - self.origin_x) / self.dx).clamp(0.0, self.ncols as f64 - 1.0);
        let fy = ((z - self.origin_y) / self.dy).clamp(0.0, self.nrows as f64 - 1.0);
        let (j, i) = (fx.floor() as usize, fy.floor() as usize);
        let (rx, ry) = (fx - j as f64, fy - i as f64);
        let j1 = (j + 1).min(self.ncols - 1);
        let i1 = (i + 1).min(self.nrows - 1);
        let f = |i: usize, j: usize| self.heights[i * self.ncols + j];
        (1.0 - ry) * ((1.0 - rx) * f(i, j) + rx * f(i, j1))
            + ry * ((1.0 - rx) * f(i1, j) + rx * f(i1, j1))
    }

    /// A grid from an ASCII ESRI `.asc` raster (projected, metres): header
    /// `ncols/nrows/xllcorner/yllcorner/cellsize`, then row-major heights
    /// (north-to-south, west-to-east).
    pub fn from_ascii(text: &str) -> Result<Self> {
        let mut ncols = None;
        let mut nrows = None;
        let mut xll = None;
        let mut yll = None;
        let mut cellsize = None;
        let mut heights = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split_whitespace();
            match parts.next().unwrap_or("").to_ascii_lowercase().as_str() {
                "ncols" => ncols = parse_asc_num(parts.next()),
                "nrows" => nrows = parse_asc_num(parts.next()),
                "xllcorner" | "xllcenter" => xll = parse_asc_num(parts.next()),
                "yllcorner" | "yllcenter" => yll = parse_asc_num(parts.next()),
                "cellsize" => cellsize = parse_asc_num(parts.next()),
                "nodata_value" => continue,
                _ => {
                    for token in line.split_whitespace() {
                        if let Ok(v) = token.parse::<f64>() {
                            heights.push(v);
                        }
                    }
                }
            }
        }
        let ncols = ncols.ok_or_else(|| Error::invalid("ascii dem", "missing ncols"))? as usize;
        let nrows = nrows.ok_or_else(|| Error::invalid("ascii dem", "missing nrows"))? as usize;
        let cellsize = cellsize.ok_or_else(|| Error::invalid("ascii dem", "missing cellsize"))?;
        let (xll, yll) = (
            xll.ok_or_else(|| Error::invalid("ascii dem", "missing xllcorner"))?,
            yll.ok_or_else(|| Error::invalid("ascii dem", "missing yllcorner"))?,
        );
        if heights.len() != ncols * nrows {
            return Err(Error::invalid(
                "ascii dem",
                "height count does not match header",
            ));
        }
        // Rows are north→south; our rows run +z (north), so reverse them.
        let mut north_up = Vec::with_capacity(heights.len());
        for row in heights.chunks_exact(ncols).rev() {
            north_up.extend_from_slice(row);
        }
        Self::new(xll, yll, cellsize, cellsize, ncols, nrows, north_up)
    }
}

fn parse_asc_num(token: Option<&str>) -> Option<f64> {
    token.and_then(|s| s.parse().ok())
}

/// A diagnostic terrain-adjusted wind field.
#[derive(Debug, Clone)]
pub struct TerrainFlowField {
    grid: HeightGrid,
    /// Fixed above-ground sampling levels, ascending.
    levels_agl: Vec<f64>,
    /// Adjusted horizontal wind per level, row-major `u` and `v` (frame x
    /// and z), `levels * ncols * nrows` each.
    u: Vec<f64>,
    v: Vec<f64>,
}

impl TerrainFlowField {
    /// Build a terrain-adjusted field. `base` seeds the first guess;
    /// `monitors` are surface wind observations `(x, z, speed_mps,
    /// from_deg)` blended near the ground and relaxed toward the base field
    /// with height. `obs_length_m` sets the horizontal influence radius,
    /// `surface_decay_m` the vertical relaxation of monitor influence, and
    /// `terrain_scale_m` the depth over which the terrain forcing (flow
    /// rising over the ground) decays.
    pub fn adjust(
        grid: HeightGrid,
        base: &dyn WindField,
        monitors: &[Monitor],
        levels_agl: &[f64],
        obs_length_m: f64,
        surface_decay_m: f64,
        terrain_scale_m: f64,
    ) -> Result<Self> {
        if levels_agl.is_empty() || levels_agl.iter().any(|h| !h.is_finite() || *h < 0.0) {
            return Err(Error::invalid(
                "levels_agl",
                "must be finite and nonnegative",
            ));
        }
        if obs_length_m <= 0.0 || surface_decay_m <= 0.0 || terrain_scale_m <= 0.0 {
            return Err(Error::invalid("terrain", "length scales must be positive"));
        }
        if levels_agl.is_empty() || levels_agl.iter().any(|h| !h.is_finite() || *h < 0.0) {
            return Err(Error::invalid(
                "levels_agl",
                "must be finite and nonnegative",
            ));
        }
        if obs_length_m <= 0.0 || surface_decay_m <= 0.0 {
            return Err(Error::invalid("terrain", "length scales must be positive"));
        }
        for m in monitors {
            if !m.x_m.is_finite() || !m.z_m.is_finite() || !m.direction_from_deg.is_finite() {
                return Err(Error::invalid("monitor", "must be finite"));
            }
            if !m.speed_mps.is_finite() || m.speed_mps < 0.0 {
                return Err(Error::invalid("monitor", "speed must be nonnegative"));
            }
        }
        let (ncols, nrows) = (grid.ncols, grid.nrows);
        let n = ncols * nrows;
        let inv_dx2 = 1.0 / (grid.dx * grid.dx);
        let inv_dy2 = 1.0 / (grid.dy * grid.dy);
        let diag = -2.0 * (inv_dx2 + inv_dy2);
        let omega = 1.6;
        let mut u = Vec::with_capacity(levels_agl.len() * n);
        let mut v = Vec::with_capacity(levels_agl.len() * n);
        let mut phi = vec![0.0; n];
        for &level in levels_agl {
            // First guess: base field + monitor correction, relaxed with height.
            let mut u0 = vec![0.0; n];
            let mut v0 = vec![0.0; n];
            for i in 0..nrows {
                for j in 0..ncols {
                    let x = grid.origin_x + j as f64 * grid.dx;
                    let z = grid.origin_y + i as f64 * grid.dy;
                    let height = grid.heights[i * ncols + j] + level;
                    let bg = base.wind(0.0, [x, height, z]);
                    let (mut cx, mut cz, mut wsum) = (0.0, 0.0, 0.0);
                    for m in monitors {
                        let mw = Wind::from_speed_direction(m.speed_mps, m.direction_from_deg);
                        let bg_m = base.wind(0.0, [m.x_m, grid.height_at(m.x_m, m.z_m), m.z_m]);
                        let dx = x - m.x_m;
                        let dz = z - m.z_m;
                        let weight = (-(dx * dx + dz * dz) / (obs_length_m * obs_length_m)).exp()
                            * (-level / surface_decay_m).exp();
                        cx += weight * (mw.x_mps - bg_m.x_mps);
                        cz += weight * (mw.z_mps - bg_m.z_mps);
                        wsum += weight;
                    }
                    if wsum > 0.0 {
                        // A regularizer keeps a lone distant monitor from
                        // overriding the background everywhere.
                        let w_floor = 0.1;
                        cx /= wsum + w_floor;
                        cz /= wsum + w_floor;
                    }
                    u0[i * ncols + j] = bg.x_mps + cx;
                    v0[i * ncols + j] = bg.z_mps + cz;
                }
            }
            // Solve div(u) = 0, u = u0 - grad(phi), including the
            // terrain-induced vertical divergence: air forced up a slope
            // converges horizontally, air descending into a lee diverges.
            // laplacian(phi) = div_h(u0) - (u0 . grad(h)) * exp(-z / H) / H.
            let level_decay = (-level / terrain_scale_m).exp() / terrain_scale_m;
            let slopes = |u0x: f64, u0z: f64, i: usize, j: usize| -> f64 {
                let dhdx = (grid.heights[i * ncols + (j + 1).min(ncols - 1)]
                    - grid.heights[i * ncols + j.saturating_sub(1)])
                    / (2.0 * grid.dx);
                let dhdz = (grid.heights[(i + 1).min(nrows - 1) * ncols + j]
                    - grid.heights[i.saturating_sub(1) * ncols + j])
                    / (2.0 * grid.dy);
                u0x * dhdx + u0z * dhdz
            };
            phi.fill(0.0);
            for _ in 0..4000 {
                let mut max_step = 0.0_f64;
                for i in 1..nrows - 1 {
                    for j in 1..ncols - 1 {
                        let ip = (i + 1).min(nrows - 1);
                        let im = i.saturating_sub(1);
                        let jp = (j + 1).min(ncols - 1);
                        let jm = j.saturating_sub(1);
                        let lap = (phi[i * ncols + jp] + phi[i * ncols + jm]) * inv_dx2
                            + (phi[ip * ncols + j] + phi[im * ncols + j]) * inv_dy2
                            + diag * phi[i * ncols + j];
                        let div = (u0[i * ncols + jp] - u0[i * ncols + jm]) / (2.0 * grid.dx)
                            + (v0[ip * ncols + j] - v0[im * ncols + j]) / (2.0 * grid.dy);
                        let rhs =
                            div - slopes(u0[i * ncols + j], v0[i * ncols + j], i, j) * level_decay;
                        let updated = phi[i * ncols + j] + omega * (rhs - lap) / diag;
                        max_step = max_step.max((updated - phi[i * ncols + j]).abs());
                        phi[i * ncols + j] = updated;
                    }
                }
                if max_step < 1e-6 {
                    break;
                }
            }
            for i in 0..nrows {
                for j in 0..ncols {
                    let edge = i == 0 || i == nrows - 1 || j == 0 || j == ncols - 1;
                    if edge {
                        // The domain is embedded in a larger uniform field: the
                        // boundary holds the background wind.
                        u.push(u0[i * ncols + j]);
                        v.push(v0[i * ncols + j]);
                        continue;
                    }
                    let jp = (j + 1).min(ncols - 1);
                    let jm = j.saturating_sub(1);
                    let ip = (i + 1).min(nrows - 1);
                    let im = i.saturating_sub(1);
                    let gx = (phi[i * ncols + jp] - phi[i * ncols + jm]) / (2.0 * grid.dx);
                    let gz = (phi[ip * ncols + j] - phi[im * ncols + j]) / (2.0 * grid.dy);
                    u.push(u0[i * ncols + j] - gx);
                    v.push(v0[i * ncols + j] - gz);
                }
            }
        }
        Ok(Self {
            grid,
            levels_agl: levels_agl.to_vec(),
            u,
            v,
        })
    }

    fn sample_level(&self, level: usize, x: f64, z: f64) -> (f64, f64) {
        let (ncols, nrows) = (self.grid.ncols, self.grid.nrows);
        let fx = ((x - self.grid.origin_x) / self.grid.dx).clamp(0.0, ncols as f64 - 1.0);
        let fy = ((z - self.grid.origin_y) / self.grid.dy).clamp(0.0, nrows as f64 - 1.0);
        let (j, i) = (fx.floor() as usize, fy.floor() as usize);
        let (rx, ry) = (fx - j as f64, fy - i as f64);
        let j1 = (j + 1).min(ncols - 1);
        let i1 = (i + 1).min(nrows - 1);
        let base = level * ncols * nrows;
        let g = |arr: &[f64], i: usize, j: usize| arr[base + i * ncols + j];
        let sample = |arr: &[f64]| {
            (1.0 - ry) * ((1.0 - rx) * g(arr, i, j) + rx * g(arr, i, j1))
                + ry * ((1.0 - rx) * g(arr, i1, j) + rx * g(arr, i1, j1))
        };
        (sample(&self.u), sample(&self.v))
    }

    /// Horizontal wind at a world position: nearest above-ground level.
    pub fn horizontal_at(&self, x: f64, y: f64, z: f64) -> (f64, f64) {
        let agl = y - self.grid.height_at(x, z);
        let mut level = 0;
        for (k, &h) in self.levels_agl.iter().enumerate() {
            if agl <= h {
                level = k;
                break;
            }
            level = k;
        }
        self.sample_level(level, x, z)
    }
}

impl WindField for TerrainFlowField {
    fn wind(&self, _time_s: f64, position_m: [f64; 3]) -> Wind {
        let (u, v) = self.horizontal_at(position_m[0], position_m[1], position_m[2]);
        Wind {
            x_mps: u,
            y_mps: 0.0,
            z_mps: v,
        }
    }
}
