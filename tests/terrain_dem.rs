//! Dem loading: SRTM `.hgt` and (behind `tiff`) GeoTIFF.

use winds_engine::{hgt_grid, SrtmResolution};

fn synthetic_hgt(resolution: SrtmResolution, heights: &[i16]) -> Vec<u8> {
    // Byte-exact SRTM layout: big-endian i16, north→south rows.
    let n = resolution.samples();
    assert_eq!(heights.len(), n * n);
    let mut out = Vec::with_capacity(heights.len() * 2);
    for &h in heights {
        out.extend_from_slice(&h.to_be_bytes());
    }
    out
}

#[test]
fn hgt_round_trips_a_known_tile() {
    // All zeros except a peak of 500 m at the north-west corner — which is
    // the first sample in the file (northernmost row) and, after the loader
    // reverses rows to be north-up, row 0.
    let n = SrtmResolution::ThreeArcSeconds.samples();
    assert_eq!(n, 1201); // guards the resolution enum
    let mut h = vec![0i16; 1201 * 1201];
    h[0] = 500; // north-west corner: first file sample
    let data = synthetic_hgt(SrtmResolution::ThreeArcSeconds, &h);
    let grid = hgt_grid(
        &data,
        SrtmResolution::ThreeArcSeconds,
        45.0,
        -93.0,
        45.0,
        -93.0,
    )
    .unwrap();
    assert_eq!((grid.ncols, grid.nrows), (1201, 1201));
    // Origin is the north-west corner when the reference is the tile corner.
    assert!((grid.origin_x - 0.0).abs() < 1e-9);
    assert!((grid.origin_y - 0.0).abs() < 1e-9);
    // The peak is at the north-west corner: file sample 0 = grid index 0.
    assert!(
        (grid.heights[0] - 500.0).abs() < 1e-9,
        "{}",
        grid.heights[0]
    );
    // SRTM-3 spacing at 45°N: ~65.5 m east-west, ~92.7 m north-south.
    assert!((grid.dx - 65.5).abs() < 1.0, "dx {}", grid.dx);
    assert!((grid.dy - 92.7).abs() < 1.0, "dy {}", grid.dy);
}

#[test]
fn hgt_validates_size() {
    assert!(hgt_grid(
        &[0u8; 10],
        SrtmResolution::ThreeArcSeconds,
        45.0,
        -93.0,
        45.0,
        -93.0
    )
    .is_err());
    // Correct byte count for 1201² samples: 1201*1201*2.
    let n = 1201usize;
    assert!(hgt_grid(
        &vec![0u8; n * n * 2 - 2],
        SrtmResolution::ThreeArcSeconds,
        45.0,
        -93.0,
        45.0,
        -93.0
    )
    .is_err());
}

#[cfg(feature = "tiff")]
#[test]
fn geotiff_decodes_a_small_raster() {
    use std::io::Cursor;
    use tiff::encoder::{colortype, TiffEncoder};

    // Build a tiny 2x3 16-bit DEM in memory, top-left origin.
    let mut buffer = Cursor::new(Vec::new());
    {
        let mut encoder = TiffEncoder::new(&mut buffer).unwrap();
        let data = vec![100u16, 200, 300, 400, 500, 600]; // 2 rows x 3 cols, row-major
        encoder
            .write_image::<colortype::Gray16>(3, 2, &data)
            .unwrap();
    }
    let bytes = buffer.into_inner();
    let grid = winds_engine::geotiff_grid(&bytes).unwrap();
    assert_eq!((grid.ncols, grid.nrows), (3, 2));
    // Rows reversed to north-up: the file's first row (100,200,300) ends up
    // as the last row, second (400,500,600) as the first.
    assert_eq!(&grid.heights[3..6], &[100.0, 200.0, 300.0]);
    assert_eq!(&grid.heights[0..3], &[400.0, 500.0, 600.0]);
    // Indexed from (0,0): no geotransform in a bare raster.
    assert_eq!(
        (grid.origin_x, grid.origin_y, grid.dx, grid.dy),
        (0.0, 0.0, 1.0, 1.0)
    );
}

#[cfg(feature = "tiff")]
#[test]
fn geotiff_reads_the_geotransform() {
    use std::io::Cursor;
    use tiff::encoder::{colortype, TiffEncoder};

    let mut buffer = Cursor::new(Vec::new());
    {
        let mut tiff = TiffEncoder::new(&mut buffer).unwrap();
        let mut image = tiff.new_image::<colortype::Gray16>(3, 2).unwrap();
        // Tie raster point (0,0) to model (500000, 4200000), 30 m pixels,
        // north-up (negative pixel-scale y).
        let tiepoint: &[f64] = &[0.0, 0.0, 0.0, 500000.0, 4200000.0, 0.0];
        let scale: &[f64] = &[30.0, -30.0, 0.0];
        image
            .encoder()
            .write_tag(tiff::tags::Tag::ModelTiepointTag, tiepoint)
            .unwrap();
        image
            .encoder()
            .write_tag(tiff::tags::Tag::ModelPixelScaleTag, scale)
            .unwrap();
        let data = vec![100u16, 200, 300, 400, 500, 600]; // 3 cols x 2 rows, row-major
        image.write_strip(&data).unwrap();
        image.finish().unwrap();
    }
    let grid = winds_engine::geotiff_grid(&buffer.into_inner()).unwrap();
    assert_eq!((grid.ncols, grid.nrows), (3, 2));
    assert!((grid.origin_x - 500000.0).abs() < 1e-6, "{}", grid.origin_x);
    assert!(
        (grid.origin_y - 4200000.0).abs() < 1e-6,
        "{}",
        grid.origin_y
    );
    assert!((grid.dx - 30.0).abs() < 1e-9, "dx {}", grid.dx);
    assert!((grid.dy - 30.0).abs() < 1e-9, "dy {}", grid.dy);
}

#[test]
fn hgt_projects_into_local_metres() {
    // Recentring a reading: the same tile but reference at the tile's
    // centre shifts the origin by -half a degree in metres.
    let n = SrtmResolution::ThreeArcSeconds.samples();
    let data = synthetic_hgt(SrtmResolution::ThreeArcSeconds, &vec![0i16; n * n]);
    let centred = hgt_grid(
        &data,
        SrtmResolution::ThreeArcSeconds,
        45.0,
        -93.0,
        44.5,
        -92.5,
    )
    .unwrap();
    // North-west corner is 0.5° north and 0.5° west of the reference:
    // ~39 km west, ~56 km north at 45° latitude.
    assert!(
        centred.origin_x < -38_000.0,
        "origin_x {}",
        centred.origin_x
    );
    assert!(centred.origin_y > 54_000.0, "origin_y {}", centred.origin_y);
}
