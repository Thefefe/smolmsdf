const ASCII: &'static str = "!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

fn noto_sans() -> ttf_parser::Face<'static> {
    const NOTO_SANS: &[u8] = include_bytes!("../assets/noto/NotoSans.ttf");
    ttf_parser::Face::parse(NOTO_SANS, 0).unwrap()
}

fn ascii_glyph_ids(face: &ttf_parser::Face) -> Vec<ttf_parser::GlyphId> {
    ASCII
        .chars()
        .map(|c| face.glyph_index(c).unwrap())
        .collect()
}

fn ascii_shapes(face: &ttf_parser::Face) -> Vec<smolmsdf::Shape> {
    ASCII
        .chars()
        .map(|c| {
            let glyph_id = face.glyph_index(c).unwrap();
            let mut shape = smolmsdf::Shape::new();
            face.outline_glyph(
                glyph_id,
                &mut shape.writer(smolmsdf::ShapeConfig::default()),
            )
            .unwrap();
            shape
        })
        .collect()
}

fn bench_shape_parse(c: &mut criterion::Criterion) {
    let face = noto_sans();
    let glyph_ids = ascii_glyph_ids(&face);
    let mut shape = smolmsdf::Shape::new();

    c.bench_function("parse_ascii_shapes", |b| {
        b.iter(|| {
            for glyph_id in glyph_ids.iter() {
                shape.clear();
                face.outline_glyph(
                    *glyph_id,
                    &mut shape.writer(smolmsdf::ShapeConfig::default()),
                );
            }
        });
    });
}

fn glyph_raster_bench(c: &mut criterion::Criterion) {
    let face = noto_sans();
    let shapes = ascii_shapes(&face);
    let mut bitmap = smolmsdf::Bitmap::with_size(0, 0);

    for (px_per_em, sdf_radius_px) in [(64.0, 4.0), (32.0, 2.0), (16.0, 1.0)] {
        let config = smolmsdf::SdfConfig::from_face(&face, px_per_em, sdf_radius_px);

        c.bench_function(&format!("rasterize_ascii_glyphs {px_per_em}px"), |b| {
            b.iter(|| {
                shapes.iter().for_each(|shape| {
                    smolmsdf::rasterize_mtsdf(criterion::black_box(shape), &mut bitmap, &config);
                })
            })
        });
    }
}

criterion::criterion_group!(benches, bench_shape_parse, glyph_raster_bench);
criterion::criterion_main!(benches);
