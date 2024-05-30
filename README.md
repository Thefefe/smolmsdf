# smolmsdf

## About

A very minimal implementation of [Chlumsky's msdf](https://github.com/Chlumsky/msdfgen/files/3050967/thesis.pdf) font rasterizer. For a full implementation see: [msdfgen](https://github.com/Chlumsky/msdfgen).

## Example
```rs
// Load a typeface
let face: ttf_parser::Face = /* ... */;

let shape_config = smolmsdf::ShapeConfig::default();
let mut shape = smolmsdf::Shape::new();

let glyph_id = face.glyph_index('a')?;

// Write the glyph into the shape
face.outline_glyph(glyph_id, &mut shape.writer(shape_config))?;

let mut bitmap = smolmsdf::Bitmap::new();

let px_per_em = 32.0;
let sdf_radius_px = 2.0;
let sdf_config = smolmsdf::SdfConfig::from_face(&face, px_per_em, sdf_radius_px);

// The returned rect can be used directly for rendering
let mesh_rect = smolmsdf::rasterize_mtsdf(&shape, &mut bitmap, &config);

// Simple error correction
smolmsdf::correct_errors_interp_sign(&mut bitmap);

```