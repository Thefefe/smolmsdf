use super::{
    math::{vec2, Rect, Vec2},
    BezierCurve, Curve, MsdfBitmap, MtsdfBitmap, SdfBitmap, Shape, SignedDistance,
};

#[derive(Debug, Clone, Copy)]
pub struct SdfConfig {
    pub unit_scale: f32,
    pub size_px: f32,
    pub sdf_radius_px: f32,
    pub min_edge_size_px: f32,
}

impl SdfConfig {
    #[cfg(feature = "ttf-parser")]
    pub fn from_face(face: &ttf_parser::Face, px_per_em: f32, sdf_radius_px: f32) -> Self {
        Self {
            unit_scale: (face.units_per_em() as f32).recip(),
            size_px: px_per_em,
            sdf_radius_px,
            min_edge_size_px: 2.0,
        }
    }
}

pub fn rasterize_mtsdf(shape: &Shape, bitmap: &mut MtsdfBitmap, config: &SdfConfig) -> Rect {
    let unit_to_px = config.size_px * config.unit_scale;
    let px_to_unit = unit_to_px.recip();

    let rect_px = Rect::new(
        (shape.shape_aabb.min.x * unit_to_px).floor(),
        (shape.shape_aabb.min.y * unit_to_px).floor(),
        (shape.shape_aabb.max.x * unit_to_px).ceil(),
        (shape.shape_aabb.max.y * unit_to_px).ceil(),
    )
    .expanded(config.sdf_radius_px);

    bitmap.resize(rect_px.size().to_array().map(|x| x as u32));

    let unit_to_norm = unit_to_px / config.sdf_radius_px;

    for (x, y, pixel) in bitmap.enumerate_pixels() {
        let offset_px = vec2(x as f32, y as f32) + vec2(0.5, 0.5);
        let p_unit = (rect_px.min + offset_px) * px_to_unit;

        *pixel = mtsd_shape(shape, p_unit, config.min_edge_size_px * px_to_unit).map(|sd| {
            let min_sd_norm = (sd.distance * unit_to_norm).clamp(-1.0, 1.0);
            ((min_sd_norm * 0.5 + 0.5) * 255.0) as u8
        });
    }

    let size_recip = config.size_px.recip();
    Rect {
        min: rect_px.min * size_recip,
        max: rect_px.max * size_recip,
    }
}

pub fn rasterize_msdf(shape: &Shape, bitmap: &mut MsdfBitmap, config: &SdfConfig) -> Rect {
    let unit_to_px = config.size_px * config.unit_scale;
    let px_to_unit = unit_to_px.recip();

    let rect_px = Rect::new(
        (shape.shape_aabb.min.x * unit_to_px).floor(),
        (shape.shape_aabb.min.y * unit_to_px).floor(),
        (shape.shape_aabb.max.x * unit_to_px).ceil(),
        (shape.shape_aabb.max.y * unit_to_px).ceil(),
    )
    .expanded(config.sdf_radius_px);

    bitmap.resize(rect_px.size().to_array().map(|x| x as u32));

    let unit_to_norm = unit_to_px / config.sdf_radius_px;

    for (x, y, pixel) in bitmap.enumerate_pixels() {
        let offset_px = vec2(x as f32, y as f32) + vec2(0.5, 0.5);
        let p_unit = (rect_px.min + offset_px) * px_to_unit;

        *pixel = msd_shape(shape, p_unit, config.min_edge_size_px * px_to_unit).map(|sd| {
            let min_sd_norm = (sd.distance * unit_to_norm).clamp(-1.0, 1.0);
            ((min_sd_norm * 0.5 + 0.5) * 255.0) as u8
        });
    }

    let size_recip = config.size_px.recip();
    Rect {
        min: rect_px.min * size_recip,
        max: rect_px.max * size_recip,
    }
}

pub fn rasterize_sdf(shape: &Shape, bitmap: &mut SdfBitmap, config: &SdfConfig) -> Rect {
    let unit_to_px = config.size_px * config.unit_scale;
    let px_to_unit = unit_to_px.recip();

    let rect_px = Rect::new(
        (shape.shape_aabb.min.x * unit_to_px).floor(),
        (shape.shape_aabb.min.y * unit_to_px).floor(),
        (shape.shape_aabb.max.x * unit_to_px).ceil(),
        (shape.shape_aabb.max.y * unit_to_px).ceil(),
    )
    .expanded(config.sdf_radius_px);

    bitmap.resize(rect_px.size().to_array().map(|x| x as u32));

    let unit_to_norm = unit_to_px / config.sdf_radius_px;

    for (x, y, pixel) in bitmap.enumerate_pixels() {
        let offset_px = vec2(x as f32, y as f32) + vec2(0.5, 0.5);
        let p_unit = (rect_px.min + offset_px) * px_to_unit;

        *pixel = {
            let sd = sd_curves(&shape.curves, p_unit);
            let min_sd_norm = (sd.distance * unit_to_norm).clamp(-1.0, 1.0);
            ((min_sd_norm * 0.5 + 0.5) * 255.0) as u8
        };
    }

    let size_recip = config.size_px.recip();
    Rect {
        min: rect_px.min * size_recip,
        max: rect_px.max * size_recip,
    }
}

fn msd_shape(shape: &Shape, p: Vec2, size_threshold: f32) -> [SignedDistance; 3] {
    let mut msd = [SignedDistance::INFINITE; 3];
    for contour in shape.contours.iter() {
        let last_edge_index = contour.edge_range.len() - 1;
        let mut color_counter = 0;
        for (edge_index, edge) in shape.edges[contour.edge_range.clone()].iter().enumerate() {
            let short_edge = edge.aabb.size().length() <= size_threshold;

            let skip_channel_index = if contour.edge_range.len() == 1 || short_edge {
                4
            } else if edge_index == 0
                || (edge_index == last_edge_index && contour.start_end_same_edge)
            {
                0
            } else {
                color_counter += 1;
                ((color_counter) % 2) + 1
            };

            let edge_sd = sd_curves(&shape.curves[edge.curve_range.clone()], p);

            for (channel_index, channel) in msd.iter_mut().enumerate() {
                if edge_sd < *channel && skip_channel_index != channel_index {
                    *channel = edge_sd;
                }
            }
        }
    }

    msd
}

fn mtsd_shape(shape: &Shape, p: Vec2, size_threshold: f32) -> [SignedDistance; 4] {
    let mut mtsd = [SignedDistance::INFINITE; 4];
    for contour in shape.contours.iter() {
        let last_edge_index = contour.edge_range.len() - 1;
        let mut color_counter = 0;
        for (edge_index, edge) in shape.edges[contour.edge_range.clone()].iter().enumerate() {
            let short_edge = edge.aabb.size().length() <= size_threshold;

            let skip_channel_index = if contour.edge_range.len() == 1 || short_edge {
                4
            } else if edge_index == 0
                || (edge_index == last_edge_index && contour.start_end_same_edge)
            {
                0
            } else {
                color_counter += 1;
                ((color_counter) % 2) + 1
            };

            let edge_sd = sd_curves(&shape.curves[edge.curve_range.clone()], p);

            for (channel_index, channel) in mtsd.iter_mut().enumerate() {
                if edge_sd <= *channel && skip_channel_index != channel_index {
                    *channel = edge_sd;
                }
            }
        }
    }

    mtsd
}

fn sd_curves(curves: &[Curve], p: Vec2) -> SignedDistance {
    let mut min_sd = SignedDistance::INFINITE;
    for curve in curves {
        let sd = curve.signed_distance(p);
        if sd < min_sd {
            min_sd = sd;
        }
    }
    min_sd
}
