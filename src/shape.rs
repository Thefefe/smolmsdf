use super::{
    curve::{BezierCurve, Curve, Linear, Quadratic, Cubic},
    math::{vec2, Rect, Vec2},
};
use std::ops::Range;

#[derive(Debug, Clone)]
pub struct Contour {
    pub edge_range: Range<usize>,
    pub aabb: Rect,
    pub start_end_same_edge: bool,
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub curve_range: Range<usize>,
    pub aabb: Rect,
}

/// Configuration for a shape writer.
#[derive(Debug, Clone, Copy)]
pub struct ShapeConfig {
    /// The sin of the minimum angle(in radians) between edges.
    pub angle_threshold_sin: f32,
    // pub size_threshold: f32,
}

impl Default for ShapeConfig {
    fn default() -> Self {
        Self {
            angle_threshold_sin: 15.0f32.to_radians().sin(),
            // size_threshold: 0.0,
        }
    }
}

/// A shape that can be rasterized as a signed distance field.
#[derive(Debug, Clone, Default)]
pub struct Shape {
    pub curves: Vec<Curve>,
    pub contours: Vec<Contour>,
    pub edges: Vec<Edge>,
    pub shape_aabb: Rect,
}

impl Shape {
    pub fn new() -> Self {
        Self {
            curves: Vec::new(),
            contours: Vec::new(),
            edges: Vec::new(),
            shape_aabb: Rect::MIN,
        }
    }

    pub fn clear(&mut self) {
        self.curves.clear();
        self.contours.clear();
        self.edges.clear();
        self.shape_aabb = Rect::MIN;
    }

    /// Constructs a writer for this shape with the given config.
    pub fn writer(&mut self, config: ShapeConfig) -> ShapeWriter {
        ShapeWriter {
            shape: self,
            config,
            current_point: Vec2::ZERO,
            last_dir: None,
            pending_move_start_point: None,
            pending_move_end_point: None,
            contour_aabb: Rect::MIN,
            edge_aabb: Rect::MIN,
        }
    }
}

#[derive(Debug)]
pub struct ShapeWriter<'a> {
    shape: &'a mut Shape,
    pub config: ShapeConfig,

    current_point: Vec2,
    last_dir: Option<Vec2>,

    pending_move_start_point: Option<Vec2>,
    pending_move_end_point: Option<Vec2>,

    contour_aabb: Rect,
    edge_aabb: Rect,
}

impl ShapeWriter<'_> {
    pub fn add_curve(&mut self, mut curve: Curve) {
        if let Some(last_dir) = self.last_dir {
            let start_dir = curve.dir_at_t(0.0);
            let angle = last_dir.normalized().cross(start_dir.normalized()).abs();
            if angle > self.config.angle_threshold_sin || last_dir.dot(start_dir) <= 0.0 {
                self.flush_edge();
            }
        }

        if let Some(new_start) = self.pending_move_start_point.take() {
            curve.move_start_point(new_start);
        }

        let curve_aabb = curve.bounding_box();
        let final_point = curve.end_point();

        self.shape.curves.push(curve);
        self.current_point = final_point;
        self.last_dir = Some(curve.dir_at_t(1.0));

        self.shape.shape_aabb = self.shape.shape_aabb.fitted_to_rect(curve_aabb);
        self.contour_aabb = self.contour_aabb.fitted_to_rect(curve_aabb);
        self.edge_aabb = self.edge_aabb.fitted_to_rect(curve_aabb);
    }

    pub fn flush_edge(&mut self) {
        let start = self.shape.edges.last().map_or(0, |e| e.curve_range.end);

        if start == self.shape.curves.len() {
            return;
        }

        let curve_range = start..self.shape.curves.len();

        // if self.edge_aabb.size().length() < self.config.size_threshold {
        //     let start_point = self.curves[curve_range.start].start_point();
        //     let end_point = self.curves[curve_range.end - 1].end_point();
        //     let center_point = (start_point + end_point) * 0.5;

        //     let first_edge_in_contour =
        //         self.contours.last().map_or(0, |c| c.edge_range.end) == self.edges.len();

        //     if !first_edge_in_contour {
        //         self.curves[curve_range.start - 1].move_end_point(center_point);
        //     } else {
        //         self.pending_move_end_point = Some(center_point);
        //     }

        //     self.pending_move_start_point = Some(center_point);

        //     self.curves.drain(curve_range.clone());
        //     self.edge_aabb = Rect::MIN;
        //     self.last_dir = None;

        //     return;
        // }

        self.shape.edges.push(Edge {
            curve_range,
            aabb: self.edge_aabb,
        });
        self.edge_aabb = Rect::MIN;
        self.last_dir = None;
    }

    pub fn flush_contour(&mut self) {
        let start = self.shape.contours.last().map_or(0, |c| c.edge_range.end);
        if start == self.shape.edges.len() {
            return;
        }

        let mut edge_range = start..self.shape.edges.len();

        if let Some(new_start) = self.pending_move_start_point.take() {
            self.shape.curves[self.shape.edges[edge_range.start].curve_range.start].move_start_point(new_start);
        }

        if let Some(new_end) = self.pending_move_end_point.take() {
            self.shape.curves[self.shape.edges[edge_range.end - 1].curve_range.end - 1].move_end_point(new_end);
        }

        if edge_range.len() == 1 && self.shape.edges[start].curve_range.len() == 1 {
            let last_curve = self.shape.curves.pop().unwrap();
            let split = last_curve.split(0.5);
            self.shape.curves.extend_from_slice(&split);

            self.shape.edges[start].aabb = split[0].bounding_box();
            self.shape.edges.push(Edge {
                curve_range: self.shape.curves.len() - 1..self.shape.curves.len(),
                aabb: split[1].bounding_box(),
            });

            edge_range.end += 1;
        }

        let end_dir = self.shape.curves.last().unwrap().dir_at_t(1.0);
        let start_dir = self.shape.curves[self.shape.edges[start].curve_range.start].dir_at_t(0.0);

        let angle = end_dir.normalized().cross(start_dir.normalized()).abs();
        let start_end_same_edge =
            angle <= self.config.angle_threshold_sin && end_dir.dot(start_dir) > 0.0;

        self.shape.contours.push(Contour {
            edge_range,
            aabb: self.contour_aabb,
            start_end_same_edge,
        });

        self.contour_aabb = Rect::MIN;
        self.last_dir = None;

        assert!(self.pending_move_end_point.is_none());
        assert!(self.pending_move_start_point.is_none());
    }
}

#[cfg(feature = "ttf-parser")]
impl ttf_parser::OutlineBuilder for ShapeWriter<'_> {
    fn move_to(&mut self, x: f32, y: f32) {
        self.current_point = vec2(x, y);
    }

    fn line_to(&mut self, x: f32, y: f32) {
        let a = self.current_point;
        let b = vec2(x, y);
        self.add_curve(Curve::Linear(Linear(a, b)));
    }

    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        let a = self.current_point;
        let b = vec2(x1, y1);
        let c = vec2(x, y);
        self.add_curve(Curve::Quadratic(Quadratic(a, b, c)));
    }

    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        let a = self.current_point;
        let b = vec2(x1, y1);
        let c = vec2(x2, y2);
        let d = vec2(x, y);
        self.add_curve(Curve::Cubic(Cubic(a, b, c, d)));
    }

    fn close(&mut self) {
        self.flush_edge();  
        self.flush_contour();
    }
}
