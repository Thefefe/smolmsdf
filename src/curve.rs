use super::math::{lerp, solve_cubic, solve_quadratic, Rect, Vec2};

#[derive(Debug, Clone, Copy)]
pub struct Line(pub Vec2, pub Vec2);

#[derive(Debug, Clone, Copy)]
pub struct Quadratic(pub Vec2, pub Vec2, pub Vec2);

#[derive(Debug, Clone, Copy)]
pub struct Cubic(pub Vec2, pub Vec2, pub Vec2, pub Vec2);

#[derive(Debug, Clone, Copy)]
pub enum Curve {
    Line(Line),
    Quad(Quadratic),
    Cubic(Cubic),
}

impl std::fmt::Display for Curve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Line(Line(p0, p1)) => write!(f, "line {p0} {p1}"),
            Self::Quad(Quadratic(p0, p1, p2)) => write!(f, "quad {p0} {p1} {p2}"),
            Self::Cubic(Cubic(p0, p1, p2, p3)) => write!(f, "cubic {p0} {p1} {p2} {p3}"),
        }
    }
}

impl From<Line> for Curve {
    fn from(value: Line) -> Self {
        Self::Line(value)
    }
}

impl From<Quadratic> for Curve {
    fn from(value: Quadratic) -> Self {
        Self::Quad(value)
    }
}

impl From<Cubic> for Curve {
    fn from(value: Cubic) -> Self {
        Self::Cubic(value)
    }
}

impl Curve {
    pub fn start_point(&self) -> Vec2 {
        match self {
            Curve::Line(line) => line.0,
            Curve::Quad(quad) => quad.0,
            Curve::Cubic(cubic) => cubic.0,
        }
    }

    pub fn end_point(&self) -> Vec2 {
        match self {
            Curve::Line(line) => line.1,
            Curve::Quad(quad) => quad.2,
            Curve::Cubic(cubic) => cubic.3,
        }
    }

    pub fn move_end_points(&mut self, new_start: Vec2, new_end: Vec2) {
        match self {
            Curve::Line(line) => *line = Line(new_start, new_end),
            Curve::Quad(quad) => {
                let Quadratic(p0, p1, p2) = *quad;

                let p10_pr_num = (p0 - p1).cross(new_start - p0);
                let p10_pr_denom = (p0 - p1).cross(p2 - p1);
                let p10_pr = p1 + (p10_pr_num / p10_pr_denom) * (p2 - p1);

                let p1_pr_num = (p2 - p10_pr).cross(new_end - p2);
                let p1_pr_denom = (p2 - p10_pr).cross(new_start - p10_pr);
                let p1_pr = p10_pr + (p1_pr_num / p1_pr_denom) * (new_start - p10_pr);

                quad.1 = p1_pr;
                quad.0 = new_start;
                quad.2 = new_end;
            }
            Curve::Cubic(cubic) => {
                cubic.1 += new_start - cubic.0;
                cubic.0 = new_start;
                cubic.2 += new_end - cubic.3;
                cubic.3 = new_end;
            }
        }
    }

    pub fn move_start_point(&mut self, new_start: Vec2) {
        self.move_end_points(new_start, self.end_point());
    }

    pub fn move_end_point(&mut self, new_end: Vec2) {
        self.move_end_points(self.start_point(), new_end);
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SignedDistance {
    pub distance: f32,
    pub orthogonality: f32,
}

impl SignedDistance {
    pub const INFINITE: Self = Self {
        distance: f32::NEG_INFINITY,
        orthogonality: 0.0,
    };
}

impl PartialOrd for SignedDistance {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if (self.distance.abs() - other.distance.abs()).abs() <= 1e-4 {
            other.orthogonality.partial_cmp(&self.orthogonality)
        } else {
            self.distance.abs().partial_cmp(&other.distance.abs())
        }
    }
}

pub trait BezierCurve: Sized {
    fn pos_at_t(&self, t: f32) -> Vec2;
    fn dir_at_t(&self, t: f32) -> Vec2;
    fn signed_distance(&self, p: Vec2) -> SignedDistance;
    fn bounding_box(&self) -> Rect;
    fn split(&self, t: f32) -> [Self; 2];
}

impl BezierCurve for Line {
    fn pos_at_t(&self, t: f32) -> Vec2 {
        self.0 + t * (self.1 - self.0)
    }

    fn dir_at_t(&self, _: f32) -> Vec2 {
        self.1 - self.0
    }

    fn signed_distance(&self, p: Vec2) -> SignedDistance {
        let t = ((p - self.0).dot(self.1 - self.0) / (self.1 - self.0).dot(self.1 - self.0))
            .clamp(0.0, 1.0);
        let v = self.pos_at_t(t) - p;
        let d = self.dir_at_t(t);

        let distance = v.length().copysign(d.cross(v));
        let orthogonality = d.normalized_or_zero().cross(-v.normalized_or_zero()).abs();

        SignedDistance {
            distance,
            orthogonality,
        }
    }

    fn bounding_box(&self) -> Rect {
        let min = self.0.min(self.1);
        let max = self.0.max(self.1);
        Rect { min, max }
    }

    fn split(&self, t: f32) -> [Self; 2] {
        let c = lerp(self.0, self.1, t);
        [Line(self.0, c), Line(c, self.1)]
    }
}

impl BezierCurve for Quadratic {
    fn pos_at_t(&self, t: f32) -> Vec2 {
        self.0 + 2.0 * t * (self.1 - self.0) + t * t * (self.2 - 2.0 * self.1 + self.0)
    }

    fn dir_at_t(&self, t: f32) -> Vec2 {
        2.0 * t * (self.2 - 2.0 * self.1 + self.0) + 2.0 * (self.1 - self.0)
    }

    fn signed_distance(&self, p: Vec2) -> SignedDistance {
        let qa = self.0 - p;
        let ab = self.1 - self.0;
        let br = self.2 - self.1 - ab;
        let a = br.dot(br);
        let b = 3.0 * ab.dot(br);
        let c = 2.0 * ab.dot(ab) + qa.dot(br);
        let d = qa.dot(ab);

        let mut closest_t = 0.0;
        let mut min_dist_sqr = f32::INFINITY;
        let mut closest_pos = Vec2::ZERO;
        for t in solve_cubic(a, b, c, d) {
            let t = t.clamp(0.0, 1.0);
            let pos = self.pos_at_t(t);
            let dist_sqr = pos.distance_squared(p);
            if dist_sqr < min_dist_sqr {
                closest_t = t;
                min_dist_sqr = dist_sqr;
                closest_pos = pos;
            }
        }

        let v = closest_pos - p;
        let d = self.dir_at_t(closest_t);

        let distance = min_dist_sqr.sqrt().copysign(d.cross(v));
        let orthogonality = d.normalized_or_zero().cross(-v.normalized_or_zero()).abs();

        SignedDistance {
            distance,
            orthogonality,
        }
    }

    fn bounding_box(&self) -> Rect {
        let t = -2.0 * (self.1 - self.0) / (2.0 * (self.2 - 2.0 * self.1 + self.0));
        let p_tx = self.pos_at_t(t.x.clamp(0.0, 1.0));
        let p_ty = self.pos_at_t(t.y.clamp(0.0, 1.0));

        let min = self.0.min(self.2).min(p_tx).min(p_ty);
        let max = self.0.max(self.2).max(p_tx).max(p_ty);
        Rect { min, max }
    }

    fn split(&self, t: f32) -> [Self; 2] {
        let a0 = lerp(self.0, self.1, t);
        let a1 = lerp(self.1, self.2, t);
        let c = lerp(a0, a1, t);

        [Quadratic(self.0, a0, c), Quadratic(c, a1, self.2)]
    }
}

impl Quadratic {
    pub fn as_cubic(&self) -> Cubic {
        Cubic(
            self.0,
            lerp(self.0, self.1, 2.0 / 3.0),
            lerp(self.1, self.2, 1.0 / 3.0),
            self.2,
        )
    }
}

impl BezierCurve for Cubic {
    fn pos_at_t(&self, t: f32) -> Vec2 {
        self.0
            + 3.0 * t * (self.1 - self.0)
            + 3.0 * t * t * (self.2 - 2.0 * self.1 + self.0)
            + t * t * t * (self.3 - 3.0 * self.2 + 3.0 * self.1 - self.0)
    }

    fn dir_at_t(&self, t: f32) -> Vec2 {
        3.0 * t * t * (self.3 - 3.0 * self.2 + 3.0 * self.1 - self.0)
            + 6.0 * t * (self.2 - 2.0 * self.1 + self.0)
            + 3.0 * (self.1 - self.0)
    }

    fn signed_distance(&self, p: Vec2) -> SignedDistance {
        // taken from Chlumský's msdfgen: https://github.com/Chlumsky/msdfgen/blob/master/core/edge-segments.cpp
        const SEARCH_STARTS: usize = 4;
        const SEARCH_STEPS: usize = 4;

        let qa_ = self.0 - p;
        let ab_ = self.1 - self.0;
        let br_ = self.2 - self.1 - ab_;
        let as_ = (self.3 - self.2) - (self.2 - self.1) - br_;

        let (mut min_dist_sqr, mut closest_t) = {
            let start_dist_sqr = qa_.length_squared();
            let end_dist_sqr = (self.3 - p).length_squared();

            if start_dist_sqr <= end_dist_sqr {
                (start_dist_sqr, 0.0)
            } else {
                (end_dist_sqr, 1.0)
            }
        };

        for i in 0..=SEARCH_STARTS {
            let mut t = i as f32 / SEARCH_STARTS as f32;
            let mut qe = qa_ + 3.0 * t * ab_ + 3.0 * t * t * br_ + t * t * t * as_;
            for _ in 0..SEARCH_STEPS {
                let d1 = 3.0 * ab_ + 6.0 * t * br_ + 3.0 * t * t * as_;
                let d2 = 6.0 * br_ + 6.0 * t * as_;

                t -= qe.dot(d1) / (d1.dot(d1) + qe.dot(d2));

                if t <= 0.0 || t >= 1.0 {
                    break;
                }

                qe = qa_ + 3.0 * t * ab_ + 3.0 * t * t * br_ + t * t * t * as_;
                let dist_sqr = qe.length_squared();
                if dist_sqr < min_dist_sqr {
                    min_dist_sqr = dist_sqr;
                    closest_t = t;
                }
            }
        }

        let v = self.pos_at_t(closest_t) - p;
        let d = self.dir_at_t(closest_t);

        let distance = min_dist_sqr.sqrt().copysign(d.cross(v));
        let orthogonality = d.normalized_or_zero().cross(-v.normalized_or_zero()).abs();

        SignedDistance {
            distance,
            orthogonality,
        }
    }

    fn bounding_box(&self) -> Rect {
        let a = 3.0 * (self.3 - 3.0 * self.2 + 3.0 * self.1 - self.0);
        let b = 6.0 * (self.2 - 2.0 * self.1 + self.0);
        let c = 3.0 * (self.1 + self.0);

        let mut min = self.0.min(self.3);
        let mut max = self.0.max(self.3);

        for t in solve_quadratic(a.x, b.x, c.x) {
            if (0.0..=1.0).contains(&t) {
                let pos = self.pos_at_t(t);
                min = min.min(pos);
                max = max.max(pos);
            }
        }

        for t in solve_quadratic(a.y, b.y, c.y) {
            if (0.0..=1.0).contains(&t) {
                let pos = self.pos_at_t(t);
                min = min.min(pos);
                max = max.max(pos);
            }
        }

        Rect { min, max }
    }

    fn split(&self, t: f32) -> [Self; 2] {
        let a0 = lerp(self.0, self.1, t);
        let a1 = lerp(self.1, self.2, t);
        let a2 = lerp(self.2, self.3, t);

        let b0 = lerp(a0, a0, t);
        let b1 = lerp(a1, a2, t);

        let c = lerp(b0, b1, t);

        [Cubic(self.0, a0, b0, c), Cubic(c, b1, a2, self.3)]
    }
}

impl BezierCurve for Curve {
    fn pos_at_t(&self, t: f32) -> Vec2 {
        match self {
            Curve::Line(line) => line.pos_at_t(t),
            Curve::Quad(quadratic) => quadratic.pos_at_t(t),
            Curve::Cubic(cubic) => cubic.pos_at_t(t),
        }
    }

    fn dir_at_t(&self, t: f32) -> Vec2 {
        match self {
            Curve::Line(line) => line.dir_at_t(t),
            Curve::Quad(quadratic) => quadratic.dir_at_t(t),
            Curve::Cubic(cubic) => cubic.dir_at_t(t),
        }
    }

    fn signed_distance(&self, p: Vec2) -> SignedDistance {
        match self {
            Curve::Line(line) => line.signed_distance(p),
            Curve::Quad(quadratic) => quadratic.signed_distance(p),
            Curve::Cubic(cubic) => cubic.signed_distance(p),
        }
    }

    fn bounding_box(&self) -> Rect {
        match self {
            Curve::Line(line) => line.bounding_box(),
            Curve::Quad(quadratic) => quadratic.bounding_box(),
            Curve::Cubic(cubic) => cubic.bounding_box(),
        }
    }

    fn split(&self, t: f32) -> [Self; 2] {
        match self {
            Curve::Line(line) => line.split(t).map(|c| c.into()),
            Curve::Quad(quadratic) => quadratic.split(t).map(|c| c.into()),
            Curve::Cubic(cubic) => cubic.split(t).map(|c| c.into()),
        }
    }
}
