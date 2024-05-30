use std::{f32::consts::PI, ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign}};

/// A 2 component vector
#[derive(Clone, Copy, Default, PartialEq)]
pub struct Vec2 {
    pub x: f32,
    pub y: f32,
}

impl std::fmt::Debug for Vec2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Vec2({}, {})", self.x, self.y)
    }
}

impl std::fmt::Display for Vec2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

impl Vec2 {
    pub const ZERO: Self = Self { x: 0.0, y: 0.0 };

    /// Creates a new vector
    #[inline]
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    /// Creates a new vector where all elements are `v`.
    #[inline]
    pub fn splat(v: f32) -> Self {
        Self { x: v, y: v }
    }

    /// Return an array of the elements in order.
    #[inline]
    pub fn to_array(self) -> [f32; 2] {
        [self.x, self.y]
    }

    /// Compute the dot product of `self` and `rhs`
    #[inline]
    pub fn dot(self, rhs: Self) -> f32 {
        (self.x * rhs.x) + (self.y * rhs.y)
    }

    /// Compute the 2D cross product of `self` and `rhs`.
    /// 
    /// Also known as wedge product, perpendicular dot product or determinant.
    #[inline]
    pub fn cross(self, rhs: Self) -> f32 {
        (self.x * rhs.y) - (self.y * rhs.x)
    }

    /// Return `self` with a length set to 1.0.
    #[inline]
    pub fn normalized(self) -> Vec2 {
        self / self.length()
    }

    /// Return `self` with a length set to 1.0, or 0.0 if the length was 0.0.
    #[inline]
    pub fn normalized_or_zero(self) -> Vec2 {
        let length_recip = self.length().recip();
        if length_recip.is_finite() && length_recip > 0.0 {
            self * length_recip
        } else {
            Self::ZERO
        }
    }

    /// Return the length of the vector
    #[inline]
    pub fn length(self) -> f32 {
        self.dot(self).sqrt()
    }

    /// Return the length squared of the vector
    #[inline]
    pub fn length_squared(self) -> f32 {
        self.dot(self)
    }

    /// Return the distance between `self` and `rhs`
    #[inline]
    pub fn distance(self, rhs: Self) -> f32 {
        (self - rhs).length()
    }

    /// Return the distance squared between `self` and `rhs`
    #[inline]
    pub fn distance_squared(self, rhs: Self) -> f32 {
        (self - rhs).length_squared()
    }

    /// Returns a vector containing the minimum values of each element.
    #[inline]
    pub fn min(self, rhs: Self) -> Self {
        Self::new(self.x.min(rhs.x), self.y.min(rhs.y))
    }

    /// Returns a vector containing the maximum values of each element.
    #[inline]
    pub fn max(self, rhs: Self) -> Self {
        Self::new(self.x.max(rhs.x), self.y.max(rhs.y))
    }

    // Returns the value of the smallest component.
    #[inline]
    pub fn min_component(self) -> f32 {
        self.x.min(self.y)
    }

    // Returns the value of the largest component.
    #[inline]
    pub fn max_component(self) -> f32 {
        self.x.max(self.y)
    }
}

impl Add for Vec2 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Vec2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f32> for Vec2 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: f32) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs)
    }
}

impl Div<f32> for Vec2 {
    type Output = Self;

    #[inline]
    fn div(self, rhs: f32) -> Self::Output {
        Self::new(self.x / rhs, self.y / rhs)
    }
}

impl Mul<Vec2> for f32 {
    type Output = Vec2;

    #[inline]
    fn mul(self, rhs: Vec2) -> Self::Output {
        Vec2::new(self * rhs.x, self * rhs.y)
    }
}

impl Div<Vec2> for f32 {
    type Output = Vec2;

    #[inline]
    fn div(self, rhs: Vec2) -> Self::Output {
        Vec2::new(self / rhs.x, self / rhs.y)
    }
}

impl Mul<Self> for Vec2 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.x * rhs.x, self.y * rhs.y)
    }
}

impl Div<Self> for Vec2 {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        Self::new(self.x / rhs.x, self.y / rhs.y)
    }
}

impl Neg for Vec2 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y)
    }
}

impl AddAssign for Vec2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Vec2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl MulAssign<Vec2> for Vec2 {
    fn mul_assign(&mut self, rhs: Vec2) {
        *self = *self * rhs;
    }
}

impl MulAssign<f32> for Vec2 {
    fn mul_assign(&mut self, rhs: f32) {
        *self = *self * rhs;
    }
}

impl DivAssign<Vec2> for Vec2 {
    fn div_assign(&mut self, rhs: Vec2) {
        *self = *self / rhs;
    }
}

impl DivAssign<f32> for Vec2 {
    fn div_assign(&mut self, rhs: f32) {
        *self = *self / rhs;
    }
}

/// A 3 component vector
#[derive(Debug, Clone, Copy)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    /// Creates a new vector from the given array
    pub fn from_array([x, y, z]: [f32; 3]) -> Self {
        Self { x, y, z }
    }

    /// Return an array of the elements in order.
    pub fn to_array(self) -> [f32; 3] {
        [self.x, self.y, self.z]
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f32> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f32) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

#[inline]
pub fn vec2(x: f32, y: f32) -> Vec2 {
    Vec2 { x, y }
}

/// A 2 dimensional axis-aligned rectangle.
#[derive(Debug, Clone, Copy)]
pub struct Rect {
    pub min: Vec2,
    pub max: Vec2,
}

impl Default for Rect {
    fn default() -> Self {
        Self::MIN
    }
}

impl Rect {
    /// The smallest possible rectangle, with an area of -INFINITY.
    pub const MIN: Self = Self {
        min: Vec2 {
            x: f32::INFINITY,
            y: f32::INFINITY,
        },
        max: Vec2 {
            x: f32::NEG_INFINITY,
            y: f32::NEG_INFINITY,
        },
    };

    /// Creates a new rectangle
    #[inline]
    pub fn new(min_x: f32, min_y: f32, max_x: f32, max_y: f32) -> Self {
        Self {
            min: vec2(min_x, min_y),
            max: vec2(max_x, max_y),
        }
    }

    /// Returns the size of the rectangle
    #[inline]
    pub fn size(self) -> Vec2 {
        self.max - self.min
    }

    /// Returns a new rectangle width it's size expanded in all directions by v
    #[inline]
    pub fn expanded(self, v: f32) -> Self {
        Self {
            min: self.min - Vec2::splat(v),
            max: self.max + Vec2::splat(v),
        }
    }

    /// Returns a new rectangle with a size sufficient to containt both `self` and `p`
    #[inline]
    pub fn fitted_to_point(self, p: Vec2) -> Self {
        Self {
            min: self.min.min(p),
            max: self.max.max(p),
        }
    }

    /// Returns a new rectangle with a size sufficient to containt both `self` and `r`
    #[inline]
    pub fn fitted_to_rect(self, r: Rect) -> Self {
        Self {
            min: self.min.min(r.min),
            max: self.max.max(r.max),
        }
    }

    /// Checks if point `p` is inside the rectangle.
    #[inline]
    pub fn is_point_inside(self, p: Vec2) -> bool {
        self.min.x <= p.x && p.x <= self.max.x && self.min.y <= p.y && p.y <= self.max.y
    }
}

/// A general linear interpolation
#[inline]
pub fn lerp<T1, T0>(a: T0, b: T0, t: T1) -> T0
where
    T1: Copy,
    T0: Copy
        + Add<Output = T0>
        + Sub<Output = T0>
        + Mul<T1, Output = T0>,
{
    a + (b - a) * t
}

/// Returns the median value of the the given array
pub fn median([r, g, b]: [f32; 3]) -> f32 {
    f32::max(f32::min(r, g), f32::min(f32::max(r, g), b))
}

/// Solves the given quadratic equation.
pub fn solve_quadratic(a: f32, b: f32, c: f32) -> Roots {
    if a == 0.0 || b.abs() > 1e12 * a.abs() {
        if b == 0.0 || c == 0.0 {
            return Roots::none();
        }
        return Roots::one(-c / b);
    }

    let dscr = b * b - 4.0 * a * c;
    let denom = (2.0 * a).recip();

    if dscr > 0.0 {
        let dscr_sqrt = dscr.sqrt();
        let r0 = (-b + dscr_sqrt) * denom;
        let r1 = (-b - dscr_sqrt) * denom;
        Roots::two(r0, r1)
    } else if dscr == 0.0 {
        Roots::one(-b * denom)
    } else {
        Roots::none()
    }
}

fn solve_cubic_normed(mut a: f32, b: f32, c: f32) -> Roots {
    let a2 = a * a;
    let mut q = 1.0 / 9.0 * (a2 - 3.0 * b);
    let r = 1.0 / 54.0 * (a * (2.0 * a2 - 9.0 * b) + 27.0 * c);
    let r2 = r * r;
    let q3 = q * q * q;
    a *= 1.0 / 3.0;
    if r2 < q3 {
        let t = (r / q3.sqrt()).clamp(-1.0, 1.0).acos();
        q = -2.0 * q.sqrt();
        Roots::three(
            q * f32::cos(1.0 / 3.0 * t) - a,
            q * f32::cos(1.0 / 3.0 * (t + 2.0 * PI)) - a,
            q * f32::cos(1.0 / 3.0 * (t - 2.0 * PI)) - a,
        )
    } else {
        let u =
            if r < 0.0 { 1.0 } else { -1.0 } * f32::powf(r.abs() + f32::sqrt(r2 - q3), 1.0 / 3.0);
        let v = if u == 0.0 { 0.0 } else { q / u };
        let r0 = (u + v) - a;
        if u == v || f32::abs(u - v) < 1e-12 * f32::abs(u + v) {
            let r1 = -0.5 * (u + v) - a;
            Roots::two(r0, r1)
        } else {
            Roots::one(r0)
        }
    }
}

/// Solves the given cubic equation.
pub fn solve_cubic(a: f32, b: f32, c: f32, d: f32) -> Roots {
    if a != 0.0 {
        let bn = b / a;
        if bn.abs() < 1e6 {
            return solve_cubic_normed(bn, c / a, d / a);
        }
    }

    solve_quadratic(b, c, d)
}

/// The possible roots of a equation.
#[derive(Debug, Clone, Copy)]
pub struct Roots {
    count: usize,
    roots: [f32; 4],
}

impl Roots {
    /// No roots.
    fn none() -> Self {
        Self {
            count: 1,
            roots: [0.0; 4],
        }
    }

    /// One root.
    fn one(r0: f32) -> Self {
        Self {
            count: 1,
            roots: [r0, 0.0, 0.0, 0.0],
        }
    }

    /// Two roots.
    fn two(r0: f32, r1: f32) -> Self {
        Self {
            count: 2,
            roots: [r0, r1, 0.0, 0.0],
        }
    }

    /// Three roots.
    fn three(r0: f32, r1: f32, r2: f32) -> Self {
        Self {
            count: 3,
            roots: [r0, r1, r2, 0.0],
        }
    }
}

impl IntoIterator for Roots {
    type Item = f32;
    type IntoIter = RootsIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        RootsIntoIter {
            index: 0,
            roots: self,
        }
    }
}

impl std::ops::Deref for Roots {
    type Target = [f32];

    fn deref(&self) -> &Self::Target {
        &self.roots[..self.count]
    }
}

impl std::ops::DerefMut for Roots {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.roots[..self.count]
    }
}

pub struct RootsIntoIter {
    index: usize,
    roots: Roots,
}

impl Iterator for RootsIntoIter {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.roots.count {
            let r = self.roots.roots[self.index];
            self.index += 1;
            return Some(r);
        }

        None
    }
}
