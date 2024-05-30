use super::{
    bitmap::Bitmap,
    math::{lerp, median, Vec3},
};

/// The simplest error correction algorithm. While it's very fast and doesn't use any additional allocations,
/// it doesn't correct all of the errors, however in most cases it removes the most noticable artifacts.
/// 
/// It works by detecting when in the interpolation between two pixels the sign of the median value changes,
/// but the signs of the medians in the ends match.
/// When this is detected both pixels' msdf components are set to contain the median values in all channels. 
pub fn correct_errors_interp_sign<P: Msd>(image: &mut Bitmap<P>) {
    let [width, height] = image.dimensions();

    for y in 0..height - 1 {
        for x in 0..width - 1 {
            let p0 = image.get_pixel([x, y]).signed_distances();
            let p1 = image.get_pixel([x + 1, y]).signed_distances();
            let p2 = image.get_pixel([x, y + 1]).signed_distances();

            let hcol = interpolation_error_simple(p0, p1);
            let vcol = interpolation_error_simple(p0, p2);

            if hcol {
                image.get_pixel_mut([x, y]).set_median();
                image.get_pixel_mut([x + 1, y]).set_median();
            }

            if vcol {
                image.get_pixel_mut([x, y]).set_median();
                image.get_pixel_mut([x, y + 1]).set_median();
            }
        }
    }
}

fn interpolation_error_simple(a: Vec3, b: Vec3) -> bool {
    let am = median(a.to_array());
    let bm = median(b.to_array());

    if am.is_sign_negative() != bm.is_sign_negative() {
        return false;
    }

    interpolation_intersections(a, b).any(|t| {
        let true_median = median(lerp(a, b, t).to_array());
        let interp_median = lerp(am, bm, t);
        true_median.is_sign_negative() != interp_median.is_sign_negative()
    })
}

fn interpolation_intersections(a: Vec3, b: Vec3) -> impl Iterator<Item = f32> {
    #[inline]
    fn lerp_intersection(ax: f32, bx: f32, ay: f32, by: f32) -> f32 {
        (ax - ay) / ((by - ay) - (bx - ax))
    }

    let xy = lerp_intersection(a.x, b.x, a.y, b.y);
    let yz = lerp_intersection(a.y, b.y, a.z, b.z);
    let xz = lerp_intersection(a.x, b.x, a.z, b.z);

    [xy, yz, xz]
        .into_iter()
        .filter(move |&t| 0.0 < t && t < 1.0)
}

/// Something that can store a multi-channel distance field.
pub trait Msd: Copy + Default {
    /// Returns the stored multi-channel distance.
    fn signed_distances(self) -> Vec3;

    // Sets the stored multi-channel distance to the median distance.
    fn set_median(&mut self);
}

impl Msd for [u8; 4] {
    fn signed_distances(self) -> Vec3 {
        Vec3::from_array([self[0], self[1], self[2]].map(|u| u as f32 / 255.0 * 2.0 - 1.0))
    }

    fn set_median(&mut self) {
        let med = u8::max(
            u8::min(self[0], self[1]),
            u8::min(u8::max(self[0], self[1]), self[2]),
        );
        *self = [med, med, med, self[3]]
    }
}

impl Msd for [u8; 3] {
    fn signed_distances(self) -> Vec3 {
        Vec3::from_array([self[0], self[1], self[2]].map(|u| u as f32 / 255.0 * 2.0 - 1.0))
    }

    fn set_median(&mut self) {
        let med = u8::max(
            u8::min(self[0], self[1]),
            u8::min(u8::max(self[0], self[1]), self[2]),
        );
        *self = [med; 3]
    }
}
