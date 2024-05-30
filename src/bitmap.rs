/// A bitmap for storing single-channel true distance fields
pub type SdfBitmap = Bitmap<u8>;
/// A bitmap for storing a multi channel signed distance fields in the rgb components
pub type MsdfBitmap = Bitmap<[u8; 3]>;
/// A bitmap for storing a multi channel signed distance fields in the rgb components,
/// and a single-channel true distance field in the alpha component
pub type MtsdfBitmap = Bitmap<[u8; 4]>;

/// A 2 dimensional bitmap, to rasterize shapes to.
pub struct Bitmap<T> {
    data: Vec<T>,
    dimensions: [u32; 2],
}

impl<T: Copy + Default> Bitmap<T> {
    /// Creates a new bitmap with a resolution of 0x0.
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            dimensions: [0; 2],
        }
    }

    /// Creates a new bitmap at the given resolution, and all pixels set to 0.
    pub fn with_size(width: u32, height: u32) -> Self {
        Self {
            data: vec![T::default(); (width * height) as usize],
            dimensions: [width, height],
        }
    }

    /// Resizes the bitmap to the provided dimensions. Additional pixels will be set to the default of the pixel type.
    pub fn resize(&mut self, new_dimensions: [u32; 2]) {
        self.dimensions = new_dimensions;
        self.data.resize(
            (new_dimensions[0] * new_dimensions[1]) as usize,
            T::default(),
        );
    }

    /// Returns the width and height of the bitmap.
    pub fn dimensions(&self) -> [u32; 2] {
        self.dimensions
    }

    /// Gets the pixel at the given coordinate.
    /// 
    /// Panics if the coordinate is outside of the bitmap.
    pub fn get_pixel(&self, [x, y]: [u32; 2]) -> T {
        self.data[(x + self.dimensions[0] * y) as usize]
    }
    
    /// Gets a mutable reference to the pixel at the given coordinate.
    /// 
    /// Panics if the coordinate is outside of the bitmap.
    pub fn get_pixel_mut(&mut self, [x, y]: [u32; 2]) -> &mut T {
        &mut self.data[(x + self.dimensions[0] * y) as usize]
    }

    /// Sets the pixel at the given coordinate.
    /// 
    /// Panics if the coordinate is outside of the bitmap.
    pub fn set_pixel(&mut self, [x, y]: [u32; 2], v: T) {
        self.data[(x + self.dimensions[0] * y) as usize] = v;
    }

    /// Returns an iterator that goes through the bitmap in a linear fashion
    /// and yields the coordinate and a mutable reference to the pixels.
    pub fn enumerate_pixels(&mut self) -> BitmapPixelEnumerateMut<T> {
        BitmapPixelEnumerateMut {
            pixels: self.data.iter_mut(),
            width: self.dimensions[0],
            x: 0,
            y: 0,
        }
    }

    /// Returns a slice to the underlying pixel data.
    pub fn as_slice(&self) -> &[T] {
        self.data.as_slice()
    }
    
    /// Returns a mutable slice to the underlying pixel data.
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        self.data.as_mut_slice()
    }
}

pub struct BitmapPixelEnumerateMut<'a, T> {
    pixels: std::slice::IterMut<'a, T>,
    width: u32,
    x: u32,
    y: u32,
}

impl<'a, T> Iterator for BitmapPixelEnumerateMut<'a, T> {
    type Item = (u32, u32, &'a mut T);

    fn next(&mut self) -> Option<Self::Item> {
        let pixel = self.pixels.next()?;

        if self.x >= self.width {
            self.x = 0;
            self.y += 1;
        }

        let (x, y) = (self.x, self.y);

        self.x += 1;

        Some((x, y, pixel))
    }
}

#[cfg(test)]
mod tests {
    use super::Bitmap;

    #[test]
    fn bitmap() {
        let mut bitmap: Bitmap<u32> = Bitmap::with_size(4, 2);

        bitmap.set_pixel([0, 0], 1);
        bitmap.set_pixel([1, 0], 2);
        bitmap.set_pixel([2, 0], 3);
        bitmap.set_pixel([3, 0], 4);

        bitmap.set_pixel([0, 1], 5);
        bitmap.set_pixel([1, 1], 6);
        bitmap.set_pixel([2, 1], 7);
        bitmap.set_pixel([3, 1], 8);

        let mut iter = bitmap.enumerate_pixels();

        assert_eq!(iter.next(), Some((0, 0, &mut 1)));
        assert_eq!(iter.next(), Some((1, 0, &mut 2)));
        assert_eq!(iter.next(), Some((2, 0, &mut 3)));
        assert_eq!(iter.next(), Some((3, 0, &mut 4)));

        assert_eq!(iter.next(), Some((0, 1, &mut 5)));
        assert_eq!(iter.next(), Some((1, 1, &mut 6)));
        assert_eq!(iter.next(), Some((2, 1, &mut 7)));
        assert_eq!(iter.next(), Some((3, 1, &mut 8)));

        assert_eq!(bitmap.get_pixel([0, 0]), 1);
        assert_eq!(bitmap.get_pixel([1, 1]), 6);
    }
}
