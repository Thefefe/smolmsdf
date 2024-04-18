pub type SdfBitmap = Bitmap<u8>;
pub type MsdfBitmap = Bitmap<[u8; 3]>;
pub type MtsdfBitmap = Bitmap<[u8; 4]>;

pub struct Bitmap<T> {
    data: Vec<T>,
    dimensions: [u32; 2],
}

impl<T: Copy + Default> Bitmap<T> {
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            data: vec![T::default(); (width * height) as usize],
            dimensions: [width, height],
        }
    }

    pub fn resize(&mut self, new_dimensions: [u32; 2]) {
        self.dimensions = new_dimensions;
        self.data.resize(
            (new_dimensions[0] * new_dimensions[1]) as usize,
            T::default(),
        );
    }

    pub fn dimensions(&self) -> [u32; 2] {
        self.dimensions
    }

    pub fn get_pixel(&self, [x, y]: [u32; 2]) -> T {
        self.data[(x + self.dimensions[0] * y) as usize]
    }

    pub fn get_pixel_mut(&mut self, [x, y]: [u32; 2]) -> &mut T {
        &mut self.data[(x + self.dimensions[0] * y) as usize]
    }

    pub fn set_pixel(&mut self, [x, y]: [u32; 2], v: T) {
        self.data[(x + self.dimensions[0] * y) as usize] = v;
    }

    pub fn enumerate_pixels(&mut self) -> BitmapPixelEnumerateMut<T> {
        BitmapPixelEnumerateMut {
            pixels: self.data.iter_mut(),
            width: self.dimensions[0],
            x: 0,
            y: 0,
        }
    }

    pub fn as_slice(&self) -> &[T] {
        self.data.as_slice()
    }

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
        let mut bitmap: Bitmap<u32> = Bitmap::new(4, 2);

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
