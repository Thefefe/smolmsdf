#![doc = include_str!("../README.md")]

mod math;
pub use math::{Rect, Vec2};

mod curve;
pub use curve::*;

mod shape;
pub use shape::*;

mod bitmap;
pub use bitmap::*;

mod rasterize;
pub use rasterize::*;

mod error_correction;
pub use error_correction::*;
