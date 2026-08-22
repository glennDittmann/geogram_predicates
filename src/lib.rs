#![forbid(unsafe_code)]

//! # Geogram Predicates (full Rust port)
//!
//! Dependency-free port of [Geogram](https://github.com/alicevision/geogram)'s robust predicates (PCK).
//! Provides exact orientation, in-circle, in-sphere, and determinant predicates using
//! multi-precision expansion arithmetic.
//!
//! Predicates are ready to use immediately. [`initialize`] and [`terminate`]
//! remain as compatibility no-ops.
//!
//! ## Example
//! ```ignore
//! use geogram_predicates::{orient_2d, Sign};
//!
//! let a = [0.0, 0.0];
//! let b = [2.0, 0.0];
//! let c = [1.0, 1.0];
//! assert_eq!(orient_2d(&a, &b, &c), Sign::Positive);
//! ```

mod exact;
mod expansion;
mod filter;
mod lifted;
mod pck;
mod side;
mod sign;
mod types;

pub use lifted::{orient_2dlifted, orient_2dlifted_sos, orient_3dlifted, orient_3dlifted_sos};
pub use pck::{
    aligned_3d, det_3d, det_4d, det_compare_4d, dot_3d, dot_compare_3d, in_circle_2d,
    in_circle_2d_sos, in_circle_3d, in_circle_3d_sos, in_circle_3dlifted, in_circle_3dlifted_sos,
    in_sphere_3d, in_sphere_3d_sos, orient_2d, orient_3d, orient_3d_inexact,
    points_are_colinear_3d, points_are_identical_2d, points_are_identical_3d, Point2d, Point3d,
    Point4d,
};
pub use side::{
    side1_sos, side2_sos, side3_3dlifted_sos, side3_sos, side4_3d, side4_3d_sos, side4_sos,
};
pub use sign::geo_sgn;
pub use types::{Sign, SosPoint};

/// Compatibility no-op; predicates need no initialization.
pub fn initialize() {}

/// Compatibility no-op.
pub fn terminate() {}

/// Compatibility no-op. Statistics are not collected by the Rust port.
pub fn show_stats() {}

/// Always returns `true`; predicates are statically initialized.
pub fn is_initialized() -> bool {
    true
}

#[cfg(test)]
mod integration_tests;
