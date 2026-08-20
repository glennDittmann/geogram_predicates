//! # Geogram Predicates (full Rust port)
//!
//! Dependency-free port of [Geogram](https://github.com/alicevision/geogram)'s robust predicates (PCK).
//! Provides exact orientation, in-circle, in-sphere, and determinant predicates using
//! multi-precision expansion arithmetic.
//!
//! ## Initialization
//! Call [`initialize()`] before using any predicate. Call [`terminate()`] when done (e.g. to print stats).
//!
//! ## Example
//! ```ignore
//! use geogram_predicates::{initialize, terminate, orient_2d, orient_3d, in_circle_2d_sos, Sign};
//!
//! initialize();
//!
//! let a = [0.0, 0.0];
//! let b = [2.0, 0.0];
//! let c = [1.0, 1.0];
//! assert_eq!(orient_2d(&a, &b, &c), Sign::Positive);
//!
//! terminate();
//! ```

mod exact;
mod expansion;
mod filter;
mod orient_2dlifted;
mod orient_3dlifted;
mod pck;
mod sign;
mod types;

pub use orient_2dlifted::orient_2dlifted_sos;
pub use orient_3dlifted::orient_3dlifted_sos;
pub use pck::{
    aligned_3d, det_3d, det_4d, dot_3d, in_circle_2d_sos, in_sphere_3d_sos, orient_2d, orient_3d,
    orient_3d_inexact, points_are_colinear_3d, points_are_identical_2d, points_are_identical_3d,
    Point2d, Point3d, Point4d,
};
pub use sign::geo_sgn;
pub use types::Sign;

static mut INITIALIZED: bool = false;

/// Must be called before using any predicate. Initializes expansion arithmetic constants.
pub fn initialize() {
    crate::expansion::expansion_initialize();
    unsafe {
        INITIALIZED = true;
    }
}

/// Call when done (e.g. to allow printing stats). No-op in this port.
pub fn terminate() {
    unsafe {
        INITIALIZED = false;
    }
}

/// Returns whether [`initialize()`] has been called.
pub fn is_initialized() -> bool {
    unsafe { INITIALIZED }
}

#[cfg(test)]
mod integration_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        if !is_initialized() {
            initialize();
        }
    }

    #[test]
    fn orient_2d_simple() {
        init();
        let a = [0.0, 0.0];
        let b = [2.0, 0.0];
        let c = [1.0, 1.0];
        assert_eq!(orient_2d(&a, &b, &c), Sign::Positive);
        assert_eq!(orient_2d(&a, &c, &b), Sign::Negative);
    }

    #[test]
    fn orient_3d_simple() {
        init();
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let d = [0.75, 0.75, 1.0];
        assert_eq!(orient_3d(&a, &b, &c, &d), Sign::Positive);
    }

    #[test]
    fn det_3d_coplanar() {
        init();
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = [7.0, 8.0, 9.0]; // a + (b-a) + (c-a) = collinear
        assert_eq!(det_3d(&a, &b, &c), Sign::Zero);
    }

    #[test]
    fn in_circle_2d_sos_simple() {
        init();
        let a = [0.0, 0.0];
        let b = [2.0, 0.0];
        let c = [1.0, 1.0];
        let p_in = [1.0, -0.4];
        let p_out = [1.0, -1.2];
        assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_in), Sign::Positive);
        assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_out), Sign::Negative);
    }

    #[test]
    fn in_sphere_3d_sos_simple() {
        init();
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let d = [0.75, 0.75, 1.0];
        let p_in = [0.75, 0.75, 0.5];
        let p_out = [0.75, 0.75, 1.5];
        assert_eq!(in_sphere_3d_sos(&a, &b, &c, &d, &p_in), Sign::Negative);
        assert_eq!(in_sphere_3d_sos(&a, &b, &c, &d, &p_out), Sign::Positive);
    }

    #[test]
    fn dot_3d_simple() {
        init();
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0]; // orthogonal
        assert_eq!(dot_3d(&a, &b, &c), Sign::Zero);
        let d = [1.0, 1.0, 0.0]; // acute with b
        assert_eq!(dot_3d(&a, &b, &d), Sign::Positive);
    }

    #[test]
    fn aligned_3d_simple() {
        init();
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0]; // collinear
        assert!(aligned_3d(&a, &b, &c));
        let d = [0.0, 1.0, 0.0]; // not collinear
        assert!(!aligned_3d(&a, &b, &d));
    }
}
