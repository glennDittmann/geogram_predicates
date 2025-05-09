//! # Geogram Predicates
//!
//! A Rust port of `geogram`s _robust predicates_.
//! An interoperability ffi with `geogram`s _robust predicates_; via `cxx` also available in the `legacy` feature.
//!
//! All functions are in written in Rust, except for `det_3d` and `det_4d` which are only available in `legacy`.
#![no_std]

#![warn(clippy::all, unused, clippy::missing_const_for_fn)]

#[cfg(test)]
mod tests;

#[cfg(feature = "legacy")]
pub use geogram_ffi::*;

use nalgebra::Point3;
use robust::{Coord, Coord3D};

mod types;
pub use types::Sign;

mod expansion;
pub use expansion::expansion::Expansion;

mod orient_2dlifted;
pub use orient_2dlifted::orient_2dlifted_sos;

mod orient_3dlifted;
pub use orient_3dlifted::orient_3dlifted_sos;

pub type Point3d = [f64; 3];
pub type Point2d = [f64; 2];

pub(crate) const FPG_UNCERTAIN_VALUE: Sign = Sign::Zero;

/// Gets the sign of a value.
///
/// ### Parameters
/// - `x` value to test
///
/// # Example
/// ```
/// use geogram_predicates::{Sign, geo_sign};
///
/// let a = 42.0;
/// let b = -42.0;
/// let c = 0.0;
///
/// assert_eq!(geo_sign(a), Sign::Positive);
/// assert_eq!(geo_sign(b), Sign::Negative);
/// assert_eq!(geo_sign(c), Sign::Zero);
/// ```
#[inline(always)]
#[must_use]
pub const fn geo_sign(value: f64) -> Sign {
    if value > 0.0 {
        Sign::Positive
    } else if value < 0.0 {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

/// Computes the orientation predicate in 3d.
///
/// Computes the sign of the signed volume of the tetrahedron `a`, `b`, `c`, `d`.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
///
/// ### Return values
/// * `Positive` - if the tetrahedron is oriented positively
/// * `Zero` - if the tetrahedron is flat
/// * `Negative` - if the tetrahedron is oriented negatively
///
/// # Example
/// ```
/// use geogram_predicates::orient_3d;
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// assert_eq!(1, orient_3d(&a, &b, &c, &d));
///```
#[must_use]
pub fn orient_3d(a: &Point3d, b: &Point3d, c: &Point3d, d: &Point3d) -> Sign {
    let orientation = robust::orient3d(
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    );

    geo_sign(-orientation)
}

/// Computes the orientation predicate in 2d.
///
/// Computes the sign of the signed area of the triangle `a`, `b`, `c`.
///
/// ### Parameters
/// - `a`, `b`, `c` vertices of the triangle
///
/// ### Return values
/// * `Positive` - if the triangle is oriented counter-clockwise
/// * `Zero` - if the triangle is flat
/// * `Negative` - if the triangle is oriented clockwise
///
/// # Example
/// ```
/// use geogram_predicates::{Sign, orient_2d};
///
/// // Define three points that form a triangle
/// let a = [0.0, 0.0];
/// let b = [2.0, 0.0];
/// let c = [1.0, 1.0];
///
/// let orientation = orient_2d(&a, &b, &c);
/// assert_eq!(Sign::Positive, orientation);
/// ```
#[must_use]
pub fn orient_2d(a: &Point2d, b: &Point2d, c: &Point2d) -> Sign {
    let orientation = robust::orient2d(
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*a) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*b) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*c) },
    );

    geo_sign(orientation)
}

/// Computes the sign of the dot product between two vectors.
///
/// ### Parameters
/// - `a`, `b`, `c`, three 3d points
///
/// ### Returns
/// - the sign of the dot product between the vectors `ab` and `ac`, true is positive
///
/// # Example
/// ```
/// use geogram_predicates::dot_3d;
///
/// // Define four points that form a matrix
/// let a = [0.0, 0.0, 0.0];
/// let b = [1.0, 0.0, 0.0];
/// let c = [0.0, 1.0, 0.0];
///
/// let dot_sign = dot_3d(&a, &b, &c); // should be orthogonal
/// assert_eq!(dot_sign, true);
/// ```
#[must_use]
pub fn dot_3d(a: &Point3d, b: &Point3d, c: &Point3d) -> bool {
    let ab_diff = Into::<Point3<_>>::into(*a) - Into::<Point3<_>>::into(*b);
    let ab_distance = ab_diff.x.hypot(ab_diff.y).hypot(ab_diff.z);
    let ac_diff = Into::<Point3<_>>::into(*a) - Into::<Point3<_>>::into(*c);
    let ac_distance = ac_diff.x.hypot(ac_diff.y).hypot(ac_diff.z);

    nalgebra::vector![ab_distance].dot(&nalgebra::vector![ac_distance]) >= 0.0
}

/// Tests whether a point is in the circum-circle of a triangle.
///
/// If the triangle `a` , `b` , `c` is oriented clockwise instead of counter-clockwise, then the result is inversed.
///
/// ### Parameters
/// - `a`, `b`, `c` vertices of the triangle
/// - `p` point to test
/// - `PERTURB` (const) - true is `Positive` false is `Negative` consistent perturbation
///
/// ### Return values
/// * `Positive` - if `p` is inside the circum-circle of `a`, `b`, `c`
/// * `Negative` - if `p` is outside the circum-circle of `a`, `b`, `c`
/// * `PERTURB` - if `p` is exactly on the circum-circle of the triangle `a`, `b`, `c`, where `perturb()` denotes a consistent perturbation, that returns either `Positive` or `Negative`
///
/// # Example
/// ```
/// use geogram_predicates::{Sign, in_circle_2d_sos};
///
/// // Define three points that form a triangle
/// let a = [0.0, 0.0];
/// let b = [2.0, 0.0];
/// let c = [1.0, 1.0];
///
/// // Define two points, to test against the triangles circum-circle
/// let p_in = [1.0, -0.4];
/// let p_out = [1.0, -1.2];
///
/// let is_in_circle_p_in = in_circle_2d_sos::<false>(&a, &b, &c, &p_in);
/// assert_eq!(Sign::Positive, is_in_circle_p_in);
/// # let is_in_circle_p_in = in_circle_2d_sos::<true>(&a, &b, &c, &p_in);
/// # assert_eq!(Sign::Positive, is_in_circle_p_in);
///
/// let is_in_circle_p_out = in_circle_2d_sos::<true>(&a, &b, &c, &p_out);
/// assert_eq!(Sign::Negative, is_in_circle_p_out);
/// # let is_in_circle_p_out = in_circle_2d_sos::<false>(&a, &b, &c, &p_out);
/// # assert_eq!(Sign::Negative, is_in_circle_p_out);
/// ```
#[must_use]
pub fn in_circle_2d_sos<const PERTURB: bool>(a: &Point2d, b: &Point2d, c: &Point2d, p: &Point2d) -> Sign {
    let incircle = robust::incircle(
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*a) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*b) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*c) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*p) },
    );

    if incircle > 0.0 {
        Sign::Positive
    } else if incircle < 0.0 {
        Sign::Negative
    } else {
        const { if PERTURB { Sign::Positive } else { Sign::Negative } }
    }
}

/// Tests whether a point is in the circum-sphere of a tetrahedron.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
/// - `p` point to test
/// - `PERTURB` (const) - what it should be if `p` is exactly on the circum-sphere, true is `Positive`, false is `Negative`
///
/// ### Return values
/// * `Negative` - if `p` is inside the circum-sphere of `a`, `b`, `c`, `d`
/// * `Positive` - if `p` is outside the circum-sphere of `a`, `b`, `c`, `d`
/// * `PERTURB` - if `p` is exactly on the circum-sphere of the tetrahedron `a`, `b`, `c`, `d`, where `perturb` denotes a consistent perturbation, that returns either `Positive` or `Negative`
///
/// # Example
/// ```
/// use geogram_predicates::{Sign, in_sphere_3d_sos};
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// // Define two points, to test against the tetrahedrons circum-sphere
/// let p_in = [0.75, 0.75, 0.5];
/// assert_eq!(Sign::Negative, in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_in));
/// # assert_eq!(Sign::Negative, in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_in));
///
/// let p_out = [0.75, 0.75, 1.5];
/// assert_eq!(1, in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_out));
/// # assert_eq!(1, in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_out));
/// ```
#[must_use]
pub fn in_sphere_3d_sos<const PERTURB: bool>(
    a: &Point3d,
    b: &Point3d,
    c: &Point3d,
    d: &Point3d,
    p: &Point3d,
) -> Sign {
    // this keeps the orientation in order as per robust docs
    // let insphere = if orient_3d(a, b, c, d) == 1 {
    //     robust::insphere(
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    //     )
    // } else {
    //     robust::insphere(
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    //     )
    // };

    let insphere = robust::insphere(
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    );

    if insphere > 0.0 {
        Sign::Positive
    } else if insphere < 0.0 {
        Sign::Negative
    } else {
        const { if PERTURB { Sign::Positive } else { Sign::Negative } }
    }
}

// /// Computes the sign of the determinant of a 3x3 matrix formed by three 3D points.
// ///
// /// ### Parameters
// /// - `a`, `b`, `c` the three points that form the matrix
// ///
// /// ### Returns
// /// - the sign of the determinant of the matrix
// ///
// /// # Example
// /// ```
// /// use geogram_predicates as gp;
// ///
// /// // Define three points that form a matrix
// /// let a = [1.0, 2.0, 3.0];
// /// let b = [4.0, 5.0, 6.0];
// /// let c = [7.0, 8.0, 9.0];
// ///
// /// let det = gp::det_3d(&a, &b, &c);
// /// assert_eq!(det, 0);
// /// ```
// #[must_use]
// pub fn det_3d(a: &Point3d, b: &Point3d, c: &Point3d) -> i8 {
//     let res = det_3d_filter(a, b, c);
//     // FIXME: this breaks tests but it's what geogram does
//     if res == 0 { det_3d_exact(a, b, c) } else { res }
// }

// #[inline]
// const fn det_3d_filter(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> i8 {
//     let delta = ((p0[0] * ((p1[1] * p2[2]) - (p1[2] * p2[1]))) - (p1[0] * ((p0[1] * p2[2]) - (p0[2] * p2[1]))))
//         + (p2[0] * ((p0[1] * p1[2]) - (p0[2] * p1[1])));

//     let max1 = p0[0].abs().max(p1[0].abs()).max(p2[0].abs());
//     let max2 = p0[1].abs().max(p0[2].abs()).max(p1[1].abs()).max(p1[2].abs());
//     let max3 = p1[1].abs().max(p1[2].abs()).max(p2[1].abs()).max(p2[2].abs());

//     let mut lower_bound_1 = max1;
//     let mut upper_bound_1 = max1;

//     if max2 < lower_bound_1 {
//         lower_bound_1 = max2;
//     } else if max2 > upper_bound_1 {
//         upper_bound_1 = max2;
//     }
//     if max3 < lower_bound_1 {
//         lower_bound_1 = max3;
//     } else if max3 > upper_bound_1 {
//         upper_bound_1 = max3;
//     }

//     if lower_bound_1 < 1.92663387981871579179e-98 || upper_bound_1 > 1.11987237108890185662e+102 {
//         FPG_UNCERTAIN_VALUE as i8
//     } else {
//         let eps = 3.11133555671680765034e-15 * ((max2 * max3) * max1);
//         if delta > eps {
//             1
//         } else if delta < -eps {
//             -1
//         } else {
//             FPG_UNCERTAIN_VALUE as i8
//         }
//     }
// }

// /// Computes the sign of the determinant of a 3x3
// /// matrix formed by three 3d points using exact arithmetics.
// ///
// /// ### Parameters
// /// p0 , p1 , p2 the three points
// ///
// /// ### Returns
// /// The sign of the determinant of the matrix.
// fn det_3d_exact(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> i8 {
//     let p0_0 = Expansion::from(p0[0]);
//     let p0_1 = Expansion::from(p0[1]);
//     let p0_2 = Expansion::from(p0[2]);

//     let p1_0 = Expansion::from(p1[0]);
//     let p1_1 = Expansion::from(p1[1]);
//     let p1_2 = Expansion::from(p1[2]);

//     let p2_0 = Expansion::from(p2[0]);
//     let p2_1 = Expansion::from(p2[1]);
//     let p2_2 = Expansion::from(p2[2]);

//     let delta: Expansion<24> = crate::expansion::expansion_det3x3!(
//         [p0_0, p0_1, p0_2],
//         [p1_0, p1_1, p1_2],
//         [p2_0, p2_1, p2_2],
//     );

//     delta.sign() as i8
// }

#[cfg(feature = "legacy")]
#[cxx::bridge(namespace = "GEOGRAM")]
mod geogram_ffi {
    // C++ types and signatures exposed to Rust.
    unsafe extern "C++" {
        include!("geogram_predicates/include/geogram_ffi.h");

        /// Computes the sign of the determinant of a 4x4 matrix formed by four 4D points.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c`, `d` the four points that form the matrix
        ///
        /// ### Returns
        /// - the sign of the determinant of the matrix
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define four points that form a matrix
        /// let a = [1.0, 2.0, 3.0, 4.0];
        /// let b = [5.0, 6.0, 7.0, 8.0];
        /// let c = [9.0, 10.0, 11.0, 12.0];
        /// let d = [13.0, 14.0, 15.0, 16.0];
        ///
        /// let det = gp::det_4d(&a, &b, &c, &d);
        /// assert_eq!(det, 0);
        /// ```
        fn det_4d(a: &[f64; 4], b: &[f64; 4], c: &[f64; 4], d: &[f64; 4]) -> i16;

        /// Computes the sign of the determinant of a 3x3 matrix formed by three 3D points.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c` the three points that form the matrix
        ///
        /// ### Returns
        /// - the sign of the determinant of the matrix
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points that form a matrix
        /// let a = [1.0, 2.0, 3.0];
        /// let b = [4.0, 5.0, 6.0];
        /// let c = [7.0, 8.0, 9.0];
        ///
        /// let det = gp::det_3d(&a, &b, &c);
        /// assert_eq!(det, 0);
        /// ```
        fn det_3d(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> i16;

        /// Needs to be called before using any predicate.
        fn initialize();

        /// Displays some statistics about predicates, including the number of calls, the number of exact arithmetics calls, and the number of Simulation of Simplicity calls.
        fn show_stats();

        /// Needs to be called at the end of the program.
        fn terminate();
    }
}

/// Tests whether three 3D points are colinear.
///
/// ### Parameters
/// - `p1` first point
/// - `p2` second point
/// - `p3` third point
///
/// ### Return values
/// - `true` - if `p1`, `p2` and `p3` are colinear
/// - `false` - otherwise
///
/// # Example
/// ```
/// use geogram_predicates::points_are_colinear_3d;
///
/// // Define three points on a line
/// let p1 = [0.0, 0.0, 0.0];
/// let p2 = [0.0, 0.0, 1.0];
/// let p3 = [0.0, 0.0, 2.0];
///
/// assert!(points_are_colinear_3d(&p1, &p2, &p3));
/// ```
#[must_use]
pub fn points_are_colinear_3d(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> bool {
    // Colinearity is tested by using four coplanarity
    // tests with four points that are not coplanar.
    const Q000: [f64; 3] = [0.0, 0.0, 0.0];
    const Q001: [f64; 3] = [0.0, 0.0, 1.0];
    const Q010: [f64; 3] = [0.0, 1.0, 0.0];
    const Q100: [f64; 3] = [1.0, 0.0, 0.0];

    orient_3d(p1, p2, p3, &Q000) == 0
    && orient_3d(p1, p2, p3, &Q001) == 0
    && orient_3d(p1, p2, p3, &Q010) == 0
    && orient_3d(p1, p2, p3, &Q100) == 0
}

/// Tests whether two 2d points are identical.
///
/// ### Parameters
/// - `p1` first point
/// - `p2` second point
///
/// ### Return values
/// - `true` - if `p1` and `p2` have exactly the same coordinates
/// - `false` - otherwise
///
/// # Example
/// ```
/// use geogram_predicates::points_are_identical_2d;
///
/// let p1 = [4.0, 2.0];
/// let p2 = [4.0, 2.0];
///
/// assert!(points_are_identical_2d(&p1, &p2));
/// ```
pub const fn points_are_identical_2d(p1: &Point2d, p2: &Point2d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1]
}

/// Tests whether two 3d points are identical.
///
/// ### Parameters
/// - `p1` first point
/// - `p2` second point
///
/// ### Return values
/// - `true` - if `p1` and `p2` have exactly the same coordinates
/// - `false` - otherwise
///
/// # Example
/// ```
/// use geogram_predicates::points_are_identical_3d;
///
/// let p1 = [4.0, 2.0, 0.42];
/// let p2 = [4.0, 2.0, 0.42];
///
/// assert!(points_are_identical_3d(&p1, &p2));
/// ```
pub const fn points_are_identical_3d(p1: &Point3d, p2: &Point3d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]
}

/// Computes the (approximate) orientation predicate in 3d.
///
/// Computes the sign of the (approximate) signed volume of the tetrahedron `a`, `b`, `c`, `d`.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
///
/// ### Return values
/// * `Positive` - if the tetrahedron is oriented positively
/// * `Zero` - if the tetrahedron is flat
/// * `Negative` - if the tetrahedron is oriented negatively
///
/// # Example
/// ```
/// use geogram_predicates as gp;
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// assert_eq!(1, gp::orient_3d_inexact(&a, &b, &c, &d));
///```
#[must_use]
pub const fn orient_3d_inexact(a: &Point3d, b: &Point3d, c: &Point3d, d: &Point3d) -> Sign {
    let a11 = b[0] - a[0];
    let a12 = b[1] - a[1];
    let a13 = b[2] - a[2];

    let a21 = c[0] - a[0];
    let a22 = c[1] - a[1];
    let a23 = c[2] - a[2];

    let a31 = d[0] - a[0];
    let a32 = d[1] - a[1];
    let a33 = d[2] - a[2];

    let delta =  det3x3(
        [a11, a12, a13],
        [a21, a22, a23],
        [a31, a32, a33],
    );

    geo_sign(delta)
}

/// Computes the determinant of a 3x3 matrix given by coefficients.
#[inline]
const fn det3x3(
    [a00, a01, a02]: Point3d,
    [a10, a11, a12]: Point3d,
    [a20, a21, a22]: Point3d,
) -> f64 {
    let m01 = a00 * a11 - a10 * a01;
    let m02 = a00 * a21 - a20 * a01;
    let m12 = a10 * a21 - a20 * a11;

    m01 * a22 - m02 * a12 + m12 * a02
}
