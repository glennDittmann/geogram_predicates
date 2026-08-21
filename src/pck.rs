//! Public API: Geogram PCK predicates (filter + exact fallback).

use crate::exact::{
    aligned_3d_exact, det_3d_exact, det_4d_exact, det_compare_4d_exact, dot_compare_3d_exact,
    orient_2d_exact, orient_3d_exact,
};
use crate::filter::{
    det_3d_filter, det_4d_filter, dot_3d_filter, orient_2d_filter, orient_3d_filter,
    side3_2d_filter, FPG_UNCERTAIN_VALUE,
};
use crate::side::{
    assert_unique_keys, side3_3dlifted, side3_3dlifted_sos, side4_3d, side4_3d_sos, side_exact,
};
use crate::{Sign, SosPoint};

/// Point in 2d (x, y).
pub type Point2d = [f64; 2];
/// Point in 3d (x, y, z).
pub type Point3d = [f64; 3];
/// Point in 4d.
pub type Point4d = [f64; 4];

/// orient_2d: sign of 2d orientation (p1-p0) x (p2-p0). Positive = counter-clockwise.
#[inline]
pub fn orient_2d(p0: &Point2d, p1: &Point2d, p2: &Point2d) -> Sign {
    let r = orient_2d_filter(p0, p1, p2);
    if r != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(r);
    }
    orient_2d_exact(p0, p1, p2)
}

/// orient_3d: sign of 3d orientation (tetrahedron volume). Positive = p3 below plane (p0,p1,p2).
#[inline]
pub fn orient_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d, p3: &Point3d) -> Sign {
    let r = orient_3d_filter(p0, p1, p2, p3);
    if r != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(r);
    }
    orient_3d_exact(p0, p1, p2, p3)
}

/// det_3d: sign of 3x3 determinant with rows p0, p1, p2.
#[inline]
pub fn det_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> Sign {
    let r = det_3d_filter(p0, p1, p2);
    if r != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(r);
    }
    det_3d_exact(p0, p1, p2)
}

/// det_4d: sign of 4x4 determinant with rows p0, p1, p2, p3.
#[inline]
pub fn det_4d(p0: &Point4d, p1: &Point4d, p2: &Point4d, p3: &Point4d) -> Sign {
    let r = det_4d_filter(p0, p1, p2, p3);
    if r != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(r);
    }
    det_4d_exact(p0, p1, p2, p3)
}

/// Exact non-SOS in-circle test. Positive is inside for a positively
/// oriented `(p0,p1,p2)`; an exactly cocircular point returns zero.
#[inline]
pub fn in_circle_2d(p0: &Point2d, p1: &Point2d, p2: &Point2d, p3: &Point2d) -> Sign {
    let filtered = side3_2d_filter(p0, p1, p2, p3, p0, p1, p2);
    if filtered != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(-filtered);
    }
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
    ];
    -side_exact(
        &[&points[0], &points[1], &points[2], &points[3]],
        &[p0, p1, p2],
        false,
    )
}

/// Keyed-SOS in-circle test. Exact cocircularity is resolved by point keys.
#[inline]
pub fn in_circle_2d_sos(
    p0: &SosPoint<2>,
    p1: &SosPoint<2>,
    p2: &SosPoint<2>,
    p3: &SosPoint<2>,
) -> Sign {
    assert_unique_keys(&[p0, p1, p2, p3]);
    let filtered = side3_2d_filter(
        &p0.coords, &p1.coords, &p2.coords, &p3.coords, &p0.coords, &p1.coords, &p2.coords,
    );
    if filtered != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(-filtered);
    }
    -side_exact(
        &[p0, p1, p2, p3],
        &[&p0.coords, &p1.coords, &p2.coords],
        true,
    )
}

/// Exact non-SOS in-circle test for a triangle embedded in 3D.
pub fn in_circle_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d, p3: &Point3d) -> Sign {
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
    ];
    -side_exact(
        &[&points[0], &points[1], &points[2], &points[3]],
        &[p0, p1, p2],
        false,
    )
}

/// Keyed-SOS in-circle test for a triangle embedded in 3D.
pub fn in_circle_3d_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
) -> Sign {
    -side_exact(
        &[p0, p1, p2, p3],
        &[&p0.coords, &p1.coords, &p2.coords],
        true,
    )
}

/// Exact non-SOS lifted in-circle test in 3D.
pub fn in_circle_3dlifted(
    p0: &Point3d,
    p1: &Point3d,
    p2: &Point3d,
    p3: &Point3d,
    heights: [f64; 4],
) -> Sign {
    -side3_3dlifted(p0, p1, p2, p3, heights, p0, p1, p2)
}

/// Keyed-SOS lifted in-circle test in 3D.
pub fn in_circle_3dlifted_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
    heights: [f64; 4],
) -> Sign {
    -side3_3dlifted_sos(p0, p1, p2, p3, heights, &p0.coords, &p1.coords, &p2.coords)
}

/// Exact non-SOS in-sphere test. Geogram requires the base tetrahedron to
/// have positive orientation. Positive means inside and zero means cospherical.
pub fn in_sphere_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d, p3: &Point3d, p4: &Point3d) -> Sign {
    -side4_3d(p0, p1, p2, p3, p4)
}

/// Keyed-SOS in-sphere test. Positive means inside for a positively oriented
/// base tetrahedron.
pub fn in_sphere_3d_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
    p4: &SosPoint<3>,
) -> Sign {
    -side4_3d_sos(p0, p1, p2, p3, p4)
}

/// dot_3d: sign of dot product (p1-p0)·(p2-p0). Positive = acute angle.
#[inline]
pub fn dot_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> Sign {
    let r = dot_3d_filter(p0, p1, p2);
    if r != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(r);
    }
    crate::exact::dot_3d_exact(p0, p1, p2)
}

/// aligned_3d: true if (p1-p0) and (p2-p0) are collinear (cross product zero).
#[inline]
pub fn aligned_3d(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> bool {
    aligned_3d_exact(p0, p1, p2)
}

/// Compare `det(p0,p1,p2,p4)` with `det(p0,p1,p2,p3)` exactly.
pub fn det_compare_4d(
    p0: &Point4d,
    p1: &Point4d,
    p2: &Point4d,
    p3: &Point4d,
    p4: &Point4d,
) -> Sign {
    det_compare_4d_exact(p0, p1, p2, p3, p4)
}

/// Exact sign of `v0.v1 - v0.v2`.
pub fn dot_compare_3d(v0: &Point3d, v1: &Point3d, v2: &Point3d) -> Sign {
    dot_compare_3d_exact(v0, v1, v2)
}

/// points_are_identical_2d: exact equality of two 2D points.
#[inline]
pub fn points_are_identical_2d(p1: &Point2d, p2: &Point2d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1]
}

/// points_are_identical_3d: exact equality of two 3D points.
#[inline]
pub fn points_are_identical_3d(p1: &Point3d, p2: &Point3d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]
}

/// points_are_colinear_3d: true if the three points lie on a line (cross product (p2-p1)×(p3-p1) is zero).
#[inline]
pub fn points_are_colinear_3d(p1: &Point3d, p2: &Point3d, p3: &Point3d) -> bool {
    aligned_3d(p1, p2, p3)
}

/// orient_3d_inexact: sign of 3D orientation using floating-point only (no exact fallback).
#[inline]
pub fn orient_3d_inexact(p0: &Point3d, p1: &Point3d, p2: &Point3d, p3: &Point3d) -> Sign {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a13 = p1[2] - p0[2];
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let a23 = p2[2] - p0[2];
    let a31 = p3[0] - p0[0];
    let a32 = p3[1] - p0[1];
    let a33 = p3[2] - p0[2];
    let det = (a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))
        + (a31 * ((a12 * a23) - (a13 * a22)));
    crate::sign::geo_sgn(det)
}
