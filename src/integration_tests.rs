//! Integration tests for the full API (included in module tree).

use crate::{self as gp, Sign};

#[test]
fn test_geo_sgn() {
    assert_eq!(gp::geo_sgn(42.0), Sign::Positive);
    assert_eq!(gp::geo_sgn(-3.14), Sign::Negative);
    assert_eq!(gp::geo_sgn(0.0), Sign::Zero);
}

#[test]
fn test_orient_2d() {
    gp::initialize();
    let a = [0.0, 0.0];
    let b = [1.0, 0.0];
    let c = [0.0, 1.0];
    assert_eq!(gp::orient_2d(&a, &b, &c), Sign::Positive);
    assert_eq!(gp::orient_2d(&a, &c, &b), Sign::Negative);
    let d = [2.0, 0.0];
    assert_eq!(gp::orient_2d(&a, &b, &d), Sign::Zero);
}

#[test]
fn test_orient_3d() {
    gp::initialize();
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    assert_eq!(gp::orient_3d(&a, &b, &c, &d), Sign::Positive);
    assert_eq!(gp::orient_3d(&a, &c, &b, &d), Sign::Negative);
    let e = [0.5, 0.5, 0.0];
    assert_eq!(gp::orient_3d(&a, &b, &c, &e), Sign::Zero);
}

#[test]
fn test_dot_3d() {
    gp::initialize();
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    assert_eq!(gp::dot_3d(&a, &b, &c), Sign::Zero);
    let d = [2.0, 0.0, 0.0];
    assert_eq!(gp::dot_3d(&a, &b, &d), Sign::Positive);
}

#[test]
fn test_points_identical() {
    let p1 = [4.0, 2.0];
    let p2 = [4.0, 2.0];
    let p3 = [4.0, 2.1];
    assert!(gp::points_are_identical_2d(&p1, &p2));
    assert!(!gp::points_are_identical_2d(&p1, &p3));

    let q1 = [1.0, 2.0, 3.0];
    let q2 = [1.0, 2.0, 3.0];
    let q3 = [1.0, 2.0, 3.1];
    assert!(gp::points_are_identical_3d(&q1, &q2));
    assert!(!gp::points_are_identical_3d(&q1, &q3));
}

#[test]
fn test_points_colinear_3d() {
    gp::initialize();
    let p1 = [0.0, 0.0, 0.0];
    let p2 = [1.0, 1.0, 1.0];
    let p3 = [2.0, 2.0, 2.0];
    assert!(gp::points_are_colinear_3d(&p1, &p2, &p3));
    let p4 = [1.0, 0.0, 0.0];
    assert!(!gp::points_are_colinear_3d(&p1, &p2, &p4));
}

#[test]
fn test_in_circle_2d() {
    gp::initialize();
    let a = [0.0, 0.0];
    let b = [1.0, 0.0];
    let c = [0.0, 1.0];
    let p_in = [0.1, 0.1];
    let p_out = [2.0, 2.0];
    assert_eq!(gp::in_circle_2d_sos(&a, &b, &c, &p_in), Sign::Positive);
    assert_eq!(gp::in_circle_2d_sos(&a, &b, &c, &p_out), Sign::Negative);
}

#[test]
fn test_in_sphere_3d() {
    gp::initialize();
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    let p_in = [0.1, 0.1, 0.1];
    let p_out = [2.0, 2.0, 2.0];
    assert_eq!(gp::in_sphere_3d_sos(&a, &b, &c, &d, &p_in), Sign::Negative);
    assert_eq!(gp::in_sphere_3d_sos(&a, &b, &c, &d, &p_out), Sign::Positive);
}

#[test]
fn test_det_4d() {
    gp::initialize();
    let a4 = [1.0, 2.0, 3.0, 4.0];
    let b4 = [5.0, 6.0, 7.0, 8.0];
    let c4 = [9.0, 10.0, 11.0, 12.0];
    let d4 = [13.0, 14.0, 15.0, 16.0];
    assert_eq!(gp::det_4d(&a4, &b4, &c4, &d4), Sign::Zero);
}

#[test]
fn test_det_3d() {
    gp::initialize();
    let a3 = [1.0, 2.0, 3.0];
    let b3 = [4.0, 5.0, 6.0];
    let c3 = [7.0, 8.0, 9.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), Sign::Zero);

    let a3 = [1.0, 0.0, 0.0];
    let b3 = [0.0, 1.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), Sign::Positive);

    let a3 = [1.0, 0.0, 0.0];
    let b3 = [0.0, -1.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), Sign::Negative);
}

#[test]
fn test_orient_2dlifted_and_3dlifted() {
    gp::initialize();
    let a2 = [0.0, 0.0];
    let b2 = [1.0, 0.0];
    let c2 = [0.0, 1.0];
    let p2 = [0.1, 0.1];
    let h = 0.0;
    let res = gp::orient_2dlifted_sos(&a2, &b2, &c2, &p2, [h, h, h, h]);
    assert_eq!(res, gp::in_circle_2d_sos(&a2, &b2, &c2, &p2));

    let a3 = [0.0, 0.0, 0.0];
    let b3 = [1.0, 0.0, 0.0];
    let c3 = [0.0, 1.0, 0.0];
    let d3 = [0.0, 0.0, 1.0];
    let p3 = [0.1, 0.1, 0.1];
    let h3 = 0.0;
    let res3 = gp::orient_3dlifted_sos(&a3, &b3, &c3, &d3, &p3, [h3, h3, h3, h3, h3]);
    assert_eq!(res3, gp::in_sphere_3d_sos(&a3, &b3, &c3, &d3, &p3));
}

#[test]
fn test_orient_3d_inexact() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    assert_eq!(gp::orient_3d_inexact(&a, &b, &c, &d), Sign::Positive);
}
