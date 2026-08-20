//! Exact predicate implementations using expansion arithmetic.

use crate::expansion::{
    expansion_diff_2, expansion_det2x2, expansion_det3x3, sign_of_expansion_det3x3,
    sign_of_expansion_det4x4,
};
use crate::sign::Sign;

/// orient_2d exact: sign of 2x2 determinant (p1-p0, p2-p0).
#[inline]
pub fn orient_2d_exact(p0: &[f64], p1: &[f64], p2: &[f64]) -> Sign {
    let a11 = expansion_diff_2(p1[0], p0[0]);
    let a12 = expansion_diff_2(p1[1], p0[1]);
    let a21 = expansion_diff_2(p2[0], p0[0]);
    let a22 = expansion_diff_2(p2[1], p0[1]);
    let delta = expansion_det2x2(&a11, &a12, &a21, &a22);
    delta.sign()
}

/// orient_3d exact: sign of 3x3 determinant (rows p1-p0, p2-p0, p3-p0).
#[inline]
pub fn orient_3d_exact(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64]) -> Sign {
    let a11 = expansion_diff_2(p1[0], p0[0]);
    let a12 = expansion_diff_2(p1[1], p0[1]);
    let a13 = expansion_diff_2(p1[2], p0[2]);
    let a21 = expansion_diff_2(p2[0], p0[0]);
    let a22 = expansion_diff_2(p2[1], p0[1]);
    let a23 = expansion_diff_2(p2[2], p0[2]);
    let a31 = expansion_diff_2(p3[0], p0[0]);
    let a32 = expansion_diff_2(p3[1], p0[1]);
    let a33 = expansion_diff_2(p3[2], p0[2]);
    let delta = expansion_det3x3(
        &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33,
    );
    delta.sign()
}

/// dot_3d exact: sign of (p1-p0)·(p2-p0).
#[inline]
pub fn dot_3d_exact(p0: &[f64], p1: &[f64], p2: &[f64]) -> Sign {
    use crate::expansion::{expansion_product, expansion_sum3};
    let u0 = expansion_diff_2(p1[0], p0[0]);
    let u1 = expansion_diff_2(p1[1], p0[1]);
    let u2 = expansion_diff_2(p1[2], p0[2]);
    let v0 = expansion_diff_2(p2[0], p0[0]);
    let v1 = expansion_diff_2(p2[1], p0[1]);
    let v2 = expansion_diff_2(p2[2], p0[2]);
    let u0v0 = expansion_product(&u0, &v0);
    let u1v1 = expansion_product(&u1, &v1);
    let u2v2 = expansion_product(&u2, &v2);
    let dot = expansion_sum3(&u0v0, &u1v1, &u2v2);
    dot.sign()
}

/// det_3d exact: sign of 3x3 determinant of rows p0, p1, p2.
#[inline]
pub fn det_3d_exact(p0: &[f64], p1: &[f64], p2: &[f64]) -> Sign {
    let a11 = crate::expansion::expansion_create(p0[0]);
    let a12 = crate::expansion::expansion_create(p0[1]);
    let a13 = crate::expansion::expansion_create(p0[2]);
    let a21 = crate::expansion::expansion_create(p1[0]);
    let a22 = crate::expansion::expansion_create(p1[1]);
    let a23 = crate::expansion::expansion_create(p1[2]);
    let a31 = crate::expansion::expansion_create(p2[0]);
    let a32 = crate::expansion::expansion_create(p2[1]);
    let a33 = crate::expansion::expansion_create(p2[2]);
    sign_of_expansion_det3x3(
        &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33,
    )
}

/// in_sphere_3d exact: sign of 4x4 in-sphere determinant (lifted). Result is -side4_3d.
#[inline]
pub fn in_sphere_3d_exact(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64], p4: &[f64]) -> Sign {
    use crate::expansion::{
        expansion_create, expansion_det2x2, expansion_product, expansion_sum3, expansion_sum4,
    };
    let a11 = expansion_diff_2(p1[0], p0[0]);
    let a12 = expansion_diff_2(p1[1], p0[1]);
    let a13 = expansion_diff_2(p1[2], p0[2]);
    let p1_0 = p1[0] - p0[0];
    let p1_1 = p1[1] - p0[1];
    let p1_2 = p1[2] - p0[2];
    let a14 = expansion_create(-(p1_0 * p1_0 + p1_1 * p1_1 + p1_2 * p1_2));
    let a21 = expansion_diff_2(p2[0], p0[0]);
    let a22 = expansion_diff_2(p2[1], p0[1]);
    let a23 = expansion_diff_2(p2[2], p0[2]);
    let p2_0 = p2[0] - p0[0];
    let p2_1 = p2[1] - p0[1];
    let p2_2 = p2[2] - p0[2];
    let a24 = expansion_create(-(p2_0 * p2_0 + p2_1 * p2_1 + p2_2 * p2_2));
    let a31 = expansion_diff_2(p3[0], p0[0]);
    let a32 = expansion_diff_2(p3[1], p0[1]);
    let a33 = expansion_diff_2(p3[2], p0[2]);
    let p3_0 = p3[0] - p0[0];
    let p3_1 = p3[1] - p0[1];
    let p3_2 = p3[2] - p0[2];
    let a34 = expansion_create(-(p3_0 * p3_0 + p3_1 * p3_1 + p3_2 * p3_2));
    let a41 = expansion_diff_2(p4[0], p0[0]);
    let a42 = expansion_diff_2(p4[1], p0[1]);
    let a43 = expansion_diff_2(p4[2], p0[2]);
    let p4_0 = p4[0] - p0[0];
    let p4_1 = p4[1] - p0[1];
    let p4_2 = p4[2] - p0[2];
    let a44 = expansion_create(-(p4_0 * p4_0 + p4_1 * p4_1 + p4_2 * p4_2));
    let m12 = expansion_det2x2(&a12, &a13, &a22, &a23);
    let m13 = expansion_det2x2(&a12, &a13, &a32, &a33);
    let m14 = expansion_det2x2(&a12, &a13, &a42, &a43);
    let m23 = expansion_det2x2(&a22, &a23, &a32, &a33);
    let m24 = expansion_det2x2(&a22, &a23, &a42, &a43);
    let m34 = expansion_det2x2(&a32, &a33, &a42, &a43);
    let z11 = expansion_product(&a21, &m34);
    let mut z12 = expansion_product(&a31, &m24);
    z12.negate();
    let z13 = expansion_product(&a41, &m23);
    let delta1 = expansion_sum3(&z11, &z12, &z13);
    let z21 = expansion_product(&a11, &m34);
    let mut z22 = expansion_product(&a31, &m14);
    z22.negate();
    let z23 = expansion_product(&a41, &m13);
    let delta2 = expansion_sum3(&z21, &z22, &z23);
    let z31 = expansion_product(&a11, &m24);
    let mut z32 = expansion_product(&a21, &m14);
    z32.negate();
    let z33 = expansion_product(&a41, &m12);
    let delta3 = expansion_sum3(&z31, &z32, &z33);
    let z41 = expansion_product(&a11, &m23);
    let mut z42 = expansion_product(&a21, &m13);
    z42.negate();
    let z43 = expansion_product(&a31, &m12);
    let delta4 = expansion_sum3(&z41, &z42, &z43);
    let r_1 = expansion_product(&delta1, &a14);
    let mut r_2 = expansion_product(&delta2, &a24);
    r_2.negate();
    let r_3 = expansion_product(&delta3, &a34);
    let mut r_4 = expansion_product(&delta4, &a44);
    r_4.negate();
    let r = expansion_sum4(&r_1, &r_2, &r_3, &r_4);
    // Our 4x4 det r matches in_sphere sign: Negative = inside, Positive = outside.
    r.sign()
}

/// det_4d exact: sign of 4x4 determinant of rows p0, p1, p2, p3.
#[inline]
pub fn det_4d_exact(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64]) -> Sign {
    let a00 = crate::expansion::expansion_create(p0[0]);
    let a01 = crate::expansion::expansion_create(p0[1]);
    let a02 = crate::expansion::expansion_create(p0[2]);
    let a03 = crate::expansion::expansion_create(p0[3]);
    let a10 = crate::expansion::expansion_create(p1[0]);
    let a11 = crate::expansion::expansion_create(p1[1]);
    let a12 = crate::expansion::expansion_create(p1[2]);
    let a13 = crate::expansion::expansion_create(p1[3]);
    let a20 = crate::expansion::expansion_create(p2[0]);
    let a21 = crate::expansion::expansion_create(p2[1]);
    let a22 = crate::expansion::expansion_create(p2[2]);
    let a23 = crate::expansion::expansion_create(p2[3]);
    let a30 = crate::expansion::expansion_create(p3[0]);
    let a31 = crate::expansion::expansion_create(p3[1]);
    let a32 = crate::expansion::expansion_create(p3[2]);
    let a33 = crate::expansion::expansion_create(p3[3]);
    sign_of_expansion_det4x4(
        &a00, &a01, &a02, &a03, &a10, &a11, &a12, &a13, &a20, &a21, &a22, &a23,
        &a30, &a31, &a32, &a33,
    )
}

/// side3_2d exact: used for in_circle_2d. in_circle_2d = -side3_2d(p0,p1,p2,p3, p0,p1,p2).
#[inline]
pub fn side3_2d_exact(
    p0: &[f64],
    p1: &[f64],
    p2: &[f64],
    p3: &[f64],
    q0: &[f64],
    q1: &[f64],
    q2: &[f64],
) -> Sign {
    use crate::expansion::{
        expansion_create, expansion_diff, expansion_product, expansion_sum, expansion_sum3,
    };

    let p1_0_p0_0 = expansion_diff_2(p1[0], p0[0]);
    let p1_1_p0_1 = expansion_diff_2(p1[1], p0[1]);
    let p1_0_sq = expansion_product(&p1_0_p0_0, &p1_0_p0_0);
    let p1_1_sq = expansion_product(&p1_1_p0_1, &p1_1_p0_1);
    let l1 = expansion_sum(&p1_0_sq, &p1_1_sq);

    let p2_0_p0_0 = expansion_diff_2(p2[0], p0[0]);
    let p2_1_p0_1 = expansion_diff_2(p2[1], p0[1]);
    let p2_0_sq = expansion_product(&p2_0_p0_0, &p2_0_p0_0);
    let p2_1_sq = expansion_product(&p2_1_p0_1, &p2_1_p0_1);
    let l2 = expansion_sum(&p2_0_sq, &p2_1_sq);

    let p3_0_p0_0 = expansion_diff_2(p3[0], p0[0]);
    let p3_1_p0_1 = expansion_diff_2(p3[1], p0[1]);
    let p3_0_sq = expansion_product(&p3_0_p0_0, &p3_0_p0_0);
    let p3_1_sq = expansion_product(&p3_1_p0_1, &p3_1_p0_1);
    let l3 = expansion_sum(&p3_0_sq, &p3_1_sq);

    let two = expansion_create(2.0);
    let q0_0_p0_0 = expansion_diff_2(q0[0], p0[0]);
    let q0_1_p0_1 = expansion_diff_2(q0[1], p0[1]);
    let a10_base = expansion_sum(
        &expansion_product(&p1_0_p0_0, &q0_0_p0_0),
        &expansion_product(&p1_1_p0_1, &q0_1_p0_1),
    );
    let a10 = expansion_product(&two, &a10_base);

    let q1_0_p0_0 = expansion_diff_2(q1[0], p0[0]);
    let q1_1_p0_1 = expansion_diff_2(q1[1], p0[1]);
    let a11_base = expansion_sum(
        &expansion_product(&p1_0_p0_0, &q1_0_p0_0),
        &expansion_product(&p1_1_p0_1, &q1_1_p0_1),
    );
    let a11 = expansion_product(&two, &a11_base);

    let q2_0_p0_0 = expansion_diff_2(q2[0], p0[0]);
    let q2_1_p0_1 = expansion_diff_2(q2[1], p0[1]);
    let a12_base = expansion_sum(
        &expansion_product(&p1_0_p0_0, &q2_0_p0_0),
        &expansion_product(&p1_1_p0_1, &q2_1_p0_1),
    );
    let a12 = expansion_product(&two, &a12_base);

    let a20_base = expansion_sum(
        &expansion_product(&p2_0_p0_0, &q0_0_p0_0),
        &expansion_product(&p2_1_p0_1, &q0_1_p0_1),
    );
    let a20 = expansion_product(&two, &a20_base);
    let a21_base = expansion_sum(
        &expansion_product(&p2_0_p0_0, &q1_0_p0_0),
        &expansion_product(&p2_1_p0_1, &q1_1_p0_1),
    );
    let a21 = expansion_product(&two, &a21_base);
    let a22_base = expansion_sum(
        &expansion_product(&p2_0_p0_0, &q2_0_p0_0),
        &expansion_product(&p2_1_p0_1, &q2_1_p0_1),
    );
    let a22 = expansion_product(&two, &a22_base);

    let a30_base = expansion_sum(
        &expansion_product(&p3_0_p0_0, &q0_0_p0_0),
        &expansion_product(&p3_1_p0_1, &q0_1_p0_1),
    );
    let a30 = expansion_product(&two, &a30_base);
    let a31_base = expansion_sum(
        &expansion_product(&p3_0_p0_0, &q1_0_p0_0),
        &expansion_product(&p3_1_p0_1, &q1_1_p0_1),
    );
    let a31 = expansion_product(&two, &a31_base);
    let a32_base = expansion_sum(
        &expansion_product(&p3_0_p0_0, &q2_0_p0_0),
        &expansion_product(&p3_1_p0_1, &q2_1_p0_1),
    );
    let a32 = expansion_product(&two, &a32_base);

    let a11a22 = expansion_product(&a11, &a22);
    let a12a21 = expansion_product(&a12, &a21);
    let b00 = expansion_diff(&a11a22, &a12a21);
    let b01 = expansion_diff(&a21, &a22);
    let b02 = expansion_diff(&a12, &a11);

    let a12a20 = expansion_product(&a12, &a20);
    let a10a22 = expansion_product(&a10, &a22);
    let b10 = expansion_diff(&a12a20, &a10a22);
    let b11 = expansion_diff(&a22, &a20);
    let b12 = expansion_diff(&a10, &a12);

    let a10a21 = expansion_product(&a10, &a21);
    let a11a20 = expansion_product(&a11, &a20);
    let b20 = expansion_diff(&a10a21, &a11a20);
    let b21 = expansion_diff(&a20, &a21);
    let b22 = expansion_diff(&a11, &a10);

    let delta = expansion_sum3(&b00, &b10, &b20);
    let b01l1 = expansion_product(&b01, &l1);
    let b02l2 = expansion_product(&b02, &l2);
    let delta_lambda0 = expansion_sum3(&b01l1, &b02l2, &b00);
    let b11l1 = expansion_product(&b11, &l1);
    let b12l2 = expansion_product(&b12, &l2);
    let delta_lambda1 = expansion_sum3(&b11l1, &b12l2, &b10);
    let b21l1 = expansion_product(&b21, &l1);
    let b22l2 = expansion_product(&b22, &l2);
    let delta_lambda2 = expansion_sum3(&b21l1, &b22l2, &b20);

    let r0 = expansion_product(&delta, &l3);
    let r1 = expansion_product(&a30, &delta_lambda0);
    let r2 = expansion_product(&a31, &delta_lambda1);
    let r3 = expansion_product(&a32, &delta_lambda2);
    let r12 = expansion_sum(&r1, &r2);
    let r123 = expansion_sum(&r12, &r3);
    let r = expansion_diff(&r0, &r123);
    r.sign()
}
