//! Floating-point filter functions. Return sign or 0 for "uncertain" (need exact path).
#![allow(clippy::excessive_precision)] // Constants are copied verbatim from Geogram's generated filters.

/// Value returned when the filter cannot determine the sign (exact arithmetic required).
pub const FPG_UNCERTAIN_VALUE: i8 = 0;

#[inline]
fn fabs(x: f64) -> f64 {
    x.abs()
}

/// orient_2d filter: sign of 2x2 determinant (p1-p0, p2-p0).
#[inline]
pub fn orient_2d_filter(p0: &[f64], p1: &[f64], p2: &[f64]) -> i8 {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let delta = (a11 * a22) - (a12 * a21);
    let mut max1 = fabs(a11);
    if max1 < fabs(a12) {
        max1 = fabs(a12);
    }
    let mut max2 = fabs(a21);
    if max2 < fabs(a22) {
        max2 = fabs(a22);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if lower_bound_1 < 5.00368081960964635413e-147 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.67597599124282407923e+153 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 8.88720573725927976811e-16 * (max1 * max2);
    if delta > eps {
        1
    } else if delta < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// orient_3d filter: sign of 3x3 determinant (rows p1-p0, p2-p0, p3-p0).
#[inline]
pub fn orient_3d_filter(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64]) -> i8 {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a13 = p1[2] - p0[2];
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let a23 = p2[2] - p0[2];
    let a31 = p3[0] - p0[0];
    let a32 = p3[1] - p0[1];
    let a33 = p3[2] - p0[2];
    let delta = (a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))
        + (a31 * ((a12 * a23) - (a13 * a22)));
    let mut max1 = fabs(a11);
    if max1 < fabs(a21) {
        max1 = fabs(a21);
    }
    if max1 < fabs(a31) {
        max1 = fabs(a31);
    }
    let mut max2 = fabs(a12);
    if max2 < fabs(a13) {
        max2 = fabs(a13);
    }
    if max2 < fabs(a22) {
        max2 = fabs(a22);
    }
    if max2 < fabs(a23) {
        max2 = fabs(a23);
    }
    let mut max3 = fabs(a32);
    if max3 < fabs(a33) {
        max3 = fabs(a33);
    }
    if max3 < fabs(a31) {
        max3 = fabs(a31);
    }
    if max3 < fabs(a21) {
        max3 = fabs(a21);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }
    if lower_bound_1 < 6.00259835109572228086e-98 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.67597599124282407923e+153 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 1.77635683940025046481e-15 * ((max2 * max3) * max1);
    if delta > eps {
        1
    } else if delta < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// det_3d filter: sign of 3x3 determinant of rows p0, p1, p2.
#[inline]
pub fn det_3d_filter(p0: &[f64], p1: &[f64], p2: &[f64]) -> i8 {
    let delta = (p0[0] * ((p1[1] * p2[2]) - (p1[2] * p2[1])))
        - (p1[0] * ((p0[1] * p2[2]) - (p0[2] * p2[1])))
        + (p2[0] * ((p0[1] * p1[2]) - (p0[2] * p1[1])));
    let mut max1 = fabs(p0[0]);
    if max1 < fabs(p1[0]) {
        max1 = fabs(p1[0]);
    }
    if max1 < fabs(p2[0]) {
        max1 = fabs(p2[0]);
    }
    let mut max2 = fabs(p0[1]);
    if max2 < fabs(p0[2]) {
        max2 = fabs(p0[2]);
    }
    if max2 < fabs(p1[1]) {
        max2 = fabs(p1[1]);
    }
    if max2 < fabs(p1[2]) {
        max2 = fabs(p1[2]);
    }
    let mut max3 = fabs(p1[1]);
    if max3 < fabs(p1[2]) {
        max3 = fabs(p1[2]);
    }
    if max3 < fabs(p2[1]) {
        max3 = fabs(p2[1]);
    }
    if max3 < fabs(p2[2]) {
        max3 = fabs(p2[2]);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }
    if lower_bound_1 < 1.92663387981871579179e-98 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.11987237108890185662e+102 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 3.11133555671680765034e-15 * ((max2 * max3) * max1);
    if delta > eps {
        1
    } else if delta < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// det_4d filter: sign of 4x4 determinant of rows p0, p1, p2, p3.
#[inline]
pub fn det_4d_filter(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64]) -> i8 {
    let m12 = (p1[0] * p0[1]) - (p0[0] * p1[1]);
    let m13 = (p2[0] * p0[1]) - (p0[0] * p2[1]);
    let m14 = (p3[0] * p0[1]) - (p0[0] * p3[1]);
    let m23 = (p2[0] * p1[1]) - (p1[0] * p2[1]);
    let m24 = (p3[0] * p1[1]) - (p1[0] * p3[1]);
    let m34 = (p3[0] * p2[1]) - (p2[0] * p3[1]);
    let m123 = (m23 * p0[2]) - (m13 * p1[2]) + (m12 * p2[2]);
    let m124 = (m24 * p0[2]) - (m14 * p1[2]) + (m12 * p3[2]);
    let m134 = (m34 * p0[2]) - (m14 * p2[2]) + (m13 * p3[2]);
    let m234 = (m34 * p1[2]) - (m24 * p2[2]) + (m23 * p3[2]);
    let r = (m234 * p0[3]) - (m134 * p1[3]) + (m124 * p2[3]) - (m123 * p3[3]);
    let mut max1 = fabs(p0[0]);
    if max1 < fabs(p1[0]) {
        max1 = fabs(p1[0]);
    }
    if max1 < fabs(p2[0]) {
        max1 = fabs(p2[0]);
    }
    if max1 < fabs(p3[0]) {
        max1 = fabs(p3[0]);
    }
    let mut max2 = fabs(p0[1]);
    if max2 < fabs(p0[2]) {
        max2 = fabs(p0[2]);
    }
    if max2 < fabs(p1[1]) {
        max2 = fabs(p1[1]);
    }
    if max2 < fabs(p1[2]) {
        max2 = fabs(p1[2]);
    }
    if max2 < fabs(p2[1]) {
        max2 = fabs(p2[1]);
    }
    if max2 < fabs(p2[2]) {
        max2 = fabs(p2[2]);
    }
    if max2 < fabs(p3[1]) {
        max2 = fabs(p3[1]);
    }
    if max2 < fabs(p3[2]) {
        max2 = fabs(p3[2]);
    }
    let mut max3 = fabs(m12);
    if max3 < fabs(m13) {
        max3 = fabs(m13);
    }
    if max3 < fabs(m14) {
        max3 = fabs(m14);
    }
    if max3 < fabs(m23) {
        max3 = fabs(m23);
    }
    if max3 < fabs(m24) {
        max3 = fabs(m24);
    }
    if max3 < fabs(m34) {
        max3 = fabs(m34);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }
    if lower_bound_1 < 2.90064983848602906132e-74 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.67597599124282407923e+153 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 5.22161699646551609832e-14 * ((max2 * max3) * max1);
    if r > eps {
        1
    } else if r < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// dot_3d filter: sign of dot product (p1-p0)*(p2-p0).
#[inline]
pub fn dot_3d_filter(p0: &[f64], p1: &[f64], p2: &[f64]) -> i8 {
    let ux = p1[0] - p0[0];
    let uy = p1[1] - p0[1];
    let uz = p1[2] - p0[2];
    let vx = p2[0] - p0[0];
    let vy = p2[1] - p0[1];
    let vz = p2[2] - p0[2];
    let dot = (ux * vx) + (uy * vy) + (uz * vz);
    let mut max1 = fabs(ux);
    if max1 < fabs(uy) {
        max1 = fabs(uy);
    }
    if max1 < fabs(uz) {
        max1 = fabs(uz);
    }
    let mut max2 = fabs(vx);
    if max2 < fabs(vy) {
        max2 = fabs(vy);
    }
    if max2 < fabs(vz) {
        max2 = fabs(vz);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if lower_bound_1 < 2.00491534964401497902e-98 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.11987237108890185662e+102 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 2.44649357902858013632e-15 * (max1 * max2);
    if dot > eps {
        1
    } else if dot < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// side4_3d filter: orient_4d in 3d lifted (4x4 matrix with last column -||pi-p0||^2).
/// Returns sign; 0 = uncertain.
#[inline]
pub fn side4_3d_filter(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64], p4: &[f64]) -> i8 {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a13 = p1[2] - p0[2];
    let p1_0_p0_0 = p1[0] - p0[0];
    let p1_1_p0_1 = p1[1] - p0[1];
    let p1_2_p0_2 = p1[2] - p0[2];
    let a14 = -((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1) + (p1_2_p0_2 * p1_2_p0_2));
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let a23 = p2[2] - p0[2];
    let p2_0_p0_0 = p2[0] - p0[0];
    let p2_1_p0_1 = p2[1] - p0[1];
    let p2_2_p0_2 = p2[2] - p0[2];
    let a24 = -((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1) + (p2_2_p0_2 * p2_2_p0_2));
    let a31 = p3[0] - p0[0];
    let a32 = p3[1] - p0[1];
    let a33 = p3[2] - p0[2];
    let p3_0_p0_0 = p3[0] - p0[0];
    let p3_1_p0_1 = p3[1] - p0[1];
    let p3_2_p0_2 = p3[2] - p0[2];
    let a34 = -((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1) + (p3_2_p0_2 * p3_2_p0_2));
    let a41 = p4[0] - p0[0];
    let a42 = p4[1] - p0[1];
    let a43 = p4[2] - p0[2];
    let p4_0_p0_0 = p4[0] - p0[0];
    let p4_1_p0_1 = p4[1] - p0[1];
    let p4_2_p0_2 = p4[2] - p0[2];
    let a44 = -((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1) + (p4_2_p0_2 * p4_2_p0_2));
    let delta1 = (a21 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a22 * a43) - (a23 * a42)))
        + (a41 * ((a22 * a33) - (a23 * a32)));
    let delta2 = (a11 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a12 * a43) - (a13 * a42)))
        + (a41 * ((a12 * a33) - (a13 * a32)));
    let delta3 = (a11 * ((a22 * a43) - (a23 * a42))) - (a21 * ((a12 * a43) - (a13 * a42)))
        + (a41 * ((a12 * a23) - (a13 * a22)));
    let delta4 = (a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))
        + (a31 * ((a12 * a23) - (a13 * a22)));
    let r = (((delta1 * a14) - (delta2 * a24)) + (delta3 * a34)) - (delta4 * a44);
    let mut max1 = fabs(a11);
    if max1 < fabs(a21) {
        max1 = fabs(a21);
    }
    if max1 < fabs(a31) {
        max1 = fabs(a31);
    }
    let mut max2 = fabs(a12);
    if max2 < fabs(a13) {
        max2 = fabs(a13);
    }
    if max2 < fabs(a22) {
        max2 = fabs(a22);
    }
    if max2 < fabs(a23) {
        max2 = fabs(a23);
    }
    let mut max3 = fabs(a22);
    if max3 < fabs(a23) {
        max3 = fabs(a23);
    }
    if max3 < fabs(a32) {
        max3 = fabs(a32);
    }
    if max3 < fabs(a33) {
        max3 = fabs(a33);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if lower_bound_1 < 2.40065289904160175849e-74 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 1.67597599124282407923e+153 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 4.76837158203125448441e-14 * ((max2 * max3) * max1);
    if r > eps {
        1
    } else if r < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}

/// side3_2d filter (2d in-circle / side predicate). in_circle_2d = -side3_2d(p0,p1,p2,p3, p0,p1,p2).
#[inline]
pub fn side3_2d_filter(
    p0: &[f64],
    p1: &[f64],
    p2: &[f64],
    p3: &[f64],
    q0: &[f64],
    q1: &[f64],
    q2: &[f64],
) -> i8 {
    let p1_0_p0_0 = p1[0] - p0[0];
    let p1_1_p0_1 = p1[1] - p0[1];
    let l1 = (p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1);
    let p2_0_p0_0 = p2[0] - p0[0];
    let p2_1_p0_1 = p2[1] - p0[1];
    let l2 = (p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1);
    let p3_0_p0_0 = p3[0] - p0[0];
    let p3_1_p0_1 = p3[1] - p0[1];
    let l3 = (p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1);
    let q0_0_p0_0 = q0[0] - p0[0];
    let q0_1_p0_1 = q0[1] - p0[1];
    let a10 = 2.0 * ((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1));
    let q1_0_p0_0 = q1[0] - p0[0];
    let q1_1_p0_1 = q1[1] - p0[1];
    let a11 = 2.0 * ((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1));
    let q2_0_p0_0 = q2[0] - p0[0];
    let q2_1_p0_1 = q2[1] - p0[1];
    let a12 = 2.0 * ((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1));
    let a20 = 2.0 * ((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1));
    let a21 = 2.0 * ((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1));
    let a22 = 2.0 * ((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1));
    let a30 = 2.0 * ((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1));
    let a31 = 2.0 * ((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1));
    let a32 = 2.0 * ((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1));
    let b00 = (a11 * a22) - (a12 * a21);
    let b01 = a21 - a22;
    let b02 = a12 - a11;
    let b10 = (a12 * a20) - (a10 * a22);
    let b11 = a22 - a20;
    let b12 = a10 - a12;
    let b20 = (a10 * a21) - (a11 * a20);
    let b21 = a20 - a21;
    let b22 = a11 - a10;
    let delta = (b00 + b10) + b20;
    let delta_lambda0 = ((b01 * l1) + (b02 * l2)) + b00;
    let delta_lambda1 = ((b11 * l1) + (b12 * l2)) + b10;
    let delta_lambda2 = ((b21 * l1) + (b22 * l2)) + b20;
    let r = (delta * l3) - ((a30 * delta_lambda0) + (a31 * delta_lambda1) + (a32 * delta_lambda2));
    let mut max1 = fabs(p2_0_p0_0);
    if max1 < fabs(p2_1_p0_1) {
        max1 = fabs(p2_1_p0_1);
    }
    let mut max2 = fabs(q0_0_p0_0);
    if max2 < fabs(q0_1_p0_1) {
        max2 = fabs(q0_1_p0_1);
    }
    if max2 < fabs(q1_0_p0_0) {
        max2 = fabs(q1_0_p0_0);
    }
    if max2 < fabs(q1_1_p0_1) {
        max2 = fabs(q1_1_p0_1);
    }
    if max2 < fabs(q2_0_p0_0) {
        max2 = fabs(q2_0_p0_0);
    }
    if max2 < fabs(q2_1_p0_1) {
        max2 = fabs(q2_1_p0_1);
    }
    let mut max3 = fabs(p1_0_p0_0);
    if max3 < fabs(p1_1_p0_1) {
        max3 = fabs(p1_1_p0_1);
    }
    let (mut lower_bound_1, mut upper_bound_1) = (max2, max2);
    if max1 < lower_bound_1 {
        lower_bound_1 = max1;
    } else if max1 > upper_bound_1 {
        upper_bound_1 = max1;
    }
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }
    if lower_bound_1 < 2.79532528033945620759e-74 {
        return FPG_UNCERTAIN_VALUE;
    }
    if upper_bound_1 > 2.59614842926741294957e+33 {
        return FPG_UNCERTAIN_VALUE;
    }
    let eps = 3.64430756537603111258e-14 * (((max3 * max2) * max1) * max2);
    if r > eps {
        1
    } else if r < -eps {
        -1
    } else {
        FPG_UNCERTAIN_VALUE
    }
}
