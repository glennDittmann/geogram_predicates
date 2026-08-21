//! Exact predicate implementations using expansion arithmetic.

use crate::expansion::{
    expansion_det2x2, expansion_det3x3, expansion_diff_2, sign_of_expansion_det3x3,
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
    let delta = expansion_det3x3(&a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33);
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

/// Exact sign of `v0.v1 - v0.v2`.
pub fn dot_compare_3d_exact(v0: &[f64], v1: &[f64], v2: &[f64]) -> Sign {
    use crate::expansion::{expansion_create, expansion_diff, expansion_product, expansion_sum3};
    let v0 = [
        expansion_create(v0[0]),
        expansion_create(v0[1]),
        expansion_create(v0[2]),
    ];
    let v1 = [
        expansion_create(v1[0]),
        expansion_create(v1[1]),
        expansion_create(v1[2]),
    ];
    let v2 = [
        expansion_create(v2[0]),
        expansion_create(v2[1]),
        expansion_create(v2[2]),
    ];
    let d01 = expansion_sum3(
        &expansion_product(&v0[0], &v1[0]),
        &expansion_product(&v0[1], &v1[1]),
        &expansion_product(&v0[2], &v1[2]),
    );
    let d02 = expansion_sum3(
        &expansion_product(&v0[0], &v2[0]),
        &expansion_product(&v0[1], &v2[1]),
        &expansion_product(&v0[2], &v2[2]),
    );
    expansion_diff(&d01, &d02).sign()
}

/// Exact collinearity test for three 3D points.
pub fn aligned_3d_exact(p0: &[f64], p1: &[f64], p2: &[f64]) -> bool {
    let u = [
        expansion_diff_2(p1[0], p0[0]),
        expansion_diff_2(p1[1], p0[1]),
        expansion_diff_2(p1[2], p0[2]),
    ];
    let v = [
        expansion_diff_2(p2[0], p0[0]),
        expansion_diff_2(p2[1], p0[1]),
        expansion_diff_2(p2[2], p0[2]),
    ];
    expansion_det2x2(&u[1], &u[2], &v[1], &v[2]).sign() == Sign::Zero
        && expansion_det2x2(&u[2], &u[0], &v[2], &v[0]).sign() == Sign::Zero
        && expansion_det2x2(&u[0], &u[1], &v[0], &v[1]).sign() == Sign::Zero
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
    sign_of_expansion_det3x3(&a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33)
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
        &a00, &a01, &a02, &a03, &a10, &a11, &a12, &a13, &a20, &a21, &a22, &a23, &a30, &a31, &a32,
        &a33,
    )
}

/// Exact sign of the 4D determinant with final row `p4 - p3`.
pub fn det_compare_4d_exact(p0: &[f64], p1: &[f64], p2: &[f64], p3: &[f64], p4: &[f64]) -> Sign {
    use crate::expansion::{expansion_create, expansion_determinant};
    let matrix = vec![
        p0.iter().map(|&x| expansion_create(x)).collect(),
        p1.iter().map(|&x| expansion_create(x)).collect(),
        p2.iter().map(|&x| expansion_create(x)).collect(),
        (0..4).map(|i| expansion_diff_2(p4[i], p3[i])).collect(),
    ];
    expansion_determinant(&matrix).sign()
}

#[cfg(test)]
mod tests {
    #![allow(clippy::needless_range_loop)]
    use super::*;
    use crate::filter::{det_3d_filter, det_4d_filter, orient_2d_filter, orient_3d_filter};

    fn sign(value: i128) -> Sign {
        if value < 0 {
            Sign::Negative
        } else if value > 0 {
            Sign::Positive
        } else {
            Sign::Zero
        }
    }

    fn det3(m: [[i128; 3]; 3]) -> i128 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    fn det4(m: [[i128; 4]; 4]) -> i128 {
        let mut result = 0;
        for column in 0..4 {
            let mut minor = [[0; 3]; 3];
            for row in 1..4 {
                let mut target = 0;
                for source in 0..4 {
                    if source != column {
                        minor[row - 1][target] = m[row][source];
                        target += 1;
                    }
                }
            }
            let term = m[0][column] * det3(minor);
            result += if column % 2 == 0 { term } else { -term };
        }
        result
    }

    fn next(seed: &mut u64) -> i128 {
        *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        ((*seed >> 32) % 33) as i128 - 16
    }

    #[test]
    fn exact_determinants_match_integer_arithmetic() {
        let mut seed = 0x5eed_u64;
        for _ in 0..500 {
            let mut m3 = [[0_i128; 3]; 3];
            for row in &mut m3 {
                for value in row {
                    *value = next(&mut seed);
                }
            }
            let f3 = m3.map(|row| row.map(|value| value as f64));
            let exact3 = det_3d_exact(&f3[0], &f3[1], &f3[2]);
            assert_eq!(exact3, sign(det3(m3)));
            let filtered3 = det_3d_filter(&f3[0], &f3[1], &f3[2]);
            if filtered3 != 0 {
                assert_eq!(Sign::from_i8(filtered3), exact3);
            }

            let mut m4 = [[0_i128; 4]; 4];
            for row in &mut m4 {
                for value in row {
                    *value = next(&mut seed);
                }
            }
            let f4 = m4.map(|row| row.map(|value| value as f64));
            let exact4 = det_4d_exact(&f4[0], &f4[1], &f4[2], &f4[3]);
            assert_eq!(exact4, sign(det4(m4)));
            let filtered4 = det_4d_filter(&f4[0], &f4[1], &f4[2], &f4[3]);
            if filtered4 != 0 {
                assert_eq!(Sign::from_i8(filtered4), exact4);
            }
        }
    }

    #[test]
    fn orientation_filters_agree_with_exact_fallbacks() {
        let mut seed = 0xc0ffee_u64;
        for _ in 0..500 {
            let points2 = core::array::from_fn::<_, 3, _>(|_| {
                [next(&mut seed) as f64, next(&mut seed) as f64]
            });
            let exact2 = orient_2d_exact(&points2[0], &points2[1], &points2[2]);
            let filtered2 = orient_2d_filter(&points2[0], &points2[1], &points2[2]);
            if filtered2 != 0 {
                assert_eq!(Sign::from_i8(filtered2), exact2);
            }

            let points3 = core::array::from_fn::<_, 4, _>(|_| {
                [
                    next(&mut seed) as f64,
                    next(&mut seed) as f64,
                    next(&mut seed) as f64,
                ]
            });
            let exact3 = orient_3d_exact(&points3[0], &points3[1], &points3[2], &points3[3]);
            let filtered3 = orient_3d_filter(&points3[0], &points3[1], &points3[2], &points3[3]);
            if filtered3 != 0 {
                assert_eq!(Sign::from_i8(filtered3), exact3);
            }
        }
    }
}
