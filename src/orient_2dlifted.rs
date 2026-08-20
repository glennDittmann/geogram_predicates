//! 2D lifted orientation (regularity test for 2D).

use crate::expansion::{expansion_det2x2, expansion_diff_2, expansion_product, expansion_sum, expansion_sum3};
use crate::pck::Point2d;
use crate::sign::Sign;
use core::cmp::Ordering;

/// Computes the 3D orientation test with lifted points (regularity test for 2D).
/// Positive = p' below the plane, Negative = above, Zero = on.
#[inline]
pub fn orient_2dlifted_sos(
    a: &Point2d,
    b: &Point2d,
    c: &Point2d,
    p: &Point2d,
    [h_a, h_b, h_c, h_p]: [f64; 4],
) -> Sign {
    let mut result = side3_2dlifted_2d_filter(a, b, c, p, [h_a, h_b, h_c, h_p]);
    if result == Sign::Zero {
        result = side3h_2d_exact_sos(a, b, c, p, [h_a, h_b, h_c, h_p], None);
    }
    result
}

#[inline]
fn side3_2dlifted_2d_filter(
    p0: &Point2d,
    p1: &Point2d,
    p2: &Point2d,
    p3: &Point2d,
    [h0, h1, h2, h3]: [f64; 4],
) -> Sign {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a13 = h0 - h1;
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let a23 = h0 - h2;
    let a31 = p3[0] - p0[0];
    let a32 = p3[1] - p0[1];
    let a33 = h0 - h3;

    let delta1 = (a21 * a32) - (a22 * a31);
    let delta2 = (a11 * a32) - (a12 * a31);
    let delta3 = (a11 * a22) - (a12 * a21);
    let r = ((delta1 * a13) - (delta2 * a23)) + (delta3 * a33);

    let max1 = (a11.abs()).max(a12.abs());
    let max2 = (a21.abs()).max(a22.abs());
    let (mut lower_bound_1, mut upper_bound_1) = (max1, max1);
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if lower_bound_1 < 5.00368081960964635413e-147 {
        return Sign::Zero;
    }
    if upper_bound_1 > 5.59936185544450928309e+101 {
        return Sign::Zero;
    }
    let eps = 8.88720573725927976811e-16 * (max1 * max2);
    let int_tmp_result = if delta3 > eps {
        Sign::Positive
    } else if delta3 < -eps {
        Sign::Negative
    } else {
        return Sign::Zero;
    };
    let delta3_sign = int_tmp_result;

    let max3 = max1.max(max2);
    let max4 = (a13.abs()).max(a23.abs()).max(a33.abs());
    let max5 = max2.max(a31.abs()).max(a32.abs());
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if max5 < lower_bound_1 {
        lower_bound_1 = max5;
    } else if max5 > upper_bound_1 {
        upper_bound_1 = max5;
    }
    if max4 < lower_bound_1 {
        lower_bound_1 = max4;
    } else if max4 > upper_bound_1 {
        upper_bound_1 = max4;
    }
    if lower_bound_1 < 1.63288018496748314939e-98 {
        return Sign::Zero;
    }
    if upper_bound_1 > 5.59936185544450928309e+101 {
        return Sign::Zero;
    }
    let eps2 = 5.11071278299732992696e-15 * ((max3 * max5) * max4);
    let int_tmp_result_ffwkcaa = if r > eps2 {
        Sign::Positive
    } else if r < -eps2 {
        Sign::Negative
    } else {
        return Sign::Zero;
    };
    delta3_sign * int_tmp_result_ffwkcaa
}

fn side3h_2d_exact_sos(
    p0: &Point2d,
    p1: &Point2d,
    p2: &Point2d,
    p3: &Point2d,
    [h0, h1, h2, h3]: [f64; 4],
    sos: Option<bool>,
) -> Sign {
    let sos = sos.unwrap_or(true);

    let a11 = expansion_diff_2(p1[0], p0[0]);
    let a12 = expansion_diff_2(p1[1], p0[1]);
    let a13 = expansion_diff_2(h0, h1);
    let a21 = expansion_diff_2(p2[0], p0[0]);
    let a22 = expansion_diff_2(p2[1], p0[1]);
    let a23 = expansion_diff_2(h0, h2);
    let a31 = expansion_diff_2(p3[0], p0[0]);
    let a32 = expansion_diff_2(p3[1], p0[1]);
    let a33 = expansion_diff_2(h0, h3);

    let delta1 = expansion_det2x2(&a21, &a22, &a31, &a32);
    let delta2 = expansion_det2x2(&a11, &a12, &a31, &a32);
    let delta3 = expansion_det2x2(&a11, &a12, &a21, &a22);

    let delta3_sign = delta3.sign();
    if delta3_sign == Sign::Zero {
        return Sign::Zero;
    }

    let r_1 = expansion_product(&delta1, &a13);
    let mut r_2 = expansion_product(&delta2, &a23);
    r_2.negate();
    let r_3 = expansion_product(&delta3, &a33);
    let r = expansion_sum3(&r_1, &r_2, &r_3);
    let r_sign = r.sign();

    if sos && r_sign == Sign::Zero {
        let mut p_sort = [p0, p1, p2, p3];
        p_sort.sort_by(lexico_compare_2d);

        for i in 0..3 {
            if std::ptr::eq(p_sort[i], p0) {
                let z1 = crate::expansion::expansion_diff(&delta2, &delta1);
                let z = expansion_sum(&z1, &delta3);
                let z_sign = z.sign();
                if z_sign != Sign::Zero {
                    return delta3_sign * z_sign;
                }
            } else if std::ptr::eq(p_sort[i], p1) {
                let delta1_sign = delta1.sign();
                if delta1_sign != Sign::Zero {
                    return delta3_sign * delta1_sign;
                }
            } else if std::ptr::eq(p_sort[i], p2) {
                let delta2_sign = delta2.sign();
                if delta2_sign != Sign::Zero {
                    return (-delta3_sign) * delta2_sign;
                }
            } else if std::ptr::eq(p_sort[i], p3) {
                return Sign::Negative;
            }
        }
    }

    delta3_sign * r_sign
}

fn lexico_compare_2d(x: &&Point2d, y: &&Point2d) -> Ordering {
    if x[0] < y[0] {
        Ordering::Less
    } else if x[0] > y[0] {
        Ordering::Greater
    } else if x[1] < y[1] {
        Ordering::Less
    } else if x[1] > y[1] {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}
