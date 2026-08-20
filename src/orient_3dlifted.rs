//! 3D lifted orientation (regularity test for 3D).

use crate::expansion::{
    expansion_det3x3, expansion_diff, expansion_diff_2, expansion_product, expansion_sum4,
};
use crate::pck::Point3d;
use crate::sign::Sign;
use core::cmp::Ordering;

/// Computes the 4D orientation test with lifted points (regularity test for 3D).
/// orient_4d is opposite to side4h: so we return -result. Negative = p' below, Positive = above.
pub fn orient_3dlifted_sos(
    a: &Point3d,
    b: &Point3d,
    c: &Point3d,
    d: &Point3d,
    p: &Point3d,
    [h_a, h_b, h_c, h_d, h_p]: [f64; 5],
) -> Sign {
    let mut result = side4h_3d_filter(a, b, c, d, p, [h_a, h_b, h_c, h_d, h_p]);
    if result == Sign::Zero {
        result = side4h_3d_exact_sos(a, b, c, d, p, [h_a, h_b, h_c, h_d, h_p]);
    }
    -result
}

#[inline]
fn side4h_3d_filter(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    p4: &[f64; 3],
    [h0, h1, h2, h3, h4]: [f64; 5],
) -> Sign {
    let a11 = p1[0] - p0[0];
    let a12 = p1[1] - p0[1];
    let a13 = p1[2] - p0[2];
    let a14 = h0 - h1;
    let a21 = p2[0] - p0[0];
    let a22 = p2[1] - p0[1];
    let a23 = p2[2] - p0[2];
    let a24 = h0 - h2;
    let a31 = p3[0] - p0[0];
    let a32 = p3[1] - p0[1];
    let a33 = p3[2] - p0[2];
    let a34 = h0 - h3;
    let a41 = p4[0] - p0[0];
    let a42 = p4[1] - p0[1];
    let a43 = p4[2] - p0[2];
    let a44 = h0 - h4;

    let delta1 = (a21 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a22 * a43) - (a23 * a42)))
        + (a41 * ((a22 * a33) - (a23 * a32)));
    let delta2 = (a11 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a12 * a43) - (a13 * a42)))
        + (a41 * ((a12 * a33) - (a13 * a32)));
    let delta3 = (a11 * ((a22 * a43) - (a23 * a42))) - (a21 * ((a12 * a43) - (a13 * a42)))
        + (a41 * ((a12 * a23) - (a13 * a22)));
    let delta4 = (a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))
        + (a31 * ((a12 * a23) - (a13 * a22)));
    let r = (((delta1 * a14) - (delta2 * a24)) + (delta3 * a34)) - (delta4 * a44);

    let max1 = (a11.abs()).max(a21.abs()).max(a31.abs());
    let max2 = (a12.abs()).max(a13.abs()).max(a22.abs()).max(a23.abs());
    let max3 = (a22.abs()).max(a23.abs()).max(a32.abs()).max(a33.abs());
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
    if lower_bound_1 < 1.63288018496748314939e-98 {
        return Sign::Zero;
    }
    if upper_bound_1 > 7.23700557733225980357e+75 {
        return Sign::Zero;
    }
    let eps = 5.11071278299732992696e-15 * ((max2 * max3) * max1);
    let int_tmp_result = if delta4 > eps {
        Sign::Positive
    } else if delta4 < -eps {
        Sign::Negative
    } else {
        return Sign::Zero;
    };
    let delta4_sign = int_tmp_result;

    let max4 = max1.max(a41.abs());
    let max5 = max2.max(max3);
    let max6 = (a14.abs()).max(a24.abs()).max(a34.abs()).max(a44.abs());
    let max7 = max3.max(a42.abs()).max(a43.abs());
    lower_bound_1 = max4;
    upper_bound_1 = max4;
    if max5 < lower_bound_1 {
        lower_bound_1 = max5;
    } else if max5 > upper_bound_1 {
        upper_bound_1 = max5;
    }
    if max6 < lower_bound_1 {
        lower_bound_1 = max6;
    } else if max6 > upper_bound_1 {
        upper_bound_1 = max6;
    }
    if max7 < lower_bound_1 {
        lower_bound_1 = max7;
    } else if max7 > upper_bound_1 {
        upper_bound_1 = max7;
    }
    if lower_bound_1 < 2.89273249588395194294e-74 {
        return Sign::Zero;
    }
    if upper_bound_1 > 7.23700557733225980357e+75 {
        return Sign::Zero;
    }
    let eps2 = 3.17768858673611390687e-14 * (((max5 * max7) * max4) * max6);
    let int_tmp_result_ffwkcaa = if r > eps2 {
        Sign::Positive
    } else if r < -eps2 {
        Sign::Negative
    } else {
        return Sign::Zero;
    };
    delta4_sign * int_tmp_result_ffwkcaa
}

fn side4h_3d_exact_sos(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    p4: &[f64; 3],
    [h0, h1, h2, h3, h4]: [f64; 5],
) -> Sign {
    let a11 = expansion_diff_2(p1[0], p0[0]);
    let a12 = expansion_diff_2(p1[1], p0[1]);
    let a13 = expansion_diff_2(p1[2], p0[2]);
    let a14 = expansion_diff_2(h0, h1);
    let a21 = expansion_diff_2(p2[0], p0[0]);
    let a22 = expansion_diff_2(p2[1], p0[1]);
    let a23 = expansion_diff_2(p2[2], p0[2]);
    let a24 = expansion_diff_2(h0, h2);
    let a31 = expansion_diff_2(p3[0], p0[0]);
    let a32 = expansion_diff_2(p3[1], p0[1]);
    let a33 = expansion_diff_2(p3[2], p0[2]);
    let a34 = expansion_diff_2(h0, h3);
    let a41 = expansion_diff_2(p4[0], p0[0]);
    let a42 = expansion_diff_2(p4[1], p0[1]);
    let a43 = expansion_diff_2(p4[2], p0[2]);
    let a44 = expansion_diff_2(h0, h4);

    let delta1 = expansion_det3x3(
        &a21, &a22, &a23, &a31, &a32, &a33, &a41, &a42, &a43,
    );
    let delta2 = expansion_det3x3(
        &a11, &a12, &a13, &a31, &a32, &a33, &a41, &a42, &a43,
    );
    let delta3 = expansion_det3x3(
        &a11, &a12, &a13, &a21, &a22, &a23, &a41, &a42, &a43,
    );
    let delta4 = expansion_det3x3(
        &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33,
    );

    let delta4_sign = delta4.sign();
    if delta4_sign == Sign::Zero {
        return Sign::Zero;
    }

    let r_1 = expansion_product(&delta1, &a14);
    let mut r_2 = expansion_product(&delta2, &a24);
    r_2.negate();
    let r_3 = expansion_product(&delta3, &a34);
    let mut r_4 = expansion_product(&delta4, &a44);
    r_4.negate();
    let r = expansion_sum4(&r_1, &r_2, &r_3, &r_4);
    let r_sign = r.sign();

    if r_sign == Sign::Zero {
        let mut p_sort = [p0, p1, p2, p3, p4];
        p_sort.sort_by(lexico_compare_3d);
        for i in 0..5 {
            if std::ptr::eq(p_sort[i], p0) {
                let z1 = expansion_diff(&delta2, &delta1);
                let z2 = expansion_diff(&delta4, &delta3);
                let z = crate::expansion::expansion_sum(&z1, &z2);
                let z_sign = z.sign();
                if z_sign != Sign::Zero {
                    return delta4_sign * z_sign;
                }
            } else if std::ptr::eq(p_sort[i], p1) {
                let delta1_sign = delta1.sign();
                if delta1_sign != Sign::Zero {
                    return delta4_sign * delta1_sign;
                }
            } else if std::ptr::eq(p_sort[i], p2) {
                let delta2_sign = delta2.sign();
                if delta2_sign != Sign::Zero {
                    return (-delta4_sign) * delta2_sign;
                }
            } else if std::ptr::eq(p_sort[i], p3) {
                let delta3_sign = delta3.sign();
                if delta3_sign != Sign::Zero {
                    return delta4_sign * delta3_sign;
                }
            } else if std::ptr::eq(p_sort[i], p4) {
                return Sign::Negative;
            }
        }
    }

    delta4_sign * r_sign
}

fn lexico_compare_3d(x: &&Point3d, y: &&Point3d) -> Ordering {
    if x[0] < y[0] {
        Ordering::Less
    } else if x[0] > y[0] {
        Ordering::Greater
    } else if x[1] < y[1] {
        Ordering::Less
    } else if x[1] > y[1] {
        Ordering::Greater
    } else {
        x[2].partial_cmp(&y[2]).unwrap_or(Ordering::Equal)
    }
}
