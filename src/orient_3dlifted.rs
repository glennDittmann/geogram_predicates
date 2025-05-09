use crate::{Expansion, Point3d, Sign, FPG_UNCERTAIN_VALUE};
use core::cmp::Ordering;

/// Computes the 4d orientation test with lifted points, i.e the regularity test for 3d.
///
/// Given four lifted points `a'`, `b'`, `c'`, `d'` in R^4, tests if the lifted point `p'` in R^4 lies below or above the hyperplane passing through the four points `a'`, `b'`, `c'`, `d'`.
///
/// Symbolic perturbation is applied, whenever the 5 vertices are not linearly independent.
///
/// The coordinates and the heights are specified in separate arguments for each vertex.
///
/// Note: if `w_i` = 0 this is equal to the in-sphere test for the tetrahedron `a`, `b`, `c`, `d` w.r.t `p`.
///
/// ### Parameters
/// - `a` ,`b`, `c`, `d` vertices of the tetrahedron
/// - `p` point to test
/// - `h_a` ,`h_b` ,`h_c`, `h_d` the heights of the lifted points, e.g. `h_a' = a.x**2 + a.y**2 - a.w`
/// - `h_p` the height of the lifted point `p'`
///
/// ### Return values
/// - `-1` - if `p'` lies below the plane
/// - `+1` - if `p'` lies above the plane
/// - `0`	- if `p'` lies exactly on the hyperplane
///
/// # Example
/// ```
/// use geogram_predicates::orient_3dlifted_sos;
/// // Define four points that form a tetrahedron
/// let a: [f64; 3] = [0.0, 0.0, 0.0];
/// let b: [f64; 3] = [2.0, 0.0, 0.0];
/// let c: [f64; 3] = [0.0, 2.0, 0.0];
/// let d: [f64; 3] = [0.75, 0.75, 1.0];
///
/// // Additionally in this scenario, each point is associated with a weight w_i
/// // And the height of a point is defined as h_i = x_i**2 + y_i**2 + z_i**2 - w_i
/// // One can interpret the height as the 4th-coordinate of a point lifted to R^4
/// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + a[2].powf(2.0) + 2.0;  // i.e. w_a = -2.0
/// let h_b = b[0].powf(2.0) + b[1].powf(2.0) + b[2].powf(2.0) - 1.0;  // i.e. w_b = 1.0
/// let h_c = c[0].powf(2.0) + c[1].powf(2.0) + c[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
/// let h_d = d[0].powf(2.0) + d[1].powf(2.0) + d[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
///
/// // Define weighted points, to test against the hyperplane, that contains the lifted tetrahedron
/// let p_below: [f64; 3] = [0.6, 0.6, 0.6];
/// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 0.28;  // i.e. w_p_below = -0.28
///
/// let p_above: [f64; 3] = [0.6, 0.6, 0.6];
/// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
///
/// let orientation_below = orient_3dlifted_sos(&a, &b, &c, &d, &p_below, [h_a, h_b, h_c, h_d, h_p_below]);
/// assert_eq!(1, orientation_below);
///
/// let orientation_above = orient_3dlifted_sos(&a, &b, &c, &d, &p_above, [h_a, h_b, h_c, h_d, h_p_above]);
/// assert_eq!(-1, orientation_above);
///
/// ```
pub fn orient_3dlifted_sos(
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
    d: &[f64; 3],
    p: &[f64; 3],
    [h_a, h_b, h_c, h_d, h_p]: [f64; 5],
) -> Sign {
    let mut result = side4h_3d_filter(a, b, c, d, p, [h_a, h_b, h_c, h_d, h_p]);
    if result == 0 {
        result = side4h_3d_exact_sos::<true>(a, b, c, d, p, [h_a, h_b, h_c, h_d, h_p]);
    }

    // orient_4d() is opposite to side4h()
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

    let mut eps;

    let max1 = (a11.abs()).max(a21.abs()).max(a31.abs());
    let max2 = (a12.abs()).max(a13.abs()).max(a22.abs()).max(a23.abs());
    let max3 = (a22.abs()).max(a23.abs()).max(a32.abs()).max(a33.abs());

    let int_tmp_result;
    let mut lower_bound_1 = max1;
    let mut upper_bound_1 = max1;
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
        return FPG_UNCERTAIN_VALUE;
    } else {
        if upper_bound_1 > 7.23700557733225980357e+75 {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = 5.11071278299732992696e-15 * ((max2 * max3) * max1);
        if delta4 > eps {
            int_tmp_result = Sign::Positive;
        } else if delta4 < -eps {
            int_tmp_result = Sign::Negative;
        } else {
            return FPG_UNCERTAIN_VALUE;
        }
    }

    let delta4_sign = int_tmp_result;
    let int_tmp_result_ffwkcaa: Sign;

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
        return FPG_UNCERTAIN_VALUE;
    } else {
        if upper_bound_1 > 7.23700557733225980357e+75 {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = 3.17768858673611390687e-14 * (((max5 * max7) * max4) * max6);
        if r > eps {
            int_tmp_result_ffwkcaa = Sign::Positive;
        } else if r < -eps {
            int_tmp_result_ffwkcaa = Sign::Negative;
        } else {
            return FPG_UNCERTAIN_VALUE;
        }
    }

    delta4_sign * int_tmp_result_ffwkcaa
}

fn side4h_3d_exact_sos<const SOS: bool>(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    p4: &[f64; 3],
    [h0, h1, h2, h3, h4]: [f64; 5],
) -> Sign {
    use crate::expansion::{
        expansion_det3x3, expansion_diff, expansion_product, expansion_sum, expansion_sum4,
    };

    let a11 = expansion_diff!(p1[0], p0[0]);
    let a12 = expansion_diff!(p1[1], p0[1]);
    let a13 = expansion_diff!(p1[2], p0[2]);
    let a14 = expansion_diff!(h0, h1);

    let a21 = expansion_diff!(p2[0], p0[0]);
    let a22 = expansion_diff!(p2[1], p0[1]);
    let a23 = expansion_diff!(p2[2], p0[2]);
    let a24 = expansion_diff!(h0, h2);

    let a31 = expansion_diff!(p3[0], p0[0]);
    let a32 = expansion_diff!(p3[1], p0[1]);
    let a33 = expansion_diff!(p3[2], p0[2]);
    let a34 = expansion_diff!(h0, h3);

    let a41 = expansion_diff!(p4[0], p0[0]);
    let a42 = expansion_diff!(p4[1], p0[1]);
    let a43 = expansion_diff!(p4[2], p0[2]);
    let a44 = expansion_diff!(h0, h4);

    // This could probably reuse some of the 2x2 co-factors, but this is easier to read
    let delta1: Expansion<6> = expansion_det3x3!(a21, a22, a23, a31, a32, a33, a41, a42, a43);
    let delta2: Expansion<6> = expansion_det3x3!(a11, a12, a13, a31, a32, a33, a41, a42, a43);
    let delta3: Expansion<6> = expansion_det3x3!(a11, a12, a13, a21, a22, a23, a41, a42, a43);
    let delta4: Expansion<6> = expansion_det3x3!(a11, a12, a13, a21, a22, a23, a31, a32, a33);

    let delta4_sign = delta4.sign();
    debug_assert!(delta4_sign != 0);

    let r_1: Expansion<2> = expansion_product!(delta1, a14);
    let mut r_2: Expansion<2> = expansion_product!(delta2, a24);
    r_2.negate();
    let r_3: Expansion<2> = expansion_product!(delta3, a34);
    let mut r_4: Expansion<2> = expansion_product!(delta4, a44);
    r_4.negate();
    let r: Expansion<4> = expansion_sum4!(r_1, r_2, r_3, r_4, 4, 4);

    let r_sign = r.sign();

    // Simulation of Simplicity (symbolic perturbation)
    if SOS && r_sign == 0 {
        let mut p_sort = [p0, p1, p2, p3, p4];
        p_sort.sort_by(lexico_compare_3d);

        for i in 0..4 {
            if p_sort[i] == p0 {
                let z1: Expansion<2> = expansion_diff!(Expansions: delta2, delta1);
                let z2: Expansion<2> = expansion_diff!(Expansions: delta4, delta3);
                let z: Expansion<1> = expansion_sum!(z1, z2);
                let z_sign = z.sign();
                if z_sign != 0 {
                    return delta4_sign * z_sign;
                }
            } else if p_sort[i] == p1 {
                let delta1_sign = delta1.sign();
                if delta1_sign != 0 {
                    return delta4_sign * delta1_sign;
                }
            } else if p_sort[i] == p2 {
                let delta2_sign = delta2.sign();
                if delta2_sign != 0 {
                    return (-delta4_sign) * delta2_sign;
                }
            } else if p_sort[i] == p3 {
                let delta3_sign = delta3.sign();
                if delta3_sign != 0 {
                    return delta4_sign * delta3_sign;
                }
            } else if p_sort[i] == p4 {
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
        x[2].total_cmp(&y[2])
    }
}
