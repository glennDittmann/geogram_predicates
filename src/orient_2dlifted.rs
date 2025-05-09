pub use crate::Expansion;
use crate::{Point2d, Sign, FPG_UNCERTAIN_VALUE};
use core::cmp::Ordering;

/// Computes the 3d orientation test with lifted points, i.e the regularity test for 2d.
///
/// Given three lifted points `a'`, `b'`, `c'` in R^3, tests if the lifted point `p'` in R^3 lies below or above the plane passing through the three points `a'`, `b'`, `c'`.
///
/// The coordinates and the heights are specified in separate arguments for each vertex.
///
/// Note: if `w_i` = 0 this is equal to the in-circle test for the triangle `a`, `b`, `c` w.r.t `p`.
///
/// ### Parameters
/// - `a` ,`b`, `c`	vertices of the triangle
/// - `p` point to test
/// - `h_a` ,`h_b` ,`h_c` the heights of the lifted points, e.g. `a' = a.x**2 + a.y**2 - a.w`
/// - `h_p` the height of the lifted point `p`
///
/// ### Return values
/// - `+1` - if p3' lies below the plane
/// - `-1` - if p3' lies above the plane
/// - `0`	- if `p'` lies exactly on the hyperplane
///
/// # Example
/// For a graphical representation see this [geogebra example](https://www.geogebra.org/m/etyzj96t) of the code below.
/// ```
/// use geogram_predicates::orient_2dlifted_sos;
///
/// // Define three points that form a triangle
/// let a: [f64; 2] = [0.0, 0.0];
/// let b: [f64; 2] = [2.0, 0.0];
/// let c: [f64; 2] = [0.0, 2.0];
///
/// // Additionally in this scenario, each point is associated with a weight w_i
/// // And the height of a point is defined as h_i = x_i**2 + y_i**2 - w_i
/// // One can interpret the height as the z-coordinate of a point lifted to R^3
/// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + 2.0;  // i.e. w_a = -2.0
/// let h_b = b[0].powf(2.0) + b[1].powf(2.0) - 1.0;  // i.e. w_b = 1.0
/// let h_c = c[0].powf(2.0) + c[1].powf(2.0) - 0.5;  // i.e. w_c = 0.5
///
/// // Define weighted points, to test against the plane, that contains the lifted triangle
/// let p_below: [f64; 2] = [0.6, 0.6];
/// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 1.28;  // i.e. w_p_below = -1.28
///
/// let p_above: [f64; 2] = [0.6, 0.6];
/// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
///
/// let orientation_below = orient_2dlifted_sos(&a, &b, &c, &p_below, [h_a, h_b, h_c, h_p_below]);
/// assert_eq!(1, orientation_below);
///
/// let orientation_above = orient_2dlifted_sos(&a, &b, &c, &p_above, [h_a, h_b, h_c, h_p_above]);
/// assert_eq!(-1, orientation_above);
/// ```
pub fn orient_2dlifted_sos(
    a: &Point2d,
    b: &Point2d,
    c: &Point2d,
    p: &Point2d,
    [h_a, h_b, h_c, h_p]: [f64; 4],
) -> Sign {
    let mut result = side3_2dlifted_2d_filter(a, b, c, p, [h_a, h_b, h_c, h_p]);
    if result == 0 {
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

    let mut eps;
    let max1 = (a11.abs()).max(a12.abs());
    let max2 = (a21.abs()).max(a22.abs());

    let int_tmp_result;
    let mut lower_bound_1 = max1;
    let mut upper_bound_1 = max1;
    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if lower_bound_1 < 5.00368081960964635413e-147 {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if upper_bound_1 > 5.59936185544450928309e+101 {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = 8.88720573725927976811e-16 * (max1 * max2);
        if delta3 > eps {
            int_tmp_result = Sign::Positive;
        } else if delta3 < -eps {
            int_tmp_result = Sign::Negative;
        } else {
            return FPG_UNCERTAIN_VALUE;
        }
    }

    let delta3_sign = int_tmp_result;
    let int_tmp_result_ffwkcaa;
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

    // I think that is a bit out of the f64 range
    if lower_bound_1 < 1.63288018496748314939e-98 {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if upper_bound_1 > 5.59936185544450928309e+101 {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = 5.11071278299732992696e-15 * ((max3 * max5) * max4);
        if r > eps {
            int_tmp_result_ffwkcaa = Sign::Positive;
        } else if r < -eps {
            int_tmp_result_ffwkcaa = Sign::Negative;
        } else {
            return FPG_UNCERTAIN_VALUE;
        }
    }

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
    use crate::expansion::{expansion_det2x2, expansion_diff, expansion_product, expansion_sum};

    let sos = sos.unwrap_or(true);

    let a11: Expansion<2> = expansion_diff!(p1[0], p0[0], 1, 1);
    let a12: Expansion<2> = expansion_diff!(p1[1], p0[1], 1, 1);
    let a13: Expansion<2> = expansion_diff!(h0, h1, 1, 1);

    let a21: Expansion<2> = expansion_diff!(p2[0], p0[0], 1, 1);
    let a22: Expansion<2> = expansion_diff!(p2[1], p0[1], 1, 1);
    let a23: Expansion<2> = expansion_diff!(h0, h2, 1, 1);

    let a31: Expansion<2> = expansion_diff!(p3[0], p0[0], 1, 1);
    let a32: Expansion<2> = expansion_diff!(p3[1], p0[1], 1, 1);
    let a33: Expansion<2> = expansion_diff!(h0, h3, 1, 1);

    let delta1: Expansion<3> = expansion_det2x2!([a21, a22], [a31, a32]);
    let delta2: Expansion<3> = expansion_det2x2!([a11, a12], [a31, a32]);
    let delta3: Expansion<3> = expansion_det2x2!([a11, a12], [a21, a22]);

    let delta3_sign = delta3.sign();
    debug_assert!(delta3_sign != 0);

    let r_1: Expansion<4> = expansion_product!(delta1, a13, 4);
    let mut r_2: Expansion<4> = expansion_product!(delta2, a23, 4);
    r_2.negate();
    let r_3: Expansion<4> = expansion_product!(delta3, a33, 4);
    let r: Expansion<6> = {
        // capacity is `r_1.length() + r_2.length()` which is 4
        let mut ab: Expansion<4> = Expansion::new();
        ab.assign_sum(&r_1, &r_2);
        let mut expansion: Expansion<6> = Expansion::new();
        expansion.assign_sum(&ab, &r_3);
        expansion
    };

    let r_sign = r.sign();

    // Simulation of Simplicity (symbolic perturbation)
    if sos && r_sign == 0 {
        let mut p_sort = [p0, p1, p2, p3];
        p_sort.sort_by(lexico_compare_2d);

        for i in 0..3 {
            if p_sort[i] == p0 {
                let z1 = {
                    let mut expansion: Expansion<2> = Expansion::with_capacity(2);
                    expansion.assign_diff(&delta2, &delta1);
                    expansion
                };
                let z: Expansion<2> = expansion_sum!(z1, delta3);
                let z_sign = z.sign();
                if z_sign != 0 {
                    return delta3_sign * z_sign;
                }
            } else if p_sort[i] == p1 {
                let delta1_sign = delta1.sign();
                if delta1_sign != 0 {
                    return delta3_sign * delta1_sign;
                }
            } else if p_sort[i] == p2 {
                let delta2_sign = delta2.sign();
                if delta2_sign != 0 {
                    return (-delta3_sign) * delta2_sign;
                }
            } else if p_sort[i] == p3 {
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
        x[0].total_cmp(&y[0])
    }
}
