//! Generic power-side predicates and their keyed symbolic perturbation.
#![allow(clippy::too_many_arguments)] // Public signatures intentionally mirror the PCK predicates.

use crate::expansion::{
    expansion_cofactor, expansion_create, expansion_determinant, expansion_diff, expansion_diff_2,
    expansion_dot_at, expansion_product, expansion_scale, expansion_sq_dist, expansion_sum,
    Expansion,
};
use crate::filter::{side4_3d_filter, FPG_UNCERTAIN_VALUE};
use crate::{Sign, SosPoint};

#[inline]
fn assert_supported_dimension<const D: usize>() {
    assert!(
        matches!(D, 3 | 4 | 6 | 7 | 8),
        "Geogram side predicates support dimensions 3, 4, 6, 7, and 8 (got {D})"
    );
}

pub(crate) fn assert_unique_keys<const D: usize>(points: &[&SosPoint<D>]) {
    for i in 0..points.len() {
        for j in (i + 1)..points.len() {
            assert_ne!(
                points[i].key, points[j].key,
                "SOS keys must be unique within a predicate call (duplicate key {})",
                points[i].key
            );
        }
    }
}

fn side_matrices<const D: usize>(
    points: &[&[f64; D]],
    queries: &[&[f64; D]],
    lifts: &[Expansion],
) -> (Vec<Vec<Expansion>>, Vec<Vec<Expansion>>) {
    let k = points.len() - 1;
    assert_eq!(queries.len(), k);
    assert_eq!(lifts.len(), k);

    let one = expansion_create(1.0);
    let mut coefficients = Vec::with_capacity(k);
    coefficients.push(vec![one.clone(); k]);
    for point in points.iter().take(k).skip(1) {
        let mut row = Vec::with_capacity(k);
        for query in queries {
            row.push(expansion_scale(
                &expansion_dot_at(point, query, points[0]),
                2.0,
            ));
        }
        coefficients.push(row);
    }

    let mut augmented = Vec::with_capacity(k + 1);
    augmented.push(vec![one; k + 1]);
    for (point_index, point) in points.iter().enumerate().skip(1) {
        let mut row = Vec::with_capacity(k + 1);
        for query in queries {
            row.push(expansion_scale(
                &expansion_dot_at(point, query, points[0]),
                2.0,
            ));
        }
        row.push(lifts[point_index - 1].clone());
        augmented.push(row);
    }
    (coefficients, augmented)
}

fn evaluate_side<const D: usize>(
    points: &[&SosPoint<D>],
    queries: &[&[f64; D]],
    lifts: Vec<Expansion>,
    sos: bool,
) -> Sign {
    if sos {
        assert_unique_keys(points);
    }
    let coordinates: Vec<&[f64; D]> = points.iter().map(|point| &point.coords).collect();
    let (coefficient_matrix, augmented) = side_matrices(&coordinates, queries, &lifts);
    let delta = expansion_determinant(&coefficient_matrix);
    let delta_sign = delta.sign();
    assert_ne!(
        delta_sign,
        Sign::Zero,
        "side predicate query simplex is degenerate"
    );

    let result = expansion_determinant(&augmented);
    if result.sign() != Sign::Zero || !sos {
        return delta_sign * result.sign();
    }

    // Reproduce Geogram's side1..side4 SOS decision using the adjugate of the
    // interpolation matrix. For side4 the PSM intentionally evaluates the SOS
    // coefficients with p3 (rather than p4), so retain that exact convention.
    let k = coefficient_matrix.len();
    let adjugate: Vec<Vec<Expansion>> = (0..k)
        .map(|row| {
            (0..k)
                .map(|column| expansion_cofactor(&coefficient_matrix, column, row))
                .collect()
        })
        .collect();
    let evaluation_row = if k == 4 { k - 1 } else { k };
    let evaluation = &augmented[evaluation_row];
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_unstable_by_key(|&index| points[index].key);
    for index in order {
        if index == k {
            return Sign::Negative;
        }
        let coefficient = if index == 0 {
            let mut correction = expansion_create(0.0);
            for row in 0..k {
                let mut row_sum = expansion_create(0.0);
                for value in adjugate[row].iter().skip(1) {
                    row_sum = expansion_sum(&row_sum, value);
                }
                correction =
                    expansion_sum(&correction, &expansion_product(&evaluation[row], &row_sum));
            }
            expansion_diff(&delta, &correction)
        } else {
            let mut coefficient = expansion_create(0.0);
            for row in 0..k {
                coefficient = expansion_sum(
                    &coefficient,
                    &expansion_product(&evaluation[row], &adjugate[row][index]),
                );
            }
            coefficient
        };
        if coefficient.sign() != Sign::Zero {
            return delta_sign * coefficient.sign();
        }
    }
    panic!("SOS perturbation could not resolve a degenerate side predicate");
}

pub(crate) fn side_exact<const D: usize>(
    points: &[&SosPoint<D>],
    queries: &[&[f64; D]],
    sos: bool,
) -> Sign {
    let lifts = points
        .iter()
        .skip(1)
        .map(|point| expansion_sq_dist(&point.coords, &points[0].coords))
        .collect();
    evaluate_side(points, queries, lifts, sos)
}

/// Power-side predicate with one-dimensional query simplex.
pub fn side1_sos<const D: usize>(p0: &SosPoint<D>, p1: &SosPoint<D>, q0: &[f64; D]) -> Sign {
    assert_supported_dimension::<D>();
    side_exact(&[p0, p1], &[q0], true)
}

/// Power-side predicate with two-dimensional query simplex.
pub fn side2_sos<const D: usize>(
    p0: &SosPoint<D>,
    p1: &SosPoint<D>,
    p2: &SosPoint<D>,
    q0: &[f64; D],
    q1: &[f64; D],
) -> Sign {
    assert_supported_dimension::<D>();
    side_exact(&[p0, p1, p2], &[q0, q1], true)
}

/// Power-side predicate with three-dimensional query simplex.
pub fn side3_sos<const D: usize>(
    p0: &SosPoint<D>,
    p1: &SosPoint<D>,
    p2: &SosPoint<D>,
    p3: &SosPoint<D>,
    q0: &[f64; D],
    q1: &[f64; D],
    q2: &[f64; D],
) -> Sign {
    assert_supported_dimension::<D>();
    side_exact(&[p0, p1, p2, p3], &[q0, q1, q2], true)
}

/// Lifted variant of [`side3_sos`] in three dimensions.
pub fn side3_3dlifted_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
    heights: [f64; 4],
    q0: &[f64; 3],
    q1: &[f64; 3],
    q2: &[f64; 3],
) -> Sign {
    let lifts = (1..4)
        .map(|i| expansion_diff_2(heights[i], heights[0]))
        .collect();
    evaluate_side(&[p0, p1, p2, p3], &[q0, q1, q2], lifts, true)
}

pub(crate) fn side3_3dlifted(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    heights: [f64; 4],
    q0: &[f64; 3],
    q1: &[f64; 3],
    q2: &[f64; 3],
) -> Sign {
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
    ];
    let lifts = (1..4)
        .map(|i| expansion_diff_2(heights[i], heights[0]))
        .collect();
    evaluate_side(
        &[&points[0], &points[1], &points[2], &points[3]],
        &[q0, q1, q2],
        lifts,
        false,
    )
}

/// Power-side predicate with four-dimensional query simplex.
pub fn side4_sos<const D: usize>(
    p0: &SosPoint<D>,
    p1: &SosPoint<D>,
    p2: &SosPoint<D>,
    p3: &SosPoint<D>,
    p4: &SosPoint<D>,
    q0: &[f64; D],
    q1: &[f64; D],
    q2: &[f64; D],
    q3: &[f64; D],
) -> Sign {
    assert_supported_dimension::<D>();
    if D == 3 {
        // As in Geogram, intrinsic dimension equals ambient dimension here;
        // the embedding tetrahedron is unnecessary and intentionally ignored.
        let points = [
            SosPoint::new([p0.coords[0], p0.coords[1], p0.coords[2]], p0.key),
            SosPoint::new([p1.coords[0], p1.coords[1], p1.coords[2]], p1.key),
            SosPoint::new([p2.coords[0], p2.coords[1], p2.coords[2]], p2.key),
            SosPoint::new([p3.coords[0], p3.coords[1], p3.coords[2]], p3.key),
            SosPoint::new([p4.coords[0], p4.coords[1], p4.coords[2]], p4.key),
        ];
        return side4_3d_sos(&points[0], &points[1], &points[2], &points[3], &points[4]);
    }
    side_exact(&[p0, p1, p2, p3, p4], &[q0, q1, q2, q3], true)
}

/// Non-SOS three-dimensional side4 predicate. Exact degeneracies return zero.
pub fn side4_3d(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> Sign {
    let filtered = side4_3d_filter(p0, p1, p2, p3, p4);
    if filtered != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(filtered);
    }
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
        SosPoint::new(*p4, 4),
    ];
    side_exact(
        &[&points[0], &points[1], &points[2], &points[3], &points[4]],
        &[p0, p1, p2, p3],
        false,
    )
}

/// Keyed SOS three-dimensional side4 predicate.
pub fn side4_3d_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
    p4: &SosPoint<3>,
) -> Sign {
    assert_unique_keys(&[p0, p1, p2, p3, p4]);
    let filtered = side4_3d_filter(&p0.coords, &p1.coords, &p2.coords, &p3.coords, &p4.coords);
    if filtered != FPG_UNCERTAIN_VALUE {
        return Sign::from_i8(filtered);
    }
    let exact = side_exact(
        &[p0, p1, p2, p3, p4],
        &[&p0.coords, &p1.coords, &p2.coords, &p3.coords],
        false,
    );
    if exact != Sign::Zero {
        return exact;
    }
    side4_3d_sos_decision([p0, p1, p2, p3, p4])
}

fn side4_3d_sos_decision(points: [&SosPoint<3>; 5]) -> Sign {
    let differences: Vec<Vec<Expansion>> = points
        .iter()
        .skip(1)
        .map(|point| {
            (0..3)
                .map(|coordinate| {
                    expansion_diff_2(point.coords[coordinate], points[0].coords[coordinate])
                })
                .collect()
        })
        .collect();
    let minor = |excluded: usize| {
        expansion_determinant(
            &differences
                .iter()
                .enumerate()
                .filter(|(row, _)| *row != excluded)
                .map(|(_, row)| row.clone())
                .collect::<Vec<_>>(),
        )
    };
    let delta1 = minor(0);
    let delta2 = minor(1);
    let delta3 = minor(2);
    let delta4 = minor(3);
    let delta4_sign = delta4.sign();
    assert_ne!(
        delta4_sign,
        Sign::Zero,
        "side4_3d base tetrahedron is degenerate"
    );

    let mut order: Vec<usize> = (0..5).collect();
    order.sort_unstable_by_key(|&index| points[index].key);
    for index in order {
        let coefficient = match index {
            0 => expansion_sum(
                &expansion_diff(&delta2, &delta1),
                &expansion_diff(&delta4, &delta3),
            ),
            1 => delta1.clone(),
            2 => {
                let mut value = delta2.clone();
                value.negate();
                value
            }
            3 => delta3.clone(),
            4 => return Sign::Negative,
            _ => unreachable!(),
        };
        if coefficient.sign() != Sign::Zero {
            return delta4_sign * coefficient.sign();
        }
    }
    panic!("SOS perturbation could not resolve a degenerate side4_3d predicate");
}
