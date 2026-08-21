//! Exact lifted orientation predicates.

use crate::expansion::{
    expansion_cofactor, expansion_create, expansion_determinant, expansion_diff_2, expansion_sum,
    Expansion,
};
use crate::{Sign, SosPoint};

fn assert_unique_keys<const D: usize>(points: &[&SosPoint<D>]) {
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

fn lifted_exact<const D: usize>(
    points: &[&SosPoint<D>],
    heights: &[f64],
    sos: bool,
    result_factor: Sign,
) -> Sign {
    assert_eq!(points.len(), D + 2);
    assert_eq!(heights.len(), points.len());
    if sos {
        assert_unique_keys(points);
    }

    let mut spatial = Vec::with_capacity(D);
    for point in points.iter().take(D + 1).skip(1) {
        spatial.push(
            (0..D)
                .map(|coordinate| {
                    expansion_diff_2(point.coords[coordinate], points[0].coords[coordinate])
                })
                .collect::<Vec<Expansion>>(),
        );
    }
    let delta_sign = expansion_determinant(&spatial).sign();
    assert_ne!(
        delta_sign,
        Sign::Zero,
        "lifted orientation requires a nondegenerate base simplex"
    );

    let mut matrix = Vec::with_capacity(D + 1);
    for (index, point) in points.iter().enumerate().skip(1) {
        let mut row = (0..D)
            .map(|coordinate| {
                expansion_diff_2(point.coords[coordinate], points[0].coords[coordinate])
            })
            .collect::<Vec<Expansion>>();
        row.push(expansion_diff_2(heights[0], heights[index]));
        matrix.push(row);
    }

    let determinant = expansion_determinant(&matrix);
    if determinant.sign() != Sign::Zero || !sos {
        return result_factor * delta_sign * determinant.sign();
    }

    let last_column = D;
    let mut p0_coefficient = expansion_create(0.0);
    let mut coefficients = Vec::with_capacity(points.len());
    for row in 0..matrix.len() {
        let cofactor = expansion_cofactor(&matrix, row, last_column);
        let last_2d_point = D == 2 && row + 1 == matrix.len();
        if last_2d_point {
            p0_coefficient = expansion_sum(&p0_coefficient, &cofactor);
            let mut point_coefficient = cofactor;
            point_coefficient.negate();
            coefficients.push(point_coefficient);
        } else {
            let mut negative_cofactor = cofactor.clone();
            negative_cofactor.negate();
            p0_coefficient = expansion_sum(&p0_coefficient, &negative_cofactor);
            coefficients.push(cofactor);
        }
    }
    // These signs mirror Geogram's specialized side3h/side4h SOS branches.
    // The final point in the 2D variant uses the terminal `NEGATIVE` term from
    // side3h_2d_exact_SOS; the 3D variant is the direct cofactor ordering.
    coefficients.insert(0, p0_coefficient);

    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_unstable_by_key(|&index| points[index].key);
    for index in order {
        let coefficient_sign = coefficients[index].sign();
        if coefficient_sign != Sign::Zero {
            return result_factor * delta_sign * coefficient_sign;
        }
    }
    panic!("SOS perturbation could not resolve a lifted orientation");
}

/// Exact non-SOS lifted orientation in 2D. Degeneracies return zero.
pub fn orient_2dlifted(
    p0: &[f64; 2],
    p1: &[f64; 2],
    p2: &[f64; 2],
    p3: &[f64; 2],
    heights: [f64; 4],
) -> Sign {
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
    ];
    lifted_exact(
        &[&points[0], &points[1], &points[2], &points[3]],
        &heights,
        false,
        Sign::Positive,
    )
}

/// Exact keyed-SOS lifted orientation in 2D.
pub fn orient_2dlifted_sos(
    p0: &SosPoint<2>,
    p1: &SosPoint<2>,
    p2: &SosPoint<2>,
    p3: &SosPoint<2>,
    heights: [f64; 4],
) -> Sign {
    lifted_exact(&[p0, p1, p2, p3], &heights, true, Sign::Positive)
}

/// Exact non-SOS lifted orientation in 3D. Degeneracies return zero.
pub fn orient_3dlifted(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    p4: &[f64; 3],
    heights: [f64; 5],
) -> Sign {
    let points = [
        SosPoint::new(*p0, 0),
        SosPoint::new(*p1, 1),
        SosPoint::new(*p2, 2),
        SosPoint::new(*p3, 3),
        SosPoint::new(*p4, 4),
    ];
    lifted_exact(
        &[&points[0], &points[1], &points[2], &points[3], &points[4]],
        &heights,
        false,
        Sign::Positive,
    )
}

/// Exact keyed-SOS lifted orientation in 3D.
pub fn orient_3dlifted_sos(
    p0: &SosPoint<3>,
    p1: &SosPoint<3>,
    p2: &SosPoint<3>,
    p3: &SosPoint<3>,
    p4: &SosPoint<3>,
    heights: [f64; 5],
) -> Sign {
    lifted_exact(&[p0, p1, p2, p3, p4], &heights, true, Sign::Positive)
}
