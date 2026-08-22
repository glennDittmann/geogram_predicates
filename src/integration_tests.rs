use crate::{self as gp, Sign, SosPoint};

#[test]
fn lifecycle_is_initialization_free() {
    assert!(gp::is_initialized());
    gp::initialize();
    gp::show_stats();
    gp::terminate();
    assert!(gp::is_initialized());
}

#[test]
fn orientations_and_determinants() {
    let a2 = [0.0, 0.0];
    let b2 = [1.0, 0.0];
    let c2 = [0.0, 1.0];
    assert_eq!(gp::orient_2d(&a2, &b2, &c2), Sign::Positive);
    assert_eq!(gp::orient_2d(&a2, &c2, &b2), Sign::Negative);

    let a3 = [0.0, 0.0, 0.0];
    let b3 = [1.0, 0.0, 0.0];
    let c3 = [0.0, 1.0, 0.0];
    let d3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::orient_3d(&a3, &b3, &c3, &d3), Sign::Positive);
    assert_eq!(gp::det_3d(&b3, &c3, &d3), Sign::Positive);
    assert_eq!(
        gp::det_4d(
            &[1.0, 0.0, 0.0, 0.0],
            &[0.0, 1.0, 0.0, 0.0],
            &[0.0, 0.0, 1.0, 0.0],
            &[0.0, 0.0, 0.0, 1.0],
        ),
        Sign::Positive
    );
}

#[test]
fn circle_and_sphere_exact_and_sos() {
    let points2 = [
        SosPoint::new([0.0, 0.0], 10),
        SosPoint::new([1.0, 0.0], 20),
        SosPoint::new([0.0, 1.0], 30),
        SosPoint::new([0.1, 0.1], 40),
        SosPoint::new([2.0, 2.0], 50),
        SosPoint::new([1.0, 1.0], 5),
    ];
    assert_eq!(
        gp::in_circle_2d(
            &points2[0].coords,
            &points2[1].coords,
            &points2[2].coords,
            &points2[3].coords
        ),
        Sign::Positive
    );
    assert_eq!(
        gp::in_circle_2d_sos(&points2[0], &points2[1], &points2[2], &points2[3]),
        Sign::Positive
    );
    assert_eq!(
        gp::in_circle_2d_sos(&points2[0], &points2[1], &points2[2], &points2[4]),
        Sign::Negative
    );
    assert_eq!(
        gp::in_circle_2d(
            &points2[0].coords,
            &points2[1].coords,
            &points2[2].coords,
            &points2[5].coords
        ),
        Sign::Zero
    );
    assert_ne!(
        gp::in_circle_2d_sos(&points2[0], &points2[1], &points2[2], &points2[5]),
        Sign::Zero
    );

    let points3 = [
        SosPoint::new([0.0, 0.0, 0.0], 10),
        SosPoint::new([1.0, 0.0, 0.0], 20),
        SosPoint::new([0.0, 1.0, 0.0], 30),
        SosPoint::new([0.0, 0.0, 1.0], 40),
        SosPoint::new([0.1, 0.1, 0.1], 50),
        SosPoint::new([2.0, 2.0, 2.0], 60),
        SosPoint::new([1.0, 1.0, 1.0], 5),
    ];
    assert_eq!(
        gp::in_sphere_3d_sos(
            &points3[0],
            &points3[1],
            &points3[2],
            &points3[3],
            &points3[4]
        ),
        Sign::Positive
    );
    assert_eq!(
        gp::in_sphere_3d_sos(
            &points3[0],
            &points3[1],
            &points3[2],
            &points3[3],
            &points3[5]
        ),
        Sign::Negative
    );
    assert_eq!(
        gp::in_sphere_3d(
            &points3[0].coords,
            &points3[1].coords,
            &points3[2].coords,
            &points3[3].coords,
            &points3[6].coords
        ),
        Sign::Zero
    );
    assert_ne!(
        gp::in_sphere_3d_sos(
            &points3[0],
            &points3[1],
            &points3[2],
            &points3[3],
            &points3[6]
        ),
        Sign::Zero
    );
}

#[test]
fn sos_is_deterministic_and_key_independent_off_degeneracy() {
    let make = |keys: [u64; 4]| {
        [
            SosPoint::new([0.0, 0.0], keys[0]),
            SosPoint::new([1.0, 0.0], keys[1]),
            SosPoint::new([0.0, 1.0], keys[2]),
            SosPoint::new([0.2, 0.2], keys[3]),
        ]
    };
    let a = make([10, 20, 30, 40]);
    let b = make([40, 10, 30, 20]);
    assert_eq!(
        gp::in_circle_2d_sos(&a[0], &a[1], &a[2], &a[3]),
        gp::in_circle_2d_sos(&b[0], &b[1], &b[2], &b[3])
    );

    let duplicate_coordinates = [
        SosPoint::new([0.0, 0.0], 10),
        SosPoint::new([1.0, 0.0], 20),
        SosPoint::new([0.0, 1.0], 30),
        SosPoint::new([0.0, 0.0], 5),
    ];
    let first = gp::in_circle_2d_sos(
        &duplicate_coordinates[0],
        &duplicate_coordinates[1],
        &duplicate_coordinates[2],
        &duplicate_coordinates[3],
    );
    let second = gp::in_circle_2d_sos(
        &duplicate_coordinates[0],
        &duplicate_coordinates[1],
        &duplicate_coordinates[2],
        &duplicate_coordinates[3],
    );
    assert_ne!(first, Sign::Zero);
    assert_eq!(first, second);
}

#[test]
fn dot_alignment_and_comparisons() {
    let origin = [0.0, 0.0, 0.0];
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
    let xx = [2.0, 0.0, 0.0];
    assert_eq!(gp::dot_3d(&origin, &x, &y), Sign::Zero);
    assert_eq!(gp::dot_compare_3d(&x, &xx, &y), Sign::Positive);
    assert!(gp::aligned_3d(&origin, &x, &xx));
    assert!(!gp::aligned_3d(&origin, &x, &y));

    let e0 = [1.0, 0.0, 0.0, 0.0];
    let e1 = [0.0, 1.0, 0.0, 0.0];
    let e2 = [0.0, 0.0, 1.0, 0.0];
    let zero = [0.0; 4];
    let e3 = [0.0, 0.0, 0.0, 1.0];
    assert_eq!(
        gp::det_compare_4d(&e0, &e1, &e2, &zero, &e3),
        Sign::Positive
    );
}

#[test]
fn side_predicates_cover_supported_dimensions() {
    fn sides<const D: usize>() {
        let p0 = SosPoint::new([0.0; D], 0);
        let mut c1 = [0.0; D];
        c1[0] = 1.0;
        let p1 = SosPoint::new(c1, 1);
        let mut c2 = [0.0; D];
        c2[1] = 1.0;
        let p2 = SosPoint::new(c2, 2);
        let mut c3 = [0.0; D];
        c3[2] = 1.0;
        let p3 = SosPoint::new(c3, 3);
        let mut c4 = [0.0; D];
        c4[0] = 0.2;
        c4[1] = 0.2;
        c4[2] = 0.2;
        let p4 = SosPoint::new(c4, 4);

        assert_eq!(gp::side1_sos(&p0, &p1, &[0.0; D]), Sign::Positive);
        assert_eq!(
            gp::side2_sos(&p0, &p1, &p2, &p0.coords, &p1.coords),
            Sign::Positive
        );
        assert_eq!(
            gp::side3_sos(&p0, &p1, &p2, &p3, &p0.coords, &p1.coords, &p2.coords,),
            Sign::Positive
        );
        assert_eq!(
            gp::side4_sos(&p0, &p1, &p2, &p3, &p4, &p0.coords, &p1.coords, &p2.coords, &p3.coords,),
            Sign::Negative
        );
        if D == 3 {
            assert_eq!(
                gp::side4_sos(&p0, &p1, &p2, &p3, &p4, &[0.0; D], &[0.0; D], &[0.0; D], &[0.0; D],),
                Sign::Negative
            );
        }
    }
    sides::<3>();
    sides::<4>();
    sides::<6>();
    sides::<7>();
    sides::<8>();

    let points = [
        SosPoint::new([0.0, 0.0, 0.0], 0),
        SosPoint::new([1.0, 0.0, 0.0], 1),
        SosPoint::new([0.0, 1.0, 0.0], 2),
        SosPoint::new([0.0, 0.0, 1.0], 3),
        SosPoint::new([0.2, 0.2, 0.2], 4),
    ];
    assert_eq!(
        gp::side2_sos(
            &points[0],
            &points[1],
            &points[2],
            &points[0].coords,
            &points[1].coords
        ),
        Sign::Positive
    );
    assert_eq!(
        gp::side3_sos(
            &points[0],
            &points[1],
            &points[2],
            &points[3],
            &points[0].coords,
            &points[1].coords,
            &points[2].coords
        ),
        Sign::Positive
    );
    assert_eq!(
        gp::side4_3d_sos(&points[0], &points[1], &points[2], &points[3], &points[4]),
        Sign::Negative
    );
    assert_eq!(
        gp::in_circle_3d_sos(&points[0], &points[1], &points[2], &points[4]),
        Sign::Positive
    );
    assert_eq!(
        gp::in_circle_3dlifted(
            &points[0].coords,
            &points[1].coords,
            &points[2].coords,
            &points[4].coords,
            [0.0, 1.0, 1.0, 0.1]
        ),
        Sign::Positive
    );
    assert_eq!(
        gp::orient_3dlifted(
            &points[0].coords,
            &points[1].coords,
            &points[2].coords,
            &points[3].coords,
            &points[4].coords,
            [0.0, 1.0, 1.0, 1.0, 0.1]
        ),
        Sign::Positive
    );
}

#[test]
fn lifted_predicates_resolve_degeneracy() {
    let p2 = [
        SosPoint::new([0.0, 0.0], 0),
        SosPoint::new([1.0, 0.0], 1),
        SosPoint::new([0.0, 1.0], 2),
        SosPoint::new([0.2, 0.2], 3),
    ];
    assert_eq!(
        gp::orient_2dlifted(
            &p2[0].coords,
            &p2[1].coords,
            &p2[2].coords,
            &p2[3].coords,
            [0.0; 4]
        ),
        Sign::Zero
    );
    assert_ne!(
        gp::orient_2dlifted_sos(&p2[0], &p2[1], &p2[2], &p2[3], [0.0; 4]),
        Sign::Zero
    );

    let p3 = [
        SosPoint::new([0.0, 0.0, 0.0], 0),
        SosPoint::new([1.0, 0.0, 0.0], 1),
        SosPoint::new([0.0, 1.0, 0.0], 2),
        SosPoint::new([0.0, 0.0, 1.0], 3),
        SosPoint::new([0.2, 0.2, 0.2], 4),
    ];
    assert_eq!(
        gp::orient_3dlifted(
            &p3[0].coords,
            &p3[1].coords,
            &p3[2].coords,
            &p3[3].coords,
            &p3[4].coords,
            [0.0; 5]
        ),
        Sign::Zero
    );
    assert_ne!(
        gp::orient_3dlifted_sos(&p3[0], &p3[1], &p3[2], &p3[3], &p3[4], [0.0; 5]),
        Sign::Zero
    );
}

#[test]
#[should_panic(expected = "SOS keys must be unique")]
fn duplicate_sos_keys_panic() {
    let p0 = SosPoint::new([0.0, 0.0], 7);
    let p1 = SosPoint::new([1.0, 0.0], 7);
    let p2 = SosPoint::new([0.0, 1.0], 8);
    let p3 = SosPoint::new([1.0, 1.0], 9);
    let _ = gp::in_circle_2d_sos(&p0, &p1, &p2, &p3);
}

#[test]
#[should_panic(expected = "support dimensions")]
fn unsupported_side_dimension_panics() {
    let p0 = SosPoint::new([0.0; 5], 0);
    let p1 = SosPoint::new([1.0; 5], 1);
    let _ = gp::side1_sos(&p0, &p1, &[0.0; 5]);
}
