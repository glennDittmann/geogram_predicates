use geogram_predicates_local::{in_circle_2d_sos, SosPoint};
use std::{fs, path::Path};
use test_utils::{next_after, predicate_2d_test, write_to_pgm};

fn main() {
    let args = std::env::args().collect::<Vec<_>>();

    if args.len() != 2 {
        usage()
    }

    let mode = args[1].as_str();

    let p0 = [next_after(12., f64::INFINITY), 12.];
    // let p0 = [12., 12.];
    let p1 = [-12., -12.];
    let p2 = [24., 24.];

    let predicate: Box<dyn Fn([f64; 2]) -> f64> = match mode {
        "naive" => Box::new(|p| naive_incircle_2d(&p0, &p1, &p2, &p)),
        "robust" => Box::new(|p| {
            let points = [
                SosPoint::new(p0, 0),
                SosPoint::new(p1, 1),
                SosPoint::new(p2, 2),
                SosPoint::new(p, 3),
            ];
            f64::from(in_circle_2d_sos(&points[0], &points[1], &points[2], &points[3]).as_i8())
        }),
        "clean" => {
            let _ = fs::remove_file("out_naive_in_circle_2d.pgm");
            let _ = fs::remove_file("out_robust_in_circle_2d.pgm");
            println!("example images removed");
            std::process::exit(1);
        }
        _ => usage(),
    };

    let predicate_results = predicate_2d_test(predicate, [0.5, 0.5], 256, 256);

    let out_path = format!("out_{}_in_circle_2d.pgm", mode);

    write_to_pgm(&predicate_results, Path::new(&out_path), 256, 256);
}

// Directly evaluate the incircle determinant.
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn naive_incircle_2d(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2], p: &[f64; 2]) -> f64 {
    let m11 = a[0] - p[0];
    let m12 = a[1] - p[1];
    let m13 = m11.powi(2) + m12.powi(2);

    let m21 = b[0] - p[0];
    let m22 = b[1] - p[1];
    let m23 = m21.powi(2) + m22.powi(2);

    let m31 = c[0] - p[0];
    let m32 = c[1] - p[1];
    let m33 = m31.powi(2) + m32.powi(2);

    -(m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31))
}

fn usage() -> ! {
    eprintln!(
        "
    Usage:
        in_circle_2d [option]

        MODES:
        naive - output an image showing the output of a naive in_circle_2d implementation
        robust - output an image showing the output of the robust in_circle_2d implementation
        OTHER:
        help - show this help message
        clean - remove the example output images
    "
    );
    std::process::exit(1);
}
