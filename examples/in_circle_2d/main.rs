
use float_extras::f64::nextafter;
use geogram_predicates::in_circle_2d_SOS;
use std::path::Path;
use test_utils::{predicate_2d_test, write_to_png};

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    
    if args.len() != 2 {
        usage(&args[0])
    }

    let mode = args[1].as_str();

    let p0 = [nextafter(12., std::f64::INFINITY), 12.];
    // let p0 = [12., 12.];
    let p1 = [-12., -12.];
    let p2 = [24., 24.];
    
    let predicate: Box<dyn Fn([f64; 2]) -> f64> = match mode {
        "naive" => Box::new(|p| naive_incircle_2d(&p0, &p1, &p2, &p)),
        "robust" => Box::new(|p| in_circle_2d_SOS(&p0, &p1, &p2, &p).into()),
        _ => unimplemented!(),
    };
    
    let predicate_results = predicate_2d_test(predicate, [0.5, 0.5], 256, 256);
    
    let out_path = format!("out_{}_in_circle_2d.png", mode);
    
    write_to_png(&predicate_results, Path::new(&out_path), 256, 256);}

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

    (m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31)) * -1.0
}

fn usage(name: &str) -> ! {
    eprintln!(
        "Usage: {} {{naive | robust}}",
        name
    );
    std::process::exit(1);
}