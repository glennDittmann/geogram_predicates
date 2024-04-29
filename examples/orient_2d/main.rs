
use geogram_predicates::orient_2d;
use std::path::Path;
use test_utils::{predicate_2d_test, write_to_png};


fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    
    if args.len() != 2 {
        usage(&args[0])
    }

    let mode = args[1].as_str();

    
    let p1 = [12., 12.];
    let p2 = [24., 24.];
    
    let predicate: Box<dyn Fn([f64; 2]) -> f64> = match mode {
        "naive" => Box::new(|p| naive_orient_2d(&p1, &p, &p2)),
        "robust" => Box::new(|p| orient_2d(&p1, &p, &p2).into()),
        _ => unimplemented!(),
    };
    
    let predicate_results = predicate_2d_test(predicate, [0.5, 0.5], 256, 256);
    
    let out_path = format!("out_{}_orient_2d.png", mode);
    
    write_to_png(&predicate_results, Path::new(&out_path), 256, 256);}

// Directly evaluate the orient2d determinant.
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn naive_orient_2d(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
}

fn usage(name: &str) -> ! {
    eprintln!(
        "Usage: {} {{naive | robust}}",
        name
    );
    std::process::exit(1);
}