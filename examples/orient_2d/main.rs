
use geogram_predicates::orient_2d;
use std::cmp::Ordering;
use std::path::Path;

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
    
    let predicate_results = orient_2d_test(predicate, [0.5, 0.5], 256, 256);
    
    let out_path = format!("out_{}_orient_2d.png", mode);
    
    write_to_png(&predicate_results, Path::new(&out_path), 256, 256);}

// Directly evaluate the orient2d determinant.
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn naive_orient_2d(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
}

fn orient_2d_test<F>(predicate: F, start: [f64; 2], width: usize, height: usize) -> Vec<Ordering>
where
    F: Fn([f64; 2]) -> f64,
{
    use float_extras::f64::nextafter;
    let mut yd = start[1];
    let mut data = Vec::with_capacity(width * height);

    for _ in 0..height {
        let mut xd = start[0];
        for _ in 0..width {
            let p = [xd, yd];
            data.push(predicate(p).partial_cmp(&0.).unwrap());
            xd = nextafter(xd, std::f64::INFINITY);
        }
        yd = nextafter(yd, std::f64::INFINITY);
    }

    data
}

fn write_to_png(data: &[Ordering], path: &Path, width: usize, height: usize) {
    assert_eq!(data.len(), width * height);

    use std::fs::File;
    use std::io::BufWriter;

    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width as u32, height as u32);
    encoder.set_color(png::ColorType::Grayscale);
    encoder.set_depth(png::BitDepth::Eight);

    let mut writer = encoder.write_header().unwrap();
    let data = data
        .iter()
        .map(|w| match w {
            Ordering::Less => 0u8,
            Ordering::Equal => 127,
            Ordering::Greater => 255,
        })
        .collect::<Vec<_>>();
    writer.write_image_data(&data).unwrap();
}

fn usage(name: &str) -> ! {
    eprintln!(
        "Usage: {} {{naive | robust}}",
        name
    );
    std::process::exit(1);
}