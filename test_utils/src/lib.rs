use float_extras::f64::nextafter;
use std::{cmp::Ordering, path::Path};

pub fn write_to_png(data: &[Ordering], path: &Path, width: usize, height: usize) {
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

pub fn predicate_2d_test<F>(
    predicate: F,
    start: [f64; 2],
    width: usize,
    height: usize,
) -> Vec<Ordering>
where
    F: Fn([f64; 2]) -> f64,
{
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
