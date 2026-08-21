use std::{cmp::Ordering, path::Path};

/// Return the adjacent binary64 value from `value` in the direction of `toward`.
pub fn next_after(value: f64, toward: f64) -> f64 {
    if value.is_nan() || toward.is_nan() {
        return f64::NAN;
    }
    if value == toward {
        return toward;
    }
    if value == 0.0 {
        return if toward > 0.0 {
            f64::from_bits(1)
        } else {
            f64::from_bits((1_u64 << 63) | 1)
        };
    }
    let bits = value.to_bits();
    let increment = (toward > value) == (value > 0.0);
    f64::from_bits(if increment { bits + 1 } else { bits - 1 })
}

/// Write predicate signs as an 8-bit binary portable graymap (P5).
pub fn write_to_pgm(data: &[Ordering], path: &Path, width: usize, height: usize) {
    assert_eq!(data.len(), width * height);

    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "P5\n{width} {height}\n255\n").unwrap();
    let pixels = data
        .iter()
        .map(|w| match w {
            Ordering::Less => 0u8,
            Ordering::Equal => 127,
            Ordering::Greater => 255,
        })
        .collect::<Vec<_>>();
    writer.write_all(&pixels).unwrap();
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
            xd = next_after(xd, f64::INFINITY);
        }
        yd = next_after(yd, f64::INFINITY);
    }
    data
}
