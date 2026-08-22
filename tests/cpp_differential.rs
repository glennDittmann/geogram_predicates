//! Opt-in differential test against the bundled Geogram PSM.
//!
//! Run with:
//! `cargo test --release --test cpp_differential -- --ignored --nocapture`

use geogram_predicates as gp;
use gp::SosPoint;
use std::env;
use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;

const DEFAULT_CASES: u64 = 10_000;
const DEFAULT_SEED: u64 = 0x6a09_e667_f3bc_c909;
const SIDE_DIMENSIONS: [usize; 5] = [3, 4, 6, 7, 8];

#[derive(Clone, Copy, Debug)]
#[repr(u8)]
enum Predicate {
    Orient2d = 0,
    Orient3d = 1,
    Det3d = 2,
    Det4d = 3,
    DetCompare4d = 4,
    Dot3d = 5,
    DotCompare3d = 6,
    Aligned3d = 7,
    Identical2d = 8,
    Identical3d = 9,
    Colinear3d = 10,
    InCircle2dSos = 11,
    InCircle3dSos = 12,
    InSphere3dSos = 13,
    Side1Sos = 14,
    Side2Sos = 15,
    Side3Sos = 16,
    Side4Sos = 17,
    Side4_3d = 18,
    Side4_3dSos = 19,
    Side3_3dLifted = 20,
    Side3_3dLiftedSos = 21,
    InCircle3dLifted = 22,
    InCircle3dLiftedSos = 23,
    Orient2dLiftedSos = 24,
    Orient3dLifted = 25,
    Orient3dLiftedSos = 26,
}

const PREDICATES: [Predicate; 27] = [
    Predicate::Orient2d,
    Predicate::Orient3d,
    Predicate::Det3d,
    Predicate::Det4d,
    Predicate::DetCompare4d,
    Predicate::Dot3d,
    Predicate::DotCompare3d,
    Predicate::Aligned3d,
    Predicate::Identical2d,
    Predicate::Identical3d,
    Predicate::Colinear3d,
    Predicate::InCircle2dSos,
    Predicate::InCircle3dSos,
    Predicate::InSphere3dSos,
    Predicate::Side1Sos,
    Predicate::Side2Sos,
    Predicate::Side3Sos,
    Predicate::Side4Sos,
    Predicate::Side4_3d,
    Predicate::Side4_3dSos,
    Predicate::Side3_3dLifted,
    Predicate::Side3_3dLiftedSos,
    Predicate::InCircle3dLifted,
    Predicate::InCircle3dLiftedSos,
    Predicate::Orient2dLiftedSos,
    Predicate::Orient3dLifted,
    Predicate::Orient3dLiftedSos,
];

impl Predicate {
    fn name(self) -> &'static str {
        match self {
            Self::Orient2d => "orient_2d",
            Self::Orient3d => "orient_3d",
            Self::Det3d => "det_3d",
            Self::Det4d => "det_4d",
            Self::DetCompare4d => "det_compare_4d",
            Self::Dot3d => "dot_3d",
            Self::DotCompare3d => "dot_compare_3d",
            Self::Aligned3d => "aligned_3d",
            Self::Identical2d => "points_are_identical_2d",
            Self::Identical3d => "points_are_identical_3d",
            Self::Colinear3d => "points_are_colinear_3d",
            Self::InCircle2dSos => "in_circle_2d_sos",
            Self::InCircle3dSos => "in_circle_3d_sos",
            Self::InSphere3dSos => "in_sphere_3d_sos",
            Self::Side1Sos => "side1_sos",
            Self::Side2Sos => "side2_sos",
            Self::Side3Sos => "side3_sos",
            Self::Side4Sos => "side4_sos",
            Self::Side4_3d => "side4_3d",
            Self::Side4_3dSos => "side4_3d_sos",
            Self::Side3_3dLifted => "side3_3dlifted (non-SOS)",
            Self::Side3_3dLiftedSos => "side3_3dlifted_sos",
            Self::InCircle3dLifted => "in_circle_3dlifted (non-SOS)",
            Self::InCircle3dLiftedSos => "in_circle_3dlifted_sos",
            Self::Orient2dLiftedSos => "orient_2dlifted_sos",
            Self::Orient3dLifted => "orient_3dlifted",
            Self::Orient3dLiftedSos => "orient_3dlifted_sos",
        }
    }
}

#[derive(Clone, Debug)]
struct Case {
    predicate: Predicate,
    dim: usize,
    points: Vec<Vec<f64>>,
    keys: Vec<u8>,
    queries: Vec<Vec<f64>>,
    heights: Vec<f64>,
}

impl Case {
    fn write_to(&self, output: &mut impl Write) -> io::Result<()> {
        let header = [
            self.predicate as u8,
            self.dim as u8,
            self.points.len() as u8,
            self.queries.len() as u8,
            self.heights.len() as u8,
        ];
        output.write_all(&header)?;
        output.write_all(&self.keys)?;
        for value in self
            .points
            .iter()
            .chain(&self.queries)
            .flat_map(|point| point.iter())
            .chain(self.heights.iter())
        {
            output.write_all(&value.to_bits().to_le_bytes())?;
        }
        Ok(())
    }

    fn point<const D: usize>(&self, index: usize) -> [f64; D] {
        assert_eq!(self.dim, D);
        std::array::from_fn(|coordinate| self.points[index][coordinate])
    }

    fn query<const D: usize>(&self, index: usize) -> [f64; D] {
        assert_eq!(self.dim, D);
        std::array::from_fn(|coordinate| self.queries[index][coordinate])
    }

    fn sos_point<const D: usize>(&self, index: usize) -> SosPoint<D> {
        SosPoint::new(self.point(index), u64::from(self.keys[index]))
    }

    fn rust_result(&self) -> i8 {
        match self.predicate {
            Predicate::Orient2d => {
                gp::orient_2d(&self.point(0), &self.point(1), &self.point(2)).as_i8()
            }
            Predicate::Orient3d => gp::orient_3d(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
            )
            .as_i8(),
            Predicate::Det3d => gp::det_3d(&self.point(0), &self.point(1), &self.point(2)).as_i8(),
            Predicate::Det4d => gp::det_4d(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
            )
            .as_i8(),
            Predicate::DetCompare4d => gp::det_compare_4d(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
                &self.point(4),
            )
            .as_i8(),
            Predicate::Dot3d => gp::dot_3d(&self.point(0), &self.point(1), &self.point(2)).as_i8(),
            Predicate::DotCompare3d => {
                gp::dot_compare_3d(&self.point(0), &self.point(1), &self.point(2)).as_i8()
            }
            Predicate::Aligned3d => {
                gp::aligned_3d(&self.point(0), &self.point(1), &self.point(2)) as i8
            }
            Predicate::Identical2d => {
                gp::points_are_identical_2d(&self.point(0), &self.point(1)) as i8
            }
            Predicate::Identical3d => {
                gp::points_are_identical_3d(&self.point(0), &self.point(1)) as i8
            }
            Predicate::Colinear3d => {
                gp::points_are_colinear_3d(&self.point(0), &self.point(1), &self.point(2)) as i8
            }
            Predicate::InCircle2dSos => {
                let p = self.sos_points::<2>();
                gp::in_circle_2d_sos(&p[0], &p[1], &p[2], &p[3]).as_i8()
            }
            Predicate::InCircle3dSos => {
                let p = self.sos_points::<3>();
                gp::in_circle_3d_sos(&p[0], &p[1], &p[2], &p[3]).as_i8()
            }
            Predicate::InSphere3dSos => {
                let p = self.sos_points::<3>();
                gp::in_sphere_3d_sos(&p[0], &p[1], &p[2], &p[3], &p[4]).as_i8()
            }
            Predicate::Side1Sos => match self.dim {
                3 => self.side1::<3>(),
                4 => self.side1::<4>(),
                6 => self.side1::<6>(),
                7 => self.side1::<7>(),
                8 => self.side1::<8>(),
                _ => unreachable!(),
            },
            Predicate::Side2Sos => match self.dim {
                3 => self.side2::<3>(),
                4 => self.side2::<4>(),
                6 => self.side2::<6>(),
                7 => self.side2::<7>(),
                8 => self.side2::<8>(),
                _ => unreachable!(),
            },
            Predicate::Side3Sos => match self.dim {
                3 => self.side3::<3>(),
                4 => self.side3::<4>(),
                6 => self.side3::<6>(),
                7 => self.side3::<7>(),
                8 => self.side3::<8>(),
                _ => unreachable!(),
            },
            Predicate::Side4Sos => match self.dim {
                3 => self.side4::<3>(),
                4 => self.side4::<4>(),
                6 => self.side4::<6>(),
                7 => self.side4::<7>(),
                8 => self.side4::<8>(),
                _ => unreachable!(),
            },
            Predicate::Side4_3d => gp::side4_3d(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
                &self.point(4),
            )
            .as_i8(),
            Predicate::Side4_3dSos => {
                let p = self.sos_points::<3>();
                gp::side4_3d_sos(&p[0], &p[1], &p[2], &p[3], &p[4]).as_i8()
            }
            // The public non-SOS circle wrapper delegates to side3_3dlifted
            // with q0..q2 = p0..p2 and flips its sign.
            Predicate::Side3_3dLifted => (-gp::in_circle_3dlifted(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
                self.heights4(),
            ))
            .as_i8(),
            Predicate::Side3_3dLiftedSos => gp::side3_3dlifted_sos(
                &self.sos_point(0),
                &self.sos_point(1),
                &self.sos_point(2),
                &self.sos_point(3),
                self.heights4(),
                &self.query(0),
                &self.query(1),
                &self.query(2),
            )
            .as_i8(),
            Predicate::InCircle3dLifted => gp::in_circle_3dlifted(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
                self.heights4(),
            )
            .as_i8(),
            Predicate::InCircle3dLiftedSos => {
                let p = self.sos_points::<3>();
                gp::in_circle_3dlifted_sos(&p[0], &p[1], &p[2], &p[3], self.heights4()).as_i8()
            }
            Predicate::Orient2dLiftedSos => {
                let p = self.sos_points::<2>();
                gp::orient_2dlifted_sos(&p[0], &p[1], &p[2], &p[3], self.heights4()).as_i8()
            }
            Predicate::Orient3dLifted => gp::orient_3dlifted(
                &self.point(0),
                &self.point(1),
                &self.point(2),
                &self.point(3),
                &self.point(4),
                self.heights5(),
            )
            .as_i8(),
            Predicate::Orient3dLiftedSos => {
                let p = self.sos_points::<3>();
                gp::orient_3dlifted_sos(&p[0], &p[1], &p[2], &p[3], &p[4], self.heights5()).as_i8()
            }
        }
    }

    fn sos_points<const D: usize>(&self) -> Vec<SosPoint<D>> {
        (0..self.points.len()).map(|i| self.sos_point(i)).collect()
    }

    fn heights4(&self) -> [f64; 4] {
        std::array::from_fn(|i| self.heights[i])
    }

    fn heights5(&self) -> [f64; 5] {
        std::array::from_fn(|i| self.heights[i])
    }

    fn side1<const D: usize>(&self) -> i8 {
        let p = self.sos_points::<D>();
        gp::side1_sos(&p[0], &p[1], &self.query(0)).as_i8()
    }

    fn side2<const D: usize>(&self) -> i8 {
        let p = self.sos_points::<D>();
        gp::side2_sos(&p[0], &p[1], &p[2], &self.query(0), &self.query(1)).as_i8()
    }

    fn side3<const D: usize>(&self) -> i8 {
        let p = self.sos_points::<D>();
        gp::side3_sos(
            &p[0],
            &p[1],
            &p[2],
            &p[3],
            &self.query(0),
            &self.query(1),
            &self.query(2),
        )
        .as_i8()
    }

    fn side4<const D: usize>(&self) -> i8 {
        let p = self.sos_points::<D>();
        gp::side4_sos(
            &p[0],
            &p[1],
            &p[2],
            &p[3],
            &p[4],
            &self.query(0),
            &self.query(1),
            &self.query(2),
            &self.query(3),
        )
        .as_i8()
    }

    fn description(&self) -> String {
        let mut result = format!(
            "{} (dimension {})\nkeys: {:?}\n",
            self.predicate.name(),
            self.dim,
            self.keys
        );
        for (label, values) in [("p", &self.points), ("q", &self.queries)] {
            for (index, point) in values.iter().enumerate() {
                let _ = writeln!(result, "{label}{index}: {}", format_values(point));
            }
        }
        if !self.heights.is_empty() {
            let _ = writeln!(result, "heights: {}", format_values(&self.heights));
        }
        result
    }
}

fn format_values(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:?} (0x{:016x})", value.to_bits()))
        .collect::<Vec<_>>()
        .join(", ")
}

#[derive(Clone)]
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    fn shuffle(&mut self, values: &mut [u8]) {
        for i in (1..values.len()).rev() {
            values.swap(i, self.next() as usize % (i + 1));
        }
    }
}

fn next_up(value: f64) -> f64 {
    if value == f64::INFINITY {
        value
    } else if value == 0.0 {
        f64::from_bits(1)
    } else if value > 0.0 {
        f64::from_bits(value.to_bits() + 1)
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

fn power_of_two(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

fn random_coord(rng: &mut SplitMix64, scale: f64) -> f64 {
    let integer = (rng.next() % 65_537) as i64 - 32_768;
    integer as f64 * scale / 16.0
}

fn random_points(rng: &mut SplitMix64, count: usize, dim: usize, scale: f64) -> Vec<Vec<f64>> {
    (0..count)
        .map(|_| (0..dim).map(|_| random_coord(rng, scale)).collect())
        .collect()
}

fn random_keys(rng: &mut SplitMix64, count: usize) -> Vec<u8> {
    let mut keys: Vec<u8> = (0..count as u8).collect();
    rng.shuffle(&mut keys);
    keys
}

fn plain_case(
    predicate: Predicate,
    dim: usize,
    count: usize,
    scenario: usize,
    rng: &mut SplitMix64,
    scale: f64,
) -> Case {
    let mut points = random_points(rng, count, dim, scale);
    match predicate {
        Predicate::Orient2d if scenario == 0 => points[2] = points[1].clone(),
        Predicate::Orient2d if scenario == 1 => {
            points[2] = points[1].clone();
            points[2][1] = next_up(points[2][1]);
        }
        Predicate::Orient3d if scenario == 0 => points[3] = points[2].clone(),
        Predicate::Orient3d if scenario == 1 => {
            points[3] = points[2].clone();
            points[3][2] = next_up(points[3][2]);
        }
        Predicate::Det3d if scenario == 0 => points[2] = points[1].clone(),
        Predicate::Det3d if scenario == 1 => {
            points[2] = points[1].clone();
            points[2][2] = next_up(points[2][2]);
        }
        Predicate::Det4d if scenario == 0 => points[3] = points[2].clone(),
        Predicate::Det4d if scenario == 1 => {
            points[3] = points[2].clone();
            points[3][3] = next_up(points[3][3]);
        }
        Predicate::DetCompare4d if scenario == 0 => points[4] = points[3].clone(),
        Predicate::DetCompare4d if scenario == 1 => {
            points[4] = points[3].clone();
            points[4][3] = next_up(points[4][3]);
        }
        Predicate::Dot3d => {
            // The PSM's public dot_3d accidentally uses det_3d_filter first.
            // Keeping p0 at zero makes that filter uncertain and exercises its
            // correct exact dot-product fallback instead.
            points = vec![vec![0.0; 3], vec![scale, 0.0, 0.0], vec![0.0; 3]];
            points[2][0] = match scenario {
                0 => 0.0,
                1 => scale * power_of_two(-40),
                2 => scale,
                _ => -scale,
            };
        }
        Predicate::DotCompare3d if scenario == 0 => points[2] = points[1].clone(),
        Predicate::DotCompare3d if scenario == 1 => {
            points[2] = points[1].clone();
            points[2][0] = next_up(points[2][0]);
        }
        Predicate::Aligned3d | Predicate::Colinear3d if scenario < 2 => {
            points = vec![
                vec![0.0; 3],
                vec![scale, 2.0 * scale, 3.0 * scale],
                vec![2.0 * scale, 4.0 * scale, 6.0 * scale],
            ];
            if scenario == 1 {
                points[2][2] = next_up(points[2][2]);
            }
        }
        Predicate::Identical2d | Predicate::Identical3d if scenario < 2 => {
            points[1] = points[0].clone();
            if scenario == 1 {
                points[1][0] = next_up(points[1][0]);
            }
        }
        _ => {}
    }
    Case {
        predicate,
        dim,
        keys: random_keys(rng, count),
        points,
        queries: Vec::new(),
        heights: Vec::new(),
    }
}

fn circle_case(
    predicate: Predicate,
    dim: usize,
    scenario: usize,
    rng: &mut SplitMix64,
    scale: f64,
) -> Case {
    let mut points = vec![vec![0.0; dim]; 4];
    points[1][0] = scale;
    points[2][1] = scale;
    match scenario {
        0 => {
            points[3][0] = scale;
            points[3][1] = scale;
        }
        1 => {
            points[3][0] = next_up(scale);
            points[3][1] = scale;
        }
        2 => {
            points[3][0] = 0.25 * scale;
            points[3][1] = 0.25 * scale;
        }
        _ => {
            points[3][0] = 2.0 * scale;
            points[3][1] = 2.0 * scale;
        }
    }
    Case {
        predicate,
        dim,
        keys: random_keys(rng, 4),
        points,
        queries: Vec::new(),
        heights: Vec::new(),
    }
}

fn sphere_case(predicate: Predicate, scenario: usize, rng: &mut SplitMix64, scale: f64) -> Case {
    let mut points = vec![vec![0.0; 3]; 5];
    points[1][0] = scale;
    points[2][1] = scale;
    points[3][2] = scale;
    points[4] = match scenario {
        0 => vec![scale; 3],
        1 => vec![next_up(scale), scale, scale],
        2 => vec![0.25 * scale; 3],
        _ => vec![2.0 * scale; 3],
    };
    Case {
        predicate,
        dim: 3,
        keys: random_keys(rng, 5),
        points,
        queries: Vec::new(),
        heights: Vec::new(),
    }
}

fn side_case(
    predicate: Predicate,
    dim: usize,
    scenario: usize,
    rng: &mut SplitMix64,
    scale: f64,
) -> Case {
    let side = match predicate {
        Predicate::Side1Sos => 1,
        Predicate::Side2Sos => 2,
        Predicate::Side3Sos => 3,
        Predicate::Side4Sos => 4,
        _ => unreachable!(),
    };
    let mut points = vec![vec![0.0; dim]; side + 1];
    for i in 1..=side.min(3) {
        points[i][i - 1] = scale;
    }
    if side == 4 {
        points[4][0] = 0.25 * scale;
        points[4][1] = 0.25 * scale;
        points[4][2] = 0.25 * scale;
    }
    if scenario == 0 {
        if side == 1 {
            // q0 at the midpoint makes p0 and p1 exactly equidistant.
        } else {
            points[side] = points[0].clone();
        }
    } else if scenario == 1 && side > 1 {
        points[side] = points[0].clone();
        points[side][0] = scale * power_of_two(-40);
    } else if side > 1 && side < 4 {
        for coordinate in &mut points[side] {
            *coordinate = random_coord(rng, scale) / 8.0;
        }
    }

    let mut queries: Vec<Vec<f64>> = (0..side).map(|i| points[i].clone()).collect();
    if side == 1 && scenario == 0 {
        queries[0][0] = 0.5 * scale;
    }
    Case {
        predicate,
        dim,
        keys: random_keys(rng, side + 1),
        points,
        queries,
        heights: Vec::new(),
    }
}

fn lifted_case(predicate: Predicate, scenario: usize, rng: &mut SplitMix64, scale: f64) -> Case {
    let orient_2d = matches!(predicate, Predicate::Orient2dLiftedSos);
    let dim = if orient_2d { 2 } else { 3 };
    let count = dim + 2;
    let mut points = vec![vec![0.0; dim]; count];
    for i in 1..=dim {
        points[i][i - 1] = scale;
    }
    points[count - 1].fill(0.25 * scale);

    let mut heights = match scenario {
        0 => vec![scale; count],
        1 => {
            let mut values = vec![scale; count];
            values[count - 1] = next_up(scale);
            values
        }
        _ => (0..count)
            .map(|_| random_coord(rng, scale) / 16.0)
            .collect(),
    };
    if matches!(
        predicate,
        Predicate::InCircle3dLifted | Predicate::InCircle3dLiftedSos
    ) {
        heights = points
            .iter()
            .map(|point| point.iter().map(|coordinate| coordinate * coordinate).sum())
            .collect();
        if scenario == 1 {
            heights[count - 1] = next_up(heights[count - 1]);
        }
    }
    let queries = if matches!(
        predicate,
        Predicate::Side3_3dLifted | Predicate::Side3_3dLiftedSos
    ) {
        points[..3].to_vec()
    } else {
        Vec::new()
    };
    Case {
        predicate,
        dim,
        keys: random_keys(rng, count),
        points,
        queries,
        heights,
    }
}

fn make_case(seed: u64, index: u64) -> Case {
    let predicate = PREDICATES[index as usize % PREDICATES.len()];
    let round = index as usize / PREDICATES.len();
    let scenario = round % 4;
    let mut rng = SplitMix64::new(seed ^ index.wrapping_mul(0xd134_2543_de82_ef95));
    let scale = power_of_two([-20, -5, 0, 10, 20][round % 5]);
    match predicate {
        Predicate::Orient2d => plain_case(predicate, 2, 3, scenario, &mut rng, scale),
        Predicate::Orient3d => plain_case(predicate, 3, 4, scenario, &mut rng, scale),
        Predicate::Det3d => plain_case(predicate, 3, 3, scenario, &mut rng, scale),
        Predicate::Det4d => plain_case(predicate, 4, 4, scenario, &mut rng, scale),
        Predicate::DetCompare4d => plain_case(predicate, 4, 5, scenario, &mut rng, scale),
        Predicate::Dot3d
        | Predicate::DotCompare3d
        | Predicate::Aligned3d
        | Predicate::Colinear3d => plain_case(predicate, 3, 3, scenario, &mut rng, scale),
        Predicate::Identical2d => plain_case(predicate, 2, 2, scenario, &mut rng, scale),
        Predicate::Identical3d => plain_case(predicate, 3, 2, scenario, &mut rng, scale),
        Predicate::InCircle2dSos => circle_case(predicate, 2, scenario, &mut rng, scale),
        Predicate::InCircle3dSos => circle_case(predicate, 3, scenario, &mut rng, scale),
        Predicate::InSphere3dSos | Predicate::Side4_3d | Predicate::Side4_3dSos => {
            sphere_case(predicate, scenario, &mut rng, scale)
        }
        Predicate::Side1Sos | Predicate::Side2Sos | Predicate::Side3Sos | Predicate::Side4Sos => {
            let dim = SIDE_DIMENSIONS[round % SIDE_DIMENSIONS.len()];
            side_case(predicate, dim, scenario, &mut rng, scale)
        }
        Predicate::Side3_3dLifted
        | Predicate::Side3_3dLiftedSos
        | Predicate::InCircle3dLifted
        | Predicate::InCircle3dLiftedSos
        | Predicate::Orient2dLiftedSos
        | Predicate::Orient3dLifted
        | Predicate::Orient3dLiftedSos => lifted_case(predicate, scenario, &mut rng, scale),
    }
}

fn parse_u64_env(name: &str, default: u64) -> u64 {
    let Ok(raw) = env::var(name) else {
        return default;
    };
    let parsed = raw
        .strip_prefix("0x")
        .or_else(|| raw.strip_prefix("0X"))
        .map_or_else(|| raw.parse(), |hex| u64::from_str_radix(hex, 16));
    parsed.unwrap_or_else(|_| panic!("{name} must be a decimal or 0x-prefixed u64, got {raw:?}"))
}

fn target_directory(root: &Path) -> PathBuf {
    match env::var_os("CARGO_TARGET_DIR") {
        Some(path) if Path::new(&path).is_absolute() => PathBuf::from(path),
        Some(path) => root.join(path),
        None => root.join("target"),
    }
}

fn compile_oracle(root: &Path) -> PathBuf {
    let directory = target_directory(root).join("cpp-differential");
    std::fs::create_dir_all(&directory).expect("could not create the C++ oracle target directory");
    let executable = directory.join(format!("oracle-{}", std::process::id()));
    let compiler = env::var_os("CXX").unwrap_or_else(|| OsString::from("c++"));
    let output = Command::new(&compiler)
        .args(["-std=c++17", "-O2", "-fno-fast-math", "-ffp-contract=off"])
        .arg(root.join("tests/cpp_oracle.cpp"))
        .arg(root.join("include/geogram_predicates_psm/Predicates_psm.cpp"))
        .arg("-I")
        .arg(root)
        .arg("-o")
        .arg(&executable)
        .output()
        .unwrap_or_else(|error| panic!("failed to run C++ compiler {compiler:?}: {error}"));
    assert!(
        output.status.success(),
        "failed to compile the Geogram C++ oracle with {compiler:?}:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    executable
}

#[test]
#[ignore = "requires a C++ compiler and is intentionally opt-in"]
fn cpp_and_rust_predicates_agree() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let count = parse_u64_env("GEOGRAM_DIFF_CASES", DEFAULT_CASES);
    let seed = parse_u64_env("GEOGRAM_DIFF_SEED", DEFAULT_SEED);
    assert!(count > 0, "GEOGRAM_DIFF_CASES must be greater than zero");

    eprintln!("compiling C++ oracle; cases={count}, seed=0x{seed:016x}");
    let executable = compile_oracle(&root);
    let mut child = Command::new(&executable)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .expect("failed to start the C++ oracle");
    let stdin = child.stdin.take().expect("oracle stdin was not piped");
    let stdout = child.stdout.take().expect("oracle stdout was not piped");

    let writer = thread::spawn(move || -> io::Result<()> {
        let mut output = BufWriter::with_capacity(1024 * 1024, stdin);
        for index in 0..count {
            make_case(seed, index).write_to(&mut output)?;
        }
        output.flush()
    });

    let mut output = BufReader::with_capacity(1024 * 1024, stdout);
    let mut first_mismatch = None;
    for index in 0..count {
        let case = make_case(seed, index);
        let expected = case.rust_result();
        let mut byte = [0_u8; 1];
        if let Err(error) = output.read_exact(&mut byte) {
            let status = child.wait().expect("failed to wait for the C++ oracle");
            panic!("C++ oracle stopped at case {index} with {status}: {error}");
        }
        let actual = byte[0] as i8;
        if actual != expected && first_mismatch.is_none() {
            first_mismatch = Some((index, expected, actual, case));
        }
    }

    let writer_result = writer.join().expect("differential input writer panicked");
    writer_result.expect("failed to stream cases to the C++ oracle");
    let status = child.wait().expect("failed to wait for the C++ oracle");
    let _ = std::fs::remove_file(&executable);
    assert!(status.success(), "C++ oracle exited with {status}");

    if let Some((index, rust, cpp, case)) = first_mismatch {
        panic!(
            "differential mismatch at case {index}\nseed: 0x{seed:016x}\nRust: {rust}\nC++: {cpp}\n{}",
            case.description()
        );
    }
    eprintln!(
        "matched all {count} cases across {} predicate variants",
        PREDICATES.len()
    );
}
