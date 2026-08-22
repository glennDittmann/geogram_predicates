//! Multi-precision expansion arithmetic (Shewchuk-style).
//! Exact predicates are implemented by computing in this number type and taking the sign.
#![allow(clippy::too_many_arguments)] // Matrix entries mirror Geogram's fixed determinant helpers.

use crate::sign::{geo_sgn, Sign};
use std::vec::Vec;

/// Shewchuk's splitter for IEEE-754 binary64 (`2^27 + 1`).
const EXPANSION_SPLITTER: f64 = 134_217_729.0;

/// Expansion: non-overlapping sum of f64 components, least significant first.
/// Same semantics as Geogram: component 0 is least significant, last is most significant.
///
/// <br>
///
/// More theoretical:
///
/// An `Expansion`\[1\]\[2\] is defined as _x = x_n + ... + x_1_
///
/// Each _x_i_ is a _component_ of _x_, represented by a p-bit significand floating-point value.
///
/// Two properties are assumed on `Expansions`
///
/// 1. Its components are ordered by magnitude, i.e. _x_n > ... > ... x_1_
///
/// 2. Its components are non-overlapping, that means the least significant non-zero bit of some component _a_ is more
/// significant than the most significant non-zero bit of some component _b_ (or vice versa) for all components of _x_.
/// `0` does not overlap any number.
///
/// ### Example Overlapping/non-Overlapping
/// `1100` and `-10.1` (1) are non-overlapping, whereas `101` and `10` (2) are.
///
/// One can see that easily by writing the numbers in rows above each other and fill with leading zeroes:
///
/// Ex. (1) `+1100.0` | Ex. (2) `101`
///
/// Ex. (2) `-0010.1` | Ex. (2) `010`
///
/// ### Notes
/// A number can be represented by multiple non-overlapping expansion, e.g. `x = 1100 + –10.1 = 1001 + 0.1 = 1000 + 1 + 0.1`.
///
/// A non-overlapping `Expansion` is favorable, because one can easily compute its sing (sign of _x_n_) or get a crude approximation (_x_n_).
#[derive(Clone, Debug)]
pub struct Expansion {
    x: Vec<f64>,
}

#[allow(dead_code)]
impl Expansion {
    #[inline]
    pub fn len(&self) -> usize {
        self.x.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.x.capacity()
    }

    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.x
    }

    #[inline]
    pub fn sign(&self) -> Sign {
        if self.x.is_empty() {
            Sign::Zero
        } else {
            geo_sgn(self.x[self.x.len() - 1])
        }
    }

    #[inline]
    pub fn negate(&mut self) {
        for v in &mut self.x {
            *v = -*v;
        }
    }

    /// Create expansion with given capacity (length 0).
    pub fn with_capacity(capa: usize) -> Self {
        Self {
            x: Vec::with_capacity(capa),
        }
    }

    /// Assign from a single double.
    pub fn assign(&mut self, a: f64) {
        self.x.clear();
        if a != 0.0 {
            self.x.push(a);
        }
    }

    /// Assign sum of two doubles.
    pub fn assign_sum_2(&mut self, a: f64, b: f64) {
        self.x.clear();
        let (hi, lo) = two_sum(a, b);
        if lo != 0.0 {
            self.x.push(lo);
        }
        if hi != 0.0 || self.x.is_empty() {
            self.x.push(hi);
        }
    }

    /// Assign difference a - b for two doubles.
    pub fn assign_diff_2(&mut self, a: f64, b: f64) {
        self.x.clear();
        let (hi, lo) = two_diff(a, b);
        if lo != 0.0 {
            self.x.push(lo);
        }
        if hi != 0.0 || self.x.is_empty() {
            self.x.push(hi);
        }
    }

    /// Assign product of two doubles.
    pub fn assign_product_2(&mut self, a: f64, b: f64) {
        self.x.clear();
        let (hi, lo) = two_product(a, b);
        if lo != 0.0 {
            self.x.push(lo);
        }
        if hi != 0.0 || self.x.is_empty() {
            self.x.push(hi);
        }
    }

    /// Assign sum of expansion and double: self := e + b.
    pub fn assign_sum_exp_double(&mut self, e: &Expansion, b: f64) {
        grow_expansion_zeroelim(e, b, self);
    }

    /// Assign difference of expansion and double: self := e - b.
    pub fn assign_diff_exp_double(&mut self, e: &Expansion, b: f64) {
        grow_expansion_zeroelim(e, -b, self);
    }

    /// Assign product of expansion and double.
    pub fn assign_product_exp_double(&mut self, e: &Expansion, b: f64) {
        scale_expansion_zeroelim(e, b, self);
    }

    /// Assign sum of two expansions.
    pub fn assign_sum_exp_exp(&mut self, a: &Expansion, b: &Expansion) {
        fast_expansion_sum_zeroelim(a, b, self);
    }

    /// Assign difference of two expansions: self := a - b.
    pub fn assign_diff_exp_exp(&mut self, a: &Expansion, b: &Expansion) {
        fast_expansion_diff_zeroelim(a, b, self);
    }

    /// Assign product of two expansions.
    pub fn assign_product_exp_exp(&mut self, a: &Expansion, b: &Expansion) {
        if a.is_empty() || b.is_empty() {
            self.x.clear();
            return;
        }
        if a.len() == 1 && b.len() == 1 {
            self.assign_product_2(a.x[0], b.x[0]);
            return;
        }
        if a.len() == 1 {
            self.assign_product_exp_double(b, a.x[0]);
            return;
        }
        if b.len() == 1 {
            self.assign_product_exp_double(a, b.x[0]);
            return;
        }
        // General case: scale and sum.
        let mut acc = Expansion::with_capacity(a.len() * 2 + 4);
        acc.assign_product_exp_double(b, a.x[0]);
        for i in 1..a.len() {
            let mut p = Expansion::with_capacity(b.len() * 2 + 4);
            p.assign_product_exp_double(b, a.x[i]);
            let mut next = Expansion::with_capacity(acc.len() + p.len());
            next.assign_sum_exp_exp(&acc, &p);
            acc = next;
        }
        self.x.clear();
        self.x.extend(acc.x.iter().cloned());
    }

    /// Assign det2x2(a11, a12, a21, a22) = a11*a22 - a12*a21.
    pub fn assign_det2x2(
        &mut self,
        a11: &Expansion,
        a12: &Expansion,
        a21: &Expansion,
        a22: &Expansion,
    ) {
        let a11a22 = expansion_product(a11, a22);
        let a12a21 = expansion_product(a12, a21);
        self.assign_diff_exp_exp(&a11a22, &a12a21);
    }

    /// Assign det3x3 (rows a1*, a2*, a3*).
    pub fn assign_det3x3(
        &mut self,
        a11: &Expansion,
        a12: &Expansion,
        a13: &Expansion,
        a21: &Expansion,
        a22: &Expansion,
        a23: &Expansion,
        a31: &Expansion,
        a32: &Expansion,
        a33: &Expansion,
    ) {
        let c11 = expansion_det2x2(a22, a23, a32, a33);
        let c12 = expansion_det2x2(a23, a21, a33, a31);
        let c13 = expansion_det2x2(a21, a22, a31, a32);
        let a11c11 = expansion_product(a11, &c11);
        let a12c12 = expansion_product(a12, &c12);
        let a13c13 = expansion_product(a13, &c13);
        let s1 = expansion_sum(&a11c11, &a12c12);
        self.assign_sum_exp_exp(&s1, &a13c13);
    }
}

impl std::ops::Index<usize> for Expansion {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        &self.x[i]
    }
}

// --- Primitives (in anonymous namespace style) ---

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bvirt = x - a;
    let avirt = x - bvirt;
    let bround = b - bvirt;
    let around = a - avirt;
    let y = around + bround;
    (x, y)
}

#[inline]
fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bvirt = x - a;
    let y = b - bvirt;
    (x, y)
}

#[inline]
fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    let bvirt = a - x;
    let avirt = x + bvirt;
    let bround = bvirt - b;
    let around = a - avirt;
    let y = around + bround;
    (x, y)
}

#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = EXPANSION_SPLITTER * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;
    (ahi, alo)
}

#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    let y = (alo * blo) - err3;
    (x, y)
}

#[inline]
fn two_product_presplit(a: f64, b: f64, bhi: f64, blo: f64) -> (f64, f64) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    let y = (alo * blo) - err3;
    (x, y)
}

#[inline]
#[allow(dead_code)]
fn square(a: f64) -> (f64, f64) {
    let x = a * a;
    let (ahi, alo) = split(a);
    let err1 = x - (ahi * ahi);
    let err3 = err1 - ((ahi + ahi) * alo);
    let y = (alo * alo) - err3;
    (x, y)
}

#[allow(dead_code)]
fn grow_expansion_zeroelim(e: &Expansion, b: f64, h: &mut Expansion) {
    h.x.clear();
    let mut q = b;
    for &enow in &e.x {
        let (qnew, hh) = two_sum(q, enow);
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    if q != 0.0 || h.x.is_empty() {
        h.x.push(q);
    }
}

fn scale_expansion_zeroelim(e: &Expansion, b: f64, h: &mut Expansion) {
    h.x.clear();
    if e.x.is_empty() {
        return;
    }
    let (bhi, blo) = split(b);
    let (mut q, hh) = two_product_presplit(e.x[0], b, bhi, blo);
    if hh != 0.0 {
        h.x.push(hh);
    }
    for &enow in e.x.iter().skip(1) {
        let (product1, product0) = two_product_presplit(enow, b, bhi, blo);
        let (sum, hh) = two_sum(q, product0);
        if hh != 0.0 {
            h.x.push(hh);
        }
        let (qnew, hh) = fast_two_sum(product1, sum);
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    if q != 0.0 || h.x.is_empty() {
        h.x.push(q);
    }
}

fn fast_expansion_sum_zeroelim(e: &Expansion, f: &Expansion, h: &mut Expansion) {
    h.x.clear();
    let elen = e.x.len();
    let flen = f.x.len();
    if elen == 0 {
        h.x.extend(f.x.iter().cloned());
        return;
    }
    if flen == 0 {
        h.x.extend(e.x.iter().cloned());
        return;
    }
    let mut eindex = 0usize;
    let mut findex = 0usize;
    let mut enow = e.x[0];
    let mut fnow = f.x[0];
    let mut q = if (fnow > enow) == (fnow > -enow) {
        eindex += 1;
        enow = if eindex < elen { e.x[eindex] } else { 0.0 };
        e.x[0]
    } else {
        findex += 1;
        fnow = if findex < flen { f.x[findex] } else { 0.0 };
        f.x[0]
    };
    while eindex < elen && findex < flen {
        let (qnew, hh) = if (fnow > enow) == (fnow > -enow) {
            let (qn, hh) = two_sum(q, enow);
            eindex += 1;
            enow = if eindex < elen { e.x[eindex] } else { 0.0 };
            (qn, hh)
        } else {
            let (qn, hh) = two_sum(q, fnow);
            findex += 1;
            fnow = if findex < flen { f.x[findex] } else { 0.0 };
            (qn, hh)
        };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    while eindex < elen {
        let (qnew, hh) = two_sum(q, enow);
        eindex += 1;
        enow = if eindex < elen { e.x[eindex] } else { 0.0 };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    while findex < flen {
        let (qnew, hh) = two_sum(q, fnow);
        findex += 1;
        fnow = if findex < flen { f.x[findex] } else { 0.0 };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    if q != 0.0 || h.x.is_empty() {
        h.x.push(q);
    }
}

fn fast_expansion_diff_zeroelim(e: &Expansion, f: &Expansion, h: &mut Expansion) {
    h.x.clear();
    let elen = e.x.len();
    let flen = f.x.len();
    if elen == 0 {
        for &fv in &f.x {
            h.x.push(-fv);
        }
        return;
    }
    if flen == 0 {
        h.x.extend(e.x.iter().cloned());
        return;
    }
    let mut eindex = 0usize;
    let mut findex = 0usize;
    let mut enow = e.x[0];
    let mut fnow = -f.x[0];
    let mut q = if (fnow > enow) == (fnow > -enow) {
        eindex += 1;
        enow = if eindex < elen { e.x[eindex] } else { 0.0 };
        e.x[0]
    } else {
        findex += 1;
        fnow = if findex < flen { -f.x[findex] } else { 0.0 };
        -f.x[0]
    };
    while eindex < elen && findex < flen {
        let (qnew, hh) = if (fnow > enow) == (fnow > -enow) {
            let (qn, hh) = two_sum(q, enow);
            eindex += 1;
            enow = if eindex < elen { e.x[eindex] } else { 0.0 };
            (qn, hh)
        } else {
            let (qn, hh) = two_sum(q, fnow);
            findex += 1;
            fnow = if findex < flen { -f.x[findex] } else { 0.0 };
            (qn, hh)
        };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    while eindex < elen {
        let (qnew, hh) = two_sum(q, enow);
        eindex += 1;
        enow = if eindex < elen { e.x[eindex] } else { 0.0 };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    while findex < flen {
        let (qnew, hh) = two_sum(q, fnow);
        findex += 1;
        fnow = if findex < flen { -f.x[findex] } else { 0.0 };
        q = qnew;
        if hh != 0.0 {
            h.x.push(hh);
        }
    }
    if q != 0.0 || h.x.is_empty() {
        h.x.push(q);
    }
}

// --- Helpers that allocate and return Expansion (replacing C++ macros) ---

pub fn expansion_create(a: f64) -> Expansion {
    let mut e = Expansion::with_capacity(2);
    e.assign(a);
    e
}

pub fn expansion_diff_2(a: f64, b: f64) -> Expansion {
    let mut e = Expansion::with_capacity(2);
    e.assign_diff_2(a, b);
    e
}

pub fn expansion_sum(a: &Expansion, b: &Expansion) -> Expansion {
    let capa = a.len() + b.len();
    let mut h = Expansion::with_capacity(capa);
    h.assign_sum_exp_exp(a, b);
    h
}

pub fn expansion_sum3(a: &Expansion, b: &Expansion, c: &Expansion) -> Expansion {
    let ab = expansion_sum(a, b);
    expansion_sum(&ab, c)
}

pub fn expansion_diff(a: &Expansion, b: &Expansion) -> Expansion {
    let capa = a.len() + b.len();
    let mut h = Expansion::with_capacity(capa);
    h.assign_diff_exp_exp(a, b);
    h
}

pub fn expansion_product(a: &Expansion, b: &Expansion) -> Expansion {
    let capa = if a.is_empty() || b.is_empty() {
        0
    } else {
        a.len() * b.len() * 2
    };
    let mut h = Expansion::with_capacity(capa.max(4));
    h.assign_product_exp_exp(a, b);
    h
}

/// Exact product of an expansion and one binary64 value.
pub fn expansion_scale(a: &Expansion, b: f64) -> Expansion {
    let mut result = Expansion::with_capacity(a.len().saturating_mul(2).max(2));
    result.assign_product_exp_double(a, b);
    result
}

/// Exact squared distance between two binary64 points.
pub fn expansion_sq_dist<const D: usize>(a: &[f64; D], b: &[f64; D]) -> Expansion {
    let mut result = expansion_create(0.0);
    for i in 0..D {
        let d = expansion_diff_2(a[i], b[i]);
        let square = expansion_product(&d, &d);
        result = expansion_sum(&result, &square);
    }
    result
}

/// Exact dot product `(p - origin) . (q - origin)`.
pub fn expansion_dot_at<const D: usize>(
    p: &[f64; D],
    q: &[f64; D],
    origin: &[f64; D],
) -> Expansion {
    let mut result = expansion_create(0.0);
    for i in 0..D {
        let u = expansion_diff_2(p[i], origin[i]);
        let v = expansion_diff_2(q[i], origin[i]);
        result = expansion_sum(&result, &expansion_product(&u, &v));
    }
    result
}

/// Exact determinant of a small square expansion matrix.
///
/// Predicate matrices in this crate are at most 5x5, so recursive cofactor
/// expansion keeps this helper compact without making exact fallbacks hot.
pub fn expansion_determinant(matrix: &[Vec<Expansion>]) -> Expansion {
    let n = matrix.len();
    assert!(matrix.iter().all(|row| row.len() == n));
    if n == 0 {
        return expansion_create(1.0);
    }
    if n == 1 {
        return matrix[0][0].clone();
    }

    let mut result = expansion_create(0.0);
    for column in 0..n {
        if matrix[0][column].sign() == Sign::Zero {
            continue;
        }
        let mut minor = Vec::with_capacity(n - 1);
        for source_row in matrix.iter().skip(1) {
            let mut row = Vec::with_capacity(n - 1);
            for (source_column, value) in source_row.iter().enumerate() {
                if source_column != column {
                    row.push(value.clone());
                }
            }
            minor.push(row);
        }
        let mut term = expansion_product(&matrix[0][column], &expansion_determinant(&minor));
        if column % 2 != 0 {
            term.negate();
        }
        result = expansion_sum(&result, &term);
    }
    result
}

/// Exact cofactor at `(row, column)` in a square expansion matrix.
pub fn expansion_cofactor(matrix: &[Vec<Expansion>], row: usize, column: usize) -> Expansion {
    let n = matrix.len();
    assert!(row < n && column < n && matrix.iter().all(|values| values.len() == n));
    let mut minor = Vec::with_capacity(n.saturating_sub(1));
    for (source_row, values) in matrix.iter().enumerate() {
        if source_row == row {
            continue;
        }
        let mut minor_row = Vec::with_capacity(n.saturating_sub(1));
        for (source_column, value) in values.iter().enumerate() {
            if source_column != column {
                minor_row.push(value.clone());
            }
        }
        minor.push(minor_row);
    }
    let mut result = expansion_determinant(&minor);
    if (row + column) % 2 != 0 {
        result.negate();
    }
    result
}

pub fn expansion_det2x2(
    a11: &Expansion,
    a12: &Expansion,
    a21: &Expansion,
    a22: &Expansion,
) -> Expansion {
    let capa = a11.len() * a22.len() * 2 + a21.len() * a12.len() * 2;
    let mut h = Expansion::with_capacity(capa.max(4));
    h.assign_det2x2(a11, a12, a21, a22);
    h
}

pub fn expansion_det3x3(
    a11: &Expansion,
    a12: &Expansion,
    a13: &Expansion,
    a21: &Expansion,
    a22: &Expansion,
    a23: &Expansion,
    a31: &Expansion,
    a32: &Expansion,
    a33: &Expansion,
) -> Expansion {
    let c11_capa = a22.len() * a33.len() * 2 + a32.len() * a23.len() * 2;
    let c12_capa = a21.len() * a33.len() * 2 + a31.len() * a23.len() * 2;
    let c13_capa = a21.len() * a22.len() * 2 + a31.len() * a32.len() * 2;
    let capa = 2 * (a11.len() * c11_capa + a12.len() * c12_capa + a13.len() * c13_capa);
    let mut h = Expansion::with_capacity(capa.max(24));
    h.assign_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
    h
}

// --- Sign of determinant (used by exact predicates) ---

#[allow(dead_code)]
pub fn sign_of_expansion_det2x2(
    a00: &Expansion,
    a01: &Expansion,
    a10: &Expansion,
    a11: &Expansion,
) -> Sign {
    expansion_det2x2(a00, a01, a10, a11).sign()
}

pub fn sign_of_expansion_det3x3(
    a00: &Expansion,
    a01: &Expansion,
    a02: &Expansion,
    a10: &Expansion,
    a11: &Expansion,
    a12: &Expansion,
    a20: &Expansion,
    a21: &Expansion,
    a22: &Expansion,
) -> Sign {
    let m01 = expansion_det2x2(a00, a10, a01, a11);
    let m02 = expansion_det2x2(a00, a20, a01, a21);
    let m12 = expansion_det2x2(a10, a20, a11, a21);
    let z1 = expansion_product(&m01, a22);
    let mut z2 = expansion_product(&m02, a12);
    z2.negate();
    let z3 = expansion_product(&m12, a02);
    let result = expansion_sum3(&z1, &z2, &z3);
    result.sign()
}

pub fn sign_of_expansion_det4x4(
    a00: &Expansion,
    a01: &Expansion,
    a02: &Expansion,
    a03: &Expansion,
    a10: &Expansion,
    a11: &Expansion,
    a12: &Expansion,
    a13: &Expansion,
    a20: &Expansion,
    a21: &Expansion,
    a22: &Expansion,
    a23: &Expansion,
    a30: &Expansion,
    a31: &Expansion,
    a32: &Expansion,
    a33: &Expansion,
) -> Sign {
    let m01 = expansion_det2x2(a10, a00, a11, a01);
    let m02 = expansion_det2x2(a20, a00, a21, a01);
    let m03 = expansion_det2x2(a30, a00, a31, a01);
    let m12 = expansion_det2x2(a20, a10, a21, a11);
    let m13 = expansion_det2x2(a30, a10, a31, a11);
    let m23 = expansion_det2x2(a30, a20, a31, a21);
    let m012_1 = expansion_product(&m12, a02);
    let mut m012_2 = expansion_product(&m02, a12);
    m012_2.negate();
    let m012_3 = expansion_product(&m01, a22);
    let m012 = expansion_sum3(&m012_1, &m012_2, &m012_3);
    let m013_1 = expansion_product(&m13, a02);
    let mut m013_2 = expansion_product(&m03, a12);
    m013_2.negate();
    let m013_3 = expansion_product(&m01, a32);
    let m013 = expansion_sum3(&m013_1, &m013_2, &m013_3);
    let m023_1 = expansion_product(&m23, a02);
    let mut m023_2 = expansion_product(&m03, a22);
    m023_2.negate();
    let m023_3 = expansion_product(&m02, a32);
    let m023 = expansion_sum3(&m023_1, &m023_2, &m023_3);
    let m123_1 = expansion_product(&m23, a12);
    let mut m123_2 = expansion_product(&m13, a22);
    m123_2.negate();
    let m123_3 = expansion_product(&m12, a32);
    let m123 = expansion_sum3(&m123_1, &m123_2, &m123_3);
    let m0123_1 = expansion_product(&m123, a03);
    let m0123_2 = expansion_product(&m023, a13);
    let m0123_3 = expansion_product(&m013, a23);
    let m0123_4 = expansion_product(&m012, a33);
    let z1 = expansion_sum(&m0123_1, &m0123_3);
    let z2 = expansion_sum(&m0123_2, &m0123_4);
    let result = expansion_diff(&z1, &z2);
    result.sign()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_product_retains_low_component() {
        let a = 134_217_729.0;
        let (high, low) = two_product(a, a);
        assert_eq!(high, 18_014_398_777_917_440.0);
        assert_eq!(low, 1.0);
    }

    #[test]
    fn expansion_cancellation_is_exact() {
        let large = expansion_create(1.0e16);
        let one = expansion_create(1.0);
        let sum = expansion_sum(&large, &one);
        let restored = expansion_diff(&sum, &large);
        assert_eq!(restored.sign(), Sign::Positive);
        assert_eq!(restored.as_slice(), &[1.0]);
    }
}
