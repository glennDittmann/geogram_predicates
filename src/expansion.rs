//! Multi-precision expansion arithmetic (Shewchuk-style).
//! Exact predicates are implemented by computing in this number type and taking the sign.

use crate::sign::{geo_sgn, Sign};
use std::vec::Vec;

/// Global splitter for split() - set by expansion::initialize().
static mut EXPANSION_SPLITTER: f64 = 0.0;

/// Expansion: non-overlapping sum of f64 components, least significant first.
/// Same semantics as Geogram: component 0 is least significant, last is most significant.
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
        if a.len() == 2 && b.len() == 2 {
            let out = two_two_product(&a.x[0..2], &b.x[0..2]);
            self.x.clear();
            for &v in &out {
                if v != 0.0 {
                    self.x.push(v);
                }
            }
            if self.x.is_empty() {
                self.x.push(0.0);
            }
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

fn expansion_splitter() -> f64 {
    unsafe { EXPANSION_SPLITTER }
}

pub fn expansion_initialize() {
    let half = 0.5f64;
    let mut expansion_epsilon = 1.0f64;
    let mut expansion_splitter = 1.0f64;
    let mut check = 1.0f64;
    let mut every_other = true;
    loop {
        let lastcheck = check;
        expansion_epsilon *= half;
        if every_other {
            expansion_splitter *= 2.0;
        }
        every_other = !every_other;
        check = 1.0 + expansion_epsilon;
        if check == 1.0 || check == lastcheck {
            break;
        }
    }
    expansion_splitter += 1.0;
    unsafe {
        EXPANSION_SPLITTER = expansion_splitter;
    }
}

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
    let c = expansion_splitter() * a;
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

fn two_two_product(a: &[f64], b: &[f64]) -> [f64; 8] {
    let (bhi, blo) = split(b[0]);
    let (mut _i, x0) = two_product_presplit(a[0], b[0], bhi, blo);
    let (_a1hi, _a1lo) = split(a[1]);
    let (_j, _0) = two_product_presplit(a[1], b[0], bhi, blo);
    let (sum_k, _1) = two_sum(_i, _0);
    let (l, _2) = fast_two_sum(_j, sum_k);
    let (b1hi, b1lo) = split(b[1]);
    let (_i, _0) = two_product_presplit(a[0], b[1], b1hi, b1lo);
    let (sum_k, x1) = two_sum(_1, _0);
    let (_j, _1) = two_sum(_2, sum_k);
    let (m, _2) = two_sum(l, _j);
    let (_j, _0) = two_product_presplit(a[1], b[1], b1hi, b1lo);
    let (n, _0) = two_sum(_i, _0);
    let (_i, x2) = two_sum(_1, _0);
    let (sum_k, _1) = two_sum(_2, _i);
    let (l, _2) = two_sum(m, sum_k);
    let (sum_k, _0) = two_sum(_j, n);
    let (_j, x3) = two_sum(_1, _0);
    let (_i, _1) = two_sum(_2, _j);
    let (m, _2) = two_sum(l, _i);
    let (_i, x4) = two_sum(_1, sum_k);
    let (sum_k, x5) = two_sum(_2, _i);
    let (x7, x6) = two_sum(m, sum_k);
    [x0, x1, x2, x3, x4, x5, x6, x7]
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

pub fn expansion_sum4(a: &Expansion, b: &Expansion, c: &Expansion, d: &Expansion) -> Expansion {
    let ab = expansion_sum(a, b);
    let cd = expansion_sum(c, d);
    expansion_sum(&ab, &cd)
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
