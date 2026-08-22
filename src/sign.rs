//! Sign type matching Geogram's PCK API.

use core::cmp::Ordering;

/// Result of a predicate: sign of a determinant or comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(i8)]
pub enum Sign {
    Negative = -1,
    Zero = 0,
    Positive = 1,
}

impl Sign {
    #[inline]
    pub const fn from_i8(i: i8) -> Self {
        match i {
            i if i > 0 => Sign::Positive,
            i if i < 0 => Sign::Negative,
            _ => Sign::Zero,
        }
    }

    #[inline]
    pub const fn as_i8(self) -> i8 {
        self as i8
    }
}

impl From<i8> for Sign {
    fn from(i: i8) -> Self {
        Self::from_i8(i)
    }
}

impl PartialEq<i8> for Sign {
    fn eq(&self, other: &i8) -> bool {
        self.as_i8() == *other
    }
}

impl PartialEq<Sign> for i8 {
    fn eq(&self, other: &Sign) -> bool {
        other.as_i8() == *self
    }
}

impl PartialOrd<i8> for Sign {
    fn partial_cmp(&self, other: &i8) -> Option<Ordering> {
        self.as_i8().partial_cmp(other)
    }
}

impl core::ops::Mul for Sign {
    type Output = Sign;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Sign::from_i8(self.as_i8().saturating_mul(rhs.as_i8()))
    }
}

impl core::ops::Neg for Sign {
    type Output = Sign;
    #[inline]
    fn neg(self) -> Self::Output {
        Sign::from_i8(-self.as_i8())
    }
}

/// Returns the sign of a floating-point value.
#[inline(always)]
pub fn geo_sgn(x: f64) -> Sign {
    if x > 0.0 {
        Sign::Positive
    } else if x < 0.0 {
        Sign::Negative
    } else {
        Sign::Zero
    }
}
