use core::cmp::Ordering;

/// A helper for functions that return signs
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
#[repr(i8)]
pub enum Sign {
    Positive = 1,
    Zero = 0,
    Negative = -1,
}

impl From<Sign> for i8 {
    fn from(sign: Sign) -> i8 {
        sign as i8
    }
}

impl TryFrom<i8> for Sign {
    type Error = ();

    /// Converts `+1`, `0`, `-1` to a `Sign` with the same memory layout,
    /// if input is another number an it's an error
    fn try_from(i8: i8) -> Result<Sign, ()> {
        if i8 == 1 {
            Ok(Sign::Positive)
        } else if i8 == -1 {
            Ok(Sign::Negative)
        } else if i8 == 0 {
            Ok(Sign::Zero)
        } else {
            Err(())
        }
    }
}

impl PartialEq<i8> for Sign {
    fn eq(&self, other: &i8) -> bool {
        (self.clone() as i8) == *other
    }
}

impl PartialEq<Sign> for i8 {
    fn eq(&self, other: &Sign) -> bool {
        (other.clone() as i8) == *self
    }
}

impl PartialOrd<i8> for Sign {
    fn partial_cmp(&self, other: &i8) -> Option<Ordering> {
        (self.clone() as i8).partial_cmp(other)
    }
}

impl core::ops::Mul for Sign {
    type Output = Sign;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // SAFETY: Sign is repr(i8) and will be (-1, 0, 1) so the product must be (1, 0, -1)
        unsafe { core::mem::transmute::<i8, Sign>(
            (self as i8) * (rhs as i8)
        ) }
    }
}

impl core::ops::Neg for Sign {
    type Output = Sign;

    #[inline]
    fn neg(self) -> Self::Output {
        // SAFETY: Sign is repr(i8) and will be (-1, 0, 1) so this is safe
        unsafe { core::mem::transmute::<i8, Sign>(
            -(self as i8)
        ) }
    }
}
