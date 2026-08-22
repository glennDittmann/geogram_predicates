//! Type aliases and re-exports for predicate results.

pub use crate::sign::Sign;

/// A point with a stable identity used by Simulation of Simplicity (SOS).
///
/// The key must identify the same logical vertex across predicate calls and
/// must be unique among all SOS points passed to one predicate invocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SosPoint<const D: usize> {
    pub coords: [f64; D],
    pub key: u64,
}

impl<const D: usize> SosPoint<D> {
    #[inline]
    pub const fn new(coords: [f64; D], key: u64) -> Self {
        Self { coords, key }
    }
}
