//! macros to construct and operate on `Expansion`s

/// Build an `Expansion` from a comma-separated list of `f64` literals or expressions.
///
/// ```
/// # use geogram_predicates::expansion;
/// let e = expansion![1.0, 2.5, 3.75];
/// assert_eq!(e.length(), 3);
/// assert_eq!(&e.data(), &[1.0, 2.5, 3.75]);
/// ```
#[macro_export]
macro_rules! expansion {
    // match zero or more comma-separated expressions, allow trailing comma
    ($($x:expr),* $(,)?) => {
        // create a fixed-length array `[f64; count]`, then call From<[f64;N]>
        $crate::Expansion::from([$($x),*])
    };
    ($x:expr) => {
        $crate::Expansion::from($x)
    };
}

/// Compute the 3×3 determinant of nine `Expansion`s, returning a new `Expansion`.
///  
/// ```ignore
/// # use geogram_predicates::expansion_det3x3;
/// expansion_det3x3!(
///     a, b, c,
///     d, e, f,
///     g, h, i,
/// )
/// // or
/// expansion_det3x3!(
///     [a, b, c],
///     [d, e, f],
///     [g, h, i],
/// )
/// ```
macro_rules! expansion_det3x3 {
    (
        $a11:expr, $a12:expr, $a13:expr,
        $a21:expr, $a22:expr, $a23:expr,
        $a31:expr, $a32:expr, $a33:expr$(,)?
    ) => {{
        // Compute exactly the capacity needed
        let cap = $crate::Expansion::det3x3_capacity(
            [&$a11, &$a12, &$a13],
            [&$a21, &$a22, &$a23],
            [&$a31, &$a32, &$a33],
        );
        // Allocate an Expansion with that capacity
        let mut e = $crate::Expansion::with_capacity(cap);
        // Perform the determinant assignment
        e.assign_det3x3(
            [&$a11, &$a12, &$a13],
            [&$a21, &$a22, &$a23],
            [&$a31, &$a32, &$a33],
        );
        e
    }};
    (
        [$a11:expr, $a12:expr, $a13:expr],
        [$a21:expr, $a22:expr, $a23:expr],
        [$a31:expr, $a32:expr, $a33:expr]$(,)?
    ) => {{
        // Compute exactly the capacity needed
        let cap = $crate::Expansion::det3x3_capacity(
            [&$a11, &$a12, &$a13],
            [&$a21, &$a22, &$a23],
            [&$a31, &$a32, &$a33],
        );
        // Allocate an Expansion with that capacity
        let mut e = $crate::Expansion::with_capacity(cap);
        // Perform the determinant assignment
        e.assign_det3x3(
            [&$a11, &$a12, &$a13],
            [&$a21, &$a22, &$a23],
            [&$a31, &$a32, &$a33],
        );
        e
    }};
}

/// Compute the 2×2 determinant of four `Expansion`s, returning a new `Expansion`.
///  
/// ```ignore
/// # use geogram_predicates::expansion_det2x2;
/// expansion_det2x2!(
///     a, b,
///     c, d,
/// )
/// // or
/// expansion_det2x2!(
///     [a, b],
///     [c, d],
/// )
/// ```
macro_rules! expansion_det2x2 {
    (
        [$a11:expr, $a12:expr],
        [$a21:expr, $a22:expr]$(,)?
    ) => {{
        // Compute exactly the capacity needed
        #[cfg(debug_assertions)]
        let cap = $crate::Expansion::det2x2_capacity(&$a11, &$a12, &$a21, &$a22);

        // Allocate an Expansion with that capacity
        let mut e = $crate::Expansion::new();

        #[cfg(debug_assertions)]
        debug_assert_eq!(cap, e.capacity());
        // Perform the determinant assignment
        e.assign_det2x2(&$a11, &$a12, &$a21, &$a22);
        e
    }};
}

macro_rules! expansion_diff {
    ($a:expr, $b:expr) => {{
        let mut expansion = Expansion::<2>::with_capacity(2);
        expansion.assign_diff(&$a.into(), &$b.into());
        expansion
    }};
    (Expansions: $a:expr, $b:expr) => {{
        let mut expansion = Expansion::<2>::with_capacity(2);
        expansion.assign_diff(&$a, &$b);
        expansion
    }};
    ($a:expr, $b:expr, $AN:literal, $BN:literal) => {{
        let mut expansion = Expansion::with_capacity(2);
        expansion.assign_diff::<$AN, $BN>(&$a.into(), &$b.into());
        expansion
    }};
}

macro_rules! expansion_sum {
    ($a:expr, $b:expr) => {{
        let mut expansion = Expansion::new();
        expansion.assign_sum(&$a, &$b);

        debug_assert_eq!(expansion.data_mut().len(), expansion.data_mut().capacity());
        expansion
    }};
}

macro_rules! expansion_sum3 {
    ($a:expr, $b:expr, $c:expr) => {{
        let capacity = $a.length() + $b.length();
        let mut ab = Expansion::with_capacity(capacity);
        ab.assign_sum(&$a, &$b);
        let mut expansion = Expansion::with_capacity(capacity + $c.length());
        expansion.assign_sum(&ab, &$c);
        expansion
    }};
}

macro_rules! expansion_sum4 {
    ($a:expr, $b:expr, $c:expr, $d:expr, $ab_capacity:literal, $cd_capacity:literal) => {{
        // let ab_capacity = $a.length() + $b.length();
        let mut ab = Expansion::<$ab_capacity>::new();
        ab.assign_sum(&$a, &$b);

        // let cd_capacity = $c.length() + $d.length();
        let mut cd = Expansion::<$cd_capacity>::new();
        cd.assign_sum(&$c, &$d);
        let mut expansion = Expansion::new();
        expansion.assign_sum(&ab, &cd);
        expansion
    }};
}

macro_rules! expansion_product {
    ($a:expr, $b:expr) => {{
        let mut expansion = Expansion::<2>::new();
        expansion.assign_product(&$a, &$b);
        expansion
    }};
    ($a:expr, $b:expr, $res_cap:literal) => {{
        let mut expansion = Expansion::<$res_cap>::new();
        expansion.assign_product(&$a, &$b);
        expansion
    }};
}

pub(crate) use {
    expansion_det2x2, expansion_det3x3, expansion_diff, expansion_product, expansion_sum,
    expansion_sum3, expansion_sum4,
};
