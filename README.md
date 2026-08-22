# Geogram Predicates

A **dependency-free** Rust port of [Geogram](https://github.com/BrunoLevy/geogram)'s robust predicates (PCK). Uses the same approach as the original: floating-point filters plus exact multi-precision expansion arithmetic (Shewchuk).

No C++ or external crates are required.

## Example

```rust
use geogram_predicates::{in_circle_2d_sos, orient_2d, Sign, SosPoint};

// 2D orientation (counter-clockwise = Positive)
let a = [0.0, 0.0];
let b = [2.0, 0.0];
let c = [1.0, 1.0];
assert_eq!(orient_2d(&a, &b, &c), Sign::Positive);

// In-circle test (Positive = inside circumcircle of triangle a,b,c)
let a = SosPoint::new(a, 10);
let b = SosPoint::new(b, 20);
let c = SosPoint::new(c, 30);
let p_in = SosPoint::new([1.0, -0.4], 40);
let p_out = SosPoint::new([1.0, -1.2], 50);
assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_in), Sign::Positive);
assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_out), Sign::Negative);
```

## Supported predicates

| Predicate | Description |
|-----------|-------------|
| `orient_2d` | Sign of 2D orientation (p1−p0)×(p2−p0) |
| `orient_3d` | Sign of 3D orientation (tetrahedron volume) |
| `in_circle_2d[_sos]`, `in_circle_3d[_sos]` | In-circle tests in 2D and embedded 3D |
| `in_sphere_3d[_sos]` | In-sphere test; positive means inside for a positive tetrahedron |
| `orient_2dlifted[_sos]`, `orient_3dlifted[_sos]` | Weighted/lifted orientation tests |
| `side1_sos` … `side4_sos` | Generic PCK power-side predicates in dimensions 3, 4, 6, 7, and 8 |
| `det_3d` | Sign of 3×3 determinant |
| `det_4d` | Sign of 4×4 determinant |
| `det_compare_4d`, `dot_compare_3d` | Exact determinant and dot-product comparisons |
| `dot_3d` | Sign of dot product (p1−p0)·(p2−p0) |
| `aligned_3d` | Are (p1−p0) and (p2−p0) collinear? |

## Simulation of Simplicity

SOS predicates take `SosPoint`s. A key is a stable identity for a logical
vertex and must be unique within a predicate call. Exact degeneracies are
resolved by key order, making results reproducible across runs and machines.
Duplicate keys panic with a contract error. Non-SOS variants return `Zero` for
an exact degeneracy.

No initialization is required. `initialize()`, `terminate()`, and
`show_stats()` are safe compatibility no-ops.

Predicate inputs must be finite `f64` values. As in the upstream PSM, NaN,
infinity, and computations that overflow binary64 are outside the supported
input domain. Side predicates also require a nondegenerate query simplex.

## Design

- **Dependency-free**: no `robust`, `nalgebra`, or C++.
- **Safe Rust**: the crate uses `#![forbid(unsafe_code)]`.
- **API aligned with Geogram**: predicate names and signs follow the bundled PSM; SOS ordering uses stable keys instead of addresses.
- **Exact arithmetic**: Shewchuk-style expansions (two_sum, two_product, etc.) with zero-elimination.

## Differential testing against Geogram

The opt-in differential test compiles the bundled C++ PSM as a batched oracle
and compares it with the Rust implementation. It requires a C++17 compiler but
does not add C++, FFI, dependencies, or unsafe code to the crate itself.

```console
cargo test --release --test cpp_differential -- --ignored --nocapture
```

The default run uses 10,000 deterministic cases. Increase the corpus or select
a reproducible seed with environment variables:

```console
GEOGRAM_DIFF_CASES=1000000 GEOGRAM_DIFF_SEED=0x123456789abcdef0 \
    cargo test --release --test cpp_differential -- --ignored --nocapture
```

Set `CXX=clang++` (or another GNU-compatible C++ driver) to override the default
`c++` command. A mismatch reports its case index, seed, predicate, keys, and the
exact binary64 bits of every input. Normal `cargo test` runs compile this test
but skip its C++ oracle.

The PSM's public `dot_3d()` contains an upstream typo: it first calls the
`det_3d` filter. Differential dot-product cases therefore force the PSM's
correct exact fallback; the Rust implementation retains the documented dot
product semantics.

## License

LGPL-3.0 OR MIT. The C++ reference in `include/geogram_predicates_psm/` is Inria’s Geogram PSM (see that directory for its license).

## References

- Shewchuk, “Robust adaptive floating-point geometric predicates,” SoCG 1996.
- Shewchuk, “Adaptive precision floating-point arithmetic and fast robust geometric predicates,” DCG 1997.
- [Geogram](https://github.com/BrunoLevy/geogram)
