# Geogram Predicates

A **dependency-free** Rust port of [Geogram](https://github.com/BrunoLevy/geogram)'s robust predicates (PCK). Uses the same approach as the original: floating-point filters plus exact multi-precision expansion arithmetic (Shewchuk).

No C++ or external crates are required.

## Example

```rust
use geogram_predicates::{initialize, orient_2d, in_circle_2d_sos, in_sphere_3d_sos, Sign};

initialize();

// 2D orientation (counter-clockwise = Positive)
let a = [0.0, 0.0];
let b = [2.0, 0.0];
let c = [1.0, 1.0];
assert_eq!(orient_2d(&a, &b, &c), Sign::Positive);

// In-circle test (Positive = inside circumcircle of triangle a,b,c)
let p_in = [1.0, -0.4];
let p_out = [1.0, -1.2];
assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_in), Sign::Positive);
assert_eq!(in_circle_2d_sos(&a, &b, &c, &p_out), Sign::Negative);

// In-sphere test (Negative = inside circumsphere)
let d = [0.0, 0.0, 0.0];
let e = [2.0, 0.0, 0.0];
let f = [0.0, 2.0, 0.0];
let g = [0.75, 0.75, 1.0];
let p = [0.75, 0.75, 0.5];
assert_eq!(in_sphere_3d_sos(&d, &e, &f, &g, &p), Sign::Negative);
```

## Supported predicates

| Predicate | Description |
|-----------|-------------|
| `orient_2d` | Sign of 2D orientation (p1−p0)×(p2−p0) |
| `orient_3d` | Sign of 3D orientation (tetrahedron volume) |
| `in_circle_2d_sos` | Is point inside circumcircle of triangle? (SOS) |
| `in_sphere_3d_sos` | Is point inside circumsphere of tetrahedron? |
| `det_3d` | Sign of 3×3 determinant |
| `det_4d` | Sign of 4×4 determinant |
| `dot_3d` | Sign of dot product (p1−p0)·(p2−p0) |
| `aligned_3d` | Are (p1−p0) and (p2−p0) collinear? |

## Initialization

Call `initialize()` before using any predicate (sets up expansion arithmetic constants). `terminate()` is optional (no-op in this port).

## Design

- **Dependency-free**: no `robust`, `nalgebra`, or C++.
- **API aligned with Geogram**: same names, signs (e.g. in_sphere: Negative = inside), and filter-then-exact pattern.
- **Exact arithmetic**: Shewchuk-style expansions (two_sum, two_product, etc.) with zero-elimination.

## License

LGPL-3.0 OR MIT. The C++ reference in `include/geogram_predicates_psm/` is Inria’s Geogram PSM (see that directory for its license).

## References

- Shewchuk, “Robust adaptive floating-point geometric predicates,” SoCG 1996.
- Shewchuk, “Adaptive precision floating-point arithmetic and fast robust geometric predicates,” DCG 1997.
- [Geogram](https://github.com/BrunoLevy/geogram)
