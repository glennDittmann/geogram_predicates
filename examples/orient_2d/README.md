# Usage

```sh
cargo run --example orient_2d naive
cargo run --example orient_2d robust
```

This will save the output in `out_{naive, robust}_orient_2d.png`. Or see the pre-computed examples in `images/`.

# Idea

This example computes the `orient_2d` predicate over a
2D grid of values for one of the inputs, while keeping the
rest of the inputs at fixed values.

For example, we compute the `orient2d` predicate on `(12.0,
12.0)`, `(c[i], c[j])`, `(24.0, 24.0)`; where `i, j` varies
in `0..256`, and `c[i]` is the `i`th float after `0.5`. In
other words, `c[i]` is obtained by starting at `0.5`, and
calling
[`nextafter`](https://docs.rs/float_extras/*/float_extras/f64/fn.nextafter.html)
`i` times.

The inputs are set up so that, if the predicates are
calculated exactly, the output is a `png` with gray values on
the main diagonal, black on the lower-left, and white on
the upper-right side of it. However, the naive versions show
that the predicate is not robust: it switches values on both
sides of the main diagonals, indicating the round-off errors.

# Credits
This example was adopted from the [georust/robust](https://github.com/georust/robust/tree/main) crate.