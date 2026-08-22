# Usage

```sh
cargo run --example in_circle_2d naive
cargo run --example in_circle_2d robust
```

This saves `out_{naive,robust}_in_circle_2d.pgm`, a dependency-free grayscale portable graymap. Pre-computed PNG renderings remain in [`images/`](../../images/).

# Idea

This example computes the `in_circle_2d` predicate over a
2D grid of values for one of the inputs, while keeping the
rest of the inputs at fixed values.

The same strategy as for the `orient_2d` [`example`](../orient_2d/) was used, i.e. generating minimal representable floating point values using `nextafter`.

# Credits
This example was adopted from the [georust/robust](https://github.com/georust/robust/tree/main) crate.
