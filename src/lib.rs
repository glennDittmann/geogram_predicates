//! # Geogram Predicates
//!
//! A crate for rust interoperability with `geogram`s _robust predicates_; via `cxx`.

pub use geogram_ffi::*;

#[cxx::bridge(namespace = "GEOGRAM")]
mod geogram_ffi {
    // Shared structs with fields visible to both languages.
    // ...

    // Rust types and signatures exposed to C++.
    // ...

    // C++ types and signatures exposed to Rust.
    unsafe extern "C++" {
        include!("geogram_predicates/include/geogram_ffi.h");

        /// Gets the sign of a value.
        ///
        /// ### Parameters
        /// - `x` value to test
        ///
        /// ### Return values
        /// - `+1` if `x` is positive
        /// - `0` if `x` is `0`
        /// - `-1` if `x` is negative
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// let a = 42.0;
        /// let b = -42.0;
        /// let c = 0.0;
        ///
        /// assert_eq!(1, gp::geo_sgn(a));
        /// assert_eq!(-1, gp::geo_sgn(b));
        /// assert_eq!(0, gp::geo_sgn(c));
        ///
        /// ```
        fn geo_sgn(x: f64) -> i16;

        /// Tests whether a point is in the circum-circle of a triangle.
        ///
        /// If the triangle `a` , `b` , `c` is oriented clockwise instead of counter-clockwise, then the result is inversed.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c` vertices of the triangle
        /// - `p` point to test
        ///
        /// ### Return values
        /// * `+1` - if `p` is inside the circum-circle of `a`, `b`, `c`
        /// * `-1` - if `p` is outside the circum-circle of `a`, `b`, `c`
        /// * `perturb()` - if `p` is exactly on the circum-circle of the triangle `a`, `b`, `c`, where `perturb()` denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points that form a triangle
        /// let a = [0.0, 0.0];
        /// let b = [2.0, 0.0];
        /// let c = [1.0, 1.0];
        ///
        /// // Define two points, to test against the triangles circum-circle
        /// let p_in = [1.0, -0.4];
        /// let p_out = [1.0, -1.2];
        ///
        /// let is_in_circle_p_in = gp::in_circle_2d_SOS(&a, &b, &c, &p_in);
        /// assert_eq!(1, is_in_circle_p_in);
        ///
        /// let is_in_circle_p_out = gp::in_circle_2d_SOS(&a, &b, &c, &p_out);
        /// assert_eq!(-1, is_in_circle_p_out);
        /// ```
        fn in_circle_2d_SOS(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2], p: &[f64; 2]) -> i16;

        /// Tests whether a point is in the circum-sphere of a tetrahedron.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c`, `d` vertices of the tetrahedron
        /// - `p` point to test
        ///
        /// ### Return values
        /// * `+1` - if `p` is inside the circum-sphere of `a`, `b`, `c`, `d`
        /// * `-1` - if `p` is outside the circum-sphere of `a`, `b`, `c`, `d`
        /// * `perturb()` - if `p` is exactly on the circum-sphere of the tetrahedron `a`, `b`, `c`, `d`, where `perturb()` denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define four points that form a tetrahedron
        /// let a = [0.0, 0.0, 0.0];
        /// let b = [2.0, 0.0, 0.0];
        /// let c = [0.0, 2.0, 0.0];
        /// let d = [0.75, 0.75, 1.0];
        ///
        /// // Define two points, to test against the tetrahedrons circum-sphere
        /// let p_in = [0.75, 0.75, 0.5];
        /// assert_eq!(1, gp::in_sphere_3d_SOS(&a, &b, &c, &d, &p_in));
        ///
        /// let p_out = [0.75, 0.75, 1.5];
        /// assert_eq!(-1, gp::in_sphere_3d_SOS(&a, &b, &c, &d, &p_out));
        /// ```
        fn in_sphere_3d_SOS(
            a: &[f64; 3],
            b: &[f64; 3],
            c: &[f64; 3],
            d: &[f64; 3],
            p: &[f64; 3],
        ) -> i16;

        /// Computes the orientation predicate in 2d.
        ///
        /// Computes the sign of the signed area of the triangle `a`, `b`, `c`.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c` vertices of the triangle
        ///
        /// ### Return values
        /// * `+1` - if the triangle is oriented counter-clockwise
        /// * `0` - if the triangle is flat
        /// * `-1` - if the triangle is oriented clockwise
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points that form a triangle
        /// let a = [0.0, 0.0];
        /// let b = [2.0, 0.0];
        /// let c = [1.0, 1.0];
        ///
        /// let orientation = gp::orient_2d(&a, &b, &c);
        /// assert_eq!(1, orientation);
        /// ```
        fn orient_2d(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> i16;

        /// Computes the 3d orientation test with lifted points, i.e the regularity test for 2d.
        ///
        /// Given three lifted points `a'`, `b'`, `c'` in R^3, tests if the lifted point `p'` in R^3 lies below or above the plane passing through the three points `a'`, `b'`, `c'`.
        ///
        /// The coordinates and the heights are specified in separate arguments for each vertex.
        ///
        /// Note: if `w_i` = 0 this is equal to the in-circle test for the triangle `a`, `b`, `c` w.r.t `p`.
        ///
        /// ### Parameters
        /// - `a` ,`b`, `c`	vertices of the triangle
        /// - `p` point to test
        /// - `h_a` ,`h_b` ,`h_c` the heights of the lifted points, e.g. `a' = a.x**2 + a.y**2 - a.w`
        /// - `h_d` the height of the lifted point `p`
        ///
        /// ### Return values
        /// - `+1` - if p3' lies below the plane
        /// - `-1` - if p3' lies above the plane
        /// - perturb()	- if `p'` lies exactly on the hyperplane, where perturb() denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// For a graphical representation see this [geogebra example](https://www.geogebra.org/m/etyzj96t) of the code below.
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points that form a triangle
        /// let a: [f64; 2] = [0.0, 0.0];
        /// let b: [f64; 2] = [2.0, 0.0];
        /// let c: [f64; 2] = [0.0, 2.0];
        ///
        /// // Additionally in this scenario, each point is associated with a weight w_i
        /// // And the height of a point is defined as h_i = x_i**2 + y_i**2 - w_i
        /// // One can interpret the height as the z-coordinate of a point lifted to R^3
        /// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + 2.0;  // i.e. w_a = -2.0
        /// let h_b = b[0].powf(2.0) + b[1].powf(2.0) - 1.0;  // i.e. w_b = 1.0
        /// let h_c = c[0].powf(2.0) + c[1].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        ///
        /// // Define weighted points, to test against the plane, that contains the lifted triangle
        /// let p_below: [f64; 2] = [0.6, 0.6];
        /// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 1.28;  // i.e. w_p_below = -1.28
        ///
        /// let p_above: [f64; 2] = [0.6, 0.6];
        /// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
        ///
        /// let orientation_below = gp::orient_2dlifted_SOS(&a, &b, &c, &p_below, h_a, h_b, h_c, h_p_below);
        /// assert_eq!(1, orientation_below);
        ///
        /// let orientation_above = gp::orient_2dlifted_SOS(&a, &b, &c, &p_above, h_a, h_b, h_c, h_p_above);
        /// assert_eq!(-1, orientation_above);
        ///
        /// ```
        #[allow(clippy::too_many_arguments)]
        fn orient_2dlifted_SOS(
            a: &[f64; 2],
            b: &[f64; 2],
            c: &[f64; 2],
            p: &[f64; 2],
            h_a: f64,
            h_b: f64,
            h_c: f64,
            h_d: f64,
        ) -> i16;

        /// Computes the orientation predicate in 3d.
        ///
        /// Computes the sign of the signed volume of the tetrahedron `a`, `b`, `c`, `d`.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c`, `d` vertices of the tetrahedron
        ///
        /// ### Return values
        /// * `+1` - if the tetrahedron is oriented positively
        /// * `0` - if the tetrahedron is flat
        /// * `-1` - if the tetrahedron is oriented negatively
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define four points that form a tetrahedron
        /// let a = [0.0, 0.0, 0.0];
        /// let b = [2.0, 0.0, 0.0];
        /// let c = [0.0, 2.0, 0.0];
        /// let d = [0.75, 0.75, 1.0];
        ///
        /// assert_eq!(1, gp::orient_3d(&a, &b, &c, &d));
        ///```
        fn orient_3d(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], d: &[f64; 3]) -> i16;

        /// Computes the 4d orientation test with lifted points, i.e the regularity test for 3d.
        ///
        /// Given four lifted points `a'`, `b'`, `c'`, `d'` in R^4, tests if the lifted point `p'` in R^4 lies below or above the hyperplane passing through the four points `a'`, `b'`, `c'`, `d'`.
        ///
        /// Symbolic perturbation is applied, whenever the 5 vertices are not linearly independent.
        ///
        /// The coordinates and the heights are specified in separate arguments for each vertex.
        ///
        /// Note: if `w_i` = 0 this is equal to the in-sphere test for the tetrahedron `a`, `b`, `c`, `d` w.r.t `p`.
        ///
        /// ### Parameters
        /// - `a` ,`b`, `c`, `d` vertices of the tetrahedron
        /// - `p` point to test
        /// - `h_a` ,`h_b` ,`h_c`, `h_d` the heights of the lifted points, e.g. `h_a' = a.x**2 + a.y**2 - a.w`
        /// - `h_p` the height of the lifted point `p'`
        ///
        /// ### Return values
        /// - `+1` - if `p'` lies below the plane
        /// - `-1` - if `p'` lies above the plane
        /// - perturb()	- if `p'` lies exactly on the hyperplane, where perturb() denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        /// // Define four points that form a tetrahedron
        /// let a: [f64; 3] = [0.0, 0.0, 0.0];
        /// let b: [f64; 3] = [2.0, 0.0, 0.0];
        /// let c: [f64; 3] = [0.0, 2.0, 0.0];
        /// let d: [f64; 3] = [0.75, 0.75, 1.0];
        ///
        /// // Additionally in this scenario, each point is associated with a weight w_i
        /// // And the height of a point is defined as h_i = x_i**2 + y_i**2 + z_i**2 - w_i
        /// // One can interpret the height as the 4th-coordinate of a point lifted to R^4
        /// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + a[2].powf(2.0) + 2.0;  // i.e. w_a = -2.0
        /// let h_b = b[0].powf(2.0) + b[1].powf(2.0) + b[2].powf(2.0) - 1.0;  // i.e. w_b = 1.0
        /// let h_c = c[0].powf(2.0) + c[1].powf(2.0) + c[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        /// let h_d = d[0].powf(2.0) + d[1].powf(2.0) + d[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        ///
        /// // Define weighted points, to test against the hyperplane, that contains the lifted tetrahedron
        /// let p_below: [f64; 3] = [0.6, 0.6, 0.6];
        /// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 0.28;  // i.e. w_p_below = -0.28
        ///
        /// let p_above: [f64; 3] = [0.6, 0.6, 0.6];
        /// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
        ///
        /// let orientation_below = gp::orient_3dlifted_SOS(&a, &b, &c, &d, &p_below, h_a, h_b, h_c, h_d, h_p_below);
        /// assert_eq!(1, orientation_below);
        ///
        /// let orientation_above = gp::orient_3dlifted_SOS(&a, &b, &c, &d, &p_above, h_a, h_b, h_c, h_d, h_p_above);
        /// assert_eq!(-1, orientation_above);
        ///
        /// ```
        #[allow(clippy::too_many_arguments)]
        fn orient_3dlifted_SOS(
            a: &[f64; 3],
            b: &[f64; 3],
            c: &[f64; 3],
            d: &[f64; 3],
            p: &[f64; 3],
            h_a: f64,
            h_b: f64,
            h_c: f64,
            h_d: f64,
            h_p: f64,
        ) -> i16;

        /// Tests whether two 2d points are identical.
        ///
        /// ### Parameters
        /// - `p1` first point
        /// - `p2` second point
        ///
        /// ### Return values
        /// - `true` - if `p1` and `p2` have exactly the same coordinates
        /// - `false` - otherwise
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// let p1 = [4.0, 2.0];
        /// let p2 = [4.0, 2.0];
        ///
        /// assert!(gp::points_are_identical_2d(&p1, &p2));
        /// ```
        fn points_are_identical_2d(p1: &[f64; 2], p2: &[f64; 2]) -> bool;

        /// Tests whether two 3d points are identical.
        ///
        /// ### Parameters
        /// - `p1` first point
        /// - `p2` second point
        ///
        /// ### Return values
        /// - `true` - if `p1` and `p2` have exactly the same coordinates
        /// - `false` - otherwise
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// let p1 = [4.0, 2.0, 0.42];
        /// let p2 = [4.0, 2.0, 0.42];
        ///
        /// assert!(gp::points_are_identical_3d(&p1, &p2));
        /// ```
        fn points_are_identical_3d(p1: &[f64; 3], p2: &[f64; 3]) -> bool;
    }
}
