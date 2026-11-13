# Express vectors, dyadics and scalars in a coordinate system

`express(expr, system)` re-expresses a vector, dyadic or scalar field in
`system`; with `variables = TRUE` the base scalars of other systems are
substituted too. `matrix_to_vector()` and `matrix_to_dyadic()` turn
component matrices into vectors and dyadics, `orthogonalize()` applies
Gram-Schmidt, `directional_derivative()` differentiates a field along a
vector, and `scalar_potential_difference()` evaluates the potential
difference of a field between two points.

## Usage

``` r
express(expr, system, system2 = NULL, variables = FALSE)

matrix_to_vector(matrix, system)

matrix_to_dyadic(matrix, system)

orthogonalize(..., orthonormal = FALSE)

directional_derivative(field, direction_vector)

scalar_potential_difference(field, coord_sys, point1, point2, origin = NULL)
```

## Arguments

- expr:

  A vector, dyadic or scalar expression.

- system, system2:

  Coordinate systems
  ([`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)).

- variables:

  Substitute the coordinate variables.

- matrix:

  A 3 by 1 (or 3 by 3) matrix.

- ...:

  Vectors.

- orthonormal:

  Normalize the result.

- field:

  A scalar or vector field.

- direction_vector:

  A vector.

- coord_sys:

  A coordinate system.

- point1, point2:

  Points.

- origin:

  For physics reference frames, the origin point used by
  `scalar_potential_difference()`; ignored for `CoordSys3D` systems.

## Value

A vector, dyadic, expression or list.

## Examples

``` r
N <- CoordSys3D("N")
A <- orient_new_axis(N, "A", Symbol("q"), N$k)
express(A$i, N)
#> (cos(q))*N.i + (sin(q))*N.j
matrix_to_vector(Matrix(list(1, 2, 3)), N)
#> N.i + 2*N.j + 3*N.k
orthogonalize(N$i + 2 * N$j, 2 * N$i + 3 * N$j)
#> [[1]]
#> N.i + 2*N.j
#> 
#> [[2]]
#> 2/5*N.i + (-1/5)*N.j
#> 
```
