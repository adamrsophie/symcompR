# Vector calculus in Cartesian coordinates

Vectors are column matrices of three components. `grad()`,
`divergence()`, `curl()` and `laplacian()` differentiate with respect to
`vars` (default `x`, `y`, `z`).
[`dot()`](https://adamrsophie.github.io/symcompR/reference/dot.md) and
[`cross()`](https://adamrsophie.github.io/symcompR/reference/dot.md)
come from the matrix module. `is_conservative()` and `is_solenoidal()`
test fields, and `scalar_potential()` integrates a conservative field.

## Usage

``` r
grad(f, vars = .xyz(), doit = TRUE)

divergence(v, vars = .xyz(), doit = TRUE)

curl(v, vars = .xyz(), doit = TRUE)

laplacian(f, vars = .xyz())

is_conservative(v, vars = .xyz())

is_solenoidal(v, vars = .xyz())

scalar_potential(v, vars = .xyz())
```

## Arguments

- f:

  A scalar expression.

- vars:

  A list of three symbols. When `vars` is not given and the field is a
  [`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
  vector or refers to base scalars, the functions of `sympy.vector` are
  used instead; `scalar_potential(v, coord_sys)` then takes the
  coordinate system.

- doit:

  Evaluate the derivatives (coordinate system fields only).

- v:

  A vector field: a list of three expressions or a column matrix.

## Value

An expression or matrix.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
curl(list(-y, x, 0))
#> Matrix(3, 1)
#> [0]
#> [0]
#> [2]
divergence(list(x, y, z))
#> 3
```
