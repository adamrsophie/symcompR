# B-splines

`bspline_basis(d, knots, n, x)` is the `n`-th B-spline of degree `d` on
the knots (zero based `n`, as in SymPy), a Piecewise in `x`.
`bspline_basis_set(d, knots, x)` lists all of them, and
`interpolating_spline(d, x, X, Y)` is the spline of degree `d` through
the points `(X, Y)`.

## Usage

``` r
bspline_basis(d, knots, n, x)

bspline_basis_set(d, knots, x)

interpolating_spline(d, x, X, Y)
```

## Arguments

- d:

  The degree.

- knots:

  The knots.

- n:

  The index, from 0.

- x:

  The variable.

- X, Y:

  Coordinates of the interpolation points.

## Value

A Piecewise, or a list of them.

## Examples

``` r
x <- Symbol("x")
bspline_basis(1, c(0, 1, 2), 0, x)
#> Piecewise((x, (x >= 0) & (x <= 1)), (2 - x, (x >= 1) & (x <= 2)), (0, True))
interpolating_spline(3, x, c(-2, 0, 1, 3, 4), c(4, 2, 1, 1, 3))
#> Piecewise((7*x^3/117 + 7*x^2/117 - 131*x/117 + 2, (x >= -2) & (x <= 1)), (10*x^3/117 - 2*x^2/117 - 122*x/117 + 77/39, (x >= 1) & (x <= 4)))
```
