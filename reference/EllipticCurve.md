# Elliptic curves

`EllipticCurve(a4, a6, a1, a2, a3, modulus)` is the curve
`y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6` over the rationals
(`modulus = 0`) or over the integers modulo a prime. The result holds
the fields `characteristic`, `discriminant`, `is_singular` and
`j_invariant`, and the methods `point(x, y, z = 1)` (also
`EllipticCurvePoint()`), `point_at_infinity()`,
[`points()`](https://rdrr.io/r/graphics/points.html) (a matrix of the
affine points over a finite field), `points_x(x)`, `torsion_points()`,
`minimal()`, [`order()`](https://rdrr.io/r/base/order.html) and
[`rank()`](https://rdrr.io/r/base/rank.html). Points support `+`, `-`,
`*` by an integer and `==`, and `point_order()` gives the order of a
point (`Inf` for infinite order).

## Usage

``` r
EllipticCurve(a4, a6, a1 = 0, a2 = 0, a3 = 0, modulus = 0)

EllipticCurvePoint(x, y, z, curve)

point_order(p)
```

## Arguments

- a4, a6, a1, a2, a3:

  Coefficients.

- modulus:

  0, or a prime.

- x, y, z:

  Projective coordinates.

- curve:

  An `EllipticCurve`.

- p:

  A point.

## Value

An `EllipticCurve` environment, or an `EllipticCurvePoint`.

## Details

Several SymPy defects are fixed:
[`order()`](https://rdrr.io/r/base/order.html) counts the point at
infinity, so it is the order of the group; a point `P` with `P == -P`
has order 2 even when `a1` or `a3` is nonzero; `point_order()` over the
rationals relies on Mazur's bound of 12 rather than integral
coordinates; `torsion_points()` works through an integral short model,
so it also finds the torsion of curves with rational coefficients; and
the discriminant of a curve with rational coefficients is exact rather
than truncated. Over a finite field,
[`points()`](https://rdrr.io/r/graphics/points.html) lists the affine
points as a numeric matrix, so the modulus must be below `2^53`, and
`points_x(x)` returns a list of `list(x, y)`.

## Examples

``` r
e <- EllipticCurve(-1, 9)
p <- e$point(0, 3)
q <- e$point(-1, 3)
p + q
#> (1, -3) 
2 * p + q
#> (35, 207) 
e$discriminant
#> [1] -34928
EllipticCurve(0, 1)$torsion_points()
#> [[1]]
#> O 
#> 
#> [[2]]
#> (-1, 0) 
#> 
#> [[3]]
#> (0, 1) 
#> 
#> [[4]]
#> (0, -1) 
#> 
#> [[5]]
#> (2, 3) 
#> 
#> [[6]]
#> (2, -3) 
#> 
f <- EllipticCurve(1, 0, modulus = 19)
f$order()
#> [1] 19
point_order(f$point(3, 7))
#> [1] 20
```
