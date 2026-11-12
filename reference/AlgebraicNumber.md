# Algebraic numbers

`AlgebraicNumber(expr, coeffs, alias)` is the element
`coeffs[1]*theta^(n-1) + ... + coeffs[n]` of the number field generated
by `theta = expr`, as SymPy's `AlgebraicNumber`; by default the element
is `theta` itself (`coeffs = c(1, 0)`). `expr` may also be
`list(minpoly, root)`. It prints as its value, or as a polynomial in
`alias` when one is given.
[`as_expr()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
gives the value,
[`coeffs()`](https://adamrsophie.github.io/symcompR/reference/degree.md)
the coordinates, and
[`minimal_polynomial()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
its minimal polynomial. `field_element(a, coeffs)` is another element of
the same field, `is_primitive_element(a)` tells whether `a` generates
its field, `to_primitive_element(a)` makes `a` the generator of a new
field, `minpoly_of_element(a)` and `to_root(a)` are the minimal
polynomial and an explicit value of `a`, and `isolate(a, eps)` a
rational interval that isolates the real algebraic number `a` among the
roots of its minimal polynomial (of width below `eps` when given).

## Usage

``` r
AlgebraicNumber(expr, coeffs = NULL, alias = NULL, gen = NULL)

field_element(a, coeffs)

is_primitive_element(a)

minpoly_of_element(a, x = NULL)

to_root(a)

to_primitive_element(a)

isolate(a, eps = NULL, fast = FALSE)
```

## Arguments

- expr:

  An algebraic number, an `AlgebraicNumber`, or `list(minpoly, root)`.

- coeffs:

  Coordinates in the power basis, highest power first.

- alias:

  A symbol (or name) printed for the generator.

- gen:

  The variable of the minimal polynomial.

- a:

  An `AlgebraicNumber` (or an algebraic number for `isolate()`).

- x:

  A variable.

- eps:

  Width of the isolating interval.

- fast:

  Use fast refinement steps.

## Value

An `AlgebraicNumber`, an expression, a logical, or a list.

## Examples

``` r
a <- AlgebraicNumber(sqrt(Integer(2)), c(1, 1))
a
#> 1 + sqrt(2)
coeffs(a)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 1
#> 
minimal_polynomial(a, Symbol("x"))
#> x^2 - 2*x - 1
AlgebraicNumber(sqrt(Integer(2)), alias = "theta")
#> theta
isolate(sqrt(Integer(2)), eps = Rational(1, 100))
#> [[1]]
#> 24/17
#> 
#> [[2]]
#> 17/12
#> 
```
