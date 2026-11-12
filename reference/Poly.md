# Polynomial objects

`Poly(expr, x, y, ...)` represents `expr` as a polynomial in the given
generators (by default, every symbol and non polynomial subexpression).
Coefficients are exact rationals; symbols that are not generators become
part of the coefficients when generators are given.

## Usage

``` r
Poly(expr, ...)

as_expr(p, ...)
```

## Arguments

- expr:

  An expression.

- ...:

  Generators.

- p:

  A `Poly`.

## Value

A `Poly` object.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
p <- Poly(x^2 + 2 * x * y + 1, x, y)
p
#> Poly(x^2 + 2*x*y + 1, x, y, domain='ZZ')
degree(p, x)
#> [1] 2
```
