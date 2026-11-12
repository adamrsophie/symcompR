# Polynomial composition, Horner form and interpolation

`compose(f, g, x)` computes `f(g)`. `horner(f, x)` rewrites in Horner
form. `interpolate(points, x)` returns the interpolating polynomial for
a list of `c(x, y)` pairs, or for values at `x = 1, 2, ...`.

## Usage

``` r
compose(f, g, x = NULL)

horner(f, x = NULL)

interpolate(points, x)
```

## Arguments

- f, g:

  Expressions.

- x:

  The variable.

- points:

  A list of `c(x, y)` pairs, or a vector of values.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
horner(x^3 + 2*x^2 + 3*x + 4, x)
#> x*(x*(x + 2) + 3) + 4
interpolate(list(c(1, 1), c(2, 4), c(3, 9)), x)
#> x^2
```
