# Finite differences

`finite_diff(expression, variable, increment)` is the forward difference
`f(x + h) - f(x)` of the expanded expression, and
`finite_diff_kauers(sum)` the difference of a
[`Sum()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
with respect to its upper limits, as in `sympy.series.kauers`.

## Usage

``` r
finite_diff(expression, variable, increment = 1)

finite_diff_kauers(sum)
```

## Arguments

- expression:

  An expression.

- variable:

  A symbol.

- increment:

  The step.

- sum:

  A `Sum`.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
finite_diff(x^2 + 2*x + 1, x)
#> 2*x + 3
n <- Symbol("n")
finite_diff_kauers(Sum(x^2, list(x, 1, n)))
#> (n + 1)^2
```
