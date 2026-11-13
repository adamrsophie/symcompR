# Pade approximant

`pade(e, x, m, n)` returns the `[m/n]` Pade approximant of `e` at `x0`,
computed from the extended Euclidean algorithm as in SymPy's
[`pade_approximant()`](https://adamrsophie.github.io/symcompR/reference/pade_approximants_gcdex.md).

## Usage

``` r
pade(e, x, m, n = m, x0 = 0)
```

## Arguments

- e:

  An expression.

- x:

  The variable.

- m, n:

  Degrees of numerator and denominator.

- x0:

  The expansion point.

## Value

A rational function.

## Examples

``` r
x <- Symbol("x")
pade(exp(x), x, 2, 2)
#> (x^2/4 + 3*x/2 + 3)/(x^2/4 - 3*x/2 + 3)
```
