# Riemann sums and principal values of integrals

`integral_as_sum(I, n, method)` approximates the definite integral `I`
by a sum over `n` subintervals (a positive integer, or by default a
positive integer symbol) with the `"midpoint"`, `"left"`, `"right"` or
`"trapezoid"` rule, as SymPy's `Integral.as_sum()`. `principal_value(I)`
is the Cauchy principal value of `I`, as SymPy's
`Integral.principal_value()`.

## Usage

``` r
integral_as_sum(I, n = NULL, method = "midpoint", evaluate = TRUE)

principal_value(I, ...)
```

## Arguments

- I:

  An unevaluated
  [`Integral()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md)
  over one variable with limits.

- n:

  The number of subintervals.

- method:

  The rule.

- evaluate:

  If `FALSE`, return the unevaluated `Sum`.

- ...:

  Passed to
  [`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md).

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
integral_as_sum(Integral(x^2, list(x, 0, 2)), 4)
#> 21/8
integral_as_sum(Integral(x^2, list(x, 0, 2)), 4, method = "trapezoid")
#> 11/4
principal_value(Integral(1/x, list(x, -1, 1)))
#> 0
```
