# Sums over the roots of a polynomial

`RootSum(expr, func, x)` is the sum of `func(r)` over the roots `r` of
the polynomial `expr` in `x`, counted with multiplicity, as SymPy's
`RootSum`. Linear factors are summed at once; when `func` is a rational
function and `auto` is `TRUE` the sum is computed exactly from the
coefficients (Viete's formulas), otherwise an unevaluated `RootSum`
remains. With `quadratic = TRUE` quadratic factors are summed over their
explicit roots.
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
sums over explicit roots when
[`roots()`](https://adamrsophie.github.io/symcompR/reference/roots.md)
finds them all, and
[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md) sums over
numerical roots.

## Usage

``` r
RootSum(expr, func = NULL, x = NULL, auto = TRUE, quadratic = FALSE)
```

## Arguments

- expr:

  A univariate polynomial.

- func:

  A
  [`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md)
  of one variable, an R function of one argument, or `NULL` for the
  identity.

- x:

  The variable of `expr`; by default its only free symbol.

- auto:

  Evaluate sums of rational functions.

- quadratic:

  Sum over the roots of quadratic factors.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
RootSum(x^3 + x + 1, Lambda(x, x^2))
#> -2
RootSum(x^3 + x + 1, Lambda(x, 1/(x + 2)))
#> 13/9
r <- RootSum(x^5 + x + 3, Lambda(x, log(x + 1)))
r
#> RootSum(w^5 + w + 3, Lambda(w, log(w + 1)))
N(r)
#> 3.14159265358979*I
```
