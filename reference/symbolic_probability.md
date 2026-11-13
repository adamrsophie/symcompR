# Unevaluated probability expressions

`Probability()`, `Expectation()`, `Variance()`, `Covariance()`,
`Moment()` and `CentralMoment()` build SymPy's unevaluated objects.
Evaluate them with
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md), or
rewrite them with
[`rewrite()`](https://adamrsophie.github.io/symcompR/reference/rewrite.md)
to `"Expectation"`, `"Probability"`, `"Integral"` or `"Sum"`;
[`expand()`](https://adamrsophie.github.io/symcompR/reference/expand.md)
applies linearity.

## Usage

``` r
Probability(prob, condition = NULL, ...)

Expectation(expr, condition = NULL, ...)

Variance(arg, condition = NULL, ...)

Covariance(arg1, arg2, condition = NULL, evaluate = TRUE, ...)

Moment(X, n, c = 0, condition = NULL, ...)

CentralMoment(X, n, condition = NULL, ...)
```

## Arguments

- prob, condition:

  A condition and an optional given condition.

- ...:

  Ignored.

- expr, arg, arg1, arg2, X:

  Random expressions.

- evaluate:

  For `Covariance()`, `FALSE` keeps the argument order.

- n:

  The order of a moment.

- c:

  The center of a moment.

## Value

A `sym` expression.

## Examples

``` r
X <- Normal("X", 0, 1)
e <- Expectation(X^2)
doit(e)
#> 1
rewrite(Variance(X), "Expectation")
#> -Expectation(X)^2 + Expectation(X^2)
```
