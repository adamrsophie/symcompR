# Expectation, probability and distribution functions

`E(expr, condition)` (alias `expectation()`) is the expected value of a
random expression; `P(condition, given_condition)` (alias
`probability()`) the probability of an event. `rv_density()` is SymPy's
[`density()`](https://rdrr.io/r/stats/density.html): a distribution, a
`Lambda`, or a finite density mapping; `rv_cdf()` is `cdf()`, and
`characteristic_function()`, `moment_generating_function()` and
`rv_quantile()` (SymPy's
[`quantile()`](https://rdrr.io/r/stats/quantile.html)) return `Lambda`
functions. Apply the result to a value with
[`apply_lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md).
For compatibility, a non-Boolean second argument to `rv_density()` and
`rv_cdf()` is a point at which to evaluate.
[`H()`](https://adamrsophie.github.io/symcompR/reference/rv_moments.md)
(alias
[`entropy()`](https://adamrsophie.github.io/symcompR/reference/rv_moments.md))
is the entropy.

## Usage

``` r
E(expr, condition = NULL, numsamples = NULL, evaluate = TRUE, ...)

expectation(expr, condition = NULL, numsamples = NULL, evaluate = TRUE, ...)

P(condition, given_condition = NULL, numsamples = NULL, evaluate = TRUE, ...)

probability(
  condition,
  given_condition = NULL,
  numsamples = NULL,
  evaluate = TRUE,
  ...
)

rv_density(expr, condition = NULL, evaluate = TRUE, numsamples = NULL, ...)

Density(expr, condition = NULL, ...)

rv_cdf(expr, condition = NULL, evaluate = TRUE, ...)

characteristic_function(expr, condition = NULL, evaluate = TRUE, ...)

moment_generating_function(expr, condition = NULL, evaluate = TRUE, ...)

rv_quantile(expr, evaluate = TRUE, ...)
```

## Arguments

- expr:

  A random expression.

- condition:

  An optional condition (Relational or Boolean).

- numsamples:

  Sampling is not supported; must be `NULL`.

- evaluate:

  `FALSE` keeps integrals and sums unevaluated.

- ...:

  Ignored, except by `Density()`: further `list(state, probability)`
  pairs of a quantum density operator
  ([density-operators](https://adamrsophie.github.io/symcompR/reference/density-operators.md)).

- given_condition:

  An optional condition for `P()`.

## Value

An expression, a distribution, a `Lambda` or a mapping.

## Examples

``` r
X <- Normal("X", 0, 1)
E(X^2)
#> 1
P(X > 0)
#> 1/2
rv_density(X)
#> NormalDistribution(0, 1)
apply_lambda(rv_cdf(X), 0)
#> 1/2
D <- Die("D", 6)
E(D)
#> 7/2
P(D > 4)
#> 1/3
```
