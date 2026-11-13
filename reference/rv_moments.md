# Moments and summary statistics of random expressions

`variance()`, `standard_deviation()` (alias `std()`), `covariance()`,
`correlation()`, `moment()`, `cmoment()` (central), `smoment()`
(standardized), `skewness()`, `kurtosis()`, `factorial_moment()`,
`coskewness()`, `entropy()` (alias `H()`) and `rv_median()` (SymPy's
[`median()`](https://rdrr.io/r/stats/median.html)), following
`sympy.stats`.

## Usage

``` r
moment(X, n, c = 0, condition = NULL, evaluate = TRUE, ...)

variance(X, condition = NULL, ...)

standard_deviation(X, condition = NULL, ...)

std(X, condition = NULL, ...)

entropy(X, condition = NULL, b = NULL, ...)

H(X, condition = NULL, b = NULL, ...)

covariance(X, Y, condition = NULL, ...)

correlation(X, Y, condition = NULL, ...)

cmoment(X, n, condition = NULL, evaluate = TRUE, ...)

smoment(X, n, condition = NULL, ...)

skewness(X, condition = NULL, ...)

kurtosis(X, condition = NULL, ...)

factorial_moment(X, n, condition = NULL, ...)

rv_median(X, evaluate = TRUE, ...)

coskewness(X, Y, Z, condition = NULL, ...)
```

## Arguments

- X, Y, Z:

  Random expressions.

- n:

  The order of a moment.

- c:

  The center of a raw moment.

- condition:

  An optional condition.

- evaluate:

  `FALSE` returns unevaluated integrals or sums.

- ...:

  Passed on to
  [`E()`](https://adamrsophie.github.io/symcompR/reference/E.md).

- b:

  Logarithm base for the entropy.

## Value

An expression; `rv_median()` returns a set.

## Examples

``` r
X <- Normal("X", 0, 1)
variance(X)
#> 1
moment(X, 4)
#> 3
skewness(Exponential("T", 1))
#> 2
rv_median(Die("D", 6))
#> {3, 4}
```
