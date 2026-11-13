# Conditional random expressions and events

`given(expr, condition)` re-expresses `expr` on the probability space
restricted by `condition`. `where(condition)` is the domain on which
`condition` holds. `independent()` and `dependent()` test whether two
random expressions depend on each other. `rs_swap(a, b)` pairs random
symbols with the same underlying symbol.

## Usage

``` r
given(expr, condition = NULL, ...)

rs_swap(a, b)

dependent(a, b)

independent(a, b)

where(condition, given_condition = NULL, ...)
```

## Arguments

- expr:

  A random expression.

- condition:

  A relational or Boolean in random symbols, or a random symbol (treated
  as fixed).

- ...:

  Ignored.

- a, b:

  Random expressions or lists of random symbols.

- given_condition:

  An optional further condition.

## Value

An expression, a domain, a logical, or a list of pairs.

## Examples

``` r
X <- Normal("X", 0, 1)
Y <- given(X, X > 0)
E(Y)
#> sqrt(2)/sqrt(pi)
where(X^2 <= 1)
#> SingleContinuousDomain(X, Interval(-1, 1))
D <- Die("D", 6)
rv_density(given(D, D > 3))
#> {4: 1/3, 5: 1/3, 6: 1/3}
```
