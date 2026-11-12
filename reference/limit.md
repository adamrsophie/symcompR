# Limits

`limit(e, x, x0, dir)` computes the limit of `e` as `x` approaches `x0`,
as SymPy's `limit()` (`Limit(e, x, x0, dir).doit(deep = FALSE)`). `dir`
is `"+"` (default), `"-"` or `"+-"` (two sided). `x0` may be `oo` or
`-oo`. `Limit()` builds an unevaluated limit. A limit that cannot be
computed is returned unevaluated.

## Usage

``` r
limit(e, x, x0, dir = "+")

Limit(e, x, x0, dir = "+")
```

## Arguments

- e:

  An expression.

- x:

  The variable.

- x0:

  The limit point.

- dir:

  Direction.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
limit(sin(x) / x, x, 0)
#> 1
limit((1 + 1/x)^x, x, oo)
#> E
limit(x^2 * exp(-x), x, oo)
#> 0
```
