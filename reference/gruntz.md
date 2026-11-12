# Limits by the Gruntz algorithm

`gruntz(e, z, z0, dir)` computes the limit of `e` as `z` approaches `z0`
with the Gruntz algorithm, as SymPy's `gruntz()`. `z0` may be `oo` or
`-oo`; `dir` (`"+"` or `"-"`) gives the side for finite `z0`.
[`limit()`](https://adamrsophie.github.io/symcompR/reference/limit.md)
tries this algorithm first.

## Usage

``` r
gruntz(e, z, z0, dir = "+")
```

## Arguments

- e:

  An expression.

- z:

  The variable, a symbol.

- z0:

  The limit point.

- dir:

  Direction, `"+"` or `"-"`.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
gruntz(exp(x) * (exp(1/x - exp(-x)) - exp(1/x)), x, oo)
#> -1
gruntz((3^x + 5^x)^(1/x), x, oo)
#> 5
gruntz(x^x, x, 0)
#> 1
```
