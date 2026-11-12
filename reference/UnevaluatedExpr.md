# Expression that does not evaluate

`UnevaluatedExpr(e)` holds `e` as an atom in arithmetic, so
`x * UnevaluatedExpr(1 / x)` stays as it is.
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
releases it.

## Usage

``` r
UnevaluatedExpr(e)
```

## Arguments

- e:

  An expression.

## Value

A `sym` UnevaluatedExpr.

## Examples

``` r
x <- Symbol("x")
x * UnevaluatedExpr(1 / x)
#> x*1/x
doit(x * UnevaluatedExpr(1 / x))
#> 1
```
