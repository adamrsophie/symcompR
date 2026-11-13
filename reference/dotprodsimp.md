# Simplify a sum of products

Simplification for the expression blowup that occurs while summing
products, as in matrix multiplication. Only algebraic structure is
handled; non-algebraic subexpressions are left alone.

## Usage

``` r
dotprodsimp(expr, withsimp = FALSE)
```

## Arguments

- expr:

  An expression.

- withsimp:

  If `TRUE`, return `list(expr, flag)`, where the flag tells whether
  further simplification may help.

## Value

The simplified expression, or a list as described for `withsimp`.

## Examples

``` r
x <- Symbol("x")
dotprodsimp(x/(x - 1) - 1/(x - 1))
#> 1
```
