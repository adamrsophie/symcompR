# Collect common numeric factors

`collect_const(expr, ...)` groups the terms of a sum that share the
numeric factors in `...` (by default every numeric factor that occurs),
as SymPy's `collect_const()`. With `Numbers = FALSE` only non-rational
factors such as `sqrt(3)` are collected.

## Usage

``` r
collect_const(expr, ..., Numbers = TRUE)
```

## Arguments

- expr:

  An expression.

- ...:

  Factors to collect.

- Numbers:

  Also collect rational numbers.

## Value

An expression.

## Examples

``` r
s <- Symbol("s")
collect_const(sqrt(3)*s + sqrt(7)*s + sqrt(3) + sqrt(7))
#> 4.37780211863347*(s + 1)
```
