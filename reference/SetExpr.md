# Set expressions

`SetExpr(set)` is an expression whose value ranges over `set`, like
SymPy's `SetExpr`. The arithmetic operators `+`, `-`, `*`, `/` and `^`
apply interval arithmetic and act elementwise on finite sets; the
resulting set is in `$args[[1]]`.

## Usage

``` r
SetExpr(set)
```

## Arguments

- set:

  A set.

## Value

A `sym` of type `SetExpr`.

## Examples

``` r
SetExpr(Interval(0, 1)) + 1
#> SetExpr(Interval(1, 2))
SetExpr(Interval(1, 2)) * SetExpr(Interval(10, 20))
#> SetExpr(Interval(10, 40))
SetExpr(FiniteSet(1, 2)) + SetExpr(Interval(0, 10))
#> SetExpr(Interval(1, 12))
```
