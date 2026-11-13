# Range of a function

`function_range(f, symbol, domain)` returns the range of the real
function `f` of `symbol` over the real set `domain`, as SymPy's
`function_range()`: the function is followed on each interval of its
continuous domain through its limits at open ends, its values at closed
ends and its critical points.

## Usage

``` r
function_range(f, symbol, domain)
```

## Arguments

- f:

  An expression.

- symbol:

  The variable.

- domain:

  A subset of the reals.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
function_range(sin(x), x, Interval(0, 2 * S_Pi))
#> Interval(-1, 1)
function_range(1 / x, x, Reals)
#> Union(Interval.open(-oo, 0), Interval.open(0, oo))
```
