# Check that a point is in the domain of a function

`domain_check(f, symbol, p)` returns `FALSE` when `f` or one of its
subexpressions becomes infinite at `symbol = p`, as SymPy's
`domain_check()`.

## Usage

``` r
domain_check(f, symbol, p)
```

## Arguments

- f:

  An expression.

- symbol:

  The variable.

- p:

  The point.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
x <- Symbol("x")
domain_check(1/x, x, 0)
#> [1] FALSE
domain_check(x^2, x, 3)
#> [1] TRUE
```
