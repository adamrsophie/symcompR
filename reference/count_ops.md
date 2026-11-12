# Count operations in an expression

The number of operations in an expression, as SymPy's `count_ops()`: a
subtraction, a division and a negation count once each, and every named
function or power counts once. With `visual = TRUE` the result is a sum
of operation names such as `2*ADD + SIN`.

## Usage

``` r
count_ops(e, visual = FALSE)
```

## Arguments

- e:

  A `sym` expression.

- visual:

  Return the operations as a symbolic sum.

## Value

An integer, or an expression when `visual = TRUE`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
count_ops(x^2 + sin(x)/y)
#> [1] 4
count_ops(x^2 + sin(x)/y, visual = TRUE)
#> ADD + DIV + POW + SIN
```
