# Piecewise defined functions

`Piecewise(list(e1, c1), list(e2, c2), ...)` takes the value of the
first expression whose condition holds. Conditions that are `FALSE` are
dropped and a `TRUE` condition ends the list. `piecewise_fold()` moves
arithmetic inside the pieces.

## Usage

``` r
Piecewise(..., evaluate = TRUE)

piecewise_fold(e, evaluate = TRUE)
```

## Arguments

- ...:

  Pairs `list(expr, condition)`.

- evaluate:

  Evaluate the `Piecewise`; with `FALSE`, `Piecewise()` keeps its pairs
  as given. For `piecewise_fold()`, evaluate the folded `Piecewise`.

- e:

  An expression.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
p <- Piecewise(list(-x, x < 0), list(x, TRUE))
subs(p, x, -3)
#> 3
diff(p, x)
#> Piecewise((-1, x < 0), (1, True))
```
