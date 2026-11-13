# Array comprehensions

`ArrayComprehension(fun, ...)` is the array of the values of `fun` when
the variables of the limits `list(var, lower, upper)` run over their
ranges; `ArrayComprehensionMap(f, ...)` takes an R function of one
argument applied to the product of the variable values.
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
expands an array with numeric bounds into an explicit array.

## Usage

``` r
ArrayComprehension(fun, ...)

ArrayComprehensionMap(f, ...)
```

## Arguments

- fun:

  An expression.

- ...:

  Limits, each `list(var, lower, upper)` or a
  [`Tuple()`](https://adamrsophie.github.io/symcompR/reference/Tuple.md).

- f:

  An R function of zero or one argument.

## Value

A `sym` array comprehension.

## Examples

``` r
i <- Symbol("i"); j <- Symbol("j")
a <- ArrayComprehension(i * j, list(i, 1, 3), list(j, 2, 4))
doit(a)
#> [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
```
