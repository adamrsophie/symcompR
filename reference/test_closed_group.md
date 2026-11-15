# Closed group test of the assumption handlers

SymPy's `test_closed_group(expr, assumptions, key)`: whether every
argument of `expr` satisfies the predicate `key` (fuzzy, with SymPy's
rule that at most one argument may fail), as used by the handlers of
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) for
sums and products.

## Usage

``` r
test_closed_group(expr, assumptions = TRUE, key)
```

## Arguments

- expr:

  An expression.

- assumptions:

  Assumptions, as for
  [`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md).

- key:

  A predicate such as `Q$positive`, or its name.

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
x <- Symbol("x")
test_closed_group(x + 1, Q$positive(x), Q$positive)
#> [1] TRUE
```
