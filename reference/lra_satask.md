# Decide a proposition with linear real arithmetic

`lra_satask()` decides propositions made of the relations `Q$eq`,
`Q$ne`, `Q$gt`, `Q$ge`, `Q$lt` and `Q$le` (and the sign predicates such
as `Q$positive`) between linear expressions of real symbols, under
assumptions of the same kind. It is the last step of
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md).
Input it does not handle (nonlinear terms, symbols that are not known to
be real or that are integers, other predicates) raises an
`UnhandledInput` condition.

## Usage

``` r
lra_satask(proposition, assumptions = TRUE)
```

## Arguments

- proposition:

  A boolean expression of relations.

- assumptions:

  Assumptions, or `TRUE` for none.

## Value

`TRUE`, `FALSE` or `NA`; inconsistent assumptions raise a `ValueError`.

## Examples

``` r
x <- Symbol("x", real = TRUE)
lra_satask(Q$gt(x, 0), Q$gt(x, 1))
#> [1] TRUE
lra_satask(Q$eq(x, 1), !Q$ne(x, 0))
#> [1] FALSE
```
