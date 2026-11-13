# Query predicates of an expression under assumptions

`AssumptionsWrapper(expr, assumptions)` answers `$is_positive`,
`$is_integer` and the other `is_<fact>` queries of `expr`: from the
expression's own assumptions when `assumptions` is `NULL`, and through
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) with
the assumptions otherwise. `is_extended_real(obj, assumptions)` and
`is_extended_nonnegative(obj, assumptions)` are the single queries.

## Usage

``` r
AssumptionsWrapper(expr, assumptions = NULL)

is_extended_real(obj, assumptions = NULL)

is_extended_nonnegative(obj, assumptions = NULL)
```

## Arguments

- expr, obj:

  An expression.

- assumptions:

  Assumptions such as `Q$positive(x)`, or `NULL`.

## Value

`AssumptionsWrapper()` returns an object whose `$is_<fact>` fields are
`TRUE`, `FALSE` or `NA`; the others return one of those.

## Examples

``` r
y <- Symbol("y")
AssumptionsWrapper(y)$is_positive
#> [1] NA
AssumptionsWrapper(y, Q$positive(y))$is_positive
#> [1] TRUE
is_extended_real(y, Q$extended_real(y))
#> [1] TRUE
```
