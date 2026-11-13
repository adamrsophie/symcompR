# Monotonicity on an interval

`monotonicity_helper(expression, predicate, interval, symbol)` tells
whether `predicate(derivative)` holds on all of `interval`, after
checking that `expression` has no singularity inside it, as SymPy's
helper behind
[`is_increasing()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md)
and its relatives. `predicate` maps the derivative to a relational, such
as `function(d) Ge(d, 0)`.

## Usage

``` r
monotonicity_helper(expression, predicate, interval = Reals, symbol = NULL)
```

## Arguments

- expression:

  An expression.

- predicate:

  A function of the derivative returning a relational.

- interval:

  A set.

- symbol:

  The variable, or `NULL` for the only free symbol.

## Value

`TRUE`, `FALSE` or `NA` when undecided.

## Examples

``` r
x <- Symbol("x")
monotonicity_helper(x^3, function(d) Ge(d, 0))
#> [1] TRUE
```
