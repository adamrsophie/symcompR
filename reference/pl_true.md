# Evaluate a proposition in a model

`pl_true(expr, model, deep)` returns `TRUE` or `FALSE` when the model (a
named list of logical values keyed by symbol names) decides `expr`, and
`NA` (SymPy's `None`) otherwise. With `deep = TRUE` an undecided result
is checked for validity and satisfiability, as SymPy's `pl_true()`.

## Usage

``` r
pl_true(expr, model = NULL, deep = FALSE)
```

## Arguments

- expr:

  A boolean expression.

- model:

  A named list of logical values.

- deep:

  Look further when the model does not decide.

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
A <- Symbol("A"); B <- Symbol("B")
pl_true(Or(A, B), list(A = TRUE))
#> [1] TRUE
pl_true(And(A, B), list(A = TRUE))
#> [1] NA
```
