# Satisfiability

`satisfiable(e)` returns a model (a named logical list) making `e` true,
or `FALSE`; with `all_models = TRUE` it returns every model. It uses the
DPLL algorithm on the CNF form, as SymPy does. `valid()` tests whether
`e` is a tautology.

## Usage

``` r
satisfiable(e, all_models = FALSE)

valid(e)
```

## Arguments

- e:

  A boolean expression.

- all_models:

  Return all models.

## Value

A named list, a list of them, or `FALSE`.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b")
satisfiable(And(a, Not(b)))
#> $a
#> [1] TRUE
#> 
#> $b
#> [1] FALSE
#> 
satisfiable(And(a, Not(a)))
#> [1] FALSE
```
