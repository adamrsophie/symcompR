# Satisfiability

`satisfiable(e)` returns a model (a named list of logical values) making
`e` true, or `FALSE`, as SymPy's `satisfiable()`. `algorithm` picks the
method: `"dpll2"` (the default; `"pycosat"`, `"minisat22"` and `"z3"`
fall back to it as in SymPy when those libraries are absent, except that
an explicit `"pycosat"` is an `ImportError`) or `"dpll"` (see
[`dpll_satisfiable()`](https://adamrsophie.github.io/symcompR/reference/dpll_satisfiable.md)).
With `all_models = TRUE` the result is the list of all models, or
`list(FALSE)`. A tautology without symbols has the model
`list(True = TRUE)`. `valid()` tests whether `e` is a tautology.

## Usage

``` r
satisfiable(
  e,
  algorithm = NULL,
  all_models = FALSE,
  minimal = FALSE,
  use_lra_theory = FALSE
)

valid(e)
```

## Arguments

- e:

  A boolean expression.

- algorithm:

  `NULL` or the name of an algorithm.

- all_models:

  Return all models.

- minimal:

  Ask for minimal models (used only by `"minisat22"`).

- use_lra_theory:

  Use the theory of linear real arithmetic, which is not available here.

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
satisfiable(Equivalent(a, b), all_models = TRUE)
#> [[1]]
#> [[1]]$a
#> [1] FALSE
#> 
#> [[1]]$b
#> [1] FALSE
#> 
#> 
#> [[2]]
#> [[2]]$a
#> [1] TRUE
#> 
#> [[2]]$b
#> [1] TRUE
#> 
#> 
```
