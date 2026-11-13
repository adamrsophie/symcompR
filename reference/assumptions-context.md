# Assumption contexts

`AssumptionsContext(...)` is a mutable set of assumptions with
`$add(...)`, `$remove(...)` and `$clear()`; its `$preds` field lists
them. `global_assumptions` is the context that
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) and
[`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md)
consult by default, as in SymPy. `assuming(..., code)` evaluates its
last argument with the other arguments added to `global_assumptions`,
and restores the previous state afterwards, like SymPy's
`with assuming(...):` block. Removing an assumption that is not in the
context raises a `KeyError`, as for a Python set.

## Usage

``` r
AssumptionsContext(...)

global_assumptions

assuming(...)
```

## Arguments

- ...:

  Assumptions such as `Q$positive(x)` or `x > 0`; for `assuming()`, the
  last argument is the code to evaluate.

## Value

`AssumptionsContext()` returns a context; `assuming()` returns the value
of the code.

## Examples

``` r
x <- Symbol("x")
assuming(Q$positive(x), ask(Q$real(x)))
#> [1] TRUE
assuming(Q$positive(x), refine(sqrt(x^2)))
#> x
ask(Q$real(x))
#> [1] NA
ctx <- AssumptionsContext(Q$integer(x))
ctx
#> AssumptionsContext({Q.integer(x)})
ask(Q$rational(x), context = ctx)
#> [1] TRUE
```
