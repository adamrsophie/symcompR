# Assumption contexts

`AssumptionsContext(...)` is a mutable set of predicates with
`$add(...)`, `$remove(...)` and `$clear()`. `global_assumptions` is the
context that
[`ask()`](https://adamrsophie.github.io/symcompR/reference/Q.md) and
[`refine()`](https://adamrsophie.github.io/symcompR/reference/Q.md)
always consult, as in SymPy. `assuming(..., code)` evaluates its last
argument with the other arguments added to `global_assumptions`, and
restores the previous state afterwards, like SymPy's
`with assuming(...):` block.

## Usage

``` r
AssumptionsContext(...)

global_assumptions

assuming(...)
```

## Arguments

- ...:

  Predicates such as `Q$positive(x)`; for `assuming()`, the last
  argument is the code to evaluate.

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
```
