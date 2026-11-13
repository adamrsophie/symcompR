# Ask whether a proposition holds under assumptions

`ask(proposition, assumptions)` evaluates a proposition built from the
predicates of [Q](https://adamrsophie.github.io/symcompR/reference/Q.md)
to `TRUE`, `FALSE` or `NA` (unknown), using the assumptions given, the
extra assumptions in `...` (joined with `&`) and those of `context`,
which is
[global_assumptions](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
unless another
[`AssumptionsContext()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
is given. This is SymPy's `ask()`: the predicates' handlers answer
first, then a SAT based search over the facts known about the
expressions involved, then the linear real arithmetic solver for
inequalities. Relations such as `x > 0` may appear in the proposition
and the assumptions. A proposition or assumptions that are not boolean
raise a `TypeError`, and inconsistent assumptions a `ValueError`.

## Usage

``` r
ask(proposition, assumptions = TRUE, ..., context = global_assumptions)
```

## Arguments

- proposition:

  A boolean expression such as `Q$positive(x)`, `Q$even(x) | Q$odd(x)`
  or `x > 1`.

- assumptions:

  Assumptions, such as `Q$positive(x) & Q$integer(y)`, or `TRUE` for
  none.

- ...:

  Further assumptions, joined to `assumptions` with `&`.

- context:

  An
  [`AssumptionsContext()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
  whose assumptions are also used.

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
ask(Q$rational(sym_pi))
#> [1] FALSE
ask(Q$even(x * y), Q$even(x) & Q$integer(y))
#> [1] TRUE
ask(Q$prime(4 * x), Q$integer(x))
#> [1] FALSE
ask(Q$odd(3 * x))
#> [1] NA
ask(x > 0, Q$positive(x))
#> [1] TRUE
```
