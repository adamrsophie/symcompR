# Define an assumption predicate

`Predicate(name)` makes a new predicate, as SymPy's `Predicate("name")`
does: it prints as `Q.name`, can be applied to one expression, and
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md)
answers it only from the assumptions it is given. With `dispatch = TRUE`
the predicate also has its own handler, so
[`register_predicate()`](https://adamrsophie.github.io/symcompR/reference/register_predicate.md)
can attach evaluation rules to it, as for a subclass of SymPy's
`Predicate`. `AppliedPredicate(p, ...)` applies a predicate to
arguments.

## Usage

``` r
Predicate(name, dispatch = FALSE)

AppliedPredicate(predicate, ...)
```

## Arguments

- name:

  The predicate name, a string.

- dispatch:

  Whether the predicate has a handler for
  [`register_predicate()`](https://adamrsophie.github.io/symcompR/reference/register_predicate.md).

- predicate:

  A predicate such as `Q$positive` or one made by `Predicate()`.

- ...:

  The arguments of the applied predicate.

## Value

`Predicate()` returns a callable predicate; `AppliedPredicate()` returns
the applied predicate as a `sym` expression.

## Examples

``` r
x <- Symbol("x")
P <- Predicate("P")
P(x)
#> Q.P(x)
ask(P(x), P(x) & Q$positive(x))
#> [1] TRUE
AppliedPredicate(Q$integer, x)
#> Q.integer(x)
```
