# Assumption predicates

`Q` holds predicate constructors such as `Q$positive(x)` or
`Q$integer(n)`, one for each assumption name. `ask(pred, ...)` answers a
predicate as `TRUE`, `FALSE` or `NA`, optionally under extra assumptions
about symbols, together with those in
[global_assumptions](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md).
`refine(e, ...)` simplifies `e` using such assumptions, for example
`refine(sqrt(x^2), Q$positive(x))` gives `x`.

## Usage

``` r
Q

ask(pred, ...)

refine(e, ...)
```

## Arguments

- pred:

  A predicate, such as `Q$positive(x)`.

- ...:

  Predicates on symbols.

- e:

  An expression.

## Value

`TRUE`, `FALSE`, `NA` or an expression.

## Examples

``` r
x <- Symbol("x")
ask(Q$positive(x^2 + 1), Q$real(x))
#> [1] TRUE
refine(sqrt(x^2), Q$positive(x))
#> x
refine(Abs(x), Q$negative(x))
#> -x
```
