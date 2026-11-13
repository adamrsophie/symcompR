# Decide a proposition by satisfiability

`satask()` is the search over facts that
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) uses
once the handlers of the predicates give no answer: it collects the
facts known about every expression in the proposition and the
assumptions, such as `Q$nonnegative(Abs(x))`, and asks whether the
proposition, or its negation, can hold together with them.

## Usage

``` r
satask(
  proposition,
  assumptions = TRUE,
  use_known_facts = TRUE,
  iterations = Inf,
  early_return = FALSE
)
```

## Arguments

- proposition:

  A boolean expression of predicates.

- assumptions:

  Assumptions, or `TRUE` for none.

- use_known_facts:

  Whether the facts between the predicates are used.

- iterations:

  How many times facts are collected from the expressions found; `Inf`
  until no new expression appears.

- early_return:

  Whether to answer from unit propagation alone when it decides the
  proposition, trusting the assumptions to be consistent.

## Value

`TRUE`, `FALSE` or `NA`. Inconsistent assumptions raise a `ValueError`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
satask(Q$zero(Abs(x)), Q$zero(x))
#> [1] TRUE
satask(Q$zero(x * y), Q$zero(x))
#> [1] TRUE
```
