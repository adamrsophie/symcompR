# Collect terms

`collect(e, syms)` groups the additive terms of `e` by powers of the
symbols (or expressions) in `syms`, as SymPy's `collect()`; it does not
expand `e` first. With `evaluate = FALSE` it returns a named list
mapping the collected powers (by their printed form) to their
coefficients. `rcollect(e, ...)` collects recursively inside the
arguments.

## Usage

``` r
collect(
  e,
  syms,
  func = NULL,
  evaluate = TRUE,
  exact = FALSE,
  distribute_order_term = TRUE
)

rcollect(e, ...)
```

## Arguments

- e:

  An expression.

- syms:

  A symbol, expression or list of them.

- func:

  Optional function applied to each coefficient.

- evaluate:

  Return the collected expression (default) or a list.

- exact:

  Collect only exact powers of the patterns.

- distribute_order_term:

  Distribute an order term over the coefficients.

- ...:

  Symbols for `rcollect()`.

## Value

An expression or a named list.

## Examples

``` r
x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b")
collect(a*x^2 + b*x^2 + a*x - b*x + 1, x)
#> x^2*(a + b) + x*(a - b) + 1
```
