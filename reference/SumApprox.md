# Approximations of sums and of functions

`SumApprox(bounds, reltol)` is an optimization (see
[`sym_optimize()`](https://adamrsophie.github.io/symcompR/reference/codegen-rewriting.md))
dropping the terms of a sum that are small next to the others, given
bounds on the symbols, as SymPy's `SumApprox`.
`SeriesApprox(bounds, reltol)` replaces functions of one bounded symbol
by their truncated series when that is accurate enough, as SymPy's
`SeriesApprox`.

## Usage

``` r
SumApprox(bounds, reltol, cost_function = NULL, priority = 1)

SeriesApprox(
  bounds,
  reltol,
  max_order = 4,
  n_point_checks = 4,
  cost_function = NULL,
  priority = 1
)
```

## Arguments

- bounds:

  A list of pairs `list(symbol, list(lo, hi))`.

- reltol:

  The relative tolerance.

- cost_function, priority:

  As in
  [`ReplaceOptim()`](https://adamrsophie.github.io/symcompR/reference/codegen-rewriting.md).

- max_order:

  The largest order of the series.

- n_point_checks:

  The number of points (even) where the series is checked.

## Value

An optimization (an R function of an expression).

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
sa <- SumApprox(list(list(x, list(0, 0.001)), list(y, list(100, 1000))), reltol = 0.01)
sym_optimize(x^2 + y^2 + 1, list(sa))
#> y^2
```
