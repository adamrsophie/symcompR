# Formal power series

`fps(f, x, x0, dir)` is SymPy's `fps()`: the formal power series of `f`
about `x0`, as a `FormalPowerSeries` (`f` itself when none is found).
`compute_fps()` returns the coefficient sequence, the power sequence and
the independent terms. `fps_truncate()` (also
[`truncate()`](https://rdrr.io/r/base/seek.html)) and `fps_polynomial()`
give the truncated series; `fps_term()` the term of a given order.

## Usage

``` r
fps_integrate(s, x = NULL)

fps_product(s, other, n = 6)

fps_compose(s, other, n = 6)

fps_inverse(s, n = 6)

fps(
  f,
  x = NULL,
  x0 = 0,
  dir = 1,
  hyper = TRUE,
  order = 4,
  rational = TRUE,
  full = FALSE
)

compute_fps(
  f,
  x,
  x0 = 0,
  dir = 1,
  hyper = TRUE,
  order = 4,
  rational = TRUE,
  full = FALSE
)

fps_term(s, pt)

fps_polynomial(s, n = 6)

fps_truncate(s, n = 6)
```

## Arguments

- s:

  A `FormalPowerSeries`.

- x:

  The variable.

- other:

  A second `FormalPowerSeries`.

- n:

  Number of terms.

- f:

  An expression.

- x0:

  The expansion point.

- dir:

  `1` or `"+"`, `-1` or `"-"`.

- hyper, order, rational, full:

  Options of SymPy's algorithms.

- pt:

  The index of the term.

## Value

A `FormalPowerSeries` or an expression.

## Examples

``` r
x <- Symbol("x")
fps_truncate(fps(exp(x), x))
#> 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O(x^6)
fps_truncate(fps(log(1 + x), x), 4)
#> x - x^2/2 + x^3/3 + O(x^4)
```
