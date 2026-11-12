# Mellin transforms

Ports of SymPy's `mellin_transform` and `inverse_mellin_transform`
(sympy/integrals/transforms.py). `mellin_transform(f, x, s)` returns
`list(F, list(a, b), cond)`: the transform, its fundamental strip
`a < re(s) < b` and an auxiliary condition; with `noconds = TRUE` only
`F`. `inverse_mellin_transform(F, s, x, strip)` inverts `F` over the
strip by rewriting it as a Meijer G function of `x`.

## Usage

``` r
mellin_transform(
  f,
  x,
  s,
  integrator = NULL,
  simplify = TRUE,
  needeval = FALSE,
  noconds = FALSE
)

inverse_mellin_transform(
  F,
  s,
  x,
  strip,
  as_meijerg = FALSE,
  needeval = FALSE,
  noconds = TRUE,
  simplify = TRUE
)
```

## Arguments

- f, F:

  Expressions.

- x, s:

  Variables.

- integrator:

  A function of `(f, x)` computing the integral over `(0, oo)`; by
  default [`integrate()`](https://rdrr.io/r/stats/integrate.html).

- simplify:

  Simplify the result.

- needeval:

  Signal an error instead of returning an unevaluated transform.

- noconds:

  Return only the transform.

- strip:

  `list(a, b)`; either may be `NULL`.

- as_meijerg:

  Leave the inverse as a Meijer G function.

## Value

The transform.

## Examples

``` r
x <- Symbol("x"); s <- Symbol("s")
mellin_transform(exp(-x), x, s)
#> [[1]]
#> gamma(s)
#> 
#> [[2]]
#> [[2]][[1]]
#> 0
#> 
#> [[2]][[2]]
#> oo
#> 
#> 
#> [[3]]
#> True
#> 
inverse_mellin_transform(gamma(s), s, x, list(0, oo))
#> exp(-x)
```
