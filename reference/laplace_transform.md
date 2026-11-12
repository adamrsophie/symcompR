# Laplace transforms

Ports of SymPy's `laplace_transform` and `inverse_laplace_transform`
(sympy/integrals/laplace.py). `laplace_transform(f, t, s)` returns
`list(F, a, cond)`: the transform, its half plane of convergence
`re(s) > a` and an auxiliary condition; with `noconds = TRUE` only `F`.
`inverse_laplace_transform(F, s, t)` returns the inverse (by default
without conditions). `expand_dirac_delta()` splits an expression linear
in DiracDelta terms into the remainder and the coefficients.

## Usage

``` r
expand_dirac_delta(expr)

laplace_transform(f, t, s, noconds = FALSE, simplify = FALSE)

inverse_laplace_transform(
  F,
  s,
  t,
  plane = NULL,
  noconds = TRUE,
  simplify = FALSE
)
```

## Arguments

- f, F, expr:

  Expressions.

- t, s:

  Variables.

- noconds:

  Return only the transform.

- simplify:

  Simplify the result.

- plane:

  The half plane of convergence, or `NULL`.

## Value

The transform, or a list with the conditions.

## Examples

``` r
t <- Symbol("t", positive = TRUE); s <- Symbol("s")
laplace_transform(exp(-2*t), t, s)
#> [[1]]
#> 1/(s + 2)
#> 
#> [[2]]
#> -2
#> 
#> [[3]]
#> True
#> 
inverse_laplace_transform(1/(s + 2), s, t)
#> exp(-2*t)
```
