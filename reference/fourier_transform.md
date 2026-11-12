# Fourier, sine, cosine and Hankel transforms

Ports of SymPy's `fourier_transform`, `inverse_fourier_transform`,
`sine_transform`, `inverse_sine_transform`, `cosine_transform`,
`inverse_cosine_transform`, `hankel_transform` and
`inverse_hankel_transform` (sympy/integrals/transforms.py), with SymPy's
conventions: the Fourier kernel is `exp(-2*pi*I*x*k)`, the sine and
cosine transforms carry `sqrt(2/pi)`, and the Hankel transform
integrates `f*besselj(nu, k*r)*r` over `(0, oo)`. When the integral is
not found the unevaluated transform is returned.

## Usage

``` r
fourier_transform(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE)

inverse_fourier_transform(
  F,
  k,
  x,
  simplify = TRUE,
  noconds = TRUE,
  needeval = FALSE
)

sine_transform(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE)

inverse_sine_transform(
  F,
  k,
  x,
  simplify = TRUE,
  noconds = TRUE,
  needeval = FALSE
)

cosine_transform(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE)

inverse_cosine_transform(
  F,
  k,
  x,
  simplify = TRUE,
  noconds = TRUE,
  needeval = FALSE
)

hankel_transform(f, r, k, nu, simplify = TRUE, noconds = TRUE)

inverse_hankel_transform(F, k, r, nu, simplify = TRUE, noconds = TRUE)
```

## Arguments

- f, F:

  Expressions.

- x, k, r:

  Variables.

- simplify:

  Simplify the result.

- noconds:

  Drop the convergence condition.

- needeval:

  Signal an error instead of returning an unevaluated transform.

- nu:

  Order of the Hankel transform.

## Value

The transform, or `list(F, cond)` when `noconds = FALSE`.

## Examples

``` r
x <- Symbol("x"); k <- Symbol("k")
fourier_transform(exp(-x^2), x, k)
#> sqrt(pi)*exp(-pi^2*k^2)
```
