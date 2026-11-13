# Joint random variables

Multivariate distributions of `sympy.stats`: `MultivariateNormal()`,
`MultivariateLaplace()`, `MultivariateT()`, `NormalGamma()`,
`MultivariateBeta()` (alias `Dirichlet()`), `MultivariateEwens()`,
`GeneralizedMultivariateLogGamma()` and its `Omega` form,
`Multinomial()`, `NegativeMultinomial()` and the hand-made `JointRV()`.
`marginal_distribution(rv, ...)` gives the marginal density of the
listed components as a `Lambda`. Index a joint random symbol with
[`rv_getitem()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md).

## Usage

``` r
marginal_distribution(rv, ...)

JointRV(symbol, pdf, set = NULL)

MultivariateNormal(name, mu, sigma)

MultivariateLaplace(name, mu, sigma)

MultivariateT(syms, mu, sigma, v)

NormalGamma(sym, mu, lamda, alpha, beta)

MultivariateBeta(syms, ...)

Dirichlet(syms, ...)

MultivariateEwens(syms, n, theta)

GeneralizedMultivariateLogGamma(syms, delta, v, lamda, mu)

GeneralizedMultivariateLogGammaOmega(syms, omega, v, lamda, mu)

Multinomial(syms, n, ...)

NegativeMultinomial(syms, k0, ...)
```

## Arguments

- rv:

  A joint random symbol.

- ...:

  Components (indices or indexed symbols) for `marginal_distribution()`,
  or further parameters.

- name, syms, symbol, sym:

  The name of the random symbol.

- mu, sigma:

  Mean (location) vector and covariance (shape) matrix, as lists or
  matrices.

- v:

  Degrees of freedom or shape parameters.

- lamda, alpha, beta, delta, omega, theta, n, k0, pdf, set:

  Distribution parameters, following SymPy.

## Value

A random symbol, or a `Lambda` for `marginal_distribution()`.

## Examples

``` r
A <- MultivariateNormal("A", list(1, 2), list(list(1, 0), list(0, 1)))
apply_lambda(rv_density(A), 1, 2)
#> 1/(2*pi)
M <- Multinomial("M", 3, list(Rational(1, 2), Rational(1, 2)))
rv_getitem(M, 0)
#> M[0]
```
