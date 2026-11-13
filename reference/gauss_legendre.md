# Gaussian quadrature

Nodes and weights of Gaussian quadrature rules to `n_digits` significant
digits, as the functions of `sympy.integrals.quadrature`:
`gauss_legendre()` on `[-1, 1]`, `gauss_laguerre()` with weight
`exp(-x)` on `[0, oo)`, `gauss_hermite()` with weight `exp(-x^2)`,
`gauss_gen_laguerre()` with weight `x^alpha exp(-x)`,
`gauss_chebyshev_t()` and `gauss_chebyshev_u()` of the first and second
kind, `gauss_jacobi()` with weight `(1 - x)^alpha (1 + x)^beta` and
`gauss_lobatto()`, whose nodes include the end points. Each returns
`list(nodes, weights)`. The nodes are the real roots of the orthogonal
polynomial, refined with exact rational arithmetic far below the last
digit before rounding.

## Usage

``` r
gauss_legendre(n, n_digits)

gauss_laguerre(n, n_digits)

gauss_hermite(n, n_digits)

gauss_gen_laguerre(n, alpha, n_digits)

gauss_chebyshev_t(n, n_digits)

gauss_chebyshev_u(n, n_digits)

gauss_jacobi(n, alpha, beta, n_digits)

gauss_lobatto(n, n_digits)
```

## Arguments

- n:

  The number of nodes.

- n_digits:

  Significant digits.

- alpha, beta:

  Parameters of the weight function (rational).

## Value

A list of two lists of numbers.

## Examples

``` r
gauss_legendre(3, 10)
#> [[1]]
#> [[1]][[1]]
#> -0.7745966692
#> 
#> [[1]][[2]]
#> 0
#> 
#> [[1]][[3]]
#> 0.7745966692
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 0.5555555556
#> 
#> [[2]][[2]]
#> 0.8888888889
#> 
#> [[2]][[3]]
#> 0.5555555556
#> 
#> 
gauss_chebyshev_t(2, 10)
#> [[1]]
#> [[1]][[1]]
#> 0.7071067812
#> 
#> [[1]][[2]]
#> -0.7071067812
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 1.570796327
#> 
#> [[2]][[2]]
#> 1.570796327
#> 
#> 
```
