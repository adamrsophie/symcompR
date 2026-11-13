# Special functions

Symbolic special functions with special values, derivatives and numeric
evaluation: the gamma family (`gamma`, `loggamma`, `polygamma`,
`digamma`, `trigamma`, `beta`, `lowergamma`, `uppergamma`), factorials
(`factorial`, `factorial2`, `subfactorial`, `binomial`,
`RisingFactorial`, `FallingFactorial`), combinatorial numbers
(`bernoulli`, `fibonacci`, `lucas`, `tribonacci`, `catalan`, `bell`,
`harmonic`, `euler`, `genocchi`, `partition`, `motzkin`), error and
exponential integrals (`erf`, `erfc`, `erfi`, `erfinv`, `erfcinv`, `Ei`,
`li`, `Si`, `Ci`, `Shi`, `Chi`, `expint`, `fresnels`, `fresnelc`), zeta
functions (`zeta`, `dirichlet_eta`, `polylog`, `lerchphi`), Bessel and
Airy functions, elliptic integrals and orthogonal polynomials
(`legendre`, `assoc_legendre`, `chebyshevt`, `chebyshevu`, `hermite`,
`hermite_prob`, `laguerre`, `assoc_laguerre`, `gegenbauer`, `jacobi`,
`Ynm`).

## Usage

``` r
loggamma(...)

polygamma(...)

digamma_sym(x)

lowergamma(...)

uppergamma(...)

factorial2(...)

subfactorial(...)

RisingFactorial(...)

FallingFactorial(...)

bernoulli(...)

fibonacci(...)

lucas(...)

tribonacci(...)

catalan(...)

bell(n, k_sym = NULL, symbols = NULL)

harmonic(...)

euler(...)

genocchi(...)

partition(...)

motzkin(...)

erf(...)

erfc(...)

erfi(...)

erfinv(...)

erfcinv(...)

Ei(...)

li(...)

Si(...)

Ci(...)

Shi(...)

Chi(...)

expint(...)

fresnels(...)

fresnelc(...)

zeta(...)

dirichlet_eta(...)

polylog(...)

lerchphi(...)

besselj(...)

bessely(...)

besseli(...)

besselk(...)

hankel1(...)

hankel2(...)

jn(...)

yn(...)

airyai(...)

airybi(...)

airyaiprime(...)

airybiprime(...)

elliptic_k(...)

elliptic_e(...)

elliptic_f(...)

elliptic_pi(...)

legendre(...)

assoc_legendre(...)

chebyshevt(...)

chebyshevu(...)

hermite(...)

hermite_prob(...)

laguerre(...)

assoc_laguerre(...)

gegenbauer(...)

jacobi(...)

Ynm(...)

sym_beta(a, b = a)

sym_factorial(x)

sym_binomial(n, k, ...)
```

## Arguments

- ...:

  Arguments.

- x, a, b, n, k:

  Arguments.

- k_sym, symbols:

  For `bell()`: the second index of the incomplete Bell polynomial and
  the list of its symbols.

## Value

A `sym` expression.

## Details

`sym_factorial()` and `sym_beta()` dispatch to the base R functions for
numeric input.

## Examples

``` r
sym_factorial(S(20))
#> 2432902008176640000
sym_binomial(Symbol("n"), 2)
#> binomial(n, 2)
zeta(S(2))
#> pi^2/6
legendre(3, Symbol("x"))
#> 5*x^3/2 - 3*x/2
```
