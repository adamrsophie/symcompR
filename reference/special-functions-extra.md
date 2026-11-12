# More special functions

`Rem(p, q)` is the remainder with the sign of `p`. `andre(n)` gives the
Euler zigzag numbers. `multigamma(x, p)` is the multivariate gamma
function, `stieltjes(n, a)` the Stieltjes constants, `riemann_xi(s)`
Riemann's xi function, `owens_t(h, a)` Owen's T function,
`erf2inv(x, y)` the inverse of
[`erf2()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
in its second argument, `marcumq(m, a, b)` the Marcum Q function.
`hyper(ap, bq, z)` is the generalized hypergeometric function and
`meijerg(a, b, z)` the Meijer G function (with `a` and `b` lists of two
lists), `appellf1()` Appell's F1. `hn1()` and `hn2()` are the spherical
Hankel functions and `jn_zeros(n, k)` the first `k` zeros of
[`jn()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md).
`SingularityFunction(x, a, n)` is the Macaulay bracket. `betainc()` and
`betainc_regularized()` are the generalized incomplete beta function
between two limits. `jtheta(n, z, q)` is the Jacobi theta function,
`jacobi_normalized()` the orthonormal Jacobi polynomial, `Ynm_c()` and
`Znm()` the conjugate and real spherical harmonics. `mathieus()`,
`mathieuc()`, `mathieusprime()` and `mathieucprime()` are symbolic only.
`exp_polar()`, `polar_lift()`, `unpolarify()`, `periodic_argument()`,
`unbranched_argument()` and `principal_branch()` work with numbers on
the Riemann surface of the logarithm. `piecewise_exclusive()` makes the
conditions of a
[`Piecewise()`](https://adamrsophie.github.io/symcompR/reference/Piecewise.md)
mutually exclusive.

## Usage

``` r
Rem(...)

andre(...)

multigamma(...)

stieltjes(...)

riemann_xi(...)

owens_t(...)

erf2inv(...)

marcumq(...)

appellf1(...)

hn1(...)

hn2(...)

jn_zeros(n, k)

SingularityFunction(...)

betainc(...)

betainc_regularized(...)

jtheta(...)

jacobi_normalized(n, a, b, x)

Ynm_c(n, m, theta, phi)

Znm(n, m, theta, phi)

mathieus(...)

mathieuc(...)

mathieusprime(...)

mathieucprime(...)

exp_polar(...)

polar_lift(...)

unpolarify(e, subs = NULL, exponents_only = FALSE)

unbranched_argument(z)

periodic_argument(z, period)

principal_branch(z, period)

piecewise_exclusive(e)

hyper(ap, bq, z)

meijerg(...)
```

## Arguments

- ...:

  Arguments of the symbolic function, in SymPy's order.

- n, x, a, b, m, z, k, theta, phi, period, e:

  Arguments.

- subs:

  For `unpolarify()`: a named list of values substituted for symbols.

- exponents_only:

  For `unpolarify()`: only unpolarify exponents.

- ap, bq:

  Lists of upper and lower parameters.

## Value

An expression, or a list of Floats for `jn_zeros()`.

## Examples

``` r
x <- Symbol("x")
Rem(-7, 3)
#> -1
andre(8)
#> 1385
multigamma(x, 2)
#> sqrt(pi)*gamma(x)*gamma(x - 1/2)
SingularityFunction(5, 2, 3)
#> 27
N(hyper(list(1, 2), list(3), Rational(1, 2)))
#> 1.54517744447956
jn_zeros(1, 2)
#> [[1]]
#> 4.49340945790906
#> 
#> [[2]]
#> 7.72525183693771
#> 
```
