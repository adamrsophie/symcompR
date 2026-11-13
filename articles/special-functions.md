# Special functions

symcompR includes the elementary and special functions of SymPy. Each
one evaluates exactly at special points, stays symbolic elsewhere, knows
its derivative and series expansion, can be rewritten in terms of
related functions, and evaluates numerically with
[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md).

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
a <- Symbol("a"); b <- Symbol("b"); nu <- Symbol("nu")
n <- Symbol("n", integer = TRUE)
xp <- Symbol("xp", positive = TRUE)
f <- Function("f")
```

A note on arguments: when a function also exists in base R (such as
[`exp()`](https://rdrr.io/r/base/Log.html) or
[`gamma()`](https://rdrr.io/r/base/Special.html)), a plain R number
gives a plain R number. Wrap numbers in
[`S()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md),
[`Integer()`](https://adamrsophie.github.io/symcompR/reference/Integer.md),
or
[`Rational()`](https://adamrsophie.github.io/symcompR/reference/Integer.md)
to get exact results.

``` r

gamma(5)
#> [1] 24
gamma(S(5))
#> 24
sym_atan2(S(1), S(1))
#> pi/4
```

## Elementary functions

``` r

acos(Rational(1, 2))
#> pi/3
asin(S(1))
#> pi/2
atan(oo)
#> pi/2
acosh(S(1))
#> 0
sinc(0)
#> 1
log(S(8), 2)
#> 3
exp(log(x))
#> x
log(exp(x))
#> log(exp(x))
log(exp(xp))
#> xp
root(8, 3)
#> 2
cbrt(-8)
#> 2*(-1)^(1/3)
real_root(-8, 3)
#> -2
```

`log(exp(x))` is not simplified for a general symbol, because the
identity fails for complex `x`; it holds for the positive symbol `xp`.
Note also that `cbrt(-8)` is the principal cube root, while
[`real_root()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md)
gives the real one.

Rewriting changes the form of an expression:

``` r

rewrite(sinh(x), "exp")
#> exp(x)/2 - exp(-x)/2
rewrite(cos(x), "exp")
#> exp(I*x)/2 + exp(-I*x)/2
rewrite(tan(x), "sin")
#> 2*sin(x)^2/sin(2*x)
```

### Complex parts and rounding

``` r

Abs(-3 + 4 * sym_I)
#> 5
re(exp(sym_I * x))
#> exp(-im(x))*cos(re(x))
arg(1 + sym_I)
#> pi/4
conjugate(2 + 3 * sym_I)
#> 2 - 3*I
sign(S(-3))
#> -1
floor(Rational(7, 2))
#> 3
ceiling(-Rational(7, 2))
#> -3
frac(Rational(7, 2))
#> 1/2
Max(1, x, 3, x)
#> Max(3, x)
Min(xp, -1)
#> -1
```

## Gamma and related functions

``` r

gamma(Rational(1, 2))
#> sqrt(pi)
gamma(Rational(-3, 2))
#> 4*sqrt(pi)/3
diff(gamma(x), x)
#> gamma(x)*polygamma(0, x)
rewrite(gamma(x + 1), "factorial")
#> factorial(x)
loggamma(S(3))
#> log(2)
polygamma(0, 1)
#> -EulerGamma
polygamma(1, 1)
#> pi^2/6
diff(digamma_sym(x), x)
#> polygamma(1, x)
sym_beta(S(2), S(3))
#> 1/12
rewrite(sym_beta(x, y), "gamma")
#> gamma(x)*gamma(y)/gamma(x + y)
lowergamma(1, x)
#> 1 - exp(-x)
uppergamma(1, x)
#> exp(-x)
uppergamma(Rational(1, 2), x)
#> sqrt(pi)*erfc(sqrt(x))
```

Factorials and binomial coefficients:

``` r

sym_factorial(S(5))
#> 120
factorial2(7)
#> 105
RisingFactorial(x, 3)
#> x*(x + 1)*(x + 2)
FallingFactorial(x, 3)
#> x*(x - 2)*(x - 1)
sym_binomial(10, 3)
#> 120
```

## Zeta and polylogarithm

``` r

zeta(2)
#> pi^2/6
zeta(4)
#> pi^4/90
zeta(-1)
#> -1/12
zeta(0)
#> -1/2
zeta(3)
#> zeta(3)
dirichlet_eta(1)
#> log(2)
polylog(2, Rational(1, 2))
#> -log(2)^2/2 + pi^2/12
expand_func(polylog(1, x))
#> -log(1 - x)
lerchphi(1, 2, 1)
#> lerchphi(1, 2, 1)
```

`zeta(3)` (Apéry’s constant) has no known closed form, so it stays
unevaluated.

## Error functions and integrals

`erf2(x, y)` is the two argument error function \mathrm{erf}(y) -
\mathrm{erf}(x).

``` r

erf(0)
#> 0
erf(oo)
#> 1
erf(-x)
#> -erf(x)
diff(erf(x), x)
#> 2*exp(-x^2)/sqrt(pi)
rewrite(erfc(x), "erf")
#> 1 - erf(x)
erfinv(0)
#> 0
erf2(0, x)
#> erf(x)
series(erf(x), x, 0, 6)
#> 2*x/sqrt(pi) - 2*x^3/(3*sqrt(pi)) + x^5/(5*sqrt(pi)) + O(x^6)
```

The exponential, logarithmic, sine, cosine, and Fresnel integrals.
[`Li()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
is the offset logarithmic integral, \mathrm{li}(x) - \mathrm{li}(2),
which vanishes at 2:

``` r

diff(Ei(x), x)
#> exp(x)/x
diff(li(x), x)
#> 1/log(x)
Li(2)
#> 0
diff(Si(x), x)
#> sin(x)/x
diff(Ci(x), x)
#> cos(x)/x
Si(oo)
#> pi/2
rewrite(expint(1, x), "Ei")
#> -Ei(x*exp_polar(I*pi)) + I*pi
E1(x)
#> expint(1, x)
fresnels(oo)
#> 1/2
diff(fresnels(x), x)
#> sin(pi*x^2/2)
```

## Bessel and Airy functions

``` r

diff(besselj(nu, x), x)
#> besselj(nu - 1, x)/2 - besselj(nu + 1, x)/2
besselj(0, 0)
#> 1
besseli(0, 0)
#> 1
besselj(-2, x)
#> besselj(2, x)
diff(bessely(0, x), x)
#> -bessely(1, x)
diff(besselk(nu, x), x)
#> -besselk(nu - 1, x)/2 - besselk(nu + 1, x)/2
rewrite(hankel1(nu, x), "besselj")
#> hankel1(nu, x)
```

Bessel functions of half integer order and the spherical Bessel
functions are elementary:

``` r

expand_func(besselj(Rational(1, 2), x))
#> besselj(1/2, x)
expand_func(jn(0, x))
#> sin(x)/x
expand_func(jn(1, x))
#> -cos(x)/x + sin(x)/x^2
```

The Airy function \mathrm{Ai} solves y'' = x y:

``` r

airyai(0)
#> 3^(1/3)/(3*gamma(2/3))
diff(airyai(x), x, 2)
#> x*airyai(x)
```

## Elliptic integrals

``` r

elliptic_k(0)
#> pi/2
elliptic_e(0)
#> pi/2
diff(elliptic_k(x), x)
#> (-(1 - x)*elliptic_k(x) + elliptic_e(x))/(2*x*(1 - x))
series(elliptic_k(x), x, 0, 3)
#> pi/2 + pi*x/8 + 9*pi*x^2/128 + O(x^3)
```

## Hypergeometric functions

`hyper(ap, bq, z)` is the generalized hypergeometric function {}\_pF_q,
and
[`meijerg()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
is the Meijer G function.
[`hyperexpand()`](https://adamrsophie.github.io/symcompR/reference/hyperexpand.md)
rewrites them in terms of named functions when possible.

``` r

hyper(list(1, 1), list(2), z)
#> hyper((1, 1), (2,), z)
hyperexpand(hyper(list(1, 1), list(2), z))
#> -log(1 - z)/z
hyperexpand(hyper(list(), list(), z))
#> exp(z)
hyperexpand(hyper(list(Rational(1, 2)), list(Rational(3, 2)), -z^2))
#> sqrt(pi)*erf(z)/(2*z)
G <- meijerg(list(list(1, 1), list()), list(list(1), list(0)), z)
G
#> meijerg(((1, 1), ()), ((1,), (0,)), z)
hyperexpand(G)
#> log(z + 1)
```

## Orthogonal polynomials

``` r

legendre(3, x)
#> 5*x^3/2 - 3*x/2
assoc_legendre(2, 1, x)
#> -3*x*sqrt(1 - x^2)
chebyshevt(4, x)
#> 8*x^4 - 8*x^2 + 1
chebyshevu(3, x)
#> 8*x^3 - 4*x
hermite(3, x)
#> 8*x^3 - 12*x
laguerre(3, x)
#> -x^3/6 + 3*x^2/2 - 3*x + 1
diff(legendre(n, x), x)
#> n*(x*legendre(n, x) - legendre(n - 1, x))/(x^2 - 1)
```

Orthogonality of the Legendre polynomials on \[-1, 1\]:

``` r

sym_integrate(legendre(2, x) * legendre(3, x), list(x, -1, 1))
#> 0
```

Spherical harmonics:

``` r

expand_func(Ynm(1, 0, x, y))
#> sqrt(3)*cos(x)/(2*sqrt(pi))
```

## Lambert W, Heaviside, and Dirac delta

``` r

LambertW(0)
#> 0
LambertW(-1 / sym_E)
#> -1
LambertW(sym_E)
#> 1
diff(LambertW(x), x)
#> LambertW(x)/(x*(LambertW(x) + 1))
Heaviside(0)
#> 1/2
diff(Heaviside(x), x)
#> DiracDelta(x)
sym_integrate(DiracDelta(x) * cos(x), list(x, -oo, oo))
#> 1
sym_integrate(DiracDelta(x - 2) * x^2, list(x, 0, 5))
#> 4
KroneckerDelta(1, 2)
#> 0
```

`SingularityFunction(x, a, n)` is the Macaulay bracket \langle x - a
\rangle^n used in beam deflection problems.

``` r

diff(SingularityFunction(x, 3, 2), x)
#> 2*SingularityFunction(x, 3, 1)
```

## Piecewise functions

[`Piecewise()`](https://adamrsophie.github.io/symcompR/reference/Piecewise.md)
takes pairs `list(expression, condition)`; the first condition that
holds selects the value, and `TRUE` acts as “otherwise”.

``` r

p <- Piecewise(list(x^2, Lt(x, 0)), list(x, TRUE))
p
#> Piecewise((x^2, x < 0), (x, True))
subs(p, x, -2)
#> 4
diff(p, x)
#> Piecewise((2*x, x < 0), (1, True))
sym_integrate(p, list(x, -1, 1))
#> 5/6
piecewise_fold(Piecewise(list(1, Lt(x, 0)), list(2, TRUE)) + 1)
#> Piecewise((2, x < 0), (3, True))
rewrite(Heaviside(x), "Piecewise")
#> Piecewise((0, x < 0), (1/2, Eq(x, 0)), (1, True))
```

`bspline_basis(d, knots, i, x)` is the i-th B-spline of degree `d`
(counting from zero), as a piecewise polynomial.

``` r

bspline_basis(1, list(0, 1, 2), 0, x)
#> Piecewise((x, (x >= 0) & (x <= 1)), (2 - x, (x >= 1) & (x <= 2)), (0, True))
bspline_basis(2, list(0, 1, 2, 3), 0, x)
#> Piecewise((x^2/2, (x >= 0) & (x <= 1)), (-x^2 + 3*x - 3/2, (x >= 1) & (x <= 2)), (x^2/2 - 3*x + 9/2, (x >= 2) & (x <= 3)), (0, True))
bspline_basis_set(2, list(0, 1, 2, 3), x)[[1]]
#> Piecewise((x^2/2, (x >= 0) & (x <= 1)), (-x^2 + 3*x - 3/2, (x >= 1) & (x <= 2)), (x^2/2 - 3*x + 9/2, (x >= 2) & (x <= 3)), (0, True))
interpolating_spline(1, x, list(1, 2, 4), list(3, 4, 7))
#> Piecewise((x + 2, (x >= 1) & (x <= 2)), (3*x/2 + 1, (x >= 2) & (x <= 4)))
```

## Numerical evaluation

`N(expr, n)` evaluates to `n` significant digits.

``` r

N(gamma(Rational(1, 3)), 30)
#> 2.67893853470774763365569294097
N(zeta(3), 30)
#> 1.20205690315959428539973816151
N(erf(1), 25)
#> 0.8427007929497148693412206
N(elliptic_k(Rational(1, 2)), 20)
#> 1.8540746773013719184
N(EulerGamma, 30)
#> 0.577215664901532860606512090082
N(Catalan, 30)
#> 0.915965594177219015054603514932
N(legendre(5, Rational(3, 10)), 20)
#> 0.34538625000000000000
```

Some functions evaluate in double precision only, so they are limited to
15 significant digits:

``` r

N(besselj(0, 1))
#> 0.765197686557967
N(airyai(1))
#> 0.135292416312881
N(LambertW(1))
#> 0.567143290409784
N(Si(1))
#> 0.946083070367183
N(li(10))
#> 6.16559950478730
N(hyper(list(1, 2), list(3), Rational(1, 2)))
#> 1.54517744447956
```

## More elementary functions

The reciprocal trigonometric and hyperbolic functions and their
inverses:

``` r

a <- Symbol("a"); b <- Symbol("b")
csc(sym_pi / 6)
#> 2
acsc(2)
#> pi/6
asec(2)
#> pi/3
acot(1)
#> pi/4
diff(acot(x), x)
#> -1/(x^2 + 1)
sech(0)
#> 1
asech(1)
#> 0
rewrite(coth(x), "exp")
#> (exp(x) + exp(-x))/(exp(x) - exp(-x))
rewrite(coth(x), "tanh")
#> 1/tanh(x)
rewrite(acoth(2), "log")
#> log(3/2)/2 + log(2)/2
rewrite(acsch(1), "log")
#> log(1 + sqrt(2))
diff(csch(x), x)
#> -coth(x)*csch(x)
im(3 + 4 * sym_I)
#> 4
```

## More special functions

Further Bessel, Airy, and related functions.
[`yn()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md),
[`hn1()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md),
and
[`hn2()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
are spherical Bessel and Hankel functions, and
[`jn_zeros()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
finds zeros of the spherical Bessel functions numerically.

``` r

rewrite(yn(0, x), "sin")
#> yn(0, x)
expand_func(yn(1, x))
#> -sin(x)/x - cos(x)/x^2
rewrite(hankel2(0, x), "besselj")
#> hankel2(0, x)
rewrite(hn1(0, x), "sin")
#> hn1(0, x)
expand_func(hn2(0, x))
#> sin(x)/x + I*cos(x)/x
lapply(jn_zeros(0, 2), N)
#> [[1]]
#> 3.14159265358979
#> 
#> [[2]]
#> 6.28318530717959
N(airybi(0))
#> 0.614926627446001
N(airyaiprime(0))
#> -0.258819403792807
airybiprime(0)
#> 3^(1/6)/gamma(1/3)
diff(airybi(x), x)
#> airybiprime(x)
```

Integrals and inverses: hyperbolic sine and cosine integrals, the
Fresnel cosine integral, the imaginary error function, and inverse error
functions.

``` r

N(Chi(1))
#> 0.837866940980208
N(Shi(1))
#> 1.05725087537573
diff(Shi(x), x)
#> sinh(x)/x
N(fresnelc(1))
#> 0.779893400376823
limit(fresnelc(x), x, oo)
#> 1/2
rewrite(erfi(x), "erf")
#> -I*erf(I*x)
N(erfi(1))
#> 1.65042575879754
erfcinv(1)
#> 0
N(erfcinv(Rational(1, 2)))
#> 0.476936276204470
erf2inv(0, y)
#> erfinv(y)
```

Elliptic integrals of the first and third kind:

``` r

N(elliptic_f(Rational(1, 2), Rational(1, 3)))
#> 0.506847756265431
elliptic_f(0, x)
#> 0
N(elliptic_pi(Rational(1, 3), Rational(1, 2)))
#> 2.31079499075428
```

Orthogonal polynomials with numeric parameters, and integer sequences:

``` r

assoc_laguerre(2, 1, x)
#> x^2/2 - 3*x + 3
gegenbauer(3, 2, x)
#> 32*x^3 - 12*x
jacobi(2, 1, 1, x)
#> 15*x^2/4 - 3/4
jacobi_normalized(2, 1, 1, x)
#> sqrt(42)*(15*x^2/4 - 3/4)/6
hermite_prob(3, x)
#> x^3 - 3*x
genocchi(6)
#> -3
lapply(0:7, tribonacci)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 1
#> 
#> [[4]]
#> 2
#> 
#> [[5]]
#> 4
#> 
#> [[6]]
#> 7
#> 
#> [[7]]
#> 13
#> 
#> [[8]]
#> 24
lapply(0:7, andre)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 1
#> 
#> [[4]]
#> 2
#> 
#> [[5]]
#> 5
#> 
#> [[6]]
#> 16
#> 
#> [[7]]
#> 61
#> 
#> [[8]]
#> 272
```

More functions from `sympy.functions.special`: the Riemann \xi function,
Stieltjes constants, the multivariate gamma function, Owen’s T function,
the Marcum Q function, Jacobi theta functions, the Appell F_1 function,
and the incomplete beta functions.

``` r

riemann_xi(2)
#> pi/6
stieltjes(0)
#> EulerGamma
multigamma(3, 2)
#> 3*pi/2
owens_t(0, a)
#> atan(a)/(2*pi)
N(owens_t(1, Rational(1, 2)))
#> 0.0430646911207854
marcumq(1, 0, b)
#> exp(-b^2/2)
N(jtheta(3, 1, Rational(1, 2)))
#> 0.505893885730485
N(appellf1(1, 1, 1, 2, Rational(1, 2), Rational(1, 3)))
#> 1.72609243471069
betainc(1, 1, 0, x)
#> betainc(1, 1, 0, x)
N(betainc_regularized(2, 3, 0, Rational(1, 2)))
#> 0.687500000000000
Rem(7, -3)
#> 1
Rem(x^2, x)
#> Rem(x^2, x)
```

The Mathieu functions are available as symbolic functions:

``` r

mathieus(a, b, x)
#> mathieus(a, b, x)
mathieuc(a, b, x)
#> mathieuc(a, b, x)
mathieusprime(a, b, x)
#> mathieusprime(a, b, x)
mathieucprime(a, b, x)
#> mathieucprime(a, b, x)
```

[`Ynm_c()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
is the complex conjugate of the spherical harmonic and
[`Znm()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
the real spherical harmonic:

``` r

expand_func(Ynm_c(1, 0, a, b))
#> sqrt(3)*cos(a)/(2*sqrt(pi))
expand_func(Znm(1, 0, a, b))
#> sqrt(3)*cos(a)/(2*sqrt(pi))
```

## Polar numbers and branches

[`exp_polar()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
is an exponential that is not reduced modulo 2\pi i, which keeps track
of the branch on the Riemann surface of the logarithm.
[`polar_lift()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md),
[`principal_branch()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md),
[`periodic_argument()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md),
and
[`unbranched_argument()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
work with such numbers, and
[`unpolarify()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
and
[`polarify()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
convert expressions.

``` r

exp_polar(0)
#> exp_polar(0)
unpolarify(exp_polar(x))
#> exp(x)
polar_lift(-1)
#> exp_polar(I*pi)
polar_lift(2 * sym_I)
#> 2*exp_polar(I*pi/2)
principal_branch(-2, 2 * sym_pi)
#> 2*exp_polar(I*pi)
periodic_argument(exp_polar(3 * sym_I * sym_pi), 2 * sym_pi)
#> pi
unbranched_argument(exp_polar(3 * sym_I * sym_pi))
#> 3*pi
polarify(exp(x) * y, subs = FALSE)
#> y*exp(x)
```

[`piecewise_exclusive()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
rewrites a `Piecewise` so that its conditions do not overlap:

``` r

piecewise_exclusive(Piecewise(list(1, x < 0), list(2, x < 1), list(3, TRUE)))
#> Piecewise((1, x < 0), (2, (x >= 0) & (x < 1)), (3, x >= 1))
```
