# Integration and integral transforms

symcompR ports SymPy’s integration machinery: a front end,
[`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md),
that tries several algorithms in turn (table lookup, the rule based
integrator that mimics integration by hand, the Risch algorithm, the
Risch-Norman heuristic, and integration with Meijer G functions), and
the integral transforms built on top of it.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); t <- Symbol("t")
s <- Symbol("s"); k <- Symbol("k")
a <- Symbol("a", positive = TRUE); p <- Symbol("p", positive = TRUE)
f <- Function("f")
```

The function is named
[`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md)
because [`integrate()`](https://rdrr.io/r/stats/integrate.html) already
exists in base R (it does numerical integration of R functions).

## Indefinite integrals

`sym_integrate(e, x)` returns an antiderivative, without the constant of
integration.

``` r

sym_integrate(x^2, x)
#> x^3/3
sym_integrate(1 / x, x)
#> log(x)
sym_integrate(x * cos(x), x)
#> x*sin(x) + cos(x)
sym_integrate(exp(x) * sin(x), x)
#> exp(x)*sin(x)/2 - exp(x)*cos(x)/2
sym_integrate(log(x), x)
#> x*log(x) - x
sym_integrate(tan(x), x)
#> -log(cos(x))
sym_integrate(sin(x)^3, x)
#> cos(x)^3/3 - cos(x)
```

Rational functions are integrated completely, and algebraic and
transcendental integrands are handled by substitution and the Risch
machinery.

``` r

sym_integrate(1 / (x^2 + 1), x)
#> atan(x)
sym_integrate(1 / (x^3 - 1), x)
#> log(x - 1)/3 - log(x^2 + x + 1)/6 - sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
sym_integrate(sqrt(1 - x^2), x)
#> x*sqrt(1 - x^2)/2 + asin(x)/2
sym_integrate(1 / sqrt(x^2 + 1), x)
#> asinh(x)
sym_integrate(x * exp(-x^2), x)
#> -exp(-x^2)/2
```

When no elementary antiderivative exists, the answer may involve special
functions such as the error function:

``` r

sym_integrate(exp(-x^2), x)
#> sqrt(pi)*erf(x)/2
```

Several variables integrate one after the other:

``` r

sym_integrate(x^2 * y, x, y)
#> x^3*y^2/6
```

## Definite and improper integrals

Limits are given as `list(x, a, b)`. Infinite limits use `oo`.

``` r

sym_integrate(x^2, list(x, 0, 1))
#> 1/3
sym_integrate(sin(x), list(x, 0, sym_pi))
#> 2
sym_integrate(abs(x), list(x, -1, 1))
#> 1
sym_integrate(log(x), list(x, 0, 1))
#> -1
sym_integrate(exp(-x), list(x, 0, oo))
#> 1
sym_integrate(1 / x^2, list(x, 1, oo))
#> 1
sym_integrate(1 / (1 + x^2), list(x, -oo, oo))
#> pi
sym_integrate(exp(-x^2), list(x, -oo, oo))
#> sqrt(pi)
sym_integrate(sin(x) / x, list(x, 0, oo))
#> pi/2
```

Assumptions on parameters matter. With p \> 0 the integral defining the
gamma function converges, and symcompR recognizes it:

``` r

sym_integrate(x^(p - 1) * exp(-x), list(x, 0, oo))
#> gamma(p)
sym_integrate(exp(-p * x), list(x, 0, oo))
#> 1/p
```

Multiple integrals take one limit list per variable, innermost first.

``` r

sym_integrate(x * y, list(x, 0, 1), list(y, 0, 2))
#> 1
sym_integrate(exp(-x^2 - y^2), list(x, -oo, oo), list(y, -oo, oo))
#> pi
```

An exact result can be evaluated to any precision with
[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md):

``` r

I1 <- sym_integrate(exp(-x^2), list(x, 0, 1))
I1
#> sqrt(pi)*erf(1)/2
N(I1, 20)
#> 0.74682413281242702540
```

## `Integral` objects

[`Integral()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md)
builds an unevaluated integral, which prints, typesets, and evaluates on
demand with
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md).

``` r

J <- Integral(sqrt(x), x)
J
#> Integral(sqrt(x), x)
doit(J)
#> 2*x^(3/2)/3
doit(Integral(x^2, list(x, 0, 3)))
#> 9
latex(Integral(exp(-x^2), list(x, 0, 1)))
#> [1] "\\int\\limits_{0}^{1} e^{- x^{2}}\\, dx"
```

## Individual algorithms

Each integration algorithm can also be called on its own.

[`manualintegrate()`](https://adamrsophie.github.io/symcompR/reference/integral_steps.md)
integrates the way a person would, with substitution, integration by
parts, and trigonometric rules.
[`integral_steps()`](https://adamrsophie.github.io/symcompR/reference/integral_steps.md)
returns the tree of rules it applied, which can be used to explain a
computation.

``` r

manualintegrate(x * sin(x), x)
#> -x*cos(x) + sin(x)
manualintegrate(1 / (x^2 + 2 * x + 5), x)
#> atan(x/2 + 1/2)/2
manualintegrate(exp(x) / (1 + exp(x)), x)
#> log(exp(x) + 1)
integral_steps(x * sin(x), x)
#> PartsRule(integrand=x*sin(x), variable=x, u=x, dv=sin(x), v_step=SinRule(integrand=sin(x), variable=x), second_step=ConstantTimesRule(integrand=-cos(x), variable=x, constant=-1, other=cos(x), substep=CosRule(integrand=cos(x), variable=x)))
```

[`risch_integrate()`](https://adamrsophie.github.io/symcompR/reference/risch_integrate.md)
implements the Risch algorithm for elementary functions. It is a
decision procedure: when it returns an unevaluated `Integral`, no
elementary antiderivative exists.

``` r

risch_integrate(x * exp(x^2), x)
#> exp(x^2)/2
risch_integrate(x * log(x), x)
#> x^2*log(x)/2 - x^2/4
risch_integrate(exp(-x^2), x)
#> Integral(exp(-x^2), x)
risch_integrate(1 / log(x), x)
#> Integral(1/log(x), x)
```

[`heurisch()`](https://adamrsophie.github.io/symcompR/reference/heurisch.md)
is the Risch-Norman heuristic,
[`ratint()`](https://adamrsophie.github.io/symcompR/reference/ratint.md)
integrates rational functions, and
[`meijerint_indefinite()`](https://adamrsophie.github.io/symcompR/reference/meijerint_indefinite.md)
and
[`meijerint_definite()`](https://adamrsophie.github.io/symcompR/reference/meijerint_indefinite.md)
integrate by converting the integrand to Meijer G functions.

``` r

heurisch(x * exp(x), x)
#> x*exp(x) - exp(x)
ratint((x^2 + 1) / (x^2 - 1), x)
#> x + log(x - 1) - log(x + 1)
ratint(1 / (x^2 + 2), x)
#> sqrt(2)*atan(sqrt(2)*x/2)/2
meijerint_indefinite(x * exp(-x), x)
#> -x*exp(-x) - exp(-x)
meijerint_definite(exp(-x^2), x, 0, oo)
#> [[1]]
#> sqrt(pi)/2
#> 
#> [[2]]
#> True
```

[`meijerint_definite()`](https://adamrsophie.github.io/symcompR/reference/meijerint_indefinite.md)
returns the value together with the condition under which it holds (here
`True`: no condition).

Lower level entry points:
[`heurisch_wrapper()`](https://adamrsophie.github.io/symcompR/reference/heurisch.md)
runs
[`heurisch()`](https://adamrsophie.github.io/symcompR/reference/heurisch.md)
and returns a `Piecewise` that covers the special parameter values,
[`ratint_ratpart()`](https://adamrsophie.github.io/symcompR/reference/ratint.md)
and
[`ratint_logpart()`](https://adamrsophie.github.io/symcompR/reference/ratint.md)
are the rational and logarithmic parts of Horowitz-Ostrogradsky rational
integration,
[`meijerint_inversion()`](https://adamrsophie.github.io/symcompR/reference/meijerint_indefinite.md)
inverts a Laplace transform with Meijer G functions, and
[`singularityintegrate()`](https://adamrsophie.github.io/symcompR/reference/transform-helpers.md)
integrates singularity functions.

``` r

b <- Symbol("b")
heurisch_wrapper(1 / (x^2 + b), x)
#> Piecewise((log(x - sqrt(-b))/(2*sqrt(-b)) - log(x + sqrt(-b))/(2*sqrt(-b)), Ne(b, 0)), (-1/x, True))
show <- function(l) vapply(l, sstr, "")
show(ratint_ratpart(x + 1, x^3 - 2 * x^2 + x, x))
#> [1] "-2/(x - 1)"   "-1/(x^2 - x)"
show(ratint_logpart(S(1), x^2 - 2, x, t)[[1]])
#> [1] "-4*t + x"  "1 - 8*t^2"
meijerint_inversion(1 / (s^2 + 1), s, t)
#> sin(t)*Heaviside(t)
singularityintegrate(SingularityFunction(x, b, 3), x)
#> SingularityFunction(x, b, 4)/4
```

[`log_to_atan()`](https://adamrsophie.github.io/symcompR/reference/ratint.md)
and
[`log_to_real()`](https://adamrsophie.github.io/symcompR/reference/ratint.md)
turn the complex logarithms of the logarithmic part into real
arctangents and logarithms:

``` r

log_to_atan(Poly(x, x, domain = "ZZ"), Poly(1, x, domain = "ZZ"))
#> 2*atan(x)
log_to_atan(Poly(x + S(1) / 2, x, domain = "QQ"), Poly(sqrt(S(3)) / 2, x, domain = "EX"))
#> 2*atan(2*sqrt(3)*x/3 + sqrt(3)/3)
log_to_real(Poly(x + 3 * t / 2 + Rational(1, 2), x, domain = "QQ[t]"),
            Poly(3 * t^2 + 1, t, domain = "ZZ"), x, t)
#> 2*sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
log_to_real(Poly(x^2 - 1, x, domain = "ZZ"), Poly(-2 * t + 1, t, domain = "ZZ"), x, t)
#> log(x^2 - 1)/2
```

## Line integrals

[`line_integrate()`](https://adamrsophie.github.io/symcompR/reference/line_integrate.md)
integrates along a parametric curve, given as
`list(t, a, b, x(t), y(t))`. Integrating x^2 + y^2 around the unit
circle gives its circumference:

``` r

line_integrate(x^2 + y^2, list(t, 0, 2 * sym_pi, cos(t), sin(t)), list(x, y))
#> 2*pi
line_integrate(x + y, list(t, 0, 1, t, t^2), list(x, y))
#> -1/12 - asinh(2)/64 + 67*sqrt(5)/96
```

## Laplace transforms

`laplace_transform(f, t, s)` computes F(s) = \int_0^\infty f(t)
e^{-st}\\dt. By default it returns the transform, the abscissa of
convergence, and auxiliary conditions; with `noconds = TRUE` only the
transform is returned.

``` r

L <- laplace_transform(exp(-a * t), t, s)
L[[1]]
#> 1/(a + s)
L[[2]]
#> -a
laplace_transform(t^2, t, s, noconds = TRUE)
#> 2/s^3
laplace_transform(sin(a * t), t, s, noconds = TRUE)
#> a/(a^2 + s^2)
laplace_transform(t * exp(2 * t), t, s, noconds = TRUE)
#> (s - 2)^(-2)
laplace_transform(Heaviside(t - 3), t, s, noconds = TRUE)
#> exp(-3*s)/s
laplace_transform(DiracDelta(t), t, s, noconds = TRUE)
#> 1
```

The transform of a derivative brings in the initial value, which is how
the Laplace transform turns linear ODEs into algebra:

``` r

laplace_transform(diff(f(t), t), t, s, noconds = TRUE)
#> s*LaplaceTransform(f(t), t, s) - f(0)
```

[`inverse_laplace_transform()`](https://adamrsophie.github.io/symcompR/reference/laplace_transform.md)
goes back. The `Heaviside(t)` factor records that the original function
is taken to vanish for t \< 0.

``` r

inverse_laplace_transform(1 / (s^2 + 1), s, t)
#> sin(t)*Heaviside(t)
inverse_laplace_transform(1 / (s * (s + 1)), s, t)
#> Heaviside(t) - exp(-t)*Heaviside(t)
inverse_laplace_transform(1 / (s - 2)^3, s, t)
#> t^2*exp(2*t)*Heaviside(t)/2
inverse_laplace_transform(exp(-2 * s) / s, s, t)
#> Heaviside(t - 2)
```

[`LaplaceTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md)
is the unevaluated form.

``` r

LaplaceTransform(t^2, t, s)
#> LaplaceTransform(t^2, t, s)
doit(LaplaceTransform(t^2, t, s))
#> 2/s^3
InverseLaplaceTransform(1 / s, s, t)
#> InverseLaplaceTransform(1/s, s, t, _None)
doit(InverseLaplaceTransform(1 / s, s, t))
#> Heaviside(t)
```

[`laplace_correspondence()`](https://adamrsophie.github.io/symcompR/reference/transform-helpers.md)
replaces unevaluated transforms of known functions by names of their
transforms,
[`laplace_initial_conds()`](https://adamrsophie.github.io/symcompR/reference/transform-helpers.md)
inserts initial values, and
[`expand_dirac_delta()`](https://adamrsophie.github.io/symcompR/reference/laplace_transform.md)
writes an expression as a linear combination of Dirac deltas (a constant
term and a coefficient for each delta).

``` r

Fs <- Function("F")
laplace_correspondence(LaplaceTransform(f(t), t, s), list(list(f, Fs)))
#> F(s)
laplace_initial_conds(s^2 * LaplaceTransform(f(t), t, s) - s * f(0) -
                        subs(diff(f(t), t), t, 0), t, list(list(f, c(2, 4))))
#> s^2*LaplaceTransform(f(t), t, s) - 2*s - 4
ed <- expand_dirac_delta((t^2 + 1) * DiracDelta(t - 1))
ed[[1]]
#> 0
show(ed[[2]][[1]])
#> [1] "DiracDelta(t - 1)" "t^2 + 1"
```

## Fourier, Mellin, sine, cosine, and Hankel transforms

The Fourier transform uses SymPy’s convention F(k) =
\int\_{-\infty}^\infty f(x) e^{-2\pi i x k}\\dx.

``` r

fourier_transform(exp(-x^2), x, k)
#> sqrt(pi)*exp(-pi^2*k^2)
inverse_fourier_transform(sqrt(sym_pi) * exp(-(sym_pi * k)^2), k, x)
#> exp(-x^2)
fourier_transform(exp(-abs(x)), x, k)
#> 2/(4*pi^2*k^2 + 1)
```

The Mellin transform \int_0^\infty x^{s-1} f(x)\\dx returns the
transform and its fundamental strip.

``` r

M <- mellin_transform(1 / (1 + x), x, s)
M[[1]]
#> pi/sin(pi*s)
M[[2]]
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1
mellin_transform(exp(-x), x, s, noconds = TRUE)
#> gamma(s)
inverse_mellin_transform(gamma(s), s, x, list(0, oo))
#> exp(-x)
```

Sine and cosine transforms, and the Hankel transform of order \nu:

``` r

sine_transform(exp(-a * x), x, k)
#> sqrt(2)*k/(sqrt(pi)*(a^2 + k^2))
cosine_transform(exp(-a * x), x, k)
#> sqrt(2)*a/(sqrt(pi)*(a^2 + k^2))
inverse_cosine_transform(sqrt(S(2)) * a / (sqrt(sym_pi) * (a^2 + k^2)), k, x)
#> exp(-a*x)
hankel_transform(1 / x, x, k, 0)
#> 1/k
inverse_hankel_transform(1 / k, k, x, 0)
#> 1/x
hankel_transform(exp(-a * x), x, k, 0)
#> a/(k^3*(a^2/k^2 + 1)^(3/2))
```

Each transform has an unevaluated class as well:
[`FourierTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md),
[`MellinTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md),
[`SineTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md),
[`CosineTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md),
and
[`HankelTransform()`](https://adamrsophie.github.io/symcompR/reference/integral-transforms.md),
with inverses. Calling
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md) on
one gives the same result as the matching function above.

``` r

FourierTransform(exp(-x^2), x, k)
#> FourierTransform(exp(-x^2), x, k)
InverseFourierTransform(exp(-k^2), k, x)
#> InverseFourierTransform(exp(-k^2), k, x)
doit(SineTransform(exp(-x), x, k))
#> sqrt(2)*k/(sqrt(pi)*(k^2 + 1))
CosineTransform(exp(-x), x, k)
#> CosineTransform(exp(-x), x, k)
InverseSineTransform(1 / k, k, x)
#> InverseSineTransform(1/k, k, x)
InverseCosineTransform(1 / (1 + k^2), k, x)
#> InverseCosineTransform(1/(k^2 + 1), k, x)
MellinTransform(exp(-x), x, k)
#> MellinTransform(exp(-x), x, k)
InverseMellinTransform(gamma(s), s, x, 0, oo)
#> InverseMellinTransform(gamma(s), s, x, 0, oo)
HankelTransform(exp(-x), x, k, 0)
#> HankelTransform(exp(-x), x, k, 0)
InverseHankelTransform(1 / k, k, x, 0)
#> InverseHankelTransform(1/k, k, x, 0)
inverse_sine_transform(sqrt(S(2)) / (sqrt(sym_pi) * k), k, x)
#> 1
```
