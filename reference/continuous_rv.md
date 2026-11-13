# Continuous random variables

Constructors for the continuous distributions of `sympy.stats`. Each
returns a random symbol named `name` on a `SingleContinuousPSpace`.
`ContinuousRV(symbol, density, set)` builds one from a density
expression in `symbol`. `BetaRV()` and `GammaRV()` are SymPy's `Beta()`
and [`Gamma()`](https://rdrr.io/r/stats/family.html), renamed so they do
not mask base R; `Wald()` is `GaussianInverse()` and `Erlang(k, l)` is
`GammaRV(k, 1/l)`.

## Usage

``` r
ContinuousRV(symbol, density, set = NULL, check = FALSE)

Arcsin(name, a = 0, b = 1)

Benini(name, alpha, beta, sigma)

BetaRV(name, alpha, beta)

BetaNoncentral(name, alpha, beta, lamda)

BetaPrime(name, alpha, beta)

BoundedPareto(name, alpha, left, right)

Cauchy(name, x0, gamma)

ChiRV(name, k)

ChiNoncentral(name, k, l)

ChiSquared(name, k)

Dagum(name, p, a, b)

Davis(name, b, n, mu)

Erlang(name, k, l)

ExGaussian(name, mean, std, rate)

Exponential(name, rate)

ExponentialPower(name, mu, alpha, beta)

FDistribution(name, d1, d2)

FisherZ(name, d1, d2)

Frechet(name, a, s = 1, m = 0)

GammaRV(name, k, theta)

GammaInverse(name, a, b)

Gumbel(name, beta, mu, minimum = FALSE)

Gompertz(name, b, eta)

Kumaraswamy(name, a, b)

Laplace(name, mu, b)

Levy(name, mu, c)

LogCauchy(name, mu, sigma)

Logistic(name, mu, s)

LogLogistic(name, alpha, beta)

LogitNormal(name, mu, s)

LogNormal(name, mean, std)

Lomax(name, alpha, lamda)

Maxwell(name, a)

Moyal(name, mu, sigma)

Nakagami(name, mu, omega)

Normal(name, mean, std)

GaussianInverse(name, mean, shape)

Wald(name, mean, shape)

Pareto(name, xm, alpha)

PowerFunction(name, alpha, a, b)

QuadraticU(name, a, b)

RaisedCosine(name, mu, s)

Rayleigh(name, sigma)

Reciprocal(name, a, b)

ShiftedGompertz(name, b, eta)

SkewNormal(name, xi, omega, alpha)

StudentT(name, nu)

Trapezoidal(name, a, b, c, d)

Triangular(name, a, b, c)

Uniform(name, left, right)

UniformSum(name, n)

VonMises(name, mu, k)

Weibull(name, alpha, beta)

WignerSemicircle(name, R)
```

## Arguments

- symbol:

  A symbol for `ContinuousRV()`.

- density:

  A density expression in `symbol`.

- set:

  The support, an interval (the real line by default).

- check:

  Whether to check the arguments.

- name:

  The name of the random variable (a string or symbol).

- a, b, c, d, alpha, beta, sigma, lamda, left, right, x0, gamma, k, l,
  p, n, mu, mean, std, rate, d1, d2, s, m, theta, minimum, eta, omega,
  shape, xm, xi, nu, R:

  Distribution parameters, named as in SymPy.

## Value

A random symbol.

## Examples

``` r
X <- Normal("X", 0, 1)
rv_density(X)
#> NormalDistribution(0, 1)
E(Exponential("T", 2))
#> 1/2
variance(Uniform("U", 0, 1))
#> 1/12
```
