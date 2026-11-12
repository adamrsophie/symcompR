# Random variables

Constructors for random symbols with exact symbolic distributions, after
`sympy.stats`: continuous `Normal`, `Uniform`, `Exponential`, `GammaRV`,
`BetaRV`, `ChiSquared`, `StudentT`, `Cauchy`, `Laplace`, `LogNormal`,
`Weibull`, `Pareto`; discrete `Bernoulli`, `Binomial`, `Poisson`,
`Geometric`, `Die` and `Coin`. Parameters may be symbolic; give them
assumptions such as `positive = TRUE` so that integrals converge.

## Usage

``` r
Normal(name, mu, sigma)

Uniform(name, a, b)

Exponential(name, rate)

GammaRV(name, k, theta)

BetaRV(name, alpha, beta)

ChiSquared(name, k)

StudentT(name, nu)

Cauchy(name, x0, gamma)

Laplace(name, mu, b)

LogNormal(name, mu, sigma)

Weibull(name, alpha, beta)

Pareto(name, x0, alpha)

Bernoulli(name, p)

Coin(name, p = Rational(1, 2))

Die(name, sides = 6)

Binomial(name, n, p)

Poisson(name, lamda)

Geometric(name, p)

E(expr)

expectation(expr)

moment(X, k)

cmoment(X, k)

variance(X)

std(X)

covariance(X, Y)

correlation(X, Y)

skewness(X)

kurtosis(X)

rv_density(X, z = Symbol("z"))

rv_cdf(X, z = Symbol("z"))

P(condition)

rv_sample(X, size = 1)
```

## Arguments

- name:

  Display name.

- mu, sigma, a, b, rate, k, theta, alpha, beta, nu, x0, gamma, p, n,
  sides, lamda:

  Distribution parameters.

- expr:

  An expression in random variables.

- X, Y:

  Random variables or expressions in them.

- z:

  Symbol for density and distribution functions.

- condition:

  A relational in one random variable.

- size:

  Number of draws.

## Value

A random symbol, expression or numeric vector.

## Details

`E(expr)` (alias `expectation()`) is the expected value, using closed
form moments when known and integration or summation otherwise.
`variance()`, `std()`, `moment()`, `cmoment()`, `skewness()`,
`kurtosis()`, `covariance()` and `correlation()` build on it.
`P(condition)` is the probability of an event in one random variable.
`rv_density(X)` gives the density or mass function and `rv_cdf(X)` the
distribution function, both as expressions in the variable `z`. Random
variables in a product or sum are treated as independent.
`rv_sample(X, n)` draws numeric samples.

## Examples

``` r
mu <- Symbol("mu", real = TRUE); s <- Symbol("s", positive = TRUE)
X <- Normal("X", mu, s)
E(X^2)
#> mu^2 + s^2
variance(2 * X + 1)
#> 4*s^2
Z <- Normal("Z", 0, 1)
P(Z > 0)
#> 1/2
D <- Die("D", 6)
E(D)
#> 7/2
```
