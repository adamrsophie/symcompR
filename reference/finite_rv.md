# Finite random variables

Constructors for the finite distributions of `sympy.stats`: `FiniteRV()`
from a mapping of values to probabilities, `DiscreteUniform()`, `Die()`,
`Bernoulli()`, `Coin()`, `Binomial()`, `BetaBinomial()`,
`Hypergeometric()`, `Rademacher()`, `IdealSoliton()` and
`RobustSoliton()`.

## Usage

``` r
FiniteRV(name, density, check = FALSE)

DiscreteUniform(name, items)

Die(name, sides = 6)

Bernoulli(name, p, succ = 1, fail = 0)

Coin(name, p = S_Half)

Binomial(name, n, p, succ = 1, fail = 0)

BetaBinomial(name, n, alpha, beta)

Hypergeometric(name, N, m, n)

Rademacher(name)

IdealSoliton(name, k)

RobustSoliton(name, k, delta, c)
```

## Arguments

- name:

  The name of the random variable.

- density:

  A [`Dict()`](https://adamrsophie.github.io/symcompR/reference/Dict.md)
  (or a list of values with names) mapping values to probabilities.

- check:

  Whether to check the arguments.

- items:

  A list of equally likely values.

- sides:

  Number of sides of the die.

- p, succ, fail:

  Success probability and the success and failure values.

- n, alpha, beta, N, m, k, delta, c:

  Distribution parameters, as in SymPy.

## Value

A random symbol.

## Examples

``` r
D <- Die("D", 6)
E(D)
#> 7/2
P(D > 3)
#> 1/2
rv_density(Coin("C"))
#> BernoulliDistribution(1/2, H, T)
```
