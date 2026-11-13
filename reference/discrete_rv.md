# Discrete random variables

Constructors for the countably infinite discrete distributions of
`sympy.stats`: `DiscreteRV()` from a probability mass function,
`FlorySchulz()`, `Geometric()`, `Hermite()`, `Logarithmic()`,
`NegativeBinomial()`, `Poisson()`, `Skellam()`, `YuleSimon()` and
`Zeta()`.

## Usage

``` r
DiscreteRV(symbol, density, set = Integers, check = FALSE)

FlorySchulz(name, a)

Geometric(name, p)

Hermite(name, a1, a2)

Logarithmic(name, p)

NegativeBinomial(name, r, p)

Poisson(name, lamda)

Skellam(name, mu1, mu2)

YuleSimon(name, rho)

Zeta(name, s)
```

## Arguments

- symbol, density, set:

  For `DiscreteRV()`: the symbol, the mass function in it, and its
  support (the integers by default).

- check:

  Whether to check the arguments.

- name:

  The name of the random variable.

- a, p, a1, a2, r, lamda, mu1, mu2, rho, s:

  Distribution parameters, named as in SymPy.

## Value

A random symbol.

## Examples

``` r
X <- Poisson("X", 3)
E(X)
#> 3
variance(Geometric("G", Rational(1, 4)))
#> 12
```
