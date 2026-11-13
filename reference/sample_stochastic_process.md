# Sample a stochastic process

`sample_stochastic_process(process)` of `sympy.stats`: a generator of
the states of a discrete Markov chain. It is an R function; each call
returns the next state. The first state is drawn uniformly from the
state space with R's random number generator (use
[`set.seed()`](https://rdrr.io/r/base/Random.html)). SymPy draws the
later states with scipy, so a second call raises SymPy's error when
scipy is missing.

## Usage

``` r
sample_stochastic_process(process)
```

## Arguments

- process:

  A stochastic process.

## Value

A function returning the next sampled state at each call.

## Examples

``` r
set.seed(1)
T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
g <- sample_stochastic_process(DiscreteMarkovChain("Y", list(0, 1), T))
g()
#> 0
```
