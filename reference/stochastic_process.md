# Stochastic processes

`BernoulliProcess()`, `PoissonProcess()`, `WienerProcess()` and
`GammaProcess()` of `sympy.stats`. A discrete time process is indexed
with
[`rv_getitem()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md)
(SymPy's `X[t]`); a continuous time process is called, `X(t)`. Process
attributes such as `"state_space"`, `"index_set"` or
`"joint_distribution"` are read with
[`rv_attr()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md).

## Usage

``` r
RandomIndexedSymbol(idx_obj, pspace = NULL)

BernoulliProcess(sym, p, success = 1, failure = 0)

PoissonProcess(sym, lamda)

WienerProcess(sym)

GammaProcess(sym, lamda, gamma)
```

## Arguments

- idx_obj:

  An indexed object or an applied function.

- pspace:

  A probability space.

- sym:

  The name of the process.

- p, success, failure:

  Probability and outcomes of a Bernoulli trial.

- lamda, gamma:

  Rate and shape parameters.

## Value

A process (an R function), or a random indexed symbol.

## Examples

``` r
B <- BernoulliProcess("B", p = Rational(1, 2))
E(rv_getitem(B, 1) + rv_getitem(B, 2))
#> 1
X <- PoissonProcess("X", 3)
E(X(2))
#> 6
```
