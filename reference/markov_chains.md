# Markov chains

`DiscreteMarkovChain()` and `ContinuousMarkovChain()` of `sympy.stats`:
finite, time homogeneous Markov chains given by a state space and a
transition matrix (or a generator matrix for continuous time). A
discrete chain is indexed with
[`rv_getitem()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md)
(SymPy's `Y[t]`), a continuous chain is called, `C(t)`. Queries go
through [`P()`](https://adamrsophie.github.io/symcompR/reference/E.md)
and [`E()`](https://adamrsophie.github.io/symcompR/reference/E.md); pass
a transition matrix, a generator matrix or a state space in the given
condition with `TransitionMatrixOf()`, `GeneratorMatrixOf()` and
`StochasticStateSpaceOf()`, combined with `&`. Methods such as
`"communication_classes"`, `"decompose"`, `"canonical_form"`,
`"fundamental_matrix"`, `"absorbing_probabilities"`, `"is_regular"`,
`"is_ergodic"`, `"is_absorbing_chain"`, `"stationary_distribution"`,
`"fixed_row_vector"`, `"limiting_distribution"` and
`"transition_probabilities"` are called with
[`rv_attr()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md).
Python tuples come back as lists.

## Usage

``` r
DiscreteMarkovChain(sym, state_space = NULL, trans_probs = NULL)

ContinuousMarkovChain(sym, state_space = NULL, gen_mat = NULL)

TransitionMatrixOf(process, matrix)

GeneratorMatrixOf(process, matrix)

StochasticStateSpaceOf(process, state_space)
```

## Arguments

- sym:

  The name of the chain.

- state_space:

  Optional state space: a list of states, a `Tuple`, a `FiniteSet` or a
  `Range`; by default `Range(n)`.

- trans_probs:

  Optional transition matrix (a `Matrix` or a `MatrixSymbol`); by
  default `MatrixSymbol("_T", n, n)`.

- gen_mat:

  Optional generator matrix of a continuous chain.

- process:

  A Markov chain.

- matrix:

  A transition or generator matrix.

## Value

A process (an R function), or an assumption to use in a given condition.

## Examples

``` r
T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
Y <- DiscreteMarkovChain("Y", list(0, 1), T)
P(Eq(rv_getitem(Y, 3), 1), Eq(rv_getitem(Y, 1), 0))
#> 5/8
rv_attr(Y, "stationary_distribution")
#> Matrix(1, 2)
#> [1/3, 2/3]
rv_attr(Y, "communication_classes")
#> [[1]]
#> ([0, 1], True, 1) 
#> 
G <- Matrix(list(list(-1, 1), list(1, -1)))
C <- ContinuousMarkovChain("C", list(0, 1), G)
rv_attr(C, "limiting_distribution")
#> Matrix(1, 2)
#> [1/2, 1/2]
```
