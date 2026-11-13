# Density operators, entropy and fidelity

`Density(list(state, prob), ...)` is the density operator of
`sympy.physics.quantum.density`:
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
expands it into outer products,
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
gives its matrix, and
[`Tr()`](https://adamrsophie.github.io/symcompR/reference/Tr.md) traces
it. [`Density()`](https://adamrsophie.github.io/symcompR/reference/E.md)
and
[`entropy()`](https://adamrsophie.github.io/symcompR/reference/rv_moments.md)
are shared with sympy.stats (see
[`Density()`](https://adamrsophie.github.io/symcompR/reference/E.md) and
[`entropy()`](https://adamrsophie.github.io/symcompR/reference/rv_moments.md));
they act as here when the arguments are `list(state, probability)`
pairs, a quantum density operator or a matrix.
[`entropy()`](https://adamrsophie.github.io/symcompR/reference/rv_moments.md)
is then the von Neumann entropy. `fidelity()` is the fidelity of two
density operators or matrices.

## Usage

``` r
fidelity(state1, state2)
```

## Arguments

- state1, state2:

  Density operators or matrices.

## Value

A density operator or an expression.

## Examples

``` r
d <- Density(list(JzKet(1, 0), 0.5), list(JzKet(1, 1), 0.5))
doit(d)
#> 0.5*|1,0><1,0| + 0.5*|1,1><1,1|
entropy(Density(list(JzKet(S(1)/2, S(1)/2), S(1)/2), list(JzKet(S(1)/2, -S(1)/2), S(1)/2)))
#> log(2)/2
```
