# Spin states

The eigenstates of `sympy.physics.quantum.spin`: `JzKet(j, m)`,
`JxKet()`, `JyKet()` and their bras, and the coupled states
`JzKetCoupled(j, m, jn, jcoupling)` (and the `Jx`, `Jy` and bra
versions) of several spin spaces with momenta `jn`. `couple()` writes
tensor products of uncoupled states in the coupled basis and
`uncouple()` the reverse.
[`rewrite()`](https://adamrsophie.github.io/symcompR/reference/rewrite.md)
with `"Jx"`, `"Jy"` or `"Jz"` changes the basis.

## Usage

``` r
SpinState(j, m)

JxKet(j, m)

JxBra(j, m)

JyKet(j, m)

JyBra(j, m)

JzKet(j, m)

JzBra(j, m)

CoupledSpinState(j, m, ...)

JxKetCoupled(j, m, ...)

JxBraCoupled(j, m, ...)

JyKetCoupled(j, m, ...)

JyBraCoupled(j, m, ...)

JzKetCoupled(j, m, ...)

JzBraCoupled(j, m, ...)

couple(expr, jcoupling_list = NULL)

uncouple(expr, jn = NULL, jcoupling_list = NULL)
```

## Arguments

- j, m:

  The angular momentum and its projection.

- ...:

  `jn` and the optional `jcoupling` of coupled states.

- expr:

  An expression of spin states.

- jcoupling_list:

  The order in which spaces are coupled.

- jn:

  The momenta of the spaces, for uncoupling an uncoupled state.

## Value

A spin state or an expression.

## Examples

``` r
JzKet(1, 0)
#> |1,0>
couple(TensorProduct(JzKet(S(1)/2, S(1)/2), JzKet(S(1)/2, S(1)/2)))
#> |1,1,j1=1/2,j2=1/2>
uncouple(JzKetCoupled(1, 0, list(S(1)/2, S(1)/2)))
#> sqrt(2)*|1/2,-1/2>|1/2,1/2>/2 + sqrt(2)*|1/2,1/2>|1/2,-1/2>/2
```
