# Bosonic operators and states

`BosonOp(name, annihilation = TRUE)` is the bosonic ladder operator of
`sympy.physics.quantum.boson` (its adjoint is the creation operator),
`BosonFockKet(n)` and `BosonFockBra(n)` are number states and
`BosonCoherentKet(alpha)` and `BosonCoherentBra(alpha)` coherent states.

## Usage

``` r
BosonOp(...)

BosonFockKet(n)

BosonFockBra(n)

BosonCoherentKet(alpha)

BosonCoherentBra(alpha)
```

## Arguments

- ...:

  The name and whether the operator annihilates.

- n:

  The occupation number.

- alpha:

  The eigenvalue of the annihilation operator.

## Value

A quantum operator or state.

## Examples

``` r
a <- BosonOp("a")
Dagger(a)
#> Dagger(a)
qapply(a * BosonFockKet(2))
#> sqrt(2)*|1>
doit(Commutator(a, Dagger(a)))
#> 1
```
