# Fermionic operators and states

`FermionOp(name, annihilation = TRUE)` is the fermionic ladder operator
of `sympy.physics.quantum.fermion` and `FermionFockKet(n)` and
`FermionFockBra(n)` its number states, with `n` 0 or 1.

## Usage

``` r
FermionOp(...)

FermionFockKet(n)

FermionFockBra(n)
```

## Arguments

- ...:

  The name and whether the operator annihilates.

- n:

  The occupation number, 0 or 1.

## Value

A quantum operator or state.

## Examples

``` r
c1 <- FermionOp("c")
doit(AntiCommutator(c1, Dagger(c1)))
#> 1
qapply(Dagger(c1) * FermionFockKet(0))
#> |1>
```
