# Particle in a box

The Hamiltonian `PIABHamiltonian()` and the eigenstates `PIABKet(n)` and
`PIABBra(n)` of a particle of mass `m` in a box of length `L`, as in
`sympy.physics.quantum.piab`.

## Usage

``` r
PIABHamiltonian(...)

PIABKet(...)

PIABBra(...)
```

## Arguments

- ...:

  The label.

## Value

A quantum operator or state.

## Examples

``` r
n <- Symbol("n")
qapply(PIABHamiltonian("H") * PIABKet(n))
#> hbar^2*pi^2*n^2*|n>/(2*L^2*m)
```
