# The particle on a ring

`pring_wavefunction()` is the normalized wavefunction and
`pring_energy()` the energy of a particle on a ring, after
`sympy.physics.pring`.

## Usage

``` r
pring_wavefunction(n, x)

pring_energy(n, m, r)
```

## Arguments

- n:

  The quantum number.

- x:

  The angle.

- m:

  The mass.

- r:

  The radius of the ring.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
pring_wavefunction(1, x)
#> sqrt(2)*exp(I*x)/(2*sqrt(pi))
pring_energy(1, Symbol("m"), Symbol("r"))
#> hbar^2/(2*m*r^2)
```
