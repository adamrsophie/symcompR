# The isotropic harmonic oscillator

`sho_R_nl()` is the radial wavefunction of the three dimensional
isotropic harmonic oscillator and `sho_E_nl()` its energy, after
`sympy.physics.sho`.

## Usage

``` r
sho_R_nl(n, l, nu, r)

sho_E_nl(n, l, hw)
```

## Arguments

- n:

  The nodal quantum number.

- l:

  The angular momentum quantum number.

- nu:

  The mass times the angular frequency over twice hbar.

- r:

  The radial coordinate.

- hw:

  The energy quantum, hbar times omega.

## Value

A `sym` expression.

## Examples

``` r
r <- Symbol("r"); nu <- Symbol("nu")
sho_R_nl(0, 0, 1, r)
#> 2*2^(3/4)*exp(-r^2)/pi^(1/4)
sho_E_nl(1, 2, Symbol("hw"))
#> 11*hw/2
```
