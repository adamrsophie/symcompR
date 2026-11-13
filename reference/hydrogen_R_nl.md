# The hydrogen atom

Wavefunctions and energies of the nonrelativistic and relativistic
hydrogen atom in atomic units, after `sympy.physics.hydrogen`.
`hydrogen_R_nl()` is the radial wavefunction, `hydrogen_Psi_nlm()` the
full wavefunction, `hydrogen_E_nl()` the nonrelativistic energy and
`hydrogen_E_nl_dirac()` the relativistic energy from the Dirac equation.
`hydrogen_Y_lm()` and `hydrogen_Z_lm()` are the expanded complex and
real spherical harmonics.

## Usage

``` r
hydrogen_R_nl(n, l, r, Z = 1)

hydrogen_Y_lm(l, m, phi, theta)

hydrogen_Z_lm(l, m, phi, theta)

hydrogen_Psi_nlm(n, l, m, r, phi, theta, Z = 1)

hydrogen_E_nl(n, Z = 1)

hydrogen_E_nl_dirac(n, l, spin_up = TRUE, Z = 1, c = Float("137.035999037"))
```

## Arguments

- n:

  The principal quantum number.

- l:

  The angular momentum quantum number.

- r:

  The radial coordinate.

- Z:

  The atomic number (1 for hydrogen).

- m:

  The magnetic quantum number.

- phi, theta:

  The azimuthal and polar angles.

- spin_up:

  `TRUE` for the electron spin up state.

- c:

  The speed of light in atomic units.

## Value

A `sym` expression.

## Examples

``` r
r <- Symbol("r")
hydrogen_R_nl(1, 0, r)
#> 2*exp(-r)
hydrogen_E_nl(2)
#> -1/8
hydrogen_E_nl_dirac(2, 0, Z = 1, c = 1)
#> -1 + sqrt(2)/2
```
