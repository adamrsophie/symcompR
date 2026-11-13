# The one dimensional quantum harmonic oscillator

`qho_psi_n()` is the wavefunction, `qho_E_n()` the energy and
`qho_coherent_state()` the coefficient of the coherent state, after
`sympy.physics.qho_1d`. The reduced Planck constant appears as `hbar`.

## Usage

``` r
qho_psi_n(n, x, m, omega)

qho_E_n(n, omega)

qho_coherent_state(n, alpha)
```

## Arguments

- n:

  The quantum number.

- x:

  The coordinate.

- m:

  The mass.

- omega:

  The angular frequency.

- alpha:

  The eigenvalue of the annihilation operator.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
qho_psi_n(0, x, 1, 1)
#> exp(-x^2/(2*hbar))/(hbar^(1/4)*pi^(1/4))
qho_E_n(1, Symbol("omega"))
#> 3*hbar*omega/2
```
