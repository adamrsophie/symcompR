# Discretization and stability margins

`gbt(tf, sample_per, alpha)` gives the numerator and denominator
coefficients of the discretization of a continuous transfer function by
the generalized bilinear transformation; `bilinear()`, `forward_diff()`
and `backward_diff()` are the cases `alpha = 1/2`, 0 and 1.
`phase_margin()` and `gain_margin()` are the phase margin (in degrees)
and the gain margin (in decibels) of a SISO system. After
`sympy.physics.control.lti`.

## Usage

``` r
gbt(tf, sample_per, alpha)

bilinear(tf, sample_per)

forward_diff(tf, sample_per)

backward_diff(tf, sample_per)

phase_margin(system)

gain_margin(system)
```

## Arguments

- tf:

  A continuous transfer function.

- sample_per:

  The sampling period.

- alpha:

  The parameter of the transformation.

- system:

  A SISO system.

## Value

A list of two coefficient lists, or a `sym` expression.

## Examples

``` r
s <- Symbol("s"); T <- Symbol("T")
bilinear(TransferFunction(1, s + 1, s), T)
#> ([T/(2*(T/2 + 1)), T/(2*(T/2 + 1))], [1, (T/2 - 1)/(T/2 + 1)]) 
gain_margin(TransferFunction(1, s^3 + 2*s^2 + 2*s + 1, s))
#> 20*log(3)/log(10)
```
