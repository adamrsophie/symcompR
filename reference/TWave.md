# A transverse wave

`TWave(amplitude, frequency, phase, time_period, n)` is a simple
transverse sinusoidal wave in a medium of refractive index `n`, given by
its frequency or its time period. Its properties are `$amplitude`,
`$frequency`, `$phase`, `$time_period`, `$n`, `$wavelength`, `$speed`,
`$angular_velocity` and `$wavenumber`. Waves of equal frequency add by
interference; a wave times a number scales its amplitude.
`w$rewrite(target)` writes the wave with `"sin"`, `"cos"` or `"exp"`, or
gives the wave equation with `"pde"`. After
`sympy.physics.optics.waves`.

## Usage

``` r
TWave(
  amplitude,
  frequency = NULL,
  phase = S_Zero,
  time_period = NULL,
  n = Symbol("n")
)
```

## Arguments

- amplitude:

  The amplitude.

- frequency:

  The frequency.

- phase:

  The phase angle.

- time_period:

  The time period.

- n:

  The refractive index of the medium.

## Value

A `TWave` object.

## Examples

``` r
f <- Symbol("f")
w1 <- TWave(Symbol("A1"), f, Symbol("phi1"))
w2 <- TWave(Symbol("A2"), f, Symbol("phi2"))
(w1 + w2)$amplitude
#> sqrt(A1^2 + 2*A1*A2*cos(phi1 - phi2) + A2^2)
w1$rewrite("cos")
#> A1*cos(-2*pi*f*t + phi1 + pi*second*f*n*x/(149896229*meter))
```
