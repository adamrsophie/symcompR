# The parameter of a Gaussian beam

`BeamParameter(wavelen, z, z_r = NULL, w = NULL, n = 1)` is the complex
beam parameter `q = z + I*z_r` of a Gaussian beam of wavelength
`wavelen` at distance `z` from the waist, given by the Rayleigh range
`z_r` or the waist `w`, in a medium of refractive index `n`. Its
properties are `$wavelen`, `$z`, `$z_r`, `$n`, `$q`, `$radius`, `$w`,
`$w_0`, `$divergence`, `$gouy` and `$waist_approximation_limit`.
`waist2rayleigh()` and `rayleigh2waist()` convert between the waist and
the Rayleigh range. After `sympy.physics.optics.gaussopt`.

## Usage

``` r
BeamParameter(wavelen, z, z_r = NULL, w = NULL, n = 1L)

waist2rayleigh(w, wavelen, n = 1L)

rayleigh2waist(z_r, wavelen)
```

## Arguments

- wavelen:

  The wavelength.

- z:

  The distance to the waist.

- z_r:

  The Rayleigh range.

- w:

  The waist.

- n:

  The refractive index of the medium.

## Value

A `BeamParameter` object, or a `sym` expression.

## Examples

``` r
p <- BeamParameter(530e-9, 1, w = 1e-3)
p$q
#> 1 + 1.88679245283019*I*pi
N(p$w_0)
#> 0.00100000000000000
waist2rayleigh(Symbol("w"), Symbol("wavelen"))
#> pi*w^2/wavelen
```
