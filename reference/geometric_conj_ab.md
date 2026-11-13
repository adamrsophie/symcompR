# Conjugation relations of geometric and Gaussian optics

`geometric_conj_ab(a, b)` is the focal distance that images a point at
distance `a` to distance `b`; `geometric_conj_af(a, f)` and
`geometric_conj_bf(b, f)` give the other distance from one distance and
the focal distance. `gaussian_conj(s_in, z_r_in, f)` gives the list
`(s_out, z_r_out, m)` for a Gaussian beam through a thin lens, and
`conjugate_gauss_beams(wavelen, waist_in, waist_out, f = f)` the list
`(s_in, s_out, f)` of the optical setup that conjugates two beams. After
`sympy.physics.optics.gaussopt`.

## Usage

``` r
geometric_conj_ab(a, b)

geometric_conj_af(a, f)

geometric_conj_bf(a, f)

gaussian_conj(s_in, z_r_in, f)

conjugate_gauss_beams(wavelen, waist_in, waist_out, ...)
```

## Arguments

- a, b:

  The distances of the object and the image.

- f:

  The focal distance.

- s_in:

  The distance of the waist before the lens.

- z_r_in:

  The Rayleigh range of the incoming beam.

- wavelen:

  The wavelength.

- waist_in, waist_out:

  The waists of the two beams.

- ...:

  The focal distance `f`, the only setup supported.

## Value

A `sym` expression or a list.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b"); f <- Symbol("f")
geometric_conj_ab(a, b)
#> a*b/(a + b)
gaussian_conj(Symbol("s_in"), Symbol("z_r_in"), f)[[1]]
#> 1/(-1/(s_in + z_r_in^2/(-f + s_in)) + 1/f)
```
