# Jones and Stokes vectors and polarization optics

`jones_vector(psi, chi)` and `stokes_vector(psi, chi, p, I)` are the
Jones and Stokes vectors of light with polarization ellipse angles `psi`
and `chi`, degree of polarization `p` and intensity `I`;
`jones_2_stokes(e)` converts a Jones vector to a Stokes vector.
`linear_polarizer()`, `phase_retarder()`, `half_wave_retarder()`,
`quarter_wave_retarder()`, `transmissive_filter()` and
`reflective_filter()` are Jones matrices of optical elements;
`mueller_matrix(J)` is the Mueller matrix of a Jones matrix and
`polarizing_beam_splitter()` the Jones matrix of a polarizing beam
splitter. After `sympy.physics.optics.polarization`.

## Usage

``` r
jones_vector(psi, chi)

stokes_vector(psi, chi, p = 1L, I = 1L)

jones_2_stokes(e)

linear_polarizer(theta = 0L)

phase_retarder(theta = 0L, delta = 0L)

half_wave_retarder(theta)

quarter_wave_retarder(theta)

transmissive_filter(T)

reflective_filter(R)

mueller_matrix(J)

polarizing_beam_splitter(
  Tp = 1L,
  Rs = 1L,
  Ts = 0L,
  Rp = 0L,
  phia = 0L,
  phib = 0L
)
```

## Arguments

- psi, chi:

  The angles of the polarization ellipse.

- p:

  The degree of polarization.

- I:

  The intensity.

- e:

  A Jones vector.

- theta:

  The angle of the fast axis.

- delta:

  The phase difference.

- T, R:

  The transmittance and reflectance.

- J:

  A Jones matrix.

- Tp, Rs, Ts, Rp:

  The transmissivities and reflectivities of the two polarizations.

- phia, phib:

  The phase differences of the two reflections.

## Value

A `symmatrix`.

## Examples

``` r
jones_vector(0, 0)
#> Matrix(2, 1)
#> [1]
#> [0]
linear_polarizer(Symbol("theta"))
#> Matrix(2, 2)
#> [cos(theta)^2         , sin(theta)*cos(theta)]
#> [sin(theta)*cos(theta), sin(theta)^2         ]
mueller_matrix(linear_polarizer(0))
#> Matrix(4, 4)
#> [1/2, 1/2, 0, 0]
#> [1/2, 1/2, 0, 0]
#> [0  , 0  , 0, 0]
#> [0  , 0  , 0, 0]
```
