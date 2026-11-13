# Lens and mirror formulas

`lens_makers_formula(n_lens, n_surr, r1, r2, d)` is the focal length of
a lens with radii of curvature `r1` and `r2` and thickness `d`.
`mirror_formula()` and `lens_formula()` give the third of the focal
length `focal_length`, the object distance `u` and the image distance
`v` from the other two (any may be infinite).
`hyperfocal_distance(f, N, c)` is the hyperfocal distance and
`transverse_magnification(si, so)` the transverse magnification. After
`sympy.physics.optics.utils`.

## Usage

``` r
lens_makers_formula(n_lens, n_surr, r1, r2, d = 0L)

mirror_formula(focal_length = NULL, u = NULL, v = NULL)

lens_formula(focal_length = NULL, u = NULL, v = NULL)

hyperfocal_distance(f, N, c)

transverse_magnification(si, so)
```

## Arguments

- n_lens, n_surr:

  The refractive indices of the lens and its surroundings, or
  [`Medium()`](https://adamrsophie.github.io/symcompR/reference/Medium.md)
  objects.

- r1, r2:

  The radii of curvature of the two surfaces.

- d:

  The thickness of the lens.

- focal_length, u, v:

  The focal length and the object and image distances.

- f:

  The focal length.

- N:

  The f-number.

- c:

  The circle of confusion.

- si, so:

  The image and object distances.

## Value

A `sym` expression.

## Examples

``` r
f <- Symbol("f"); u <- Symbol("u")
lens_formula(focal_length = f, u = u)
#> f*u/(f + u)
mirror_formula(u = oo, v = Symbol("v"))
#> v
transverse_magnification(30, 15)
#> -2
```
