# Refraction and reflection at an interface

`refraction_angle(incident, medium1, medium2, normal, plane)` is the
angle of refraction for an angle of incidence, or the direction (a
matrix or a `Ray3D`) of the refracted ray for an incident direction
given as a matrix, a sequence or a `Ray3D` with the normal or the plane
of the interface. `deviation()` is the angle of deviation of the
refracted ray. `fresnel_coefficients(angle, medium1, medium2)` is the
list of reflection and transmission coefficients `(R_p, R_s, T_p, T_s)`,
or `(R_p, R_s)` past the critical angle. `brewster_angle()` and
`critical_angle()` are the Brewster and critical angles. The media are
[`Medium()`](https://adamrsophie.github.io/symcompR/reference/Medium.md)
objects or refractive indices. After `sympy.physics.optics.utils`.

## Usage

``` r
refraction_angle(incident, medium1, medium2, normal = NULL, plane = NULL)

fresnel_coefficients(angle_of_incidence, medium1, medium2)

deviation(incident, medium1, medium2, normal = NULL, plane = NULL)

brewster_angle(medium1, medium2)

critical_angle(medium1, medium2)
```

## Arguments

- incident:

  The angle of incidence or the incident direction.

- medium1, medium2:

  The media on the two sides.

- normal:

  The normal of the interface.

- plane:

  The `Plane` of the interface.

- angle_of_incidence:

  The angle of incidence.

## Value

A `sym` expression, a `symmatrix`, a `Ray3D` or a list.

## Examples

``` r
refraction_angle(0.5, 1, 2)
#> 0.242069983487306
refraction_angle(list(1, 1, 1), 1, 1, list(0, 0, 1))
#> Matrix(3, 1)
#> [1 ]
#> [1 ]
#> [-1]
critical_angle(1.33, 1)
#> 0.850908514477849
```
