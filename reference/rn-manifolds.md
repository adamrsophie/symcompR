# Predefined manifolds of sympy.diffgeom.rn

`R2`, `R3`, their patches `R2_origin` and `R3_origin`, the coordinate
systems `R2_r` (rectangular), `R2_p` (polar), `R3_r`, `R3_c`
(cylindrical) and `R3_s` (spherical), and the transformation laws
`relations_2d` and `relations_3d`, as in SymPy;
[`diffgeom_rn()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
returns them all in a list. `BaseScalarField(cs, i)`,
`BaseVectorField(cs, i)` and `CoordinateSymbol(cs, i)` are SymPy's
classes for the coordinate function, the base vector and the coordinate
symbol of index `i`.

## Usage

``` r
BaseScalarField(coord_sys, index)

BaseVectorField(coord_sys, index)

CoordinateSymbol(coord_sys, index)

R2

R2_origin

relations_2d

R2_r

R2_p

R3

R3_origin

relations_3d

R3_r

R3_c

R3_s
```

## Format

Manifolds, patches, coordinate systems and lists.

## Arguments

- coord_sys:

  A coordinate system.

- index:

  A coordinate index, from 0.

## Value

`BaseScalarField()`, `BaseVectorField()` and `CoordinateSymbol()` return
`sym` expressions.

## Examples

``` r
R2_r
#> rectangular
BaseScalarField(R2_r, 0)
#> x
BaseVectorField(R2_p, 1)
#> e_theta
```
