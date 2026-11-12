# Properties of coordinate systems

[`base_vectors()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
and
[`base_scalars()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
give the basis of a system (they also accept a diffgeom coordinate
system), `lame_coefficients()` its Lame coefficients,
`transformation_to_parent()`, `transformation_from_parent()` and
`transformation_from_parent_function()` its transformation equations,
`rotation_matrix(sys, other)` the direction cosine matrix,
`change_of_basis_matrix_from(sys, other)` the change of basis matrix,
`scalar_map(sys, other)` the base scalars of `sys` in terms of those of
`other`, `position_wrt()` the position of a system or point relative to
another, and `express_coordinates()` the coordinates of a point.

## Usage

``` r
lame_coefficients(sys)

transformation_to_parent(sys)

transformation_from_parent(sys)

transformation_from_parent_function(sys)

rotation_matrix(sys, other)

change_of_basis_matrix_from(sys, other)

position_wrt(sys, other)

express_coordinates(sys, system)

scalar_map(sys, other)
```

## Arguments

- sys:

  A coordinate system (or a point for `position_wrt()` and
  `express_coordinates()`).

- other:

  Another coordinate system or point.

- system:

  A coordinate system.

## Value

A list, matrix, function or vector.

## Examples

``` r
A <- CoordSys3D("A")
B <- orient_new_axis(A, "B", Symbol("q"), A$k)
rotation_matrix(A, B)
#> Matrix(3, 3)
#> [cos(q), -sin(q), 0]
#> [sin(q), cos(q) , 0]
#> [0     , 0      , 1]
scalar_map(B, A)
#> {B.x: A.x*cos(q) + A.y*sin(q), B.y: -A.x*sin(q) + A.y*cos(q), B.z: A.z} 
```
