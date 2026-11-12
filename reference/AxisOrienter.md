# Orienters

`AxisOrienter(angle, axis)`, `BodyOrienter()`, `SpaceOrienter()` and
`QuaternionOrienter()` describe rotations for
[`orient_new()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md).
`orienter_rotation_matrix()` gives the rotation matrix (an axis orienter
needs the system the axis is expressed in).

## Usage

``` r
AxisOrienter(angle, axis)

BodyOrienter(angle1, angle2, angle3, rot_order)

SpaceOrienter(angle1, angle2, angle3, rot_order)

QuaternionOrienter(q0, q1, q2, q3)

orienter_rotation_matrix(orienter, system = NULL)
```

## Arguments

- angle, angle1, angle2, angle3:

  Angles.

- axis:

  A vector.

- rot_order:

  A rotation order such as `"123"` or `"XYZ"`.

- q0, q1, q2, q3:

  Quaternion parameters.

- orienter:

  An orienter.

- system:

  A coordinate system (axis orienters only).

## Value

An orienter, or a matrix.

## Examples

``` r
A <- CoordSys3D("A")
a <- Symbol("a")
orienter_rotation_matrix(AxisOrienter(a, A$k), A)
#> Matrix(3, 3)
#> [cos(a) , sin(a), 0]
#> [-sin(a), cos(a), 0]
#> [0      , 0     , 1]
orienter_rotation_matrix(BodyOrienter(a, 0, 0, "123"))
#> Matrix(3, 3)
#> [1, 0      , 0     ]
#> [0, cos(a) , sin(a)]
#> [0, -sin(a), cos(a)]
```
