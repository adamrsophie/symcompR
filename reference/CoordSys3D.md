# Three dimensional coordinate systems

`CoordSys3D(name)` creates a Cartesian coordinate system. Its base
scalars and base vectors are reached with `$`: `C$x`, `C$y`, `C$z`,
`C$i`, `C$j`, `C$k`, and its origin with `C$origin`. A `transformation`
of `"spherical"` or `"cylindrical"` (or `"cartesian"`) gives curvilinear
systems, whose scalars default to `r, theta, phi` and `r, theta, z`; a
list of two lists `list(list(x, y, z), list(f1, f2, f3))` or an R
function of three arguments gives a user defined orthogonal system, and
`list(rotation_matrix, location)` a rotated and translated child of
`parent`.

## Usage

``` r
CoordSys3D(
  name,
  transformation = NULL,
  parent = NULL,
  location = NULL,
  rotation_matrix = NULL,
  vector_names = NULL,
  variable_names = NULL
)

BaseScalar(index, system)

locate_new(sys, name, position, vector_names = NULL, variable_names = NULL)

orient_new(
  sys,
  name,
  orienters,
  location = NULL,
  vector_names = NULL,
  variable_names = NULL
)

orient_new_axis(
  sys,
  name,
  angle,
  axis,
  location = NULL,
  vector_names = NULL,
  variable_names = NULL
)

orient_new_body(
  sys,
  name,
  angle1,
  angle2,
  angle3,
  rotation_order,
  location = NULL,
  vector_names = NULL,
  variable_names = NULL
)

orient_new_space(
  sys,
  name,
  angle1,
  angle2,
  angle3,
  rotation_order,
  location = NULL,
  vector_names = NULL,
  variable_names = NULL
)

orient_new_quaternion(
  sys,
  name,
  q0,
  q1,
  q2,
  q3,
  location = NULL,
  vector_names = NULL,
  variable_names = NULL
)

create_new(
  sys,
  name,
  transformation,
  variable_names = NULL,
  vector_names = NULL
)
```

## Arguments

- name:

  The name of the system.

- transformation:

  See Details.

- parent:

  The parent system, or `NULL`.

- location:

  The position of the origin, a vector in `parent`.

- rotation_matrix:

  The rotation matrix with respect to `parent`.

- vector_names, variable_names:

  Three names for the base vectors and base scalars.

- index:

  An integer 0, 1 or 2.

- system:

  A coordinate system.

- sys:

  A coordinate system, or (for `locate_new()`) a point.

- position:

  The position of the new origin or point.

- orienters:

  An orienter or a list of orienters.

- angle, angle1, angle2, angle3:

  Rotation angles.

- axis:

  The rotation axis, a vector.

- rotation_order:

  A rotation order string such as `"123"` or `"XYZ"`.

- q0, q1, q2, q3:

  Quaternion parameters.

## Value

A `CoordSys3D`, or a point for `locate_new()` on a point.

## Details

`locate_new()`, `orient_new()`, `orient_new_axis()`,
`orient_new_body()`, `orient_new_space()`, `orient_new_quaternion()` and
`create_new()` create child systems, as the methods of SymPy's
`CoordSys3D` do.

## Examples

``` r
N <- CoordSys3D("N")
N$x * N$i + N$y * N$j
#> N.x*N.i + N.y*N.j
q <- Symbol("q")
A <- orient_new_axis(N, "A", q, N$k)
express(A$i, N)
#> (cos(q))*N.i + (sin(q))*N.j
S <- create_new(N, "S", transformation = "spherical")
lame_coefficients(S)
#> (1, S.r, S.r*sin(S.theta)) 
```
