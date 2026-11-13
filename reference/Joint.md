# Joints between bodies

The joints of SymPy's `sympy.physics.mechanics`: `PinJoint()` (one
rotation about `joint_axis`), `PrismaticJoint()` (one translation along
`joint_axis`), `CylindricalJoint()` (a rotation and a translation along
the same axis), `PlanarJoint()` (a rotation about the parent
interframe's x axis and two translations in its y z plane),
`SphericalJoint()` (three rotations) and `WeldJoint()` (no motion). Each
joint orients the child frame and sets the velocities of the attachment
points when it is created. `Joint()` is the abstract base and cannot be
instantiated. Attributes (`coordinates`, `speeds`, `kdes`,
`parent_point`, `child_point`, `parent_interframe`, `child_interframe`,
`joint_axis` and others) are reached with `$`.

## Usage

``` r
Joint(
  name,
  parent,
  child,
  coordinates = NULL,
  speeds = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL,
  parent_axis = NULL,
  child_axis = NULL,
  parent_joint_pos = NULL,
  child_joint_pos = NULL
)

PinJoint(
  name,
  parent,
  child,
  coordinates = NULL,
  speeds = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL,
  parent_axis = NULL,
  child_axis = NULL,
  joint_axis = NULL,
  parent_joint_pos = NULL,
  child_joint_pos = NULL
)

PrismaticJoint(
  name,
  parent,
  child,
  coordinates = NULL,
  speeds = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL,
  parent_axis = NULL,
  child_axis = NULL,
  joint_axis = NULL,
  parent_joint_pos = NULL,
  child_joint_pos = NULL
)

CylindricalJoint(
  name,
  parent,
  child,
  rotation_coordinate = NULL,
  translation_coordinate = NULL,
  rotation_speed = NULL,
  translation_speed = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL,
  joint_axis = NULL
)

PlanarJoint(
  name,
  parent,
  child,
  rotation_coordinate = NULL,
  planar_coordinates = NULL,
  rotation_speed = NULL,
  planar_speeds = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL
)

SphericalJoint(
  name,
  parent,
  child,
  coordinates = NULL,
  speeds = NULL,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL,
  rot_type = "BODY",
  amounts = NULL,
  rot_order = 123
)

WeldJoint(
  name,
  parent,
  child,
  parent_point = NULL,
  child_point = NULL,
  parent_interframe = NULL,
  child_interframe = NULL
)
```

## Arguments

- name:

  Joint name.

- parent, child:

  The parent and child bodies.

- coordinates, speeds:

  Generalized coordinates and speeds; created with
  [`dynamicsymbols()`](https://adamrsophie.github.io/symcompR/reference/dynamicsymbols.md)
  when `NULL`.

- parent_point, child_point:

  Attachment points, as points or position vectors from the mass
  centers.

- parent_interframe, child_interframe:

  Intermediate frames, as frames or vectors to align with the x axis.

- parent_axis, child_axis, parent_joint_pos, child_joint_pos:

  Deprecated arguments of SymPy.

- joint_axis:

  The axis of rotation or translation.

- rotation_coordinate, translation_coordinate, rotation_speed,
  translation_speed:

  Coordinates and speeds of a cylindrical or planar joint.

- planar_coordinates, planar_speeds:

  The two planar coordinates and speeds of a planar joint.

- rot_type, amounts, rot_order:

  Orientation of a spherical joint.

## Value

A joint object.

## Examples

``` r
P <- RigidBody("P")
C <- RigidBody("C")
J <- PinJoint("J", P, C)
J$coordinates
#> Matrix(1, 1)
#> [q_J(t)]
J$kdes
#> Matrix(1, 1)
#> [u_J(t) - Derivative(q_J(t), t)]
C$frame$dcm(P$frame)
#> Matrix(3, 3)
#> [1, 0           , 0          ]
#> [0, cos(q_J(t)) , sin(q_J(t))]
#> [0, -sin(q_J(t)), cos(q_J(t))]
```
