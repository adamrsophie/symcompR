# Inertia dyadics and the Inertia tuple

`inertia()` builds the inertia dyadic with the given moments and
products of inertia in `frame`, as SymPy's `inertia()`.
`inertia_of_point_mass()` is the inertia dyadic of a point mass at
`pos_vec`. `Inertia()` pairs an inertia dyadic with its reference point,
as SymPy's `Inertia` named tuple; `Inertia_from_inertia_scalars()` is
SymPy's `Inertia.from_inertia_scalars()`.

## Usage

``` r
inertia(frame, ixx, iyy, izz, ixy = 0, iyz = 0, izx = 0)

inertia_of_point_mass(mass, pos_vec, frame)

Inertia(dyadic, point)

Inertia_from_inertia_scalars(
  point,
  frame,
  ixx,
  iyy,
  izz,
  ixy = 0,
  iyz = 0,
  izx = 0
)
```

## Arguments

- frame:

  A
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md).

- ixx, iyy, izz:

  Moments of inertia.

- ixy, iyz, izx:

  Products of inertia.

- mass:

  Mass of the point.

- pos_vec:

  Position vector of the point from the reference point.

- dyadic:

  An inertia dyadic.

- point:

  The reference point, a
  [`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md).

## Value

A `PhysDyadic`, or an `Inertia` tuple.

## Examples

``` r
N <- ReferenceFrame("N")
inertia(N, 1, 2, 3)
#> (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)
O <- PhysPoint("O")
Inertia(inertia(N, 1, 2, 3), O)
#> Inertia(dyadic=(N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z), point=O)
```
