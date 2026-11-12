# Polyhedra

`Polyhedron(corners, faces, pgroup)` is SymPy's polyhedron: corners,
faces (stored in canonical rotation order, as a sorted set) and a group
of permutations of the corners. `poly_rotate()` applies a permutation
(or the generator of `pgroup` with the given zero based index) to the
corners in place, `poly_reset()` restores them and `poly_copy()` copies.
`poly_corners()`, `poly_faces()`, `poly_edges()`, `poly_pgroup()`,
`poly_size()`, `poly_array_form()` and `poly_cyclic_form()` are the
properties. `tetrahedron()`, `cube()`, `octahedron()`, `dodecahedron()`
and `icosahedron()` return new copies of the standard polyhedra, and
`polyhedron_faces(name)` gives their faces.

## Usage

``` r
Polyhedron(corners, faces = list(), pgroup = list())

poly_corners(P)

poly_faces(P)

poly_pgroup(P)

poly_size(P)

poly_array_form(P)

poly_cyclic_form(P)

poly_edges(P)

poly_rotate(P, perm)

poly_reset(P)

poly_copy(P)

tetrahedron()

cube()

octahedron()

dodecahedron()

icosahedron()

polyhedron_faces(name)
```

## Arguments

- corners:

  A vector or list of corners.

- faces:

  A list of faces (vectors of corners).

- pgroup:

  A list of permutations.

- P:

  A polyhedron.

- perm:

  A permutation, or a zero based index into the group.

- name:

  One of `"tetrahedron"`, `"cube"`, `"octahedron"`, `"dodecahedron"`,
  `"icosahedron"`.

## Value

A `Polyhedron` (an environment) or one of its properties.

## Examples

``` r
C <- cube()
poly_rotate(C, 0)
fp_str(poly_corners(C))
#> [1] "(1, 2, 3, 0, 5, 6, 7, 4)"
fp_str(poly_edges(tetrahedron()))
#> [1] "{(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}"
```
