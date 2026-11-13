# Coordinate variables of a reference frame

A `CoordinateSym` is a symbol tied to a
[`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
and one of its three axes, as SymPy's
`sympy.physics.vector.CoordinateSym`. Each frame creates its own three;
`N[1]` returns the first. Two coordinate symbols are equal only when
frame and index agree.

## Usage

``` r
CoordinateSym(name, frame, index)
```

## Arguments

- name:

  Display name.

- frame:

  The `ReferenceFrame` it belongs to.

- index:

  The axis, 0, 1 or 2 as in SymPy.

## Value

A `sym` symbol.

## Examples

``` r
A <- ReferenceFrame("A")
A[1]
#> A_x
CoordinateSym("A_x", A, 0)
#> A_x
```
