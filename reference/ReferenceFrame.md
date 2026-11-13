# Reference frames for classical mechanics

`ReferenceFrame()` creates a frame with three orthonormal unit vectors,
as SymPy's `sympy.physics.vector.ReferenceFrame`. Frames are mutable and
compared by identity. Unit vectors are `N$x`, `N$y` and `N$z` (or
`N["x"]`), and `N[1]` to `N[3]` are the coordinate variables. Orient a
frame with the methods `orient_axis`, `orient_dcm`, `orient_explicit`,
`orient_body_fixed`, `orient_space_fixed`, `orient_quaternion`, `orient`
and `orientnew`, and query it with `dcm`, `ang_vel_in`, `ang_acc_in`,
`variable_map` and `partial_velocity`, all reached with `$`.

## Usage

``` r
ReferenceFrame(name, indices = NULL, latexs = NULL, variables = NULL)
```

## Arguments

- name:

  Frame name.

- indices:

  Optional three strings naming the unit vectors.

- latexs:

  Optional three LaTeX strings for the unit vectors.

- variables:

  Optional three names for the coordinate variables.

## Value

A `ReferenceFrame`.

## Examples

``` r
N <- ReferenceFrame("N")
q <- dynamicsymbols("q")
A <- N$orientnew("A", "Axis", list(q, N$z))
A$dcm(N)
#> Matrix(3, 3)
#> [cos(q(t)) , sin(q(t)), 0]
#> [-sin(q(t)), cos(q(t)), 0]
#> [0         , 0        , 1]
A$ang_vel_in(N)
#> Derivative(q(t), t)*N.z
```
