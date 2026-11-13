# Points of classical mechanics

`PhysPoint()` creates a point of `sympy.physics.vector` (SymPy's
`Point`; the name differs from the geometry
[`Point()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
of this package). Points are mutable: `set_pos`, `set_vel` and `set_acc`
record position, velocity and acceleration, `pos_from`, `vel` and `acc`
query them, and `locatenew`, `v1pt_theory`, `v2pt_theory`,
`a1pt_theory`, `a2pt_theory` and `partial_velocity` follow SymPy. All
are reached with `$`.

## Usage

``` r
PhysPoint(name)
```

## Arguments

- name:

  Point name.

## Value

A `PhysPoint`.

## Examples

``` r
N <- ReferenceFrame("N")
O <- PhysPoint("O")
q <- dynamicsymbols("q")
P <- O$locatenew("P", q * N$x)
O$set_vel(N, 0)
P$vel(N)
#> Derivative(q(t), t)*N.x
```
