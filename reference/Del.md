# The del operator

`Del()` returns the vector differential operator. Called on a scalar
field it gives the gradient; `delop & v` gives the divergence,
`cross(delop, v)` the curl and `v & delop` a function computing the
directional derivative along `v`. Results stay unevaluated unless
`doit = TRUE`.

## Usage

``` r
Del()

del_gradient(delop, field, doit = FALSE)

del_dot(delop, field, doit = FALSE)

del_cross(delop, field, doit = FALSE)
```

## Arguments

- delop:

  The del operator.

- field:

  A scalar or vector field.

- doit:

  Evaluate the derivatives.

## Value

A function of class `"Del"`.

## Examples

``` r
C <- CoordSys3D("C")
delop <- Del()
delop(C$x * C$y * C$z)
#> (Derivative(C.x*C.y*C.z, C.x))*C.i + (Derivative(C.x*C.y*C.z, C.y))*C.j + (Derivative(C.x*C.y*C.z, C.z))*C.k
doit(delop & (C$x * C$i))
#> 1
```
