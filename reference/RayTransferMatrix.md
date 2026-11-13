# Ray transfer matrices of Gaussian optics

`RayTransferMatrix(A, B, C, D)` (or a 2 by 2 matrix) is the ABCD matrix
of an optical element; `FreeSpace(d)`, `FlatRefraction(n1, n2)`,
`CurvedRefraction(R, n1, n2)`, `FlatMirror()`, `CurvedMirror(R)` and
`ThinLens(f)` are the matrices of common elements.
`GeometricRay(h, angle)` (or a 2 by 1 matrix) is a ray. The product
`%*%` of a matrix with another matrix, a ray or a
[`BeamParameter()`](https://adamrsophie.github.io/symcompR/reference/BeamParameter.md)
gives the transformed object. `M$A`, `M$B`, `M$C` and `M$D` are the
entries; `ray$height` and `ray$angle` the parts of a ray. After
`sympy.physics.optics.gaussopt`.

## Usage

``` r
RayTransferMatrix(...)

FreeSpace(d)

FlatRefraction(n1, n2)

CurvedRefraction(R, n1, n2)

FlatMirror()

CurvedMirror(R)

ThinLens(f)

GeometricRay(...)
```

## Arguments

- ...:

  The four entries, or one 2 by 2 matrix.

- d:

  A distance.

- n1, n2:

  The refractive indices before and after the interface.

- R:

  The radius of curvature.

- f:

  The focal distance.

## Value

A `symmatrix` of class `RayTransferMatrix` (or `GeometricRay`).

## Examples

``` r
d <- Symbol("d"); f <- Symbol("f")
ThinLens(f) %*% FreeSpace(d)
#> Matrix(2, 2)
#> [1   , d       ]
#> [-1/f, -d/f + 1]
FreeSpace(d) %*% GeometricRay(Symbol("h"), Symbol("angle"))
#> Matrix(2, 1)
#> [angle*d + h]
#> [angle      ]
```
