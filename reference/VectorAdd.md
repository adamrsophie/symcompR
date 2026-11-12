# Vectors and dyadics of a coordinate system

`VectorAdd()`, `VectorMul()`, `DyadicAdd()` and `DyadicMul()` build sums
and scalar multiples of vectors and dyadics, as the arithmetic operators
do. `VectorZero()` and `DyadicZero()` are the zero vector and the zero
dyadic (`Vector.zero` and `Dyadic.zero` in SymPy).
`BaseVector(index, system)` is a basis vector of a
[`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
and `BaseDyadic(v1, v2)` the outer product of two basis vectors.
`measure_number()` and
[`base_vector()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
give the parts of a `VectorMul`;
[`components()`](https://adamrsophie.github.io/symcompR/reference/categories.md)
the components of a vector or dyadic as a named list.

## Usage

``` r
VectorAdd(...)

VectorMul(...)

VectorZero()

DyadicAdd(...)

DyadicMul(...)

DyadicZero()

BaseVector(index, system)

BaseDyadic(vector1, vector2)

measure_number(v)
```

## Arguments

- ...:

  Vectors or dyadics, and scalars.

- index:

  An integer 0, 1 or 2.

- system:

  A
  [`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md).

- vector1, vector2:

  Base vectors.

- v:

  A vector or dyadic.

## Value

A vector or dyadic.

## Details

Vectors support `+`, `-`, multiplication and division by scalars, `&`
for the dot product and `|` for the outer product; `^` raises an error
as SymPy's `**` does. Use
[`cross()`](https://adamrsophie.github.io/symcompR/reference/dot.md) for
the cross product.

## Examples

``` r
C <- CoordSys3D("C")
v <- 3 * C$i + C$x * C$j
v
#> 3*C.i + C.x*C.j
VectorAdd(C$i, C$j, C$i)
#> 2*C.i + C.j
VectorMul(2, C$k)
#> 2*C.k
C$i | C$j
#> (C.i|C.j)
```
