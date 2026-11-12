# Vector operations

`magnitude()`, `normalize()`,
[`projection()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md),
`to_matrix()`, `separate()` and `vector_equals()` port the methods of
SymPy's `Vector`. `outer_sym()` is the outer product (also `v1 | v2`).
[`dot()`](https://adamrsophie.github.io/symcompR/reference/dot.md) and
[`cross()`](https://adamrsophie.github.io/symcompR/reference/dot.md)
accept vectors and dyadics as well as matrices.

## Usage

``` r
magnitude(v)

normalize(v)

outer_sym(v, other)

vector_equals(v, other)

to_matrix(v, system, second_system = NULL)

separate(v)
```

## Arguments

- v, other:

  Vectors.

- system, second_system:

  Coordinate systems
  ([`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)).

## Value

An expression, vector, dyadic, matrix or list.

## Examples

``` r
C <- CoordSys3D("C")
v <- 3 * C$i + 4 * C$j
magnitude(v)
#> 5
normalize(v)
#> 3/5*C.i + 4/5*C.j
to_matrix(v, C)
#> Matrix(3, 1)
#> [3]
#> [4]
#> [0]
outer_sym(C$i, C$j)
#> (C.i|C.j)
```
