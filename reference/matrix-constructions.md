# Matrix constructions

`rot_axis1()`, `rot_axis2()` and `rot_axis3()` are the rotation matrices
used by SymPy (clockwise about the axis when looking down it);
`rot_ccw_axis1()` and friends rotate counterclockwise.
`jordan_cell(e, n)` builds a Jordan block, `casoratian(seqs, n)` the
Casoratian determinant, `randMatrix(r, c, min, max)` a random integer
matrix, `symarray(prefix, shape)` a matrix of symbols, `permanent()`
(alias `per()`) the permanent, and `matrix_multiply_elementwise()` the
Hadamard product. `ImmutableMatrix()`, `MutableMatrix()`,
`MutableDenseMatrix()` and `ImmutableDenseMatrix()` are
[`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md).

## Usage

``` r
rot_axis1(theta)

rot_axis2(theta)

rot_axis3(theta)

rot_ccw_axis1(theta)

rot_ccw_axis2(theta)

rot_ccw_axis3(theta)

jordan_cell(e, n)

casoratian(seqs, n, zero = TRUE)

randMatrix(r, c = r, min = 0, max = 99)

symarray(prefix, shape)

permanent(M)

per(M)

matrix_multiply_elementwise(A, B)

ImmutableMatrix(...)

MutableMatrix(...)

MutableDenseMatrix(...)

ImmutableDenseMatrix(...)
```

## Arguments

- theta:

  Angle.

- e:

  Diagonal value.

- n:

  Size, or the index variable for `casoratian()`.

- seqs:

  A list of sequences in `n`.

- zero:

  Evaluate at `n = 0, 1, ...` (SymPy's default) instead of
  `n, n + 1, ...`.

- r, c:

  Dimensions.

- min, max:

  Entry range.

- prefix:

  Symbol prefix.

- shape:

  Dimensions.

- M, A, B:

  Matrices.

- ...:

  Arguments for
  [`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md).

## Value

A matrix or expression.

## Examples

``` r
rot_axis3(sym_pi / 2)
#> Matrix(3, 3)
#> [0 , 1, 0]
#> [-1, 0, 0]
#> [0 , 0, 1]
jordan_cell(Symbol("a"), 3)
#> Matrix(3, 3)
#> [a, 1, 0]
#> [0, a, 1]
#> [0, 0, a]
permanent(Matrix(list(list(1, 2), list(3, 4))))
#> 10
```
