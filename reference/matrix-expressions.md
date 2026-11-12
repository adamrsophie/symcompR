# Matrix expressions

`MatrixSymbol(name, n, m)` is a symbolic `n` by `m` matrix.
`Identity(n)`, `ZeroMatrix(m, n)` and `OneMatrix(m, n)` are the special
matrices. `+`, `-`, `*` and `^` build the canonical `MatAdd`, `MatMul`
and `MatPow` (like terms combine, identities and zeros simplify, shapes
are checked and mismatches raise a `ShapeError`).
[`t()`](https://rdrr.io/r/base/t.html) transposes and distributes;
`Transpose()`, `Inverse()` and `Adjoint()` build unevaluated nodes
(`A^-1` is `Inverse(A)`). `Trace()`, `Determinant()` and `Permanent()`
are scalar functions of a matrix;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
applies their rules. `MatrixElement(A, i, j)` is the entry with SymPy's
zero based indices. `as_explicit()` turns an expression with integer
shape into an explicit
[`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md).
`mat_shape()` gives the shape and `matrix_symbols()` the matrix symbols.

## Usage

``` r
MatrixSymbol(name, n, m)

Identity(n)

ZeroMatrix(n, m)

OneMatrix(n, m)

mat_shape(e)

MatMul(..., evaluate = TRUE)

MatAdd(..., evaluate = TRUE)

MatPow(A, k)

Transpose(A)

Inverse(A)

Adjoint(A)

Trace(A)

Determinant(A)

Permanent(A)

HadamardProduct(...)

HadamardPower(A, k)

KroneckerProduct(...)

MatrixElement(A, i, j)

FunctionMatrix(n, m, f)

DiagMatrix(A)

diagonalize_vector(A)

DiagonalOf(A)

DiagonalMatrix(A)

PermutationMatrix(perm)

MatrixPermute(A, perm, axis = 0)

DotProduct(A, B)

MatrixSlice(A, rowslice, colslice)

blockcut(A, rowsizes, colsizes)

BlockMatrix(blocks)

BlockDiagMatrix(...)

block_collapse(e)

matrix_symbols(e)

as_explicit(e)
```

## Arguments

- name:

  The name.

- n, m:

  Dimensions.

- e, A, B:

  A matrix expression.

- ...:

  Matrices, or for `MatMul()` matrices and scalars.

- evaluate:

  If `FALSE`, `MatAdd()` and `MatMul()` keep their arguments as given,
  as SymPy's default; explicit matrices are held as arguments either way
  when symbolic matrices are present.

- k:

  An exponent.

- i, j:

  Zero based indices.

- f:

  A
  [`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md)
  of two variables giving the entries.

- perm:

  A
  [`Permutation()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
  or a vector of zero based images.

- axis:

  0 to permute rows, 1 for columns.

- rowslice, colslice:

  A zero based index, `c(start, stop)` or `c(start, stop, step)`; `NA`
  (or `NULL` in a list) takes the default end, and negative values count
  from the end, as Python slices.

- rowsizes, colsizes:

  Block sizes along rows and columns.

- blocks:

  A list of rows, each a list of matrices.

## Value

A matrix expression, a scalar expression or an explicit matrix.

## Examples

``` r
A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
A * B + A * B
#> 2*A*B
t(A * B)
#> B.T*A.T
A^-1 * A
#> I
as_explicit(A * B)
#> Matrix(2, 2)
#> [A[0, 0]*B[0, 0] + A[0, 1]*B[1, 0], A[0, 0]*B[0, 1] + A[0, 1]*B[1, 1]]
#> [A[1, 0]*B[0, 0] + A[1, 1]*B[1, 0], A[1, 0]*B[0, 1] + A[1, 1]*B[1, 1]]
doit(Trace(A + B))
#> Trace(A) + Trace(B)
```
