# Row and column operations

Functional versions of SymPy's `Matrix` row and column methods. Indices
are 1 based, as everywhere in R. `row_insert()` and `col_insert()`
insert the rows or columns of `other` before position `pos`; `row_del()`
and `col_del()` delete one; `row_join()` and `col_join()` concatenate;
`extract()` picks rows and columns; `permute()`, `permute_rows()` and
`permute_cols()` apply a permutation given as a list of swaps or as an
image vector (a SymPy `Permutation` is taken 0 based);
`elementary_row_op()` and `elementary_col_op()` apply `"n->kn"`,
`"n<->m"` or `"n->n+km"`; `row_swap()` and `col_swap()` swap two rows or
columns. `diagonal()`, `upper_triangular()` and `lower_triangular()`
take the `k`th diagonal and the triangles; `vech()` stacks the lower
triangle; `rot90()` rotates; `mat_row()`, `mat_col()`, `tolist()` and
`matrix_values()` read entries; `multiply_elementwise()` is the Hadamard
product.

## Usage

``` r
row_insert(M, pos, other)

col_insert(M, pos, other)

row_del(M, i)

col_del(M, j)

row_join(M, other)

col_join(M, other)

extract(M, rows, cols)

permute(M, perm, orientation = "rows", direction = "forward")

permute_rows(M, perm, direction = "forward")

permute_cols(M, perm, direction = "forward")

elementary_row_op(
  M,
  op = "n->kn",
  row = NULL,
  k = NULL,
  row1 = NULL,
  row2 = NULL
)

elementary_col_op(
  M,
  op = "n->kn",
  col = NULL,
  k = NULL,
  col1 = NULL,
  col2 = NULL
)

row_swap(M, i1, i2)

col_swap(M, i1, i2)

diagonal(M, k = 0)

upper_triangular(M, k = 0)

lower_triangular(M, k = 0)

vech(M, diagonal = TRUE, check_symmetry = TRUE)

rot90(M, k = 1)

mat_row(M, i)

mat_col(M, j)

tolist(M)

matrix_values(M)

multiply_elementwise(M, other)
```

## Arguments

- M:

  A `symmatrix`.

- pos, i, j, i1, i2:

  Positions (1 based).

- other:

  A `symmatrix`.

- rows, cols:

  Row and column indices, or logical vectors.

- perm:

  A list of swaps, an image vector or a `Permutation`.

- orientation:

  `"rows"` or `"cols"`.

- direction:

  `"forward"` or `"backward"`.

- op:

  The elementary operation.

- row, col, row1, row2, col1, col2:

  The rows or columns it acts on.

- k:

  A scale factor, or for `diagonal()` and the triangles the offset of
  the diagonal; for `rot90()` the number of quarter turns.

- diagonal:

  For `vech()`: include the diagonal.

- check_symmetry:

  For `vech()`: require a symmetric matrix.

## Value

A `symmatrix` (`tolist()` and `matrix_values()` return lists).

## Examples

``` r
M <- Matrix(list(list(1, 2), list(3, 4)))
row_insert(M, 2, Matrix(list(list(0, 0))))
#> Matrix(3, 2)
#> [1, 2]
#> [0, 0]
#> [3, 4]
col_del(M, 1)
#> Matrix(2, 1)
#> [2]
#> [4]
permute_rows(M, list(list(1, 2)))
#> Matrix(2, 2)
#> [3, 4]
#> [1, 2]
elementary_row_op(M, "n->n+km", row = 2, k = -3, row2 = 1)
#> Matrix(2, 2)
#> [1, 2 ]
#> [0, -2]
diagonal(M)
#> Matrix(1, 2)
#> [1, 4]
```
