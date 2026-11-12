# Matrix constructors

`sym_diag(...)` is SymPy's [`diag()`](https://rdrr.io/r/base/diag.html)
(block diagonal from scalars and matrices; renamed so it does not mask
[`base::diag()`](https://rdrr.io/r/base/diag.html)). `banded(...)`
builds a banded matrix from a named list of diagonals,
`list("0" = 1, "-1" = c(2, 3))`, with the size given as `banded(n, d)`
or `banded(r, c, d)` or inferred. `rot_givens(i, j, theta, dim)` is the
Givens rotation in the `(i, j)` plane (zero based). `sym_trace(M)` is
the trace (renamed so it does not mask
[`base::trace()`](https://rdrr.io/r/base/trace.html)).
`SparseMatrix(r, c, entries)` takes `list(list(i, j, value), ...)` with
zero based indices.

## Usage

``` r
sym_diag(...)

sym_trace(M)

banded(...)

rot_givens(i, j, theta, dim = 3)

SparseMatrix(r, c = NULL, entries = NULL)

MutableSparseMatrix(r, c = NULL, entries = NULL)

ImmutableSparseMatrix(r, c = NULL, entries = NULL)
```

## Arguments

- ...:

  Diagonal entries or blocks; for `banded()`, the size and the
  diagonals.

- M:

  A matrix or matrix expression.

- i, j:

  Plane indices.

- theta:

  The angle.

- dim:

  The dimension.

- r, c:

  Dimensions.

- entries:

  Nonzero entries.

## Value

A matrix.

## Examples

``` r
sym_diag(1, 2, 3)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 2, 0]
#> [0, 0, 3]
banded(3, list("0" = 1, "1" = 2))
#> Matrix(3, 3)
#> [1, 2, 0]
#> [0, 1, 2]
#> [0, 0, 1]
rot_givens(0, 1, pi/3)
#> Matrix(3, 3)
#> [0.5               , 0.866025403784439, 0]
#> [-0.866025403784439, 0.5              , 0]
#> [0                 , 0                , 1]
```
