# Symbolic matrices

`Matrix()` builds a matrix of expressions from a list of rows, from an R
matrix, from a vector (a column), or from dimensions and a function
`f(i, j)` or a data list in row-major order. `zeros()`, `ones()`,
`eye()` and `diag_matrix()` build special matrices. `+`, `-` and scalar
`*` work as usual; `*` between two matrices of equal shape is
elementwise, `%*%` is the matrix product and `%^%` the matrix power. `[`
indexes like R matrices.

## Usage

``` r
Matrix(data, ncol = NULL, f = NULL)

zeros(n, m = n)

ones(n, m = n)

eye(n, m = n)

diag_matrix(...)

A %^% n
```

## Arguments

- data:

  A list of rows, an R matrix, a vector, or a number of rows.

- ncol:

  Number of columns, when `data` is a number of rows.

- f:

  A function of `(i, j)` or a list of entries in row-major order.

- n, m:

  Dimensions.

- ...:

  Diagonal entries for `diag_matrix()`.

- A:

  A matrix.

## Value

A `symmatrix`.

## Examples

``` r
x <- Symbol("x")
M <- Matrix(list(list(1, x), list(x, 1)))
sym_det(M)
#> 1 - x^2
M %*% M
#> Matrix(2, 2)
#> [x^2 + 1, 2*x    ]
#> [2*x    , x^2 + 1]
Matrix(2, 2, function(i, j) i + j)
#> Matrix(2, 2)
#> [2, 3]
#> [3, 4]
```
