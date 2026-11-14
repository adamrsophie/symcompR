# Minor submatrix

`minor_submatrix(M, i, j)` deletes row `i` and column `j` (1 based), as
SymPy's `Matrix.minor_submatrix()`.

## Usage

``` r
minor_submatrix(M, i, j)
```

## Arguments

- M:

  A `symmatrix`.

- i, j:

  The row and column to delete.

## Value

A `symmatrix`.

## Examples

``` r
minor_submatrix(Matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9))), 2, 2)
#> Matrix(2, 2)
#> [1, 3]
#> [7, 9]
```
