# The set of matrices of a shape with entries in a set

`MatrixSet(n, m, set)` is the set of `n` by `m` matrices whose entries
lie in `set`. `matrix_set_contains(M, X)` tests membership, as SymPy's
`X in M`.

## Usage

``` r
MatrixSet(n, m, set)

matrix_set_contains(M, X)
```

## Arguments

- n, m:

  The shape.

- set:

  A set.

- M:

  A `MatrixSet`.

- X:

  A matrix or matrix expression.

## Value

A `sym` set, or `TRUE`, `FALSE` or `NA`.

## Examples

``` r
M <- MatrixSet(2, 2, Reals)
matrix_set_contains(M, Matrix(list(list(1, 2), list(3, 4))))
#> [1] TRUE
```
