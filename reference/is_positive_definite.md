# Definiteness of a matrix

Fuzzy tests of whether a matrix is positive or negative (semi)definite,
or indefinite. A matrix that is not Hermitian is tested through
`M + M^H`. Each returns `TRUE`, `FALSE` or `NA` when the answer cannot
be decided.

## Usage

``` r
is_positive_definite(M)

is_positive_semidefinite(M)

is_negative_definite(M)

is_negative_semidefinite(M)

is_indefinite(M)
```

## Arguments

- M:

  A matrix.

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
is_positive_definite(Matrix(list(list(2, -1), list(-1, 2))))
#> [1] TRUE
is_indefinite(Matrix(list(list(1, 2), list(2, 1))))
#> [1] TRUE
```
