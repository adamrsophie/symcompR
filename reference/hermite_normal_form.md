# Hermite and Smith normal forms

`hermite_normal_form(A, D)` is the Hermite normal form of an integer
matrix (SymPy's `hermite_normal_form`, optionally modulo `D`);
`smith_normal_form(A)` and `invariant_factors(A)` give the Smith normal
form and its diagonal.

## Usage

``` r
hermite_normal_form(A, D = NULL, check_rank = FALSE)

invariant_factors(A)

smith_normal_form(A)
```

## Arguments

- A:

  An integer matrix.

- D:

  Optional positive integer modulus (a multiple of the determinant).

- check_rank:

  With `D`, use the modulo method only for full row rank.

## Value

A matrix, or for `invariant_factors()` a Tuple of integers.

## Examples

``` r
A <- Matrix(list(list(2, 7, 17, 29, 41), list(3, 11, 19, 31, 43), list(5, 13, 23, 37, 47)))
hermite_normal_form(A)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 2, 1]
#> [0, 0, 1]
smith_normal_form(Matrix(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14))))
#> Matrix(3, 3)
#> [1, 0 , 0 ]
#> [0, 10, 0 ]
#> [0, 0 , 30]
```
