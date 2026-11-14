# Structural tests of a matrix

`is_upper()`, `is_lower()`, `is_diagonal()`, `is_hermitian()`,
`is_anti_symmetric()`, `is_upper_hessenberg()`, `is_lower_hessenberg()`,
`is_weakly_diagonally_dominant()`, `is_strongly_diagonally_dominant()`,
`is_Identity()`, `is_symbolic()` and `is_echelon()` follow the
properties and methods of SymPy's `MatrixBase`. The fuzzy tests return
`NA` where SymPy returns `None`.

## Usage

``` r
is_upper(M)

is_lower(M)

is_upper_hessenberg(M)

is_lower_hessenberg(M)

is_diagonal(M)

is_hermitian(M)

is_anti_symmetric(M, simplify = TRUE)

is_weakly_diagonally_dominant(M)

is_strongly_diagonally_dominant(M)

is_Identity(M)

is_symbolic(M)

is_echelon(M, iszerofunc = NULL)
```

## Arguments

- M:

  A `symmatrix`.

- simplify:

  For `is_anti_symmetric()`: simplify the entries before testing
  (`TRUE`), or a function to do so.

- iszerofunc:

  For `is_echelon()`: the zero test.

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
M <- Matrix(list(list(1, 2), list(0, 3)))
is_upper(M)
#> [1] TRUE
is_lower(M)
#> [1] FALSE
is_diagonal(M)
#> [1] FALSE
```
