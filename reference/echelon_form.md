# Echelon forms

`echelon_form(M)` is the row echelon form of SymPy's
`Matrix.echelon_form()` (with `with_pivots = TRUE` also the 1 based
pivot columns); `rref_rhs(M, rhs)` reduces `M` augmented with the
identity and `rhs` and returns `list(reduced M, reduced rhs)`.

## Usage

``` r
echelon_form(M, iszerofunc = NULL, simplify = FALSE, with_pivots = FALSE)

rref_rhs(M, rhs)
```

## Arguments

- M:

  A `symmatrix`.

- iszerofunc:

  The zero test.

- simplify:

  A simplification function, or `FALSE` for
  [`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md).

- with_pivots:

  Also return the pivot columns.

- rhs:

  The right hand side.

## Value

A `symmatrix` or a list.

## Examples

``` r
echelon_form(Matrix(list(list(1, 2), list(3, 4))))
#> Matrix(2, 2)
#> [1, 2 ]
#> [0, -2]
rref_rhs(Matrix(list(list(1, 2), list(3, 4))), Matrix(list(5, 6)))
#> [[1]]
#> Matrix(2, 2)
#> [1, 0]
#> [0, 1]
#> 
#> [[2]]
#> Matrix(2, 1)
#> [-4 ]
#> [9/2]
#> 
```
