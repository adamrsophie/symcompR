# More matrix decompositions

`LUdecomposition_Simple()` returns `list(LU, perm)`: the combined `L`
and `U` factors and the row swaps; `LUdecompositionFF()` the fraction
free `list(P, L, DD, U)` with `P*M = L*DD^-1*U`; `rank_decomposition()`
`list(C, F)` with `M = C*F`; `singular_value_decomposition()`
`list(U, S, V)` with `M = U*S*V^H`; `upper_hessenberg_decomposition()`
`list(H, P)` with `M = P*H*P^H`; `bidiagonal_decomposition()`
`list(U, B, V)` and `bidiagonalize()` the bidiagonal form, all as
SymPy's `Matrix` methods. Row swaps are 1 based.

## Usage

``` r
LUdecomposition_Simple(
  M,
  iszerofunc = NULL,
  simpfunc = NULL,
  rankcheck = FALSE
)

LUdecompositionFF(M)

rank_decomposition(M, iszerofunc = NULL, simplify = FALSE)

singular_value_decomposition(M)

upper_hessenberg_decomposition(M)

bidiagonal_decomposition(M, upper = TRUE)

bidiagonalize(M, upper = TRUE)
```

## Arguments

- M:

  A `symmatrix`.

- iszerofunc:

  The zero test (`x$is_zero` by default).

- simpfunc:

  An optional simplification function for pivot search.

- rankcheck:

  Raise when the matrix is rank deficient.

- simplify:

  For `rank_decomposition()`: simplify during row reduction.

- upper:

  For the bidiagonal forms: upper (`TRUE`) or lower.

## Value

A list of matrices, or a matrix.

## Examples

``` r
M <- Matrix(list(list(1, 2), list(3, 4)))
LUdecomposition_Simple(M)
#> [[1]]
#> Matrix(2, 2)
#> [1, 2 ]
#> [3, -2]
#> 
#> [[2]]
#> list()
#> 
rank_decomposition(Matrix(list(list(1, 2), list(2, 4))))
#> $C
#> Matrix(2, 1)
#> [1]
#> [2]
#> 
#> $F
#> Matrix(1, 2)
#> [1, 2]
#> 
upper_hessenberg_decomposition(Matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9))))$H
#> Matrix(3, 3)
#> [1                                                                                                                                , -(4 + sqrt(65))*(42/sqrt(49 + (4 + sqrt(65))^2) + 4*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 2                                                                                                                                                                                                                                                                                                                                                                                                       , -7*(42/sqrt(49 + (4 + sqrt(65))^2) + 4*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 3                                                                                                                                                                                                                                                                                                                                                                                                       ]
#> [-2*(4 + sqrt(65))*(4*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2) + 49/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 4, -2*(4 + sqrt(65))*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 5 - (4 + sqrt(65))*(2*(4 + sqrt(65))*(-2*(4 + sqrt(65))*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 5)/sqrt(49 + (4 + sqrt(65))^2) + 14*(-2*(4 + sqrt(65))*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 6)/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2), -2*(4 + sqrt(65))*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 6 - 7*(2*(4 + sqrt(65))*(-2*(4 + sqrt(65))*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 5)/sqrt(49 + (4 + sqrt(65))^2) + 14*(-2*(4 + sqrt(65))*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 6)/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2)]
#> [-14*(4*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2) + 49/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 7              , -14*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) - (4 + sqrt(65))*(14*(-14*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 9)/sqrt(49 + (4 + sqrt(65))^2) + 2*(4 + sqrt(65))*(-14*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 8)/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 8                                          , -14*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) - 7*(14*(-14*(63/sqrt(49 + (4 + sqrt(65))^2) + 6*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 9)/sqrt(49 + (4 + sqrt(65))^2) + 2*(4 + sqrt(65))*(-14*(56/sqrt(49 + (4 + sqrt(65))^2) + 5*(4 + sqrt(65))/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 8)/sqrt(49 + (4 + sqrt(65))^2))/sqrt(49 + (4 + sqrt(65))^2) + 9                                          ]
```
