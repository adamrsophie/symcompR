# Matrix decompositions

`LUdecomposition()` returns `list(L, U, perm)` with `P*M = L*U`, where
`perm` lists the row swaps (1 based) that make `P`; `QRdecomposition()`
returns `list(Q, R)` by Gram-Schmidt; `cholesky()` returns `L` with
`M = L*t(L)` for symmetric positive definite matrices;
`LDLdecomposition()` returns `list(L, D)`. `GramSchmidt()`
orthogonalizes a list of vectors.

## Usage

``` r
LUdecomposition(M, iszerofunc = NULL, simpfunc = NULL, rankcheck = FALSE)

GramSchmidt(vlist, orthonormal = FALSE)

QRdecomposition(M)

cholesky(M, hermitian = TRUE)

LDLdecomposition(M, hermitian = TRUE)
```

## Arguments

- M:

  A `symmatrix`.

- iszerofunc, simpfunc, rankcheck:

  For `LUdecomposition()`: the zero test, an optional simplification for
  the pivot search, and whether to raise for a rank deficient matrix.

- vlist:

  A list of column vectors.

- orthonormal:

  Normalize the vectors.

- hermitian:

  For `cholesky()` and `LDLdecomposition()`: require a Hermitian matrix
  and use conjugates (SymPy's default); `FALSE` requires a symmetric
  matrix.

## Value

A list of matrices.
