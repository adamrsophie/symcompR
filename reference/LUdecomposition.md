# Matrix decompositions

`LUdecomposition()` returns `list(L, U, perm)` with `P*M = L*U`;
`QRdecomposition()` returns `list(Q, R)` by Gram-Schmidt; `cholesky()`
returns `L` with `M = L*t(L)` for symmetric positive definite matrices;
`LDLdecomposition()` returns `list(L, D)`. `GramSchmidt()`
orthogonalizes a list of vectors.

## Usage

``` r
LUdecomposition(M)

GramSchmidt(vlist, orthonormal = FALSE)

QRdecomposition(M)

cholesky(M, hermitian = TRUE)

LDLdecomposition(M, hermitian = TRUE)
```

## Arguments

- M:

  A `symmatrix`.

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
