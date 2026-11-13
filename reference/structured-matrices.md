# Structured matrix expressions

`CompanionMatrix(poly)` is the companion matrix of a monic univariate
polynomial. `DFT(n)` and `IDFT(n)` are the unitary discrete Fourier
transform matrix and its inverse. `LofLU()`, `UofLU()`, `LofCholesky()`,
`UofCholesky()`, `QofQR()`, `RofQR()`, `EigenVectors()`,
`EigenValues()`, `UofSVD()`, `SofSVD()` and `VofSVD()` are the factors
of a matrix factorization, with the shape of their argument; `mat_lu()`,
`mat_qr()`, `mat_eig()` and `mat_svd()` give the factors as a list.

## Usage

``` r
CompanionMatrix(poly)

LofLU(arg)

UofLU(arg)

LofCholesky(arg)

UofCholesky(arg)

QofQR(arg)

RofQR(arg)

EigenVectors(arg)

EigenValues(arg)

UofSVD(arg)

SofSVD(arg)

VofSVD(arg)

mat_lu(expr)

mat_qr(expr)

mat_eig(expr)

mat_svd(expr)

DFT(n)

IDFT(n)
```

## Arguments

- poly:

  A monic univariate
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

- arg, expr:

  A matrix expression.

- n:

  The size, a nonnegative integer.

## Value

A `sym` matrix expression, or a list of them.

## Examples

``` r
x <- Symbol("x")
as_explicit(DFT(2))
#> Matrix(2, 2)
#> [sqrt(2)/2, sqrt(2)/2 ]
#> [sqrt(2)/2, -sqrt(2)/2]
QofQR(MatrixSymbol("X", 3, 3))
#> QofQR(X)
```
