# Inverse by a named method

`inverse_GE()`, `inverse_LU()`, `inverse_ADJ()`, `inverse_CH()`,
`inverse_LDL()`, `inverse_QR()` and `inverse_BLOCK()` are
`inverse(M, method = )` for each of SymPy's inversion methods.

## Usage

``` r
inverse_GE(M, iszerofunc = NULL)

inverse_LU(M, iszerofunc = NULL)

inverse_ADJ(M, iszerofunc = NULL)

inverse_CH(M, iszerofunc = NULL)

inverse_LDL(M, iszerofunc = NULL)

inverse_QR(M, iszerofunc = NULL)

inverse_BLOCK(M, iszerofunc = NULL)
```

## Arguments

- M:

  A square `symmatrix`.

- iszerofunc:

  The zero test.

## Value

A `symmatrix`.

## Examples

``` r
A <- Matrix(list(list(2, 1), list(1, 3)))
inverse_CH(A)
#> Matrix(2, 2)
#> [3/5 , -1/5]
#> [-1/5, 2/5 ]
inverse_ADJ(A)
#> Matrix(2, 2)
#> [3/5 , -1/5]
#> [-1/5, 2/5 ]
```
