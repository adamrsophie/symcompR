# Power of a square matrix

`matrix_pow(M, n, method)` is SymPy's `Matrix.pow()`: integer powers by
repeated squaring with intermediate simplification, symbolic and
rational powers from the Jordan normal form
(`_matrix_pow_by_jordan_blocks`), so a matrix that is not diagonalizable
works too. `M %^% n` calls it.

## Usage

``` r
matrix_pow(M, n, method = NULL)
```

## Arguments

- M:

  A square `symmatrix`.

- n:

  The exponent.

- method:

  `NULL` (choose), `"multiply"`, `"mulsimp"`, `"jordan"` or `"cayley"`.

## Value

A `symmatrix`, or an unevaluated `MatPow`.

## Examples

``` r
n <- Symbol("n")
matrix_pow(Matrix(list(list(1, 1), list(0, 1))), n)
#> Matrix(2, 2)
#> [1, n]
#> [0, 1]
Matrix(list(list(2, 1), list(0, 2))) %^% n
#> Matrix(2, 2)
#> [2^n, 2^(n - 1)*n]
#> [0  , 2^n        ]
```
