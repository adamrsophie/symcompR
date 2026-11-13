# Matrix utilities of the quantum module

The SymPy matrix helpers of `sympy.physics.quantum.matrixutils`:
`matrix_tensor_product()` is the Kronecker product of its arguments,
`matrix_eye(n)` and `matrix_zeros(m, n)` the identity and zero matrices,
`matrix_dagger(m)` the conjugate transpose, `to_sympy(m)` (one argument)
the identity on symbolic matrices and `flatten_scalar(m)` the entry of a
1 by 1 matrix. Only the `"sympy"` format exists here; NumPy and SciPy
formats raise an error.

## Usage

``` r
matrix_tensor_product(...)

matrix_eye(n, format = "sympy")

matrix_zeros(m, n, format = "sympy")

matrix_dagger(e)

matrix_to_zero(e)
```

## Arguments

- ...:

  Matrices.

- n, m:

  Dimensions.

- format:

  The matrix format, only `"sympy"`.

- e:

  A matrix or an expression.

## Value

A matrix or an expression.

## Examples

``` r
matrix_tensor_product(Matrix(list(list(1, 0), list(0, 1))), Matrix(list(list(0, 1), list(1, 0))))
#> Matrix(4, 4)
#> [0, 1, 0, 0]
#> [1, 0, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
```
