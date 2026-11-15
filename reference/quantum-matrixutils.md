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

sympy_to_numpy(m, dtype = "complex")

to_numpy(m, dtype = "complex")

numpy_to_sympy(m)
```

## Arguments

- ...:

  Matrices.

- n:

  The number of rows.

- format:

  The matrix format: `"sympy"`, or `"numpy"` for a base R complex matrix
  (SymPy's NumPy format).

- m:

  The number of columns for `matrix_zeros()`; a matrix, or a number, for
  `sympy_to_numpy()`.

- e:

  A matrix or an expression.

- dtype:

  `"complex"` or `"double"`.

## Value

A matrix or an expression.

## Details

`sympy_to_numpy()` and `to_numpy()` turn a SymPy matrix into a base R
matrix, the counterpart of a NumPy array; `numpy_to_sympy()` turns a
base R matrix back into a `Matrix`.

## Examples

``` r
matrix_tensor_product(Matrix(list(list(1, 0), list(0, 1))), Matrix(list(list(0, 1), list(1, 0))))
#> Matrix(4, 4)
#> [0, 1, 0, 0]
#> [1, 0, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
```
