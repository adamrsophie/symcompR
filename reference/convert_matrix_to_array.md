# Convert between matrix and array expressions

`convert_matrix_to_array()` rewrites a matrix expression as an array
expression: products become contractions of tensor products, transposes
become axis permutations, traces contractions and Hadamard products
diagonals. `convert_array_to_matrix()` recognizes matrix expressions in
an array expression, as SymPy's functions of the same names.

## Usage

``` r
convert_matrix_to_array(expr)

convert_array_to_matrix(expr)
```

## Arguments

- expr:

  A matrix or array expression.

## Value

An array or matrix expression.

## Examples

``` r
k <- Symbol("k")
A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
cg <- convert_matrix_to_array(A * t(B))
cg
#> ArrayContraction(ArrayTensorProduct(A, B), (1, 3))
convert_array_to_matrix(cg)
#> A*B.T
```
