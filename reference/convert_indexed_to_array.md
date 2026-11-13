# Convert between indexed and array expressions

`convert_indexed_to_array(expr, first_indices)` turns an expression in
indexed matrix or array elements, with sums over repeated indices, into
an array expression whose axes follow `first_indices`.
`convert_array_to_indexed(expr, indices)` writes an array expression as
an indexed expression with the given free indices.

## Usage

``` r
convert_indexed_to_array(expr, first_indices = NULL)

convert_array_to_indexed(expr, indices)
```

## Arguments

- expr:

  An expression.

- first_indices:

  Optional list of indices giving the order of the axes.

- indices:

  A list of index symbols, one per axis.

## Value

An array expression or an indexed expression.

## Examples

``` r
k <- Symbol("k"); i <- Symbol("i"); j <- Symbol("j"); d <- Symbol("d")
M <- MatrixSymbol("M", k, k); N <- MatrixSymbol("N", k, k)
convert_indexed_to_array(Sum(MatrixElement(M, i, d) * MatrixElement(N, d, j), list(d, 0, k - 1)), list(i, j))
#> ArrayContraction(ArrayTensorProduct(M, N), (1, 2))
convert_array_to_indexed(ArraySymbol("A", list(3, 3)), list(i, j))
#> A[i, j]
```
