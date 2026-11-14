# Block decompositions from the graph of a matrix

`connected_components_decomposition(M)` returns `list(P, B)` with a
`PermutationMatrix` `P` and a `BlockDiagMatrix` `B` such that
`M = P^T B P`; `strongly_connected_components_decomposition(M, lower)`
does the same with a block triangular `BlockMatrix`, as SymPy.

## Usage

``` r
connected_components_decomposition(M)

strongly_connected_components_decomposition(M, lower = TRUE)
```

## Arguments

- M:

  A square `symmatrix`.

- lower:

  Make the blocks lower triangular (`TRUE`) or upper.

## Value

A list of two matrix expressions.

## Examples

``` r
A <- Matrix(list(list(1, 0, 2), list(0, 3, 0), list(4, 0, 5)))
connected_components_decomposition(A)
#> [[1]]
#> PermutationMatrix((0, 2, 1))
#> 
#> [[2]]
#> Matrix([
#> [Matrix([
#> [1, 2],
#> [4, 5]]),             0],
#> [                        0, Matrix([[3]])]])
#> 
```
