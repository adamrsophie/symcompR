# Left eigenvectors and singular values

`left_eigenvects()` returns the eigenvectors of the transpose as row
vectors, and `singular_values()` returns the square roots of the
eigenvalues of `M^H M` (or `M M^H` for a wide matrix) in decreasing
order.

## Usage

``` r
left_eigenvects(M, ...)

singular_values(M)
```

## Arguments

- M:

  A matrix.

- ...:

  Flags passed to
  [`eigenvects()`](https://adamrsophie.github.io/symcompR/reference/charpoly.md).

## Value

For `left_eigenvects()`, a list of `list(value, multiplicity, vectors)`;
for `singular_values()`, a list of expressions.

## Examples

``` r
left_eigenvects(Matrix(list(list(0, 1), list(1, 0))))
#> [[1]]
#> [[1]]$value
#> -1
#> 
#> [[1]]$multiplicity
#> [1] 1
#> 
#> [[1]]$vectors
#> [[1]]$vectors[[1]]
#> Matrix(1, 2)
#> [-1, 1]
#> 
#> 
#> 
#> [[2]]
#> [[2]]$value
#> 1
#> 
#> [[2]]$multiplicity
#> [1] 1
#> 
#> [[2]]$vectors
#> [[2]]$vectors[[1]]
#> Matrix(1, 2)
#> [1, 1]
#> 
#> 
#> 
singular_values(Matrix(list(list(1, 0), list(0, 2))))
#> [[1]]
#> 2
#> 
#> [[2]]
#> 1
#> 
```
