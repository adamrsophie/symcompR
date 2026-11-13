# Determinant, inverse, trace and friends

`sym_det()` computes the determinant (Bareiss elimination by default, or
`"berkowitz"` or `"laplace"`). `inverse()` (also `solve(M)`) inverts by
Gauss-Jordan elimination; `adjugate()`, `cofactor()`, `minor()`,
`trace_matrix()`, `rank_matrix()`, `rref()`, `nullspace()`,
`columnspace()` and `rowspace()` follow SymPy. Given a characteristic
curve such as
[`TendonForceLengthDeGroote2016()`](https://adamrsophie.github.io/symcompR/reference/characteristic-curves.md),
`inverse()` returns the constructor of its inverse curve, as SymPy's
`Function.inverse()`.

## Usage

``` r
inverse(M, method = "GE")

# S3 method for class 'symmatrix'
solve(a, b, ...)

rref(M)

rank_matrix(M)

nullspace(M)

columnspace(M)

rowspace(M)

minor(M, i, j)

cofactor(M, i, j)

cofactor_matrix(M)

adjugate(M)

trace_matrix(M)

sym_det(x, ...)
```

## Arguments

- method:

  Algorithm name.

- a:

  A `symmatrix`.

- b:

  Optional right hand side.

- ...:

  Further arguments, such as `method`.

- i, j:

  Row and column.

- x, M:

  A `symmatrix`.

## Value

An expression, matrix or list.

## Examples

``` r
a <- Symbol("a")
M <- Matrix(list(list(a, 1), list(1, a)))
sym_det(M)
#> a^2 - 1
inverse(M)
#> Matrix(2, 2)
#> [a/(a^2 - 1) , -1/(a^2 - 1)]
#> [-1/(a^2 - 1), a/(a^2 - 1) ]
rref(Matrix(list(list(1, 2, 3), list(4, 5, 6))))
#> $matrix
#> Matrix(2, 3)
#> [1, 0, -1]
#> [0, 1, 2 ]
#> 
#> $pivots
#> [1] 1 2
#> 
```
