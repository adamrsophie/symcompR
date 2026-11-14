# Special matrices

`jordan_block(size, eigenvalue, band)` is a Jordan block, `companion(p)`
the companion matrix of a monic univariate `Poly`, `wilkinson(n)` the
pair of Wilkinson matrices `list(W_minus, W_plus)`, and
`get_diag_blocks(M)` the square blocks along the diagonal, as SymPy's
`Matrix` class methods.

## Usage

``` r
jordan_block(size = NULL, eigenvalue = NULL, band = "upper", eigenval = NULL)

companion(p)

wilkinson(n)

get_diag_blocks(M)
```

## Arguments

- size:

  The size of the block.

- eigenvalue:

  The diagonal entry.

- band:

  `"upper"` or `"lower"`: where the ones go.

- eigenval:

  An alias of `eigenvalue`.

- p:

  A monic univariate `Poly`.

- n:

  The Wilkinson order: the matrices are `2n + 1` square.

- M:

  A `symmatrix`.

## Value

A `symmatrix`, or a list of them.

## Examples

``` r
jordan_block(3, 2)
#> Matrix(3, 3)
#> [2, 1, 0]
#> [0, 2, 1]
#> [0, 0, 2]
x <- Symbol("x")
companion(Poly(x^2 - 3*x + 2, x))
#> Matrix(2, 2)
#> [0, -2]
#> [1, 3 ]
get_diag_blocks(diag_matrix(1, Matrix(list(list(1, 2), list(3, 4)))))
#> [[1]]
#> Matrix(1, 1)
#> [1]
#> 
#> [[2]]
#> Matrix(2, 2)
#> [1, 2]
#> [3, 4]
#> 
```
