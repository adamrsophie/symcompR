# Characteristic polynomial, eigenvalues and eigenvectors

`charpoly(M, lambda)` returns `sym_det(lambda*I - M)`. `eigenvals()`
returns the roots of the characteristic polynomial with multiplicities;
`eigenvects()` returns a list of `list(value, multiplicity, vectors)`.
`diagonalize()` returns `list(P, D)` with `M = P D P^-1`.
`is_diagonalizable()` tests whether that is possible.

## Usage

``` r
charpoly(M, lambda = Symbol("lambda"))

eigenvals(M)

eigenvects(M)

is_diagonalizable(M)

diagonalize(M)
```

## Arguments

- M:

  A square `symmatrix`.

- lambda:

  Symbol for the characteristic polynomial.

## Value

An expression or list.

## Examples

``` r
M <- Matrix(list(list(2, 1), list(1, 2)))
charpoly(M)
#> lambda^2 - 4*lambda + 3
eigenvals(M)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 1
```
