# Jordan normal form

`jordan_form(M)` returns `list(P, J)` with `M = P J P^-1` and `J` in
Jordan normal form, as SymPy's `Matrix.jordan_form()`; with
`calc_transform = FALSE` only `J`. Rational matrices are worked in the
algebraic field of each factor of the characteristic polynomial, as
SymPy does, so the columns of `P` agree with SymPy's.

## Usage

``` r
jordan_form(M, calc_transform = TRUE)
```

## Arguments

- M:

  A square `symmatrix`.

- calc_transform:

  Also compute the transform `P`.

## Value

`list(P = , J = )`, or the matrix `J`.

## Examples

``` r
jordan_form(Matrix(list(list(1, 1, 0), list(0, 1, 0), list(0, 0, 2))))
#> $P
#> Matrix(3, 3)
#> [0, 1, 0]
#> [0, 0, 1]
#> [1, 0, 0]
#> 
#> $J
#> Matrix(3, 3)
#> [2, 0, 0]
#> [0, 1, 1]
#> [0, 0, 1]
#> 
```
