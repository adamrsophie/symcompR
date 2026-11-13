# Matrix exponential of a constant matrix times t

`matrix_exp_jordan_form(A, t)` returns matrices `P` and `expJ` with
`exp(A*t) = P * expJ * P^-1`, as SymPy's `matrix_exp_jordan_form()`:
`expJ` is block diagonal, with real blocks for complex conjugate pairs
of eigenvalues of a real matrix.

## Usage

``` r
matrix_exp_jordan_form(A, t)
```

## Arguments

- A:

  A square matrix independent of `t`.

- t:

  The symbol.

## Value

A tuple `(P, expJ)`.

## Examples

``` r
t <- Symbol("t")
matrix_exp_jordan_form(Matrix(list(list(1, 1), list(0, 1))), t)
#> (Matrix([[1, 0], [0, 1]]), Matrix([[exp(t), t*exp(t)], [0, exp(t)]])) 
```
