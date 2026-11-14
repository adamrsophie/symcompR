# Condition number of a matrix

`condition_number(M)` is the ratio of the largest to the smallest
singular value, as SymPy's `Matrix.condition_number()`.

## Usage

``` r
condition_number(M)
```

## Arguments

- M:

  A `symmatrix`.

## Value

An expression.

## Examples

``` r
condition_number(Matrix(list(list(1, 0), list(0, 10))))
#> 10
```
