# Matrix unit

`MatrixUnit(rows, cols, i, j)` is the `rows` by `cols` matrix with a one
at position `(i, j)` (zero based) and zeros elsewhere: the derivative of
the element `A[i, j]` with respect to `A`.

## Usage

``` r
MatrixUnit(rows, cols, i, j)
```

## Arguments

- rows, cols:

  The shape.

- i, j:

  The position of the unit entry.

## Value

A `sym` matrix expression.

## Examples

``` r
k <- Symbol("k")
MatrixUnit(k, k, 0, 1)
#> MatrixUnit(k, k, 0, 1)
as_explicit(MatrixUnit(2, 2, 0, 1))
#> Matrix(2, 2)
#> [0, 1]
#> [0, 0]
```
