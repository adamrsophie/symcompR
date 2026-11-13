# Shape and rank of array and matrix expressions

`get_shape()` returns the shape of an array expression, a matrix
expression or an explicit array, as a list of dimensions (an empty list
for a scalar). `get_ndim()` returns the number of axes, and `get_rank()`
is its older name.

## Usage

``` r
get_shape(expr)

get_ndim(expr)

get_rank(expr)
```

## Arguments

- expr:

  An expression.

## Value

`get_shape()` returns a list of `sym` dimensions; `get_ndim()` and
`get_rank()` return an integer.

## Examples

``` r
A <- ArraySymbol("A", c(2, 3, 4))
get_shape(A)
#> [[1]]
#> 2
#> 
#> [[2]]
#> 3
#> 
#> [[3]]
#> 4
#> 
get_ndim(ArrayContraction(A, c(0, 1)))
#> Error: contracting indices of different dimensions
```
