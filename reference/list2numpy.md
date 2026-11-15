# Lists and matrices as R arrays

SymPy's `list2numpy()` and `matrix2numpy()`: a list as an R vector and a
matrix as an R matrix. With the default `dtype = "list"` (SymPy's
`object`) the entries stay expressions, in a list vector or list matrix;
`"complex"` or `"double"` converts them to numbers.

## Usage

``` r
list2numpy(l, dtype = "list")

matrix2numpy(m, dtype = "list")
```

## Arguments

- l:

  A list.

- dtype:

  `"list"`, `"complex"` or `"double"`.

- m:

  A matrix.

## Value

An R vector or matrix.

## Examples

``` r
x <- Symbol("x")
matrix2numpy(Matrix(list(list(1, x), list(x, 2))))
#>      [,1] [,2]
#> [1,] 1    x   
#> [2,] x    2   
list2numpy(list(1, 2), dtype = "double")
#> [1] 1 2
```
