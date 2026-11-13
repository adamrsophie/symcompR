# Spin quantum numbers

`m_values(j)` gives the dimension `2*j + 1` and the values
`j, j - 1, ..., -j` of the projection quantum number, as SymPy's
`spin.m_values`.

## Usage

``` r
m_values(j)
```

## Arguments

- j:

  An integer or half-integer.

## Value

A list of the size and the list of m values.

## Examples

``` r
m_values(S(3) / 2)
#> [[1]]
#> 4
#> 
#> [[2]]
#> [[2]][[1]]
#> 3/2
#> 
#> [[2]][[2]]
#> 1/2
#> 
#> [[2]][[3]]
#> -1/2
#> 
#> [[2]][[4]]
#> -3/2
#> 
#> 
```
