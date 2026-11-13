# Restricted growth strings

`RGS_generalized(m)` returns the matrix of generalized unrestricted
growth string counts, `RGS_enum(m)` the number of restricted growth
strings of length `m` (a Bell number), `RGS_unrank(rank, m)` the string
of a given rank and `RGS_rank(rgs)` the rank of a string, after
`sympy.combinatorics.partitions`.

## Usage

``` r
RGS_generalized(m)

RGS_enum(m)

RGS_unrank(rank, m)

RGS_rank(rgs)
```

## Arguments

- m:

  The superset size.

- rank:

  A rank.

- rgs:

  A restricted growth string.

## Value

A `Matrix`, an `Integer` or an integer vector.

## Examples

``` r
RGS_generalized(3)
#> Matrix(4, 4)
#> [1, 1, 1, 1]
#> [1, 2, 3, 0]
#> [2, 5, 0, 0]
#> [5, 0, 0, 0]
RGS_enum(4)
#> 15
RGS_unrank(7, 5)
#> [1] 0 0 1 0 2
RGS_rank(c(0, 1, 2, 1, 2, 3))
#> 154
```
