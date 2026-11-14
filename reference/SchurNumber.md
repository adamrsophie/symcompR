# Schur numbers and sum-free partitions

`SchurNumber(k)` is the largest `n` for which `1, ..., n` splits into
`k` sum-free sets (known for `k <= 5`, otherwise unevaluated), and
`schur_lower_bound(k)` its lower bound. `schur_partition(n)` splits
`1, ..., n` into the fewest sum-free sets that SymPy's construction
finds.

## Usage

``` r
SchurNumber(k)

schur_lower_bound(s)

schur_partition(n)
```

## Arguments

- k:

  A nonnegative integer or symbol.

- s:

  A `SchurNumber()` expression.

- n:

  A positive integer.

## Value

An expression, or a list of integer vectors.

## Examples

``` r
SchurNumber(4)
#> 44
schur_lower_bound(SchurNumber(8))
#> 5039
schur_partition(9)
#> [[1]]
#> [1] 3 2
#> 
#> [[2]]
#> [1] 6 9 5 8
#> 
#> [[3]]
#> [1] 1 4 7
#> 
```
