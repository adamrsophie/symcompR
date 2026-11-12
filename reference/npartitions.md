# Number of partitions

`npartitions(n)` is the number of ways of writing `n` as a sum of
positive integers, as
[`partition()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md)
gives symbolically.

## Usage

``` r
npartitions(n)
```

## Arguments

- n:

  A nonnegative integer.

## Value

An integer (a `gmp` big integer from `n = 300` on).

## Examples

``` r
npartitions(25)
#> [1] 1958
npartitions(1000)
#> Big Integer ('bigz') :
#> [1] 24061467864032622473692149727991
```
