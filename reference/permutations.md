# Combinatorial generators

`permutations(x, k)`, `combinations(x, k)`,
`combinations_with_replacement(x, k)`, `subsets(x)`,
`set_partitions(n)`, `derangements(x)`, `cartes(...)` (Cartesian
product) and `gray_code(n)` return lists. `nT(n, k)` counts partitions
of `n` into `k` parts and `stirling(n, k, kind)` gives Stirling numbers.

## Usage

``` r
permutations(x, k = length(x))

combinations(x, k)

combinations_with_replacement(x, k)

subsets(x, k = NULL)

set_partitions(n)

derangements(x)

cartes(...)

gray_code(n)

nT(n, k)

stirling(n, k, kind = 2, signed = FALSE)
```

## Arguments

- x:

  A vector or list.

- k:

  Size.

- n:

  An integer.

- ...:

  Vectors for `cartes()`.

- kind:

  1 or 2 for Stirling numbers of the first or second kind.

- signed:

  For first kind, return signed numbers.

## Value

A list, or a `sym` integer.

## Examples

``` r
length(permutations(1:4))
#> [1] 24
stirling(5, 2)
#> 15
```
