# Permutations, combinations and derangements

`nP(n, k)` counts the permutations of `k` items taken from `n` (a number
of distinct items, or a multiset given as a vector, list, or named list
of counts), `nC(n, k)` the combinations, with `replacement` allowing
repeats, and `nD(i)` the derangements of the multiset `i` (or
`nD(n = 5)` of `n` distinct items, or `nD(m = list(...))` of a multiset
given by multiplicities); all as SymPy's functions of the same names.
With `k` omitted, `nP()` and `nC()` sum over all sizes.

## Usage

``` r
nP(n, k = NULL, replacement = FALSE)

nC(n, k = NULL, replacement = FALSE)

nD(i = NULL, n = NULL, m = NULL)
```

## Arguments

- n:

  A nonnegative integer or a multiset.

- k:

  The number of items taken, or `NULL` for all sizes.

- replacement:

  Allow items to be taken again.

- i:

  A multiset (vector, list, or named list of counts).

- m:

  For `nD()`: a multiset given as counts of multiplicities.

## Value

An `Integer`.

## Examples

``` r
nP(3, 2)
#> 6
nP(c("a", "a", "b"), 2)
#> 3
nC(c("a", "a", "b"), 2)
#> 2
nC(5, 2, replacement = TRUE)
#> 15
nD(c(1, 2, 3))
#> 2
nD(n = 4)
#> 9
```
