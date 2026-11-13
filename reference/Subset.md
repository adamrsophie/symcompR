# Subsets

`Subset(subset, superset)` is SymPy's `Subset`: a subset of an ordered
superset, with its ranks in binary, lexicographic and Gray code order.
The fields `$subset`, `$superset`, `$size`, `$superset_size`,
`$cardinality`, `$rank_binary`, `$rank_lexicographic` and `$rank_gray`
mirror SymPy's. The `subset_next_*()`, `subset_prev_*()` and
`subset_iterate_*()` functions step through the subsets in each order;
`subset_unrank_binary()` and `subset_unrank_gray()` build the subset of
a rank, `subset_from_bitlist()` and `bitlist_from_subset()` convert to
and from bit strings, and `subset_indices()` gives the positions of a
subset in a superset. `ksubsets(superset, k)` lists the subsets of size
`k`.

## Usage

``` r
Subset(subset, superset)

subset_iterate_binary(a, k)

subset_next_binary(a)

subset_prev_binary(a)

subset_next_lexicographic(a)

subset_prev_lexicographic(a)

subset_iterate_graycode(a, k)

subset_next_gray(a)

subset_prev_gray(a)

subset_from_bitlist(super_set, bitlist)

bitlist_from_subset(subset, superset)

subset_unrank_binary(rank, superset)

subset_unrank_gray(rank, superset)

subset_indices(subset, superset)

ksubsets(superset, k)
```

## Arguments

- subset, superset:

  Vectors or lists.

- a:

  A `Subset`.

- k:

  Steps to move, or for `ksubsets()` the subset size.

- super_set:

  A vector or list.

- bitlist:

  A bit string.

- rank:

  A rank.

## Value

A `Subset`, a bit string, an integer, or a list.

## Examples

``` r
a <- Subset(c("c", "d"), c("a", "b", "c", "d"))
subset_next_binary(a)$subset
#> [1] "b"
a$rank_gray
#> [1] 2
subset_unrank_gray(10, c("a", "b", "c", "d"))$subset
#> [1] "a" "b" "c" "d"
ksubsets(1:3, 2)
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 1 3
#> 
#> [[3]]
#> [1] 2 3
#> 
```
