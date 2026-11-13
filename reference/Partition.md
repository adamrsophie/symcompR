# Set partitions

`Partition(...)` builds the partition of a set into the given disjoint
blocks, as SymPy's `Partition`. Each block is a vector, a list or a
`FiniteSet`. The fields `$members`, `$size`, `$partition` (the blocks as
sorted lists), `$RGS` (the restricted growth string) and `$rank` mirror
SymPy's attributes. `p + k` and `p - k` move `k` steps in rank order,
`<` and `<=` compare by `partition_sort_key()`, `>` and `>=` test for
(proper) supersets of blocks as SymPy's sets do, and
`partition_from_rgs(rgs, elements)` rebuilds a partition from its
restricted growth string.

## Usage

``` r
Partition(...)

partition_sort_key(p, order = NULL)

partition_from_rgs(rgs, elements)
```

## Arguments

- ...:

  Blocks.

- p:

  A `Partition`.

- order:

  Ordering for the members (`NULL` keeps their order).

- rgs:

  A restricted growth string (0 based block indices).

- elements:

  The elements, as a vector or list.

## Value

A `Partition`; `partition_sort_key()` returns a sort key.

## Examples

``` r
p <- Partition(c(1, 2, 3), 4)
p
#> Partition({4}, {1, 2, 3})
p$RGS
#> [1] 0 0 0 1
p$rank
#> 1
(p + 1)$partition
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> [[1]][[2]]
#> 2
#> 
#> [[1]][[3]]
#> 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 3
#> 
#> 
partition_from_rgs(c(0, 1, 2, 0, 1), c(1, 3, 2, 4, 5))
#> Partition({2}, {1, 4}, {3, 5})
```
