# Multiset partition enumeration

Knuth's Algorithm M as SymPy's `MultisetPartitionTraverser`:
`enum_all(mult)` lists every partition of the multiset with component
multiplicities `mult`, `enum_small(mult, ub)` those with fewer than `ub`
parts, `enum_large(mult, lb)` those with more than `lb` parts and
`enum_range(mult, lb, ub)` those in between; `count_partitions(mult)`
counts them. Each partition is a list of parts, each part a vector of
component indices (from 0), or of `components` when given.
`multiset_partitions_taocp(mult)` is `enum_all()`.

## Usage

``` r
enum_all(mult, components = NULL)

enum_small(mult, ub, components = NULL)

enum_large(mult, lb, components = NULL)

enum_range(mult, lb, ub, components = NULL)

multiset_partitions_taocp(mult, components = NULL)

count_partitions(mult)
```

## Arguments

- mult:

  Multiplicities of the components.

- components:

  Labels of the components.

- lb, ub:

  Bounds on the number of parts.

## Value

A list of partitions, or a count.

## Examples

``` r
enum_all(c(2, 2))
#> [[1]]
#> [[1]][[1]]
#> [1] 0 0 1 1
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 0 0 1
#> 
#> [[2]][[2]]
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 0 0
#> 
#> [[3]][[2]]
#> [1] 1 1
#> 
#> 
#> [[4]]
#> [[4]][[1]]
#> [1] 0 0
#> 
#> [[4]][[2]]
#> [1] 1
#> 
#> [[4]][[3]]
#> [1] 1
#> 
#> 
#> [[5]]
#> [[5]][[1]]
#> [1] 0 1 1
#> 
#> [[5]][[2]]
#> [1] 0
#> 
#> 
#> [[6]]
#> [[6]][[1]]
#> [1] 0 1
#> 
#> [[6]][[2]]
#> [1] 0 1
#> 
#> 
#> [[7]]
#> [[7]][[1]]
#> [1] 0 1
#> 
#> [[7]][[2]]
#> [1] 0
#> 
#> [[7]][[3]]
#> [1] 1
#> 
#> 
#> [[8]]
#> [[8]][[1]]
#> [1] 0
#> 
#> [[8]][[2]]
#> [1] 0
#> 
#> [[8]][[3]]
#> [1] 1 1
#> 
#> 
#> [[9]]
#> [[9]][[1]]
#> [1] 0
#> 
#> [[9]][[2]]
#> [1] 0
#> 
#> [[9]][[3]]
#> [1] 1
#> 
#> [[9]][[4]]
#> [1] 1
#> 
#> 
count_partitions(c(2, 2, 1))
#> [1] 26
```
