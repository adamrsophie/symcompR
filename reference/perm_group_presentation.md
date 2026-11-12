# Presentation of a permutation group

Returns an `FpGroup` isomorphic to the permutation group `G`, built by
SymPy's algorithm (Handbook of Computational Group Theory, 6.1). A
cyclic group gives the presentation on one generator; the trivial group
gives the free group of rank 0.

## Usage

``` r
perm_group_presentation(G, eliminate_gens = TRUE)
```

## Arguments

- G:

  A permutation group.

- eliminate_gens:

  Remove redundant generators.

## Value

An `FpGroup` or a free group.

## Examples

``` r
P <- perm_group_presentation(DihedralGroup(4))
fp_order(P)
#> 8
```
