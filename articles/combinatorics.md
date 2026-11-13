# Combinatorics

This vignette covers the combinatorial generators for permutations,
combinations, partitions, and related objects, together with the
iteration utilities of `sympy.utilities.iterables` and the multiset
partition enumerators. Permutations and groups are covered in
[`vignette("permutation-groups", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/permutation-groups.md)
and
[`vignette("fp-groups", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/fp-groups.md).

``` r

library(symcompR)
```

## Enumerating combinatorial objects

The generators return lists of vectors. A small helper prints each
element on one line.

``` r

flat <- function(l) vapply(l, paste, character(1), collapse = " ")
length(permutations(1:4))
#> [1] 24
flat(combinations(1:4, 2))
#> [1] "1 2" "1 3" "1 4" "2 3" "2 4" "3 4"
length(subsets(1:3, 2))
#> [1] 3
flat(multiset_permutations(c(1, 1, 2)))
#> [1] "1 1 2" "1 2 1" "2 1 1"
flat(generate_derangements(c(1, 2, 3)))
#> [1] "2 3 1" "3 1 2"
length(multiset_partitions(1:4))
#> [1] 15
flat(necklaces(3, 2))
#> [1] "0 0 0" "0 0 1" "0 1 1" "1 1 1"
flat(gray_code(3))
#> [1] "0 0 0" "0 0 1" "0 1 1" "0 1 0" "1 1 0" "1 1 1" "1 0 1" "1 0 0"
```

Consecutive elements of the Gray code differ in exactly one bit.

[`partitions()`](https://adamrsophie.github.io/symcompR/reference/iterables.md)
returns integer partitions as named vectors, mapping each part to its
multiplicity;
[`ordered_partitions()`](https://adamrsophie.github.io/symcompR/reference/iterables.md)
returns them as vectors.

``` r

partitions(4)[[2]]
#> 3 1 
#> 1 1
flat(ordered_partitions(5, 3))
#> [1] "1 1 3" "1 2 2"
```

Counting functions:

``` r

nT(5, 2)
#> 2
stirling(5, 2)
#> 15
stirling(5, 2, kind = 1)
#> 50
```

`stirling(5, 2)` is the number of ways to split 5 items into 2 nonempty
sets; with `kind = 1` it counts permutations of 5 items with 2 cycles.

## Prüfer sequences

A labeled tree on n nodes corresponds to a unique sequence of n - 2
labels.

``` r

prufer_from_tree(list(c(0, 1), c(0, 2), c(0, 3)), 4)
#> [1] 0 0
flat(tree_from_prufer(c(0, 0)))
#> [1] "0 1" "0 2" "0 3"
```

The sequence `0 0` encodes the star with center 0 and edges to 1, 2, and
3.

## Iterables

More generators from `sympy.utilities.iterables`. Partitions of a
sequence into consecutive pieces, sets into blocks, and integers into
powers of two:

``` r

pp <- function(l) vapply(l, function(b) paste0("(", vapply(b, paste, "", collapse = " "), ")", collapse = ""), "")
pp(kbins(0:3, 2))
#> [1] "(0)(1 2 3)" "(0 1)(2 3)" "(0 1 2)(3)"
pp(sequence_partitions(1:4, 2))
#> [1] "(1)(2 3 4)" "(1 2)(3 4)" "(1 2 3)(4)"
pp(sequence_partitions_empty(1:3, 2))
#> [1] "()(1 2 3)" "(1)(2 3)"  "(1 2)(3)"  "(1 2 3)()"
flat(binary_partitions(6))
#> [1] "4 2"         "4 1 1"       "2 2 2"       "2 2 1 1"     "2 1 1 1 1"  
#> [6] "1 1 1 1 1 1"
length(set_partitions(4))
#> [1] 15
flat(derangements(1:4))
#> [1] "2 1 4 3" "2 3 4 1" "2 4 1 3" "3 1 4 2" "3 4 1 2" "3 4 2 1" "4 1 2 3"
#> [8] "4 3 1 2" "4 3 2 1"
flat(combinations_with_replacement(1:3, 2))
#> [1] "1 1" "1 2" "1 3" "2 2" "2 3" "3 3"
flat(cartes(1:2, c("a", "b")))
#> [1] "1 a" "1 b" "2 a" "2 b"
flat(multiset_combinations(c(a = 2, b = 1), 2))
#> [1] "a a" "a b"
flat(multiset_derangements(c("a", "a", "b", "b")))
#> [1] "b b a a"
flat(signed_permutations(c(0, 1)))
#> [1] "0 1"  "0 -1" "1 0"  "-1 0"
flat(generate_involutions(4))
#>  [1] "0 1 2 3" "0 1 3 2" "0 2 1 3" "0 3 2 1" "1 0 2 3" "1 0 3 2" "2 1 0 3"
#>  [8] "2 3 0 1" "3 1 2 0" "3 2 1 0"
flat(generate_oriented_forest(4))
#> [1] "0 1 2 3" "0 1 2 2" "0 1 2 1" "0 1 2 0" "0 1 1 1" "0 1 1 0" "0 1 0 1"
#> [8] "0 1 0 0" "0 0 0 0"
flat(bracelets(5, 2))
#> [1] "0 0 0 0 0" "0 0 0 0 1" "0 0 0 1 1" "0 0 1 0 1" "0 0 1 1 1" "0 1 0 1 1"
#> [7] "0 1 1 1 1" "1 1 1 1 1"
```

Utilities for sequences:

``` r

multiset(c(1, 1, 2, 3, 1))
#> 1 2 3 
#> 3 1 1
uniq(c(1, 4, 1, 5, 4, 2, 1, 2))
#> [1] 1 4 5 2
common_prefix(c(1, 2, 3), c(1, 2, 5))
#> [1] 1 2
common_suffix(c(1, 2, 3), c(9, 2, 3))
#> [1] 2 3
rotate_left(c(0, 1, 2), 1)
#> [1] 1 2 0
rotate_right(c(0, 1, 2), 1)
#> [1] 2 0 1
least_rotation(c(3, 1, 5, 1, 2))
#> [1] 3
minlex(c(1, 2, 0))
#> [1] 0 1 2
unlist(roundrobin(c("A", "B", "C"), "D", c("E", "F")))
#> [1] "A" "D" "E" "B" "F" "C"
flat(reshape_seq(0:7, list(2)))
#> [1] "0 1" "2 3" "4 5" "6 7"
```

Graph algorithms on a vertex list and an edge list:

``` r

flat(connected_components(list(1:5, list(c(1, 2), c(4, 5)))))
#> [1] "1 2" "3"   "4 5"
flat(strongly_connected_components(list(1:5, list(c(1, 2), c(2, 1), c(3, 4),
                                                  c(4, 5), c(5, 3)))))
#> [1] "1 2"   "3 4 5"
```

A random derangement, with a seeded generator:

``` r

sympy_seed(1)
d <- random_derangement(1:5)
all(d != 1:5)
#> [1] TRUE
```

## Multiset partitions

The enumerators of `sympy.utilities.enumerative` follow Knuth’s
algorithm 7.1.2.5M. A multiset is given by its multiplicities: `c(2, 2)`
is \\0, 0, 1, 1\\.

``` r

count_partitions(c(2, 2, 1))
#> [1] 26
length(multiset_partitions_taocp(c(2, 2)))
#> [1] 9
length(enum_all(c(2, 2)))
#> [1] 9
pp(enum_small(c(2, 2), 2))
#> [1] "(0 0 1 1)"  "(0 0 1)(1)" "(0 0)(1 1)" "(0 1 1)(0)" "(0 1)(0 1)"
pp(enum_large(c(2, 2), 2))
#> [1] "(0 0)(1)(1)"  "(0 1)(0)(1)"  "(0)(0)(1 1)"  "(0)(0)(1)(1)"
pp(enum_range(c(2, 2), 1, 3))
#> [1] "(0 0 1)(1)"  "(0 0)(1 1)"  "(0 0)(1)(1)" "(0 1 1)(0)"  "(0 1)(0 1)" 
#> [6] "(0 1)(0)(1)" "(0)(0)(1 1)"
```

`enum_small(m, ub)` lists the partitions with at most `ub` parts,
`enum_large(m, lb)` those with more than `lb` parts, and
`enum_range(m, lb, ub)` those in between.
