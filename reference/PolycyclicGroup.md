# Polycyclic groups and collectors

`perm_polycyclic_group(G)` returns the polycyclic group of a solvable
permutation group `G` (SymPy's `G.polycyclic_group()`), with its
polycyclic generating sequence `$pcgs`, series `$pc_series`, relative
orders `$relative_order` and `Collector` `$collector`. The collector's
`$pc_presentation` is an ordered dictionary of power and conjugate
relations; `pc_presentation()` lists it as pairs. `collected_word()`,
`exponent_vector()`, `pc_depth()`, `leading_exponent()`,
`induced_pcgs()`, `constructive_membership_test()`,
`minimal_uncollected_subword()` and `pc_relations()` are the methods of
SymPy's `Collector`.

## Usage

``` r
PolycyclicGroup(pc_sequence, pc_series, relative_order, collector = NULL)

perm_polycyclic_group(G)

pc_is_prime_order(PC)

pc_length(PC)

Collector(
  pcgs,
  pc_series,
  relative_order,
  free_group_ = NULL,
  pc_presentation_ = NULL
)

pc_presentation(C)

minimal_uncollected_subword(C, word)

pc_relations(C)

collected_word(C, word)

exponent_vector(C, element)

pc_depth(C, element)

leading_exponent(C, element)

induced_pcgs(C, gens)

constructive_membership_test(C, ipcgs, g)
```

## Arguments

- pc_sequence, pcgs:

  A list of permutations.

- pc_series:

  A list of permutation groups.

- relative_order:

  An integer vector.

- collector:

  A collector (made when `NULL`).

- G:

  A solvable permutation group.

- PC:

  A polycyclic group.

- free_group\_, pc_presentation\_:

  Optional free group and presentation.

- C:

  A collector.

- word:

  A word of the collector's free group.

- element:

  A permutation.

- gens, ipcgs:

  Lists of permutations.

- g:

  A permutation.

## Value

A `PolycyclicGroup`, a `Collector`, a word, a list or an integer.

## Examples

``` r
PC <- perm_polycyclic_group(SymmetricGroup(4))
PC$relative_order
#> [1] 2 3 2 2
exponent_vector(PC$collector, SymmetricGroup(4)$gens[[1]])
#> [1] 1 0 0 0
```
