# Finitely presented groups

`FpGroup(F, relators)` is the group with the generators of the free
group `F` and the given relators; `fp_group()` also returns its
generators. The `fp_*` functions are the methods of SymPy's `FpGroup`:
`fp_order()` (by coset enumeration), `fp_index()`,
`fp_coset_enumeration()`, `fp_coset_table()`, `fp_subgroup()` (a
presentation of a subgroup by the Reidemeister Schreier algorithm,
optionally with the homomorphism into the group), `fp_quotient()`,
`fp_reduce()`, `fp_equals()` and `fp_make_confluent()` (rewriting
systems), `fp_most_frequent_generator()` and `fp_random()`.

## Usage

``` r
FpGroup(fr_grp, relators)

fp_group(fr_grp, relators = list())

xfp_group(fr_grp, relators = list())

fp_make_confluent(G)

fp_reduce(G, word)

fp_equals(G, word1, word2)

fp_subgroup(G, gens, C = NULL, homomorphism = FALSE)

fp_quotient(G, N)

fp_coset_enumeration(
  G,
  H,
  strategy = "relator_based",
  max_cosets = NULL,
  draft = NULL,
  incomplete = FALSE
)

fp_coset_table(
  G,
  H,
  strategy = "relator_based",
  max_cosets = NULL,
  draft = NULL,
  incomplete = FALSE
)

fp_order(G, strategy = "relator_based")

fp_most_frequent_generator(G)

fp_random(G)

fp_index(G, H, strategy = "relator_based")
```

## Arguments

- fr_grp:

  A free group.

- relators:

  A list of words.

- G:

  An `FpGroup`.

- word, word1, word2:

  Words.

- C:

  A coset table to use.

- homomorphism:

  Also return the homomorphism from the subgroup.

- N:

  A list of words, or a normal `FpSubgroup`.

- H, gens:

  A list of words generating a subgroup.

- strategy:

  `"relator_based"` (HLT) or `"coset_table_based"` (Felsch).

- max_cosets, draft, incomplete:

  Passed to the coset enumeration.

## Value

An `FpGroup` (an environment), a list, a word, an integer, `oo` or a
logical.

## Examples

``` r
g <- free_group("x, y")
x <- g[[2]]; y <- g[[3]]
f <- FpGroup(g[[1]], list(x^4, y^2, x * y * x^-1 * y))
fp_order(f)
#> 8
fp_str(fp_subgroup(f, list(x * y)))
#> [1] "<fp group on the generators (y_0,)>"
```
