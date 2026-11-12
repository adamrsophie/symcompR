# Permutation group methods of finitely presented groups

SymPy answers these questions about a finite `FpGroup` through an
isomorphic permutation group from coset enumeration:
`fp_to_perm_group()` returns that group and the isomorphism,
`fp_is_abelian()`, `fp_is_cyclic()`, `fp_is_solvable()` and
`fp_is_nilpotent()` test properties, `fp_center()`,
`fp_derived_subgroup()`, `fp_derived_series()`,
`fp_lower_central_series()`, `fp_centralizer()` and
`fp_normal_closure()` give generators of subgroups as words, and
`fp_elements()`, `fp_abelian_invariants()` and `fp_composition_series()`
are as in SymPy. `fp_is_abelian()`, `fp_center()` and `fp_elements()`
also accept free groups.

## Usage

``` r
fp_to_perm_group(G)

fp_derived_series(G)

fp_lower_central_series(G)

fp_center(G)

fp_derived_subgroup(G)

fp_centralizer(G, other)

fp_normal_closure(G, other)

fp_is_abelian(G)

fp_is_nilpotent(G)

fp_is_solvable(G)

fp_elements(G)

fp_is_cyclic(G)

fp_abelian_invariants(G)

fp_composition_series(G)
```

## Arguments

- G:

  An `FpGroup` (or a free group where noted).

- other:

  A list of words.

## Value

A list, a logical or a list of permutation groups.

## Examples

``` r
g <- free_group("x, y")
x <- g[[2]]; y <- g[[3]]
G <- FpGroup(g[[1]], list(x^2, y^8, x * y * x^-1 * y))
fp_str(fp_center(G))
#> [1] "[y^4]"
fp_is_abelian(G)
#> [1] FALSE
```
