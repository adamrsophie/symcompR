# Homomorphisms of permutation groups

`orbit_homomorphism(G, omega)` is the action of `G` on its invariant set
`omega` and `block_homomorphism(G, blocks)` the action on a block system
(as returned by
[`minimal_block()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)),
both as in SymPy. `hom_image()`, `hom_kernel()`, `hom_apply()` and
`hom_invert()` are the homomorphism's
[`image()`](https://rdrr.io/r/graphics/image.html),
[`kernel()`](https://rdrr.io/r/stats/kernel.html), call and
[`invert()`](https://adamrsophie.github.io/symcompR/reference/gcd.md).

## Usage

``` r
orbit_homomorphism(G, omega)

block_homomorphism(G, blocks)

hom_image(hom)

hom_kernel(hom)

hom_apply(hom, g)

hom_invert(hom, g)
```

## Arguments

- G:

  A permutation group.

- omega:

  A vector of points.

- blocks:

  A block system.

- hom:

  A homomorphism.

- g:

  A permutation.

## Value

A `GroupHomomorphism`, a group or a permutation.

## Examples

``` r
D <- DihedralGroup(6)
h <- block_homomorphism(D, minimal_block(D, c(0, 3)))
group_order(hom_image(h))
#> [1] 6
```
