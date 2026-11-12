# Group homomorphisms

`homomorphism(domain, codomain, gens, images)` is the homomorphism that
maps the generators `gens` of `domain` to `images` (the other generators
go to the identity); domain and codomain are permutation groups,
finitely presented groups or free groups. With `check = TRUE` the images
are checked against the relators of the domain (the presentation of a
permutation group). Between permutation groups the result is the
package's `GroupHomomorphism`; otherwise it is an `FpGroupHomomorphism`.
`hom_is_injective()`, `hom_is_surjective()`, `hom_is_isomorphism()` and
`hom_is_trivial()` test properties of either kind, `hom_images()` gives
the images of the generators (SymPy's `images` dictionary),
`hom_compose(T, S)` is `T` after `S`, `hom_factor()` gives a surjection
and an injection whose composite is `T`, `hom_restrict_to()` restricts
to a subgroup and `hom_invert_subgroup()` gives a preimage subgroup.
`group_isomorphism()` searches an isomorphism and `is_isomorphic()`
tests for one. Applying, inverting, kernels and images are
[`hom_apply()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md),
[`hom_invert()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md),
[`hom_kernel()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md)
and
[`hom_image()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md).

## Usage

``` r
GroupHomomorphism(domain, codomain, images)

hom_images(T)

hom_is_injective(T)

hom_is_surjective(T)

hom_is_isomorphism(T)

hom_is_trivial(T)

hom_factor(T)

hom_compose(T, other)

hom_restrict_to(T, H)

hom_invert_subgroup(T, H)

homomorphism(domain, codomain, gens, images = list(), check = TRUE)

group_isomorphism(G, H, isomorphism = TRUE)

is_isomorphic(G, H)
```

## Arguments

- domain, codomain, G, H:

  Groups.

- images:

  A list of elements of the codomain; for `GroupHomomorphism()` an
  ordered dictionary as built by `homomorphism()`.

- T, other:

  Homomorphisms.

- gens:

  A list of generators of the domain.

- check:

  Check that the images define a homomorphism.

- isomorphism:

  Also return the isomorphism.

## Value

A homomorphism (an environment), a group or a logical.

## Examples

``` r
g <- free_group("a, b")
a <- g[[2]]; b <- g[[3]]
G <- FpGroup(g[[1]], list(a^3, b^3, (a * b)^2))
A <- AlternatingGroup(4)
T <- homomorphism(G, A, list(a, b), list(Permutation(0, 1, 2, size = 4), Permutation(1, 2, 3)))
hom_is_isomorphism(T)
#> [1] TRUE
is_isomorphic(G, A)
#> [1] TRUE
```
