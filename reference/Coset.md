# Cosets and the symmetric group of a degree

`Coset(g, H, G, dir)` is SymPy's `Coset`: the right coset `H*g` for
`dir = "+"` and the left coset `g*H` for `dir = "-"`; `G` defaults to
the full symmetric group, `SymmetricPermutationGroup(n)`. `G * g` and
`g * G` build cosets too. `is_left_coset()`, `is_right_coset()` and
`coset_as_list()` are its methods.

## Usage

``` r
Coset(g, H, G = NULL, dir = "+")

is_left_coset(x)

is_right_coset(x)

coset_as_list(x)

SymmetricPermutationGroup(n)
```

## Arguments

- g:

  A permutation.

- H, G:

  Permutation groups.

- dir:

  `"+"` or `"-"`.

- x:

  A `Coset`.

- n:

  The degree.

## Value

A `Coset`, a `SymmetricPermutationGroup`, a logical or a list.

## Examples

``` r
G <- PermutationGroup(Permutation(1, 2), Permutation(0, 1))
cst <- G * Permutation(1, 2)
is_right_coset(cst)
#> [1] TRUE
length(coset_as_list(cst))
#> [1] 6
```
