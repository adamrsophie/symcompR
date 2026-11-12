# Galois group of a polynomial

Computes the Galois group of an irreducible polynomial over the integers
or rationals of degree at most 6, as SymPy's `galois_group()`, with the
algorithms of Cohen's *A Course in Computational Algebraic Number
Theory* (Section 6.3) and SymPy's lookup tables of resolvent
coefficients. When a resolvent has a repeated root the polynomial is
replaced by a Tschirnhausen transformation, tried in SymPy's fixed order
unless `randomize = TRUE`.

## Usage

``` r
galois_group(
  f,
  gens = NULL,
  by_name = FALSE,
  max_tries = 30,
  randomize = FALSE
)
```

## Arguments

- f:

  An irreducible polynomial expression or `Poly`.

- gens:

  Optional generator (or list of generators) for converting `f`.

- by_name:

  Return the group by name, an element of
  [S1TransitiveSubgroups](https://adamrsophie.github.io/symcompR/reference/TransitiveSubgroups.md)
  to
  [S6TransitiveSubgroups](https://adamrsophie.github.io/symcompR/reference/TransitiveSubgroups.md),
  instead of a `PermutationGroup`.

- max_tries:

  Attempts allowed in the steps that generate Tschirnhausen
  transformations.

- randomize:

  Use random coefficients (from the Python compatible generator) for the
  Tschirnhausen transformations.

## Value

A `galois_group_result`: a list `(G, alt)` printed like SymPy's tuple,
where `G` is the group or its name and `alt` says whether the group is
contained in the alternating group of the degree of `f`.

## Examples

``` r
x <- Symbol("x")
galois_group(x^4 + 1)
#> (PermutationGroup([
#>     (0 1)(2 3),
#>     (0 2)(1 3)]), True)
galois_group(x^4 - 2, by_name = TRUE)
#> (<S4TransitiveSubgroups.D4: 'D4'>, False)
group_order(galois_group(x^3 + 2)[[1]])
#> [1] 6
```
