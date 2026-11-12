# Seed SymPy's random generator

`sympy_seed(a)` is `sympy.core.random.seed(a)`: it seeds the Mersenne
Twister generator that the randomized group algorithms
([`random_pr()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md),
[`group_random()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md),
[`normal_closure()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md),
[`sylow_subgroup()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md))
and
[`perm_random()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
draw from, so that they reproduce SymPy's values.

## Usage

``` r
sympy_seed(a)
```

## Arguments

- a:

  An integer seed.

## Value

`NULL`, invisibly.

## Examples

``` r
sympy_seed(3)
perm_random(5)
#> (1 2 3 4)
```
