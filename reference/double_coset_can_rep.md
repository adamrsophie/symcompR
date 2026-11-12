# Butler-Portugal canonical representative of a double coset

Port of SymPy's `double_coset_can_rep()`, the Butler-Portugal algorithm
for tensor canonicalization with dummy indices. The tensor is the
permutation `g` (size `n + 2`, the last two points carry the sign); the
slot symmetry group `S` has minimal BSGS `b_S`, `sgens`; the dummy
symmetry group `D` is built from `dummies` and `sym`. The result is the
minimum of `D*g*S` under lexicographic order, or `0` when the tensor
vanishes.

## Usage

``` r
double_coset_can_rep(dummies, sym, b_S, sgens, S_transversals, g)
```

## Arguments

- dummies:

  List of integer vectors of dummy indices, one per index type, ordered
  `d0, -d0, d1, -d1, ...`.

- sym:

  List of metric symmetries, one per index type: `NULL`, `0` or `1`.

- b_S:

  Base of a minimal slot symmetry BSGS.

- sgens:

  Strong generators of the slot symmetry (permutations).

- S_transversals:

  Transversals of the slot BSGS, from
  [`get_transversals()`](https://adamrsophie.github.io/symcompR/reference/tensor-bsgs.md).

- g:

  The permutation representing the tensor.

## Value

An integer vector (array form), or `0`.

## Examples

``` r
gens <- list(Permutation(c(2, 1, 0, 3, 4, 5, 7, 6)), Permutation(c(4, 1, 2, 3, 0, 5, 7, 6)))
base <- c(0, 2)
tr <- get_transversals(base, gens)
double_coset_can_rep(list(0:5), list(0), base, gens, tr, Permutation(c(4, 2, 0, 1, 3, 5, 6, 7)))
#> [1] 0 1 2 3 4 5 7 6
double_coset_can_rep(list(0:5), list(0), base, gens, tr, Permutation(c(4, 1, 3, 0, 5, 2, 6, 7)))
#> [1] 0
```
