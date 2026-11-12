# Canonicalize a tensor formed by component tensors

Port of SymPy's `canonicalize()`. The tensor is the permutation `g`; the
free indices come first, then the dummy indices in the order
`d0, -d0, d1, -d1, ...`. Each argument in `...` describes a type of
component tensor as `list(base_i, gens_i, n_i, sym_i)`: a BSGS of its
slot symmetry, the number of such tensors, and their commutation
symmetry (`NULL` none, `0` commuting, `1` anticommuting). The algorithm
first minimizes the free indices with
[`canonical_free()`](https://adamrsophie.github.io/symcompR/reference/canonical_free.md),
then the dummy indices with
[`double_coset_can_rep()`](https://adamrsophie.github.io/symcompR/reference/double_coset_can_rep.md).
`riemann_bsgs` is the BSGS of the Riemann tensor symmetry.

## Usage

``` r
canonicalize(g, dummies, msym, ...)
```

## Arguments

- g:

  The permutation representing the tensor.

- dummies:

  Integer vector of dummy indices, or a list of them, one per index
  type.

- msym:

  Metric symmetry (`NULL`, `0` or `1`), or a list of them, one per index
  type.

- ...:

  Tensor types, each `list(base_i, gens_i, n_i, sym_i)`.

## Value

The canonical array form (integer vector), or `0` when the tensor
vanishes.

## Examples

``` r
b2a <- get_symmetric_group_sgs(2, 1)
g <- Permutation(c(1, 3, 0, 5, 4, 2, 6, 7))
canonicalize(g, 0:5, 0, list(b2a[[1]], b2a[[2]], 1, 0), list(b2a[[1]], b2a[[2]], 2, 0))
#> [1] 0
canonicalize(g, 0:5, 0, list(b2a[[1]], b2a[[2]], 1, 0), list(b2a[[1]], b2a[[2]], 2, 1))
#> [1] 0 2 1 4 3 5 7 6
canonicalize(Permutation(c(0, 2, 3, 1, 4, 5)), 2:3, 0,
             list(riemann_bsgs[[1]], riemann_bsgs[[2]], 1, 0))
#> [1] 0 2 1 3 5 4
```
