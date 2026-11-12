# BSGS of tensor symmetries

Ports of SymPy's BSGS constructors for tensor slot symmetries.
`get_symmetric_group_sgs(n, antisym)` is the minimal BSGS of a
(anti)symmetric tensor of rank `n`; `bsgs_direct_product()` is the
direct product of two BSGS; `get_transversals()` gives the basic
transversals of a BSGS; `get_minimal_bsgs()` recomputes a BSGS with
minimal base (or `NULL`); `tensor_gens()` is the BSGS of `n` tensors of
the same type, with fixed slots `list_free_indices` and commutation
symmetry `sym`; `gens_products()` combines tensors of different types,
each given as `list(base_i, gens_i, free_i, sym_i)`. `riemann_bsgs` is
the BSGS of the Riemann tensor.

## Usage

``` r
bsgs_direct_product(base1, gens1, base2, gens2, signed = TRUE)

get_symmetric_group_sgs(n, antisym = FALSE)

riemann_bsgs

get_transversals(base, gens)

get_minimal_bsgs(base, gens)

tensor_gens(base, gens, list_free_indices, sym = 0)

gens_products(...)
```

## Format

`riemann_bsgs`: `list(base, gens)` for the Riemann tensor.

## Arguments

- signed:

  Whether the permutations are signed.

- n:

  Rank of the tensor.

- antisym:

  `FALSE` for a symmetric, `TRUE` for an antisymmetric tensor.

- base, base1, base2:

  Bases (0 based integer vectors).

- gens, gens1, gens2:

  Strong generators (lists of permutations).

- list_free_indices:

  List of the slots occupied by fixed indices, one integer vector per
  tensor.

- sym:

  Commutation symmetry: `NULL`, `0` or `1`.

- ...:

  Tensor types, each `list(base_i, gens_i, free_i, sym_i)`.

## Value

`get_symmetric_group_sgs()`, `bsgs_direct_product()` and
`get_minimal_bsgs()` return `list(base, gens)`; `tensor_gens()` and
`gens_products()` return `list(size, base, gens)`; `get_transversals()`
returns a list of named lists (point to array form).

## Examples

``` r
get_symmetric_group_sgs(3)
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [[2]][[1]]
#> (4)(0 1)
#> 
#> [[2]][[2]]
#> (4)(1 2)
#> 
#> 
b1 <- get_symmetric_group_sgs(1); b2 <- get_symmetric_group_sgs(2)
bsgs_direct_product(b1[[1]], b1[[2]], b2[[1]], b2[[2]])
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [[2]][[1]]
#> (4)(1 2)
#> 
#> 
get_transversals(riemann_bsgs[[1]], riemann_bsgs[[2]])
#> [[1]]
#> [[1]]$`0`
#> [1] 0 1 2 3 4 5
#> 
#> [[1]]$`1`
#> [1] 1 0 2 3 5 4
#> 
#> [[1]]$`2`
#> [1] 2 3 0 1 4 5
#> 
#> [[1]]$`3`
#> [1] 3 2 0 1 5 4
#> 
#> 
#> [[2]]
#> [[2]]$`2`
#> [1] 0 1 2 3 4 5
#> 
#> [[2]]$`3`
#> [1] 0 1 3 2 5 4
#> 
#> 
b3 <- get_symmetric_group_sgs(3)
tensor_gens(b3[[1]], b3[[2]], list(1, 0))
#> [[1]]
#> [1] 8
#> 
#> [[2]]
#> [1] 0 4
#> 
#> [[3]]
#> [[3]][[1]]
#> (7)(0 2)
#> 
#> [[3]][[2]]
#> (7)(4 5)
#> 
#> 
gens_products(list(b2[[1]], b2[[2]], list(integer(), integer()), 0))
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> [1] 0 2
#> 
#> [[3]]
#> [[3]][[1]]
#> (5)(0 1)
#> 
#> [[3]][[2]]
#> (5)(2 3)
#> 
#> [[3]][[3]]
#> (5)(0 2)(1 3)
#> 
#> 
get_minimal_bsgs(c(2, 0), list(Permutation(c(1, 0, 2, 3, 5, 4)), Permutation(c(2, 3, 0, 1, 4, 5))))
#> [[1]]
#> [1] 0 2
#> 
#> [[2]]
#> [[2]][[1]]
#> (0 1)(4 5)
#> 
#> [[2]][[2]]
#> (5)(0 2)(1 3)
#> 
#> [[2]][[3]]
#> (2 3)(4 5)
#> 
#> 
```
