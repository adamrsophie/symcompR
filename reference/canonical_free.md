# Canonicalize a tensor with respect to its free indices

Port of SymPy's `canonical_free()`: uses the slot symmetry with BSGS
`base`, `gens` to bring the free indices (the first `num_free` indices)
to the lowest slots in lexicographic order.

## Usage

``` r
canonical_free(base, gens, g, num_free)
```

## Arguments

- base, gens:

  BSGS of the slot permutation group.

- g:

  The permutation representing the tensor.

- num_free:

  Number of free indices.

## Value

An integer vector (array form).

## Examples

``` r
gens <- list(Permutation(c(1, 0, 2, 3, 5, 4)), Permutation(c(2, 3, 0, 1, 4, 5)),
             Permutation(c(0, 1, 3, 2, 5, 4)))
canonical_free(c(0, 2), gens, Permutation(c(2, 1, 0, 3, 4, 5)), 4)
#> [1] 0 3 1 2 5 4
```
