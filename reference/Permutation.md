# Permutations

`Permutation()` builds a permutation of `0..n-1`, following SymPy's
constructor: an integer vector of length two or more is the array form
(`Permutation(c(2, 0, 1))`), a list gives disjoint or composed cycles
(`Permutation(list(c(0, 1), c(2, 3)))`), several integers form one cycle
(`Permutation(1, 2)`), a single integer `n` is the identity on `0..n`
(`Permutation(3)`), and no argument is the empty permutation. `size`
pads with fixed points and must be named. `*` composes (`p * q` applies
`p` first, then `q`), `^` with a number is a power, `p ^ q` with a
permutation is the conjugate `q^-1 p q`, `i ^ p` is the image of the
integer `i`, and `p + k` or `p - k` move `k` steps in lexicographic
rank.

## Usage

``` r
Permutation(..., size = NULL)

array_form(p)

perm_size(p)

perm_inverse(p)

apply_perm(p, i)

cyclic_form(p)

full_cyclic_form(p)

perm_list(p, size = NULL)

perm_support(p)

perm_rmul(...)

perm_prod(perms)

perm_mul_inv(p, q)

perm_commutes_with(p, q)

perm_call(p, ...)

perm_transpositions(p)

perm_from_sequence(x, key = NULL)

perm_atoms(p)

perm_next_lex(p)

perm_unrank_nonlex(n, r)

perm_rank_nonlex(p)

perm_next_nonlex(p)

perm_rank(p)

perm_cardinality(p)

perm_parity(p)

perm_is_even(p)

perm_is_odd(p)

perm_is_singleton(p)

perm_is_empty(p)

perm_is_identity(p)

perm_ascents(p)

perm_descents(p)

perm_max(p)

perm_min(p)

perm_inversions(p)

perm_commutator(p, q)

perm_signature(p)

perm_sign(p)

perm_order(p)

perm_length(p)

perm_cycle_structure(p)

perm_cycles(p)

perm_index(p)

perm_runs(p)

perm_inversion_vector(p)

perm_rank_trotterjohnson(p)

perm_unrank_trotterjohnson(size, r)

perm_next_trotterjohnson(p)

perm_precedence_matrix(p)

perm_precedence_distance(p, q)

perm_adjacency_matrix(p)

perm_adjacency_distance(p, q)

perm_positional_distance(p, q)

perm_josephus(m, n, s = 1)

perm_from_inversion_vector(inversion)

perm_random(n)

unrank_lex(n, r)

perm_resize(p, n)
```

## Arguments

- ...:

  For `Permutation()`, the array form, a list of cycles, integers
  forming one cycle, a `Permutation` or a `Cycle`. For `perm_call()`, an
  index, a vector to permute, or integers of a cycle to multiply by. For
  `perm_rmul()`, permutations.

- size:

  Total size; pads with fixed points.

- p, q:

  Permutations.

- i:

  An element or vector of elements.

- perms:

  A list of permutations.

- x:

  A character string (its characters) or a vector.

- key:

  An optional function mapping each element to its sort key.

- n:

  A size.

- r:

  A rank.

- m, s:

  Josephus parameters: every `m`th of `n` is removed until `s` remain.

- inversion:

  An inversion vector.

## Value

A `Permutation`, an integer, logical, vector or list.

## Details

The helpers mirror the methods of SymPy's `Permutation`: `array_form()`,
`cyclic_form()`, `full_cyclic_form()`, `perm_size()`, `perm_list()`,
`perm_support()`, `perm_inverse()`, `apply_perm()`, `perm_call()`,
`perm_rmul()`, `perm_prod()`, `perm_mul_inv()`, `perm_commutes_with()`,
`perm_commutator()`, `perm_transpositions()`, `perm_atoms()`,
`perm_order()`, `perm_parity()`, `perm_sign()`, `perm_signature()`,
`perm_is_even()`, `perm_is_odd()`, `perm_is_identity()`,
`perm_is_singleton()`, `perm_is_empty()`, `perm_cycle_structure()`,
`perm_cycles()`, `perm_length()`, `perm_max()`, `perm_min()`,
`perm_ascents()`, `perm_descents()`, `perm_runs()`, `perm_index()`,
`perm_inversions()`, `perm_inversion_vector()`, `perm_cardinality()`,
`perm_resize()`.

## Examples

``` r
p <- Permutation(c(1, 2, 0))
q <- Permutation(list(c(0, 1)), size = 3)
p * q
#> (1 2)
perm_order(p)
#> [1] 3
cyclic_form(p * q)
#> [[1]]
#> [1] 1 2
#> 
Permutation(1, 2, size = 4)
#> (3)(1 2)
```
