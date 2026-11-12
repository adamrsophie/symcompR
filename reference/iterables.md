# Iterables

Ports of the iterators of `sympy.utilities.iterables`, each returning a
list of everything the SymPy generator yields, in SymPy's order.
Sequences are R vectors or lists.

Small helpers from `sympy.utilities.iterables`: `flatten()`,
`unflatten()`, `group()` (runs of equal elements), `take()`,
`variations()`, `sift()`, `has_dups()`, `has_variety()`, `rotations()`,
`prefixes()`, `postfixes()`, `numbered_symbols()`, `dict_merge()` and
`topological_sort()`.

## Usage

``` r
reshape_seq(seq, how)

multiset(seq)

ibin(n, bits = NULL, str = FALSE)

common_prefix(...)

common_suffix(...)

strongly_connected_components(G)

connected_components(G)

rotate_left(x, y)

rotate_right(x, y)

least_rotation(x, key = identity)

minlex(seq, directed = TRUE, key = identity)

multiset_combinations(m, n, g = NULL)

multiset_permutations(m, size = NULL, g = NULL)

ordered_partitions(n, m = NULL, sort = TRUE)

binary_partitions(n)

uniq(seq)

generate_bell(n)

generate_involutions(n)

multiset_derangements(s)

generate_derangements(s)

random_derangement(t, choice = NULL, strict = TRUE)

necklaces(n, k, free = FALSE)

bracelets(n, k)

generate_oriented_forest(n)

runs(seq, op = function(a, b) a > b)

sequence_partitions(l, n)

sequence_partitions_empty(l, n)

kbins(l, k, ordered = NULL)

permute_signs(t)

signed_permutations(t)

roundrobin(...)

partitions(n, m = NULL, k = NULL, size = FALSE)

multiset_partitions(m, k = NULL)

flatten(x)

unflatten(x, n = 2)

group(x, multiple = TRUE)

take(x, n)

variations(x, n, repetition = FALSE)

sift(x, f)

has_dups(x)

has_variety(x)

rotations(x)

prefixes(x)

postfixes(x)

numbered_symbols(prefix = "x", n = 10, start = 0)

dict_merge(...)

topological_sort(vertices, edges)
```

## Arguments

- seq, x, l, s, t:

  A sequence (vector or list); for `take()` also a generator.

- how:

  A nested list of counts.

- n, y, size:

  Integers.

- bits:

  Number of bits, or `"all"` for all bit strings of length `n`.

- str:

  Return strings.

- ...:

  Lists to merge.

- G:

  A graph `list(V, E)` with `E` a list of pairs.

- key:

  A key function for comparisons.

- directed:

  Directed rotation only.

- m:

  A sequence, or a named vector of counts, for the multiset functions;
  for `partitions()` the largest number of parts.

- g:

  Internal: the grouped form of `m`, as in SymPy.

- sort:

  Sort the partitions.

- choice:

  A function choosing an element of a vector.

- strict:

  Fail when no derangement exists.

- k:

  An integer; for `partitions()` the largest part, for
  `multiset_partitions()` the number of parts.

- free:

  Count rotations and reflections as equal.

- op:

  A binary predicate for `runs()`.

- ordered:

  `NULL`, `0`, `1`, `10` or `11`, as SymPy's `kbins()`.

- multiple:

  For `group()`, return runs instead of value and count.

- repetition:

  Allow repeated elements.

- f:

  A key function.

- prefix:

  Symbol prefix.

- start:

  First index.

- vertices, edges:

  A graph, with edges as pairs.

## Value

A list (or a vector or string where SymPy returns one).

A list or logical.

## Examples

``` r
ordered_partitions(5, 3)
#> [[1]]
#> [1] 1 1 3
#> 
#> [[2]]
#> [1] 1 2 2
#> 
kbins(1:4, 2)
#> [[1]]
#> [[1]][[1]]
#> [1] 1
#> 
#> [[1]][[2]]
#> [1] 2 3 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 1 2
#> 
#> [[2]][[2]]
#> [1] 3 4
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 1 2 3
#> 
#> [[3]][[2]]
#> [1] 4
#> 
#> 
multiset_permutations(c(1, 1, 2))
#> [[1]]
#> [1] 1 1 2
#> 
#> [[2]]
#> [1] 1 2 1
#> 
#> [[3]]
#> [1] 2 1 1
#> 
binary_partitions(5)
#> [[1]]
#> [1] 4 1
#> 
#> [[2]]
#> [1] 2 2 1
#> 
#> [[3]]
#> [1] 2 1 1 1
#> 
#> [[4]]
#> [1] 1 1 1 1 1
#> 
flatten(list(1, list(2, list(3))))
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
group(c(1, 1, 2, 3, 3, 3))
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3 3 3
#> 
topological_sort(1:4, list(c(1, 2), c(2, 3), c(1, 4)))
#> [1] 1 2 4 3
```
