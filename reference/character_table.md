# Character tables of finite groups

`character_table(G)` computes the character table of a permutation
group, as SymPy's `PermutationGroup.character_table()` and
`CharacterTable.from_perm_group()`: rows are the irreducible characters,
with the trivial character first, and columns the conjugacy classes in
the order of
[`conjugacy_classes()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md).
Symmetric groups use `character_table_sn()`, the Murnaghan-Nakayama
rule, whose rows are indexed by the partitions of `n` in reverse
lexicographic order and whose columns by the partitions in lexicographic
order. Other groups use `dixon_character_table()`, Dixon's algorithm on
a list of conjugacy classes over a prime field, with the values lifted
to the smallest cyclotomic field `QQ(zeta_k)` that holds them.

## Usage

``` r
character_table(G)

dixon_character_table(conjugacy_classes)

character_table_sn(n)

ct_tolist(tbl)

ct_as_matrix(tbl)

ct_shape(tbl)

ct_zeta_order(tbl)

ct_conjugacy_class_reps(tbl)
```

## Arguments

- G:

  A permutation group.

- conjugacy_classes:

  A list of conjugacy classes, each a list of permutations; one of them
  must hold the identity.

- n:

  The degree of the symmetric group, a positive integer.

- tbl:

  A `CharacterTable`.

## Value

A `CharacterTable`; the accessors return a list of lists of expressions,
a matrix, an integer vector, an integer or a list of permutations.

## Details

The table is a `CharacterTable` object: `ct_tolist()` gives its rows as
lists of expressions in the root of unity `zeta_k` (an
[`AlgebraicNumber()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)),
`ct_as_matrix()` the matrix, `ct_shape()` its dimensions,
`ct_zeta_order()` the order `k` of the root of unity (1 for an integer
table) and `ct_conjugacy_class_reps()` the first element of each
conjugacy class.

## Examples

``` r
ct_tolist(character_table_sn(3))
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> [[1]][[2]]
#> 1
#> 
#> [[1]][[3]]
#> 1
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 2
#> 
#> [[2]][[2]]
#> 0
#> 
#> [[2]][[3]]
#> -1
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> 1
#> 
#> [[3]][[2]]
#> -1
#> 
#> [[3]][[3]]
#> 1
#> 
#> 
ct_as_matrix(character_table(AlternatingGroup(4)))
#> Matrix(4, 4)
#> [1, 1         , 1         , 1 ]
#> [1, -1 - zeta3, zeta3     , 1 ]
#> [1, zeta3     , -1 - zeta3, 1 ]
#> [3, 0         , 0         , -1]
```
