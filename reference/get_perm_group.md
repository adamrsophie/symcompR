# Permutation group of a named transitive subgroup

`get_perm_group(name)` returns the permutation group for a name from
[S1TransitiveSubgroups](https://adamrsophie.github.io/symcompR/reference/TransitiveSubgroups.md)
to
[S6TransitiveSubgroups](https://adamrsophie.github.io/symcompR/reference/TransitiveSubgroups.md),
as the method of the same name in SymPy. The constructors
`four_group()`, `M20()`, `S3_in_S6()`, `A4_in_S6()`, `S4m()`, `S4p()`,
`A4xC2()`, `S4xC2()`, `G18()`, `G36m()`, `G36p()`, `G72()`, `PSL2F5()`
and `PGL2F5()` build the groups that are not symmetric, alternating,
cyclic or dihedral, with SymPy's generators.

## Usage

``` r
get_perm_group(name)

four_group()

M20()

S3_in_S6()

A4_in_S6()

S4m()

S4p()

A4xC2()

S4xC2()

G18()

G36m()

G36p()

G72()

PSL2F5()

PGL2F5()
```

## Arguments

- name:

  A `TransitiveSubgroupName`.

## Value

A `PermutationGroup`.

## Examples

``` r
G <- get_perm_group(S4TransitiveSubgroups$V)
group_order(G)
#> [1] 4
group_order(M20())
#> [1] 20
```
