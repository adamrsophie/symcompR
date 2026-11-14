# Numbers of groups of a given order

`groups_count(n)` is the number of groups of order `n` up to isomorphism
where SymPy knows it; `is_nilpotent_number(n)`, `is_abelian_number(n)`
and `is_cyclic_number(n)` tell whether every group of order `n` is
nilpotent, abelian or cyclic.

## Usage

``` r
groups_count(n)

is_nilpotent_number(n)

is_abelian_number(n)

is_cyclic_number(n)
```

## Arguments

- n:

  A positive integer.

## Value

An integer or a logical.

## Examples

``` r
groups_count(16)
#> [1] 14
groups_count(30)
#> [1] 4
is_cyclic_number(15)
#> [1] TRUE
```
