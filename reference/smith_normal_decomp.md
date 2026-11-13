# Smith normal form decomposition and test

`smith_normal_decomp(m)` returns the Smith normal form `smf` of a
[`DomainMatrix()`](https://adamrsophie.github.io/symcompR/reference/DomainMatrix.md)
over a principal ideal domain with invertible matrices `s` and `t` such
that `smf == s * m * t`; `is_smith_normal_form(m)` tests whether a
matrix is in Smith normal form.
[`smith_normal_form()`](https://adamrsophie.github.io/symcompR/reference/hermite_normal_form.md),
[`invariant_factors()`](https://adamrsophie.github.io/symcompR/reference/hermite_normal_form.md)
and
[`hermite_normal_form()`](https://adamrsophie.github.io/symcompR/reference/hermite_normal_form.md)
also accept a `DomainMatrix`.

## Usage

``` r
smith_normal_decomp(m)

is_smith_normal_form(m)
```

## Arguments

- m:

  A `DomainMatrix`.

## Value

A tuple of three matrices, or `TRUE` or `FALSE`.

## Examples

``` r
m <- DM(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14)), ZZ)
smith_normal_decomp(m)
#> (DomainMatrix([[1, 0, 0], [0, 10, 0], [0, 0, 30]], (3, 3), ZZ), DomainMatrix([[0, 1, -1], [1, -4, 0], [0, -2, 3]], (3, 3), ZZ), DomainMatrix([[1, 1, 10], [0, -1, -2], [0, 1, 3]], (3, 3), ZZ)) 
is_smith_normal_form(DM(list(list(1, 0), list(0, 2)), ZZ))
#> [1] TRUE
```
