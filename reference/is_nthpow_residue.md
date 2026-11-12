# Power residues and n-th roots modulo m

`is_nthpow_residue(a, n, m)` tests whether `x^n = a (mod m)` has a
solution, and `nthroot_mod(a, n, p)` returns the smallest solution,
every solution with `all_roots = TRUE`, or `NULL` (an empty vector) when
there is none.

## Usage

``` r
is_nthpow_residue(a, n, m)

nthroot_mod(a, n, p, all_roots = FALSE)
```

## Arguments

- a:

  An integer.

- n:

  Exponent.

- m, p:

  Modulus.

- all_roots:

  Return all roots.

## Value

A logical, an integer, a vector or `NULL`.

## Examples

``` r
is_nthpow_residue(2, 3, 7)
#> [1] FALSE
nthroot_mod(11, 4, 19)
#> [1] 8
nthroot_mod(68, 3, 109, all_roots = TRUE)
#> [1] 23 32 54
```
