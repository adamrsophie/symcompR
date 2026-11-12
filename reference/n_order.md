# Orders and primitive roots

`n_order(a, n)` is the multiplicative order of `a` modulo `n`.
`is_primitive_root(a, p)` tests whether `a` generates the units modulo
`p`, and `primitive_root(p)` returns the smallest primitive root, or
with `smallest = FALSE` one found quickly; it returns `NULL` when none
exists.

## Usage

``` r
n_order(a, n)

is_primitive_root(a, p)

primitive_root(p, smallest = TRUE)
```

## Arguments

- a:

  An integer coprime to the modulus.

- n, p:

  Modulus, greater than 1.

- smallest:

  Return the smallest primitive root.

## Value

An integer, logical or `NULL`.

## Examples

``` r
n_order(3, 7)
#> [1] 6
is_primitive_root(3, 10)
#> [1] TRUE
primitive_root(19)
#> [1] 2
primitive_root(50)
#> [1] 3
```
