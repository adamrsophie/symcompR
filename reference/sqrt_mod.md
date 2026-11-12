# Square roots modulo n

`sqrt_mod(a, p)` returns one root of `x^2 = a (mod p)`, at most `p / 2`,
or `NULL`; `all_roots = TRUE` returns every root in increasing order.
`sqrt_mod_iter(a, p)` returns every root in the order SymPy's iterator
produces them. `quadratic_residues(p)` lists the squares modulo `p`, and
`is_quad_residue(a, p)` tests one.

## Usage

``` r
sqrt_mod(a, p, all_roots = FALSE)

sqrt_mod_iter(a, p, domain = NULL)

quadratic_residues(p)

is_quad_residue(a, p)
```

## Arguments

- a:

  An integer.

- p:

  Modulus.

- all_roots:

  Return all roots.

- domain:

  Function applied to each root, for example
  [`Integer()`](https://adamrsophie.github.io/symcompR/reference/Integer.md).

## Value

An integer, a vector, `NULL` or a logical.

## Examples

``` r
sqrt_mod(11, 43)
#> [1] 21
sqrt_mod(17, 32, all_roots = TRUE)
#> [1]  7  9 23 25
quadratic_residues(13)
#> [1]  0  1  3  4  9 10 12
is_quad_residue(21, 100)
#> [1] TRUE
```
