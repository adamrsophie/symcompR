# Chinese remainder theorem

`crt(m, v)` returns `c(x, M)` with `x %% m[i] == v[i] %% m[i]` for every
`i` and `M` the product of the moduli, or `NULL` when the congruences
are inconsistent; `symmetric = TRUE` puts `x` within half of `M`. Moduli
that are not coprime are handled by `solve_congruence()`, whose
arguments are pairs `c(remainder, modulus)`. `crt1(m)` precomputes for
repeated use with `crt2(m, v, mm, e, s)`. `symmetric_residue(a, m)` maps
`a` in `[0, m)` to the residue nearest zero.

## Usage

``` r
crt(m, v, symmetric = FALSE, check = TRUE)

crt1(m)

crt2(m, v, mm, e, s, symmetric = FALSE)

symmetric_residue(a, m)

solve_congruence(..., symmetric = FALSE, check = TRUE)
```

## Arguments

- m:

  Moduli.

- v:

  Remainders.

- symmetric:

  Return the residue nearest zero.

- check:

  Check the result, falling back to `solve_congruence()`.

- mm, e, s:

  Values from `crt1()`.

- a:

  A residue.

- ...:

  Pairs `c(remainder, modulus)`.

## Value

`c(x, M)` or `NULL`; `crt1()` returns `list(mm, e, s)`.

## Examples

``` r
crt(c(99, 97, 95), c(49, 76, 65))
#> [1] 639985 912285
solve_congruence(c(2, 3), c(3, 5), c(2, 7))
#> [1]  23 105
symmetric_residue(3, 6)
#> [1] 3
```
