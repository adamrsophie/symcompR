# Integer helpers

`igcd()` and `ilcm()` are the gcd and lcm of integers.
`mod_inverse(a, m)` is the inverse of `a` modulo `m` (for a rational
`a`, its inverse in the field). `integer_nthroot(y, n)` returns
`list(root, exact)` with `root` the floor of the `n`-th root;
`integer_log(y, x)` returns `list(e, exact)` with
`x^e <= y < x^(e + 1)`. `trailing(n)` counts trailing zero bits and
`num_digits(n, base)` counts digits. `arity(f)` is the number of
arguments of an R function.

## Usage

``` r
igcd(...)

ilcm(...)

mod_inverse(a, m)

integer_nthroot(y, n)

integer_log(y, x)

trailing(n)

num_digits(n, base = 10)

arity(f)
```

## Arguments

- ...:

  Integers.

- a, m, y, x, n:

  Integers.

- base:

  Base.

- f:

  An R function.

## Value

An integer, `sym` number or list.

## Examples

``` r
igcd(12, 18, -8)
#> [1] 2
mod_inverse(3, 11)
#> [1] 4
integer_nthroot(26, 3)
#> $root
#> [1] 2
#> 
#> $exact
#> [1] FALSE
#> 
integer_log(125, 5)
#> $e
#> [1] 3
#> 
#> $exact
#> [1] TRUE
#> 
```
