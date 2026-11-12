# Square free and other cores, and digital roots

`core(n, t)` is the `t`-free part of `n`: the product of `p^(e %% t)`
over the prime powers `p^e` of `n`. `dra(n, b)` is the additive digital
root of `n` in base `b` and `drm(n, b)` the multiplicative one.

## Usage

``` r
core(n, t = 2)

dra(n, b)

drm(n, b)
```

## Arguments

- n:

  A positive integer (any integer for `dra()` and `drm()`).

- t:

  Power, at least 2.

- b:

  Base, at least 2.

## Value

An integer.

## Examples

``` r
core(24)
#> [1] 6
core(5^6 * 7^5, 3)
#> [1] 49
dra(19, 12)
#> [1] 8
drm(277777788888899, 10)
#> [1] 0
```
