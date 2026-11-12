# Number theoretic transform

`ntt(seq, prime)` and `intt(seq, prime)` are the discrete Fourier
transform and its inverse over the integers modulo a prime of the form
`m 2^k + 1`. The sequence is padded with zeros to a power of two.

## Usage

``` r
ntt(seq, prime)

intt(seq, prime)
```

## Arguments

- seq:

  A vector of integers.

- prime:

  A prime modulus.

## Value

A vector of integers.

## Examples

``` r
ntt(c(1, 2, 3, 5), 3 * 2^8 + 1)
#> [1]  11 581 766 184
intt(ntt(c(1, 2, 3, 5), 769), 769)
#> [1] 1 2 3 5
```
