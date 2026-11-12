# Hexadecimal digits of pi

`pi_hex_digits(n, prec)` returns `prec` hexadecimal digits of pi
starting after the `n`-th, by the Bailey-Borwein-Plouffe formula, so
`pi_hex_digits(0)` begins with the leading 3.

## Usage

``` r
pi_hex_digits(n, prec = 14)
```

## Arguments

- n:

  Starting position, 0 or more.

- prec:

  Number of digits.

## Value

A string.

## Examples

``` r
pi_hex_digits(0)
#> [1] "3243f6a8885a30"
pi_hex_digits(1, 10)
#> [1] "243f6a8885"
```
