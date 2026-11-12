# Digits of an integer

`digits(n, b)` returns the base `b` digits of `n`, most significant
first, preceded by `b` (or `-b` for negative `n`) as in SymPy; `digits`
pads with zeros to that many digits. `count_digits(n, b)` counts each
digit value and `is_palindromic(n, b)` tests whether the digits read the
same backwards.

## Usage

``` r
digits(n, b = 10, digits = NULL)

count_digits(n, b = 10)

is_palindromic(n, b = 10)
```

## Arguments

- n:

  An integer.

- b:

  Base, at least 2.

- digits:

  Minimum number of digits, or `NULL`.

## Value

A numeric vector, a named count vector, or a logical.

## Examples

``` r
digits(35)
#> [1] 10  3  5
digits(27, b = 2)
#> [1] 2 1 1 0 1 1
digits(-3958, 27)
#> [1] -27   5  11  16
count_digits(1111339)
#> 1 3 9 
#> 4 2 1 
is_palindromic(88, 3)
#> [1] FALSE
```
