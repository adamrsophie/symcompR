# Binomial and multinomial coefficients

`binomial_coefficients(n)` returns the coefficients of `(x + y)^n` named
by the exponent pair, as in `"1,2"`, and `binomial_coefficients_list(n)`
the row of Pascal's triangle. `multinomial_coefficients(m, n)` returns
the coefficients of `(x1 + ... + xm)^n` named by exponent tuples.
`multinomial_coefficients_iterator(m, n)` gives the same coefficients in
the order of SymPy's iterator, which is cheaper when `m` is much larger
than `n`. Values are doubles while exact, otherwise a named list of big
integers.

## Usage

``` r
binomial_coefficients(n)

binomial_coefficients_list(n)

multinomial_coefficients(m, n)

multinomial_coefficients_iterator(m, n)
```

## Arguments

- n:

  Power.

- m:

  Number of terms.

## Value

A named vector (or list), or a vector for
`binomial_coefficients_list()`.

## Examples

``` r
binomial_coefficients(4)
#> 0,4 1,3 2,2 3,1 4,0 
#>   1   4   6   4   1 
binomial_coefficients_list(5)
#> [1]  1  5 10 10  5  1
multinomial_coefficients(3, 2)
#> 0,0,2 0,1,1 0,2,0 1,0,1 1,1,0 2,0,0 
#>     1     2     1     2     2     1 
```
