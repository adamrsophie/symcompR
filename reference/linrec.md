# Linear recurrences with constant coefficients

`linrec(coeffs, init, n)` evaluates the `n`-th term of the sequence with
`a(n) = coeffs[1] a(n - 1) + ... + coeffs[k] a(n - k)` and initial terms
`init` (padded with zeros), by SymPy's square and reduce method in
`O(k^2 log n)` operations. Coefficients and initial terms may be
symbolic.

## Usage

``` r
linrec(coeffs, init, n)
```

## Arguments

- coeffs:

  Coefficients of the recurrence.

- init:

  Initial terms.

- n:

  A nonnegative integer.

## Value

A `sym` expression.

## Examples

``` r
linrec(c(1, 1), c(0, 1), 10)
#> 55
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
linrec(list(x, y, z), c(1, 1, 1), 4)
#> x^2 + x*y + x*z + y + z
```
