# Limits of sequences

`limit_seq(expr, n)` is the limit of the sequence `expr` as the integer
`n` grows, as SymPy's `limit_seq`: it rewrites Fibonacci numbers and
factorials, tries the limit at infinity, and otherwise applies the
Stolz-Cesaro theorem to the dominant terms of numerator and denominator.
Oscillating sequences give `AccumBounds`; `NULL` means no answer was
found. `dominant(expr, n)` is the term of a sum that grows fastest in
`n`.

## Usage

``` r
dominant(expr, n)

limit_seq(expr, n = NULL, trials = 5L)
```

## Arguments

- expr:

  An expression in `n`.

- n:

  The index symbol.

- trials:

  How many times to apply Stolz-Cesaro.

## Value

An expression, `AccumBounds` or `NULL`.

## Examples

``` r
n <- Symbol("n", integer = TRUE)
limit_seq(n / (n + 1), n)
#> 1
limit_seq((-1)^n, n)
#> AccumBounds(-1, 1)
dominant(5 * n^3 + 4 * n^2 + n + 1, n)
#> 5*n^3
```
