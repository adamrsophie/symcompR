# Solve univariate inequalities

Solves polynomial and rational inequalities in one real variable by
locating the real roots of numerator and denominator and testing the
sign on each piece, like SymPy's `solve_univariate_inequality()`.
Returns a relational (default) or a set with `relational = FALSE`.

## Usage

``` r
solve_univariate_inequality(expr, x, relational = TRUE, domain = Reals)
```

## Arguments

- expr:

  An inequality such as `x^2 > 4`.

- x:

  The variable.

- relational:

  Return a relational instead of a set.

- domain:

  The domain, by default the real line.

## Value

A boolean expression or a set.

## Examples

``` r
x <- Symbol("x")
solve_univariate_inequality(x^2 > 4, x)
#> ((-oo < x) & (x < -2)) | ((2 < x) & (x < oo))
solve_univariate_inequality((x - 1) / (x + 2) <= 0, x, relational = FALSE)
#> Interval.Lopen(-2, 1)
```
