# Numerical and real roots

`nroots(f, x)` returns all roots as floating point numbers (via
[`polyroot()`](https://rdrr.io/r/base/polyroot.html) with Newton
polishing). `real_roots(f, x)` returns the real roots with multiplicity,
exact where possible. `count_roots(f, x, a, b)` counts real roots in
`[a, b]` with a Sturm sequence. `sturm()` returns the Sturm sequence.

## Usage

``` r
nroots(f, x = NULL)

real_roots(f, x = NULL)

sturm(f, x = NULL)

count_roots(f, x = NULL, a = -Inf, b = Inf)
```

## Arguments

- f:

  A polynomial.

- x:

  The variable; when `k` is missing, `x` is the index and the variable
  is the single free symbol of `f`.

- a, b:

  Interval ends (default the whole line).

## Value

A list of expressions, or an integer.

## Examples

``` r
x <- Symbol("x")
nroots(x^5 - x + 1, x)
#> [[1]]
#> -1.16730397826142
#> 
#> [[2]]
#> -0.181232444469875 - 1.08395410131771*I
#> 
#> [[3]]
#> -0.181232444469875 + 1.08395410131771*I
#> 
#> [[4]]
#> 0.764884433600585 - 0.352471546031726*I
#> 
#> [[5]]
#> 0.764884433600585 + 0.352471546031726*I
#> 
count_roots(x^3 - 2*x, x)
#> [1] 3
```
