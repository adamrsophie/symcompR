# Fourier series

`fourier_series(f, limits, finite)` is the Fourier series of `f` over
`limits = list(x, a, b)` (default `list(x, -pi, pi)`). Trigonometric
polynomials give a `FiniteFourierSeries`. Use
[`truncate()`](https://rdrr.io/r/base/seek.html) for the first `n`
nonzero terms, `sigma_approximation()` for the Lanczos smoothed sum, and
`shift()`, `shiftx()`, `scale(s, k)`, `scalex()` to transform the series
as SymPy does; `+` and `-` add series with the same period.

## Usage

``` r
fourier_series(f, limits = NULL, finite = TRUE)

# S3 method for class 'sym'
truncate(con, n, ...)

sigma_approximation(s, n = 3)

shift(s, k)

shiftx(s, k)

# S3 method for class 'sym'
scale(x, center = TRUE, scale = TRUE)

scalex(s, k)
```

## Arguments

- f:

  An expression.

- limits:

  `list(x, a, b)` or `list(a, b)`.

- finite:

  Detect trigonometric polynomials.

- n:

  Number of terms; `NULL` in
  [`truncate()`](https://rdrr.io/r/base/seek.html) returns all terms of
  a finite series.

- ...:

  Unused.

- s, x, con:

  A Fourier series (`con` for
  [`truncate()`](https://rdrr.io/r/base/seek.html), `x` for
  [`scale()`](https://rdrr.io/r/base/scale.html)).

- k:

  The shift or scale amount.

- center:

  The factor for `scale(s, k)` (the second argument of the
  [`base::scale()`](https://rdrr.io/r/base/scale.html) generic).

- scale:

  Ignored; part of the
  [`base::scale()`](https://rdrr.io/r/base/scale.html) generic.

## Value

A Fourier series or an expression.

## Examples

``` r
x <- Symbol("x")
s <- fourier_series(x^2, list(x, -pi, pi))
truncate(s, 4)
#> -39.4784176043574*cos(0.318309886183791*pi*x)/pi^2 + 9.86960440108936*cos(0.636619772367581*pi*x)/pi^2 - 4.3864908449286*cos(0.954929658551372*pi*x)/pi^2 + 3.28986813369645
# \donttest{
fourier_series(sin(x)^2)
#> FiniteFourierSeries(sin(x)^2, (x, -pi, pi), (1/2, {2: -1/2}, {}))
# }
```
