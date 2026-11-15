# Power series rings

A port of SymPy's `sympy.polys.series`:
`power_series_ring(symbol, K, prec)` returns the ring of power series in
one variable over `ZZ` or `QQ` truncated at `O(x^prec)`, together with
its generator, as `list(R, x)`. Series are added, multiplied, divided
and raised to integer powers with the usual operators; an exact
polynomial has no
[`O()`](https://adamrsophie.github.io/symcompR/reference/O.md) term
until it grows past the ring precision. The ring `R` holds the methods
of SymPy's ring as functions: `R$compose(s, t)`, `R$inverse(s)`,
`R$reversion(s)`, `R$differentiate(s)`, `R$truncate(s, n)`,
`R$square(s)`, `R$from_list(coeffs, prec)`, `R$from_expr(e)`,
`R$to_list(s)`, `R$to_expr(s)` and `R$order_term()`; over `QQ` also
`R$integrate()`, `R$sqrt()`, `R$log()`, `R$log1p()`, `R$exp()`,
`R$expm1()`, `R$sin()`, `R$cos()`, `R$tan()`, `R$sinh()`, `R$cosh()`,
`R$tanh()`, `R$asin()`, `R$asinh()`, `R$atan()`, `R$atanh()` and
`R$hypot()`. `PowerSeriesRingRing()` and `PowerSeriesRingField()` build
the rings over `ZZ` and `QQ` directly.

## Usage

``` r
power_series_ring(symbol, K, prec = 6)

PowerSeriesRingRing(domain = ZZ, symbol = "x", prec = 6)

PowerSeriesRingField(domain = QQ, symbol = "x", prec = 6)
```

## Arguments

- symbol:

  The variable name, or a symbol.

- K:

  The ground domain, `ZZ` or `QQ`.

- prec:

  The ring precision.

- domain:

  The ground domain.

## Value

`power_series_ring()` a list of the ring and its generator; the other
constructors a ring.

## Examples

``` r
r <- power_series_ring("x", QQ, 6)
R <- r[[1]]; x <- r[[2]]
(1 + x + x^2)^3
#> 1 + 3*x + 6*x^2 + 7*x^3 + 6*x^4 + 3*x^5 + O(x^6) 
R$sin(x)
#> x - 1/6*x^3 + 1/120*x^5 + O(x^6) 
R$compose(R$exp(x) - 1, x + x^2)
#> x + 3/2*x^2 + 7/6*x^3 + 25/24*x^4 + 27/40*x^5 + O(x^6) 
R$reversion(R$sin(x))
#> x + 1/6*x^3 + 3/40*x^5 + O(x^6) 
R$to_expr(R$tan(x))
#> x + x^3/3 + 2*x^5/15 + O(x^6)
```
