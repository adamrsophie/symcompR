# Content and primitive part

`as_content_primitive(e)` returns `list(c, p)` with `e = c*p`, `c` a
positive rational and `p` with integer coefficients of gcd one, as
SymPy's `Expr.as_content_primitive()`. `radical = TRUE` also takes out
common radicals of a sum; `clear = FALSE` keeps a denominator inside a
sum when some term would keep an integer coefficient.

## Usage

``` r
as_content_primitive(e, radical = FALSE, clear = TRUE)
```

## Arguments

- e:

  An expression.

- radical, clear:

  See above.

## Value

A list of the content and the primitive part.

## Examples

``` r
x <- Symbol("x")
as_content_primitive((2*x + 4)^2)
#> [[1]]
#> 4
#> 
#> [[2]]
#> (x + 2)^2
#> 
as_content_primitive(2*sqrt(2) + 4*sqrt(10), radical = TRUE)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 15.4775377654197
#> 
```
