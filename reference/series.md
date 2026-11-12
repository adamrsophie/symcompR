# Series expansion

`series(e, x, x0, n)` is SymPy's `Expr.series`: the expansion of `e`
about `x = x0` up to (but excluding) `(x - x0)^n`, with an
[`O()`](https://adamrsophie.github.io/symcompR/reference/O.md) term
(none when the expansion is exact). `dir = "-"` expands from the left;
`x0` may be `oo`. With `n = NULL` it returns SymPy's lazy series, a
generator of the terms one at a time:
[`take()`](https://adamrsophie.github.io/symcompR/reference/iterables.md)
its first terms. `nseries()` is `Expr.nseries` and `leadterm()` returns
the coefficient and exponent of the leading term at 0;
`as_leading_term()` the leading term itself.

## Usage

``` r
series(e, x = NULL, x0 = 0, n = 6, dir = "+", logx = NULL, cdir = 0)

nseries(e, x = NULL, x0 = 0, n = 6, dir = "+", logx = NULL, cdir = 0)

leadterm(e, x, logx = NULL, cdir = 0)

aseries(e, x, n = 6, bound = 0, hir = FALSE)

as_leading_term(e, x, logx = NULL, cdir = 0)
```

## Arguments

- e:

  An expression.

- x:

  The variable (default: the only free symbol).

- x0:

  The expansion point.

- n:

  Number of terms (order of truncation), or `NULL` for the lazy series.

- dir:

  `"+"` or `"-"`.

- logx:

  Optional symbol standing for `log(x)`.

- cdir:

  The direction as a complex number (overrides `dir`).

- bound, hir:

  Arguments of SymPy's `aseries()`: the recursion bound for exponential
  terms and whether to return the hierarchical form.

## Value

An expression (a generator for `n = NULL`); for `leadterm()` a list
`(coeff, exponent)`.

## Examples

``` r
x <- Symbol("x")
series(sin(x), x, 0, 8)
#> x - x^3/6 + x^5/120 - x^7/5040 + O(x^8)
series(exp(x) / x, x, 0, 3)
#> 1/x + 1 + x/2 + x^2/6 + O(x^3)
series(log(1 + x), x)
#> x - x^2/2 + x^3/3 - x^4/4 + x^5/5 + O(x^6)
take(series(cos(x), x, 1, n = NULL), 3)
#> [[1]]
#> cos(1)
#> 
#> [[2]]
#> -(x - 1)*sin(1)
#> 
#> [[3]]
#> -(x - 1)^2*cos(1)/2
#> 
```
