# Small helpers

`polarify(e, subs, lift)` replaces the symbols of `e` by polar dummies
and returns `list(expression, mapping)` (or the expression alone when
`subs = FALSE`). `sym_prod(l)` is the product of a list (SymPy's
[`prod()`](https://rdrr.io/r/base/prod.html)). `sym_reshape(seq, how)`
(SymPy's [`reshape()`](https://rdrr.io/r/stats/reshape.html), renamed so
it does not mask
[`stats::reshape()`](https://rdrr.io/r/stats/reshape.html)) groups a
sequence by the counts in `how`. `parallel_poly_from_expr(exprs)` builds
polynomials over common generators. `sqf_norm(f, extension)` returns
`list(s, g, r)` with `g = f(x - s*a)` for the algebraic number `a` and
its square free norm `r`. `gff()` is not implemented, as in SymPy.
`ilex()`, `igrlex()` and `igrevlex()` are the inverse monomial order
keys. `carmichael` collects SymPy's Carmichael helpers.

## Usage

``` r
polarify(e, subs = TRUE, lift = FALSE)

sym_prod(l)

sym_reshape(seq, how)

parallel_poly_from_expr(exprs, ...)

sqf_norm(f, extension)

gff(f, ...)

ilex(m)

igrlex(m)

igrevlex(m)

carmichael
```

## Arguments

- e, f:

  An expression.

- subs:

  Replace symbols by polar dummies.

- lift:

  Wrap numbers in
  [`polar_lift()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md).

- l, seq:

  A list.

- how:

  A vector of counts, possibly nested.

- exprs:

  Expressions.

- ...:

  Generators.

- extension:

  An algebraic number.

- m:

  An exponent vector.

## Value

Varies.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
polarify(exp(x)*y)
#> [[1]]
#> _y*exp(_x)
#> 
#> [[2]]
#> [[2]][[1]]
#> [[2]][[1]][[1]]
#> _x
#> 
#> [[2]][[1]][[2]]
#> x
#> 
#> 
#> [[2]][[2]]
#> [[2]][[2]][[1]]
#> _y
#> 
#> [[2]][[2]][[2]]
#> y
#> 
#> 
#> 
sym_reshape(1:6, 2)
#> [[1]]
#> [[1]][[1]]
#> [1] 1
#> 
#> [[1]][[2]]
#> [1] 2
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 3
#> 
#> [[2]][[2]]
#> [1] 4
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 5
#> 
#> [[3]][[2]]
#> [1] 6
#> 
#> 
sqf_norm(x^2 - 2, sqrt(Integer(3)))
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> 
#> [[2]]
#> x^2 - 2*sqrt(3)*x + 1
#> 
#> [[3]]
#> x^4 - 10*x^2 + 1
#> 
```
