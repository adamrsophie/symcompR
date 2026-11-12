# Sparse polynomial rings and rational function fields

`ring(symbols, domain)` returns a polynomial ring over `domain` followed
by its generators, as SymPy's `ring()`; `xring()` returns the ring and a
list of the generators, `vring()` assigns the generators in `envir` and
returns the ring, and `sring(exprs)` builds the ring of a list of
expressions and returns it with their elements. `field()`, `xfield()`,
`vfield()` and `sfield()` do the same for fields of rational functions.
Ring elements are polynomials kept expanded; field elements are rational
functions kept in lowest terms.

## Usage

``` r
# S3 method for class 'sym_ring'
format(x, ...)

# S3 method for class 'sym_ring'
print(x, ...)

# S3 method for class 'sym_ring_element'
format(x, ...)

# S3 method for class 'sym_ring_element'
print(x, ...)

# S3 method for class 'sym_ring_element'
Ops(e1, e2)

# S3 method for class 'sym_ring_element'
chooseOpsMethod(x, y, mx, my, cl, reverse)

ring(symbols, domain = ZZ, order = "lex")

xring(symbols, domain = ZZ, order = "lex")

vring(symbols, domain = ZZ, order = "lex", envir = parent.frame())

sring(exprs, ...)

field(symbols, domain = ZZ, order = "lex")

xfield(symbols, domain = ZZ, order = "lex")

vfield(symbols, domain = ZZ, order = "lex", envir = parent.frame())

sfield(exprs, ...)
```

## Arguments

- ...:

  Unused.

- e1, e2, x:

  Ring or field elements.

- y, mx, my, cl, reverse:

  Arguments of the
  [`chooseOpsMethod()`](https://rdrr.io/r/base/chooseOpsMethod.html)
  generic.

- symbols:

  Names as a comma separated string, or a list of symbols.

- domain:

  The coefficient domain.

- order:

  The monomial order (only `"lex"`).

- envir:

  The environment receiving the generators.

- exprs:

  An expression or a list of them.

## Value

A list with the ring or field and the generators or elements.

## Examples

``` r
r <- ring("a,b", ZZ); R <- r[[1L]]; a <- r[[2L]]; b <- r[[3L]]
R
#> Polynomial ring in a, b over ZZ with lex order
(a + b)^2
#> a^2 + 2*a*b + b^2
f <- field("u,v", QQ); u <- f[[2L]]; v <- f[[3L]]
(u^2 - v^2)/(u - v)
#> u + v
```
