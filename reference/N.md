# Numerical evaluation

`N()` (alias `evalf()`) evaluates an expression to a floating point
number, as SymPy's `evalf`: arbitrary precision binary arithmetic on gmp
integers that tracks the accuracy of the real and imaginary parts and
raises the working precision where terms cancel. Functions without an
arbitrary precision evaluator fall back to R double precision, capped at
15 digits.

## Usage

``` r
N(
  e,
  digits = 15,
  subs = list(),
  maxn = 100,
  chop = FALSE,
  strict = FALSE,
  quad = NULL,
  verbose = FALSE
)

evalf(
  e,
  digits = 15,
  subs = list(),
  maxn = 100,
  chop = FALSE,
  strict = FALSE,
  quad = NULL,
  verbose = FALSE
)
```

## Arguments

- e:

  A `sym` expression.

- digits:

  Significant digits.

- subs:

  Named list of values for symbols.

- maxn:

  Maximum working precision in digits.

- chop:

  `TRUE` to replace tiny real or imaginary parts by exact zeros, or a
  number giving the size of a tiny part.

- strict:

  `TRUE` to raise an error when a part cannot be evaluated to full
  accuracy within `maxn`.

- quad:

  Quadrature method for integrals (`NULL` for tanh-sinh).

- verbose:

  `TRUE` to print the intermediate results.

## Value

A `sym` Float, a complex `a + b*I` with Float parts, or the expression
with numeric subexpressions evaluated when free symbols remain.

## Examples

``` r
N(sym_pi)
#> 3.14159265358979
N(sqrt(S(2)), 5)
#> 1.4142
x <- Symbol("x")
N(sin(x), subs = list(x = 1))
#> 0.841470984807897
N(parse_expr("(pi + pi*I)**2"), chop = TRUE)
#> 19.7392088021787*I
```
