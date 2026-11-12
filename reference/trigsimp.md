# Trigonometric simplification

`trigsimp()` simplifies trigonometric and hyperbolic expressions as
SymPy's `trigsimp()`. The default `method = "matching"` runs `futrig()`,
a greedy search over the Fu rules. `"fu"` runs
[`fu()`](https://adamrsophie.github.io/symcompR/reference/fu.md),
`"groebner"` reduces against the ideal of trigonometric identities with
`trigsimp_groebner()`, `"combined"` runs the Groebner method with
tangent hints and then `futrig()`, and `"old"` (or `old = TRUE`) runs
the pattern based algorithm of earlier SymPy versions. `inverse = TRUE`
also simplifies inverse trigonometric functions of their inverses.
`exptrigsimp()` rewrites exponentials into trigonometric or hyperbolic
functions where shorter. `futrig()` is the "matching" method itself.

## Usage

``` r
trigsimp_groebner(
  expr,
  hints = list(),
  quick = FALSE,
  order = "grlex",
  polynomial = FALSE
)

futrig(e, hyper = TRUE)

trigsimp(e, inverse = FALSE, method = "matching", old = FALSE, ...)

trigsimp_old(
  expr,
  first = TRUE,
  method = "matching",
  deep = FALSE,
  recursive = FALSE,
  ...
)

exptrigsimp(expr)
```

## Arguments

- hints, quick, order, polynomial:

  See `trigsimp_groebner()` in SymPy: generator hints (integers and
  function names), a quick search, the monomial order and whether to
  keep the result polynomial.

- e, expr:

  An expression, or a matrix of expressions.

- hyper:

  Also simplify hyperbolic functions through their trigonometric
  counterparts.

- inverse:

  Also simplify `asin(sin(x))` and the like.

- method:

  One of `"matching"`, `"fu"`, `"groebner"`, `"combined"` or `"old"`.

- old:

  Use the old method.

- ...:

  Options of the method: `measure` for `"fu"`; `hints`, `quick`, `order`
  and `polynomial` for `"groebner"`; `deep`, `recursive` and `first` for
  `"old"`.

- first, deep, recursive:

  Options of the old method: split the expression into separable factors
  first, simplify inside every argument, and simplify common
  subexpressions one at a time.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
trigsimp(sin(x)^4 - cos(x)^4)
#> -cos(2*x)
trigsimp(sin(x) / cos(x))
#> tan(x)
trigsimp(2 * sin(x)^2 + 2 * cos(x)^2, method = "groebner")
#> 2
```
