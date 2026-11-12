# Simplify an expression

Runs SymPy's `simplify()` heuristics: rational simplification,
`powsimp`, cancellation and `together`, `factor_terms`, then the
specialized simplifications (Piecewise, hypergeometric, Bessel,
trigonometric, logarithmic, combinatorial, sums, integrals and products)
and keeps the shortest result under `measure`, returning the input when
the result is more than `ratio` times as complex.

## Usage

``` r
simplify(
  e,
  ratio = 1.7,
  measure = count_ops,
  rational = FALSE,
  inverse = FALSE,
  doit = TRUE,
  ...
)
```

## Arguments

- e:

  An expression, matrix or list.

- ratio:

  Accept a result only if it is not more than `ratio` times as complex
  as the input (`Inf` keeps the shortest).

- measure:

  The complexity measure,
  [`count_ops()`](https://adamrsophie.github.io/symcompR/reference/count_ops.md)
  by default.

- rational:

  `FALSE` keeps Floats, `TRUE` converts them to rationals, `NULL`
  converts them and back at the end.

- inverse:

  Also combine functions with their inverses
  ([`inversecombine()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)).

- doit:

  Evaluate unevaluated objects at the end.

- ...:

  Ignored.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
simplify(sin(x)^2 + cos(x)^2)
#> 1
simplify((x^3 + x^2 - x - 1) / (x^2 + 2*x + 1))
#> x - 1
```
