# Simplification helpers

`collect_abs()` joins products of absolute values, `|x|*|y| = |x*y|`.
`inversecombine()` simplifies functions of their inverses such as
`log(exp(x))` and `asin(sin(x))`. `ratsimp()` writes a rational function
as polynomial plus a proper part. `sum_simplify()` and
`product_simplify()` combine sums and products with the same limits, as
the SymPy functions of the same names.

## Usage

``` r
collect_abs(expr)

inversecombine(expr)

ratsimp(expr)

sum_simplify(e, kw = .simp_kw())

product_simplify(e, kw = .simp_kw())
```

## Arguments

- e, expr:

  An expression.

- kw:

  Options of
  [`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
  as a list (internal).

## Value

An expression.

## Examples

``` r
x <- Symbol("x", real = TRUE); y <- Symbol("y", real = TRUE)
collect_abs(abs(x) * abs(y))
#> Abs(x*y)
inversecombine(log(exp(x)))
#> x
```
