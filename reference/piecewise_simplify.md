# Simplify a Piecewise expression

`piecewise_simplify(expr)` simplifies the expressions and conditions of
a Piecewise, rewrites univariate conditions as disjoint intervals, uses
equality conditions to simplify their expressions, and merges segments
that agree, as SymPy's `piecewise_simplify()`. `simplify_univariate()`
writes a boolean in one variable as a union of intervals.
`distribute_or_over_and()` puts a boolean into disjunctive form.

## Usage

``` r
distribute_or_over_and(expr)

piecewise_simplify(expr, deep = TRUE, kw = .simp_kw(), ...)

simplify_univariate(expr)
```

## Arguments

- expr:

  A Piecewise or boolean expression.

- deep:

  Unused; kept for SymPy's signature.

- kw:

  Options of
  [`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
  (internal).

- ...:

  Ignored.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
piecewise_simplify(Piecewise(list(x^2, x < 0), list(x^2, x > 0), list(0, TRUE)))
#> Piecewise((x^2, (x > 0) | (x < 0)), (0, True))
simplify_univariate((x > 0) & (x < 2) | (x > 1))
#> (x >= 1) | ((1 <= x) & (x <= 2)) | ((x <= 2) & (0 < x))
```
