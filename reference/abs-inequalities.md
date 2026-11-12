# Rational and absolute value inequalities

`solve_rational_inequalities(eqs)` solves `eqs`, a list of lists whose
inner entries `list(list(n, d), rel)` mean `n/d rel 0`; inner lists are
joined by `And` and the outer list by `Or`.
`reduce_abs_inequality(expr, rel, x)` reduces `expr rel 0` with nested
absolute values, and `reduce_abs_inequalities(exprs, x)` a list of
`list(expr, rel)` pairs.

## Usage

``` r
solve_rational_inequalities(eqs)

reduce_abs_inequality(expr, rel, x)

reduce_abs_inequalities(exprs, x)
```

## Arguments

- eqs:

  Nested list of inequalities.

- expr:

  An expression.

- rel:

  A relation: `"<"`, `"<="`, `">"`, `">="`, `"=="` or `"!="`.

- x:

  The real variable.

- exprs:

  A list of `list(expr, rel)` pairs.

## Value

A set, or a relational condition on `x`.

## Examples

``` r
x <- Symbol("x")
solve_rational_inequalities(list(list(list(list(x - 1, x + 1), "<"))))
#> Interval.open(-1, 1)
reduce_abs_inequality(Abs(x - 5) - 3, "<", x)
#> (2 < x) & (x < 8)
```
