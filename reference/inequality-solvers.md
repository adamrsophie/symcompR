# Inequalities in one real variable

These follow SymPy's `sympy.solvers.inequalities`.
`solve_poly_inequality(poly, rel)` solves `poly rel 0` for a
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md) and
returns a list of intervals. `solve_poly_inequalities(polys)` joins the
solutions of a list of `list(poly, rel)` pairs.
`solve_rational_inequalities(eqs)` solves `eqs`, a list of lists whose
entries `list(list(numer, denom), rel)` mean `numer/denom rel 0`; inner
lists are joined by `And`, the outer list by `Or`.
`reduce_rational_inequalities(exprs, gen)` does the same for relationals
or `list(expr, rel)` pairs. `reduce_abs_inequality(expr, rel, gen)` and
`reduce_abs_inequalities(exprs, gen)` handle nested absolute values, and
`reduce_inequalities(inequalities, symbols)` reduces a system of
inequalities.

## Usage

``` r
solve_poly_inequality(poly, rel)

solve_poly_inequalities(polys)

solve_rational_inequalities(eqs)

reduce_rational_inequalities(exprs, gen, relational = TRUE)

reduce_abs_inequality(expr, rel, gen)

reduce_abs_inequalities(exprs, gen)

reduce_inequalities(inequalities, symbols = list())
```

## Arguments

- poly:

  A `Poly` (or a polynomial expression).

- rel:

  A relation: `"<"`, `"<="`, `">"`, `">="`, `"=="` or `"!="`.

- polys:

  A list of `list(poly, rel)` pairs.

- eqs:

  A nested list of inequalities.

- exprs:

  A nested list (rational) or a list of `list(expr, rel)` pairs
  (absolute values).

- gen:

  The real variable.

- relational:

  Return a relational instead of a set.

- expr:

  An expression.

- inequalities:

  A relational or a list of them.

- symbols:

  The symbols to solve for.

## Value

A list of sets, a set or a relational.

## Examples

``` r
x <- Symbol("x")
solve_poly_inequality(Poly(x^2 - 1, x), ">")
#> [[1]]
#> Interval.open(-oo, -1)
#> 
#> [[2]]
#> Interval.open(1, oo)
#> 
reduce_rational_inequalities(list(list(x^2 < 4)), x)
#> (-2 < x) & (x < 2)
reduce_abs_inequality(Abs(x - 5) - 3, "<", x)
#> (2 < x) & (x < 8)
reduce_inequalities(list(x^2 < 4, x > 0), x)
#> (0 < x) & (x < 2)
```
