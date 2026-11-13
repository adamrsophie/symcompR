# Solve a univariate inequality

`solve_univariate_inequality(expr, x)` solves a real inequality in one
variable, as SymPy's function of the same name: the critical points come
from
[`solvify()`](https://adamrsophie.github.io/symcompR/reference/solvify.md),
the singularities and the discontinuities of the domain, and a test
point decides each interval between them.

## Usage

``` r
solve_univariate_inequality(
  expr,
  x,
  relational = TRUE,
  domain = Reals,
  continuous = FALSE
)
```

## Arguments

- expr:

  A relational such as `x^2 > 4`.

- x:

  The variable.

- relational:

  Return a relational (default) instead of a set.

- domain:

  A subset of the real line, `Reals` by default.

- continuous:

  Whether `expr` is known to be continuous on `domain`.

## Value

A boolean expression, or a set when `relational = FALSE`.

## Examples

``` r
x <- Symbol("x")
solve_univariate_inequality(x^2 > 4, x)
#> ((-oo < x) & (x < -2)) | ((2 < x) & (x < oo))
solve_univariate_inequality((x - 1) / (x + 2) <= 0, x, relational = FALSE)
#> Interval.Lopen(-2, 1)
```
