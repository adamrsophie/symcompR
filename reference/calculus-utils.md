# Calculus utilities

`singularities(f, x, domain)` returns the points of `domain` where `f`
is undefined (zeros of denominators and of logarithm arguments); the
default domain is `Reals` for a real `x` and `Complexes` otherwise.
`stationary_points()` solves `f'(x) = 0` on a domain; `minimum()` and
`maximum()` return the extreme values on a closed interval or the real
line, including one-sided limits at poles. `is_increasing()`,
`is_strictly_increasing()`, `is_decreasing()`,
`is_strictly_decreasing()`, `is_monotonic()` and `is_convex()` test the
sign of derivatives on an interval. `periodicity(f, x)` returns the
period of trigonometric expressions. `euler_equations(L, funcs, vars)`
gives the list of Euler-Lagrange equations of a Lagrangian, for
derivatives of any order. `finite_diff_weights(order, x_list, x0)`
returns Fornberg weights `w[[m + 1]][[k]][[j]]`: derivative order `m`,
first `k` grid points, point `j`. `apply_finite_diff()` applies the
full-grid weights.

## Usage

``` r
singularities(f, x, domain = NULL)

stationary_points(f, x, domain = Reals)

minimum(f, x, domain = Reals)

maximum(f, x, domain = Reals)

is_increasing(f, interval = Reals, x = NULL)

is_strictly_increasing(f, interval = Reals, x = NULL)

is_decreasing(f, interval = Reals, x = NULL)

is_strictly_decreasing(f, interval = Reals, x = NULL)

is_monotonic(f, interval = Reals, x = NULL)

is_convex(f, x = NULL, domain = Reals)

periodicity(f, x)

euler_equations(L, funcs = NULL, vars = NULL)

finite_diff_weights(order, x_list, x0 = 1)

apply_finite_diff(order, x_list, y_list, x0 = 0)
```

## Arguments

- f, L:

  Expressions.

- x:

  The variable.

- domain:

  An interval, `Reals` or `Complexes`; for `singularities()`, `NULL`
  picks from the assumptions on `x`.

- interval:

  An interval, or `Reals`.

- funcs:

  Function applications such as `u(x)`, a list of them, or `NULL` for
  all undefined functions in `L`.

- vars:

  Independent variables, or `NULL` for the arguments of the first
  function.

- order:

  Derivative order.

- x_list, y_list:

  Grid points and values.

- x0:

  Evaluation point.

## Value

A set, expression, logical or list.

## Examples

``` r
x <- Symbol("x")
singularities(1 / (x^2 - 1), x)
#> {-1, 1}
minimum(x^2 - 2*x, x)
#> -1
is_increasing(exp(x), Reals, x)
#> [1] TRUE
finite_diff_weights(1, list(-1, 0, 1), 0)[[2]][[3]]
#> [[1]]
#> -1/2
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> 1/2
#> 
```
