# Linear programming

`lpmin(f, constr)` and `lpmax(f, constr)` minimize and maximize the
linear function `f` subject to the linear constraints `constr` (a list
of non-strict inequalities and equalities) by the simplex method, as
SymPy's functions of the same names; they return `list(optimum, argopt)`
with the values of the symbols in a named list. Infeasible problems
signal an `InfeasibleLPError` and unbounded ones an `UnboundedLPError`
condition. `linprog(c, A, b, A_eq, b_eq, bounds)` minimizes `c x`
subject to `A x <= b`, `A_eq x = b_eq` and the bounds (by default
`x >= 0`) in matrix form, and `show_linprog()` returns the objective and
constraints that `linprog()` would solve.

## Usage

``` r
lpmin(f, constr)

lpmax(f, constr)

linprog(c, A = NULL, b = NULL, A_eq = NULL, b_eq = NULL, bounds = NULL)

show_linprog(c, A = NULL, b = NULL, A_eq = NULL, b_eq = NULL, bounds = NULL)
```

## Arguments

- f:

  The objective, linear in its symbols.

- constr:

  A list of constraints.

- c:

  Objective coefficients.

- A, b:

  Inequality constraints `A x <= b` (lists of rows, vectors or
  matrices).

- A_eq, b_eq:

  Equality constraints.

- bounds:

  `NULL`, a pair `list(lo, hi)` for all variables, a list of pairs, or a
  named list of pairs keyed by 0 based index; `NULL` in a pair means no
  bound.

## Value

A list of the optimum and the solution.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
lpmax(x - y, list(Le(x, y + 2), Ge(x, 0), Ge(y, 0)))
#> [[1]]
#> 2
#> 
#> [[2]]
#> [[2]]$x
#> 2
#> 
#> [[2]]$y
#> 0
#> 
#> 
lpmin(x, list(Ge(y, 1), Ge(x, y)))
#> [[1]]
#> 1
#> 
#> [[2]]
#> [[2]]$x
#> 1
#> 
#> [[2]]$y
#> 1
#> 
#> 
linprog(c(1, -1), list(c(1, 1)), 2, bounds = list(`1` = list(NULL, NULL)))
#> [[1]]
#> -2
#> 
#> [[2]]
#> [[2]][[1]]
#> 0
#> 
#> [[2]][[2]]
#> 2
#> 
#> 
```
