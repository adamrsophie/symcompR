# First order partial differential equations

`classify_pde(eq, func)` lists the solution methods that apply to `eq`,
a PDE in an unknown function of two variables such as `f(x, y)`.
`pdsolve(eq, func)` solves first order linear equations with constant or
variable coefficients; the arbitrary function is `F` unless `solvefun`
is given. `checkpdesol(eq, sol)` returns `list(ok, residual)`.
`pde_separate(eq, fun, sep, strategy)` separates variables with the
product (`"mul"`) or sum (`"add"`) ansatz `sep`, returning the two
separated sides or `NULL`.

## Usage

``` r
classify_pde(eq, func = NULL, dict = FALSE)

pdsolve(eq, func = NULL, hint = "default", solvefun = NULL)

checkpdesol(eq, sol, func = NULL, solve_for_func = TRUE)

pde_separate(eq, fun, sep, strategy = "mul")

pde_separate_add(eq, fun, sep)

pde_separate_mul(eq, fun, sep)
```

## Arguments

- eq:

  A PDE, as an equation or an expression equal to zero.

- func:

  The unknown function application, found when `NULL`.

- dict:

  Return a named list with the order, the default hint and the matches.

- hint:

  A method name or `"default"`.

- solvefun:

  The arbitrary function.

- sol:

  A solution, as an equation or the right side.

- solve_for_func:

  Solve `sol` for `func` first when needed.

- fun:

  The unknown function application, such as `u(x, t)`.

- sep:

  A list of function applications, such as `list(X(x), T(t))`.

- strategy:

  `"mul"` or `"add"`.

## Value

A character vector, a list or an equation.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); f <- Function("f")
u <- f(x, y)
pdsolve(diff(u, x) + diff(u, y) + u)
#> Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
classify_pde(x*diff(u, x) + diff(u, y))
#> [1] "1st_linear_variable_coeff"
```
