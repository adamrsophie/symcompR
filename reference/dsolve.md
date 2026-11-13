# Solve ordinary differential equations

`dsolve(eq, func)` solves the ordinary differential equation `eq` (an
expression taken equal to zero, or an
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md)) for
the unknown function `func`, as SymPy's `dsolve()`. The hints of
[`classify_ode()`](https://adamrsophie.github.io/symcompR/reference/classify_ode.md)
are tried in SymPy's order; `hint` names one of them, or the meta hints
`"default"`, `"all"`, `"all_Integral"` and `"best"`. A list of equations
is solved as a system of ODEs.

## Usage

``` r
dsolve(
  eq,
  func = NULL,
  hint = "default",
  simplify = TRUE,
  ics = NULL,
  xi = NULL,
  eta = NULL,
  x0 = 0,
  n = 6,
  ...
)
```

## Arguments

- eq:

  The equation, or a list of equations for a system.

- func:

  The unknown function such as `f(x)`, or a list of functions for a
  system; detected from the equation when `NULL`.

- hint:

  The solution method.

- simplify:

  Whether to simplify the solution with `odesimp()`.

- ics:

  Initial conditions: a list of pairs `list(f(x0), value)` or
  `list(Subs(Derivative(f(x), x, n), x, x0), value)`; the triple form
  `list(x0, n, value)` is accepted too.

- xi, eta:

  Infinitesimals for the `lie_group` hint.

- x0:

  The point of expansion for the power series hints.

- n:

  The number of power series terms.

- ...:

  Further options passed on to the hints.

## Value

An equation, or a list of equations; a named list for the `"all"` hints.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
dsolve(Derivative(f(x), x, 2) - 3 * Derivative(f(x), x) + 2 * f(x), f(x))
#> Eq(f(x), (C1 + C2*exp(x))*exp(x))
dsolve(Derivative(f(x), x) - f(x), f(x), ics = list(list(f(0), 1)))
#> Eq(f(x), exp(x))
```
