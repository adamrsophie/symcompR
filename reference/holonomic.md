# Holonomic functions

`DifferentialOperators(base, generator)` returns the algebra and `Dx`;
the base is `holonomic_ring(x)` (SymPy's `QQ.old_poly_ring(x)`).
Operators combine with `+`, `-`, `*` and `^`.
`HolonomicFunction(ann, x, x0, y0)` is a holonomic function; `+`, `*`,
`^`, `integrate_holonomic()`, `diff_holonomic()`, `composition()`,
`to_sequence()`, `series_holonomic()`, `to_hyper()`,
`to_expr_holonomic()`, `to_meijerg()`, `change_ics()`, `shift_x()`,
`change_x()`, `evalf_holonomic()` and `holonomic_indicial()` are SymPy's
methods. `expr_to_holonomic()`, `from_hyper()` and `from_meijerg()`
convert expressions. `RecurrenceOperators()` and `HolonomicSequence()`
are the recurrence counterparts.

## Usage

``` r
holonomic_ring(x, dom = "QQ")

DifferentialOperators(base, generator = "Dx")

HolonomicFunction(annihilator, x, x0 = 0, y0 = NULL)

integrate_holonomic(f, limits, initcond = FALSE)

diff_holonomic(f, x = NULL, n = 1L)

composition(f, expr, ...)

holonomic_degree(f)

RecurrenceOperators(base, generator = "Sn")

HolonomicSequence(recurrence, u0 = list())

to_sequence(f, lb = TRUE)

shift_x(f, a)

change_x(f, z)

series_holonomic(f, n = 6, coefficient = FALSE, order = TRUE)

holonomic_indicial(f)

to_hyper(f, as_list = FALSE)

to_expr_holonomic(f)

change_ics(f, b, lenics = NULL)

to_meijerg(f)

evalf_holonomic(f, points, method = "RK4", step = 0.05, derivatives = FALSE)

from_hyper(func, x0 = 0, evalf = FALSE)

from_meijerg(func, x0 = 0, evalf = FALSE, initcond = TRUE, domain = "QQ")

expr_to_holonomic(
  func,
  x = NULL,
  x0 = 0,
  y0 = NULL,
  lenics = NULL,
  domain = NULL,
  initcond = TRUE
)
```

## Arguments

- x, z, n:

  Symbols.

- dom:

  The coefficient domain name.

- base:

  A base ring from `holonomic_ring()`.

- generator:

  The generator name.

- annihilator:

  A differential operator.

- x0:

  The point of the initial conditions.

- y0:

  Initial conditions: a list, or a named list keyed by exponents for
  singular conditions.

- f:

  A holonomic function.

- limits:

  `x` or `list(x, a, b)`.

- initcond:

  Compute initial conditions.

- expr, func:

  An expression.

- ...:

  Initial point and conditions for `composition()`.

- recurrence:

  A recurrence operator.

- u0:

  Initial values.

- lb:

  Return the lower bound too.

- coefficient:

  Return coefficients.

- order:

  Add an Order term.

- as_list:

  Return a list of terms.

- b, a:

  Points.

- lenics:

  The number of initial conditions.

- points:

  Points for numerical evaluation.

- method:

  `"RK4"` or `"Euler"`.

- step:

  The step size.

- derivatives:

  Return derivatives too.

- evalf:

  Evaluate initial conditions numerically.

- domain:

  The coefficient domain.

## Value

An operator, holonomic function, sequence or expression.

## Examples

``` r
x <- Symbol("x")
R <- DifferentialOperators(holonomic_ring(x), "Dx"); Dx <- R$Dx
p <- HolonomicFunction(Dx - 1, x, 0, list(1))
q <- HolonomicFunction(Dx^2 + 1, x, 0, list(0, 1))
p + q
#> HolonomicFunction((-1) + (1)*Dx + (-1)*Dx^2 + (1)*Dx^3, x, 0, [1, 2, 1]) 
p * q
#> HolonomicFunction((2) + (-2)*Dx + (1)*Dx^2, x, 0, [0, 1]) 
to_expr_holonomic(q * q)
#> sin(x)^2
```
