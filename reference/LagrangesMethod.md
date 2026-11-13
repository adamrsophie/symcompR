# Lagrange's method

`LagrangesMethod()` forms the equations of motion from a Lagrangian, as
SymPy's `LagrangesMethod`. Call `form_lagranges_equations()` (reached
with `$`), then read `mass_matrix`, `forcing`, `mass_matrix_full`,
`forcing_full`,
[`rhs()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
`solve_multipliers()` or `linearize()`.

## Usage

``` r
LagrangesMethod(
  Lagrangian,
  qs,
  forcelist = NULL,
  bodies = NULL,
  frame = NULL,
  hol_coneqs = NULL,
  nonhol_coneqs = NULL
)
```

## Arguments

- Lagrangian:

  The Lagrangian, a scalar expression.

- qs:

  Generalized coordinates.

- forcelist:

  Loads.

- bodies:

  Bodies.

- frame:

  The inertial
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md),
  needed with loads.

- hol_coneqs, nonhol_coneqs:

  Holonomic and nonholonomic constraints.

## Value

A `LagrangesMethod` object.

## Examples

``` r
q <- dynamicsymbols("q")
qd <- dynamicsymbols("q", 1)
m <- Symbol("m"); k <- Symbol("k")
lm <- LagrangesMethod(m * qd^2 / 2 - k * q^2 / 2, list(q))
lm$form_lagranges_equations()
#> Matrix(1, 1)
#> [k*q(t) + m*Derivative(q(t), (t, 2))]
```
