# Symbolic equations of motion

`SymbolicSystem()` stores equations of motion in one of SymPy's three
forms (explicit, implicit with a combined mass matrix, or implicit with
separate kinematic and dynamic parts), as SymPy's `SymbolicSystem`.
Attributes such as `coordinates`, `speeds`, `states`,
`comb_explicit_rhs`, `comb_implicit_mat` and `comb_implicit_rhs`, and
the methods `compute_explicit_form()`, `dynamic_symbols()` and
`constant_symbols()`, are reached with `$`.

## Usage

``` r
SymbolicSystem(
  coord_states,
  right_hand_side,
  speeds = NULL,
  mass_matrix = NULL,
  coordinate_derivatives = NULL,
  alg_con = NULL,
  output_eqns = list(),
  coord_idxs = NULL,
  speed_idxs = NULL,
  bodies = NULL,
  loads = NULL
)
```

## Arguments

- coord_states:

  Coordinates, or all states when `speeds` is `NULL`.

- right_hand_side:

  The right hand side of the equations.

- speeds:

  Generalized speeds.

- mass_matrix:

  Mass matrix of the implicit forms.

- coordinate_derivatives:

  Explicit kinematic equations.

- alg_con:

  Indices of algebraic constraint equations.

- output_eqns:

  A named list of output equations.

- coord_idxs, speed_idxs:

  Zero based indices of the coordinates and speeds in `coord_states`.

- bodies, loads:

  Bodies and loads of the system.

## Value

A `SymbolicSystem` object.

## Examples

``` r
x <- dynamicsymbols("x"); v <- dynamicsymbols("v")
s <- SymbolicSystem(list(x, v), Matrix(list(v, -x)))
s$states
#> Matrix(2, 1)
#> [x(t)]
#> [v(t)]
s$comb_explicit_rhs
#> Matrix(2, 1)
#> [v(t) ]
#> [-x(t)]
```
