# Kane's method

`KanesMethod()` forms the equations of motion by Kane's method, as
SymPy's `KanesMethod`. Call `kanes_equations(bodies, loads)` (reached
with `$`) to form Fr and Fr\*, then read `mass_matrix`, `forcing`,
`mass_matrix_full`, `forcing_full`,
[`rhs()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
`kindiffdict()` or `linearize()`.

## Usage

``` r
KanesMethod(
  frame,
  q_ind,
  u_ind,
  kd_eqs = NULL,
  q_dependent = NULL,
  configuration_constraints = NULL,
  u_dependent = NULL,
  velocity_constraints = NULL,
  acceleration_constraints = NULL,
  u_auxiliary = NULL,
  bodies = NULL,
  forcelist = NULL,
  explicit_kinematics = TRUE,
  kd_eqs_solver = "LU",
  constraint_solver = "LU",
  nonholonomic_constraints = NULL
)
```

## Arguments

- frame:

  The inertial
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md).

- q_ind, u_ind:

  Independent generalized coordinates and speeds.

- kd_eqs:

  Kinematic differential equations.

- q_dependent, u_dependent:

  Dependent coordinates and speeds.

- configuration_constraints, velocity_constraints,
  acceleration_constraints, nonholonomic_constraints:

  Constraint equations.

- u_auxiliary:

  Auxiliary speeds.

- bodies, forcelist:

  Bodies and loads.

- explicit_kinematics:

  Solve the kinematic equations explicitly.

- kd_eqs_solver, constraint_solver:

  Linear solvers: `"LU"` or a function of `A` and `b`.

## Value

A `KanesMethod` object.

## Examples

``` r
q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
N <- ReferenceFrame("N")
P <- PhysPoint("P")
P$set_vel(N, u * N$x)
km <- KanesMethod(N, list(q), list(u), kd_eqs = list(dynamicsymbols("q", 1) - u))
invisible(km$kanes_equations(list(Particle("p", P, Symbol("m"))),
                             list(Force(P, -Symbol("k") * q * N$x))))
km$mass_matrix
#> Matrix(1, 1)
#> [m]
km$forcing
#> Matrix(1, 1)
#> [-k*q(t)]
```
