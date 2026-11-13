# Differential geometry

`Manifold(name, dim)`, `Patch(name, manifold)` and
`CoordSystem(name, patch, symbols, relations)` build coordinate systems;
`relations` is a list of `list(c(from, to), syms, exprs)` entries.
`transformation()`, `cs_transform()` (SymPy's `CoordSystem.transform`),
`jacobian_dg()` (SymPy's `CoordSystem.jacobian`),
`jacobian_determinant()`, `point()`, `point_to_coords()`,
`base_scalar()`, `base_scalars()`, `base_vector()`, `base_vectors()`,
`base_oneform()` and `base_oneforms()` are SymPy's methods.
`Differential()`, `TensorProduct()`, `WedgeProduct()`, `Commutator()`,
`LieDerivative()`, `BaseCovarDerivativeOp()` and `CovarDerivativeOp()`
build fields and operators, applied with `rcall(field, ...)`. When no
argument holds a field, `TensorProduct()` and `Commutator()` are those
of `sympy.physics.quantum` (see
[`Operator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md)).
`contravariant_order()`, `covariant_order()`, `vectors_in_basis()`,
`twoform_to_matrix()`, `metric_to_Christoffel_1st()`,
`metric_to_Christoffel_2nd()`, `metric_to_Riemann_components()`,
`metric_to_Ricci_components()`, `intcurve_series()` and
`intcurve_diffequ()` are SymPy's functions. `R2`, `R2_r`, `R2_p`, `R3`,
`R3_r`, `R3_c` and `R3_s` (from `sympy.diffgeom.rn`) are built with
`diffgeom_rn()`.

## Usage

``` r
Manifold(name, dim)

Patch(name, manifold)

CoordSystem(name, patch, symbols = NULL, relations = list())

coord_symbols(cs)

transformation(cs, sys)

cs_transform(cs, sys, coordinates = NULL)

jacobian_dg(cs, sys, coordinates = NULL)

jacobian_determinant(cs, sys, coordinates = NULL)

point(cs, coords)

point_to_coords(cs, pt)

base_scalar(cs, coord_index)

base_scalars(cs)

base_vector(cs, coord_index)

base_vectors(cs)

base_oneform(cs, coord_index)

base_oneforms(cs)

Differential(form_field)

TensorProduct(...)

WedgeProduct(...)

Commutator(v1, v2)

LieDerivative(v_field, expr)

BaseCovarDerivativeOp(coord_sys, index, christoffel)

CovarDerivativeOp(wrt, christoffel)

rcall(expr, ...)

contravariant_order(expr, .strict = FALSE)

covariant_order(expr, .strict = FALSE)

vectors_in_basis(expr, to_sys)

twoform_to_matrix(expr)

metric_to_Christoffel_1st(expr)

metric_to_Christoffel_2nd(expr)

metric_to_Riemann_components(expr)

metric_to_Ricci_components(expr)

intcurve_series(
  vector_field,
  param,
  start_point,
  n = 6,
  coord_sys = NULL,
  coeffs = FALSE
)

intcurve_diffequ(vector_field, param, start_point, coord_sys = NULL)

diffgeom_rn()
```

## Arguments

- name:

  A name.

- dim:

  The dimension.

- manifold, patch, cs, sys, coord_sys:

  A manifold, patch or coordinate system.

- symbols:

  Coordinate symbols.

- relations:

  Transformations between systems.

- coordinates, coords:

  Coordinates.

- pt:

  A point.

- coord_index:

  An index from 0.

- form_field, expr, v_field, v1, v2, wrt, vector_field:

  Fields.

- ...:

  Fields for `rcall()`, `TensorProduct()` and `WedgeProduct()`.

- index:

  An index.

- christoffel:

  Christoffel symbols, an array.

- .strict:

  Return -1 for expressions unrelated to the module.

- to_sys:

  A coordinate system.

- param:

  The curve parameter.

- start_point:

  A point.

- n:

  The number of terms.

- coeffs:

  Return the coefficients.

## Value

A sym object, matrix, array or list.

## Examples

``` r
rn <- diffgeom_rn()
x <- base_scalar(rn$R2_r, 0); y <- base_scalar(rn$R2_r, 1)
e_x <- base_vector(rn$R2_r, 0)
rcall(e_x, x^2 * y)
#> 2*x*y
```
