# Matrix representations

`represent(expr, basis = ...)` gives the matrix (or wavefunction) of a
quantum expression in a basis, as SymPy's
`sympy.physics.quantum.represent`. The basis is an operator, a state or
a class name such as `"JzOp"`; without one the default basis of the
object is used. `rep_innerproduct()` and `rep_expectation()` build the
position representation of states and operators, `get_basis()` finds the
basis state, `enumerate_states()` lists states with consecutive labels
and `integrate_result()` integrates over the unities inserted by
`represent()`.

## Usage

``` r
represent(expr, ...)

rep_innerproduct(expr, ...)

rep_expectation(expr, ...)

integrate_result(orig_expr, result, ...)

get_basis(expr, basis = NULL, replace_none = TRUE, ...)

enumerate_states(...)
```

## Arguments

- expr:

  A quantum expression.

- ...:

  Options such as `nqubits`, `j` or `format`.

- orig_expr, result:

  The represented expression and its matrix.

- basis:

  A basis operator, state or class name.

- replace_none:

  Use the default basis when `basis` is `NULL`.

## Value

A matrix, an expression or a list of states.

## Examples

``` r
represent(IdentityOperator(2))
#> Matrix(2, 2)
#> [1, 0]
#> [0, 1]
```
