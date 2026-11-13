# Operators and their eigenstates

`operators_to_state(ops)` gives the eigenstate of an operator or of a
set (an R list) of commuting operators, and `state_to_operators(state)`
the operators of a state, as SymPy's
`sympy.physics.quantum.operatorset`. Classes are given by name, for
example `"XOp"` or `list("J2Op", "JzOp")`.

## Usage

``` r
operators_to_state(operators, ...)

state_to_operators(state, ...)
```

## Arguments

- operators:

  An operator, a class name or a list of either.

- ...:

  Options for the state or operator constructors.

- state:

  A state or a state class name.

## Value

A state or operator, a list of operators, a class name or `NULL`.

## Examples

``` r
operators_to_state(XOp())
#> |x>
state_to_operators(XKet())
#> X
```
