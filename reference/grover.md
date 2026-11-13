# Grover's search

`superposition_basis(nqubits)` is the equal superposition of all basis
states, `OracleGate(nqubits, f)` flips the sign of the states for which
the R function `f` is `TRUE`, `WGate(nqubits)` is the Grover diffusion
operator, and `grover_iteration()` and `apply_grover()` run the
algorithm, as in `sympy.physics.quantum.grover`.

## Usage

``` r
superposition_basis(nqubits)

OracleGate(...)

WGate(...)

grover_iteration(qstate, oracle)

apply_grover(oracle, nqubits, iterations = NULL)
```

## Arguments

- nqubits:

  The number of qubits.

- ...:

  For `OracleGate()`, the number of qubits and the search function; for
  `WGate()`, the number of qubits.

- qstate:

  A state.

- oracle:

  An `OracleGate()` or, for `apply_grover()`, a search function.

- iterations:

  The number of iterations.

## Value

A state or a gate.

## Examples

``` r
f <- function(q) isTRUE(q == IntQubit(2, nqubits = 2))
apply_grover(f, 2)
#> |2>
```
