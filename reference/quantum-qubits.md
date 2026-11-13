# Qubits

The qubit states of `sympy.physics.quantum.qubit`: `Qubit("0101")` or
`Qubit(0, 1, 0, 1)`, its bra `QubitBra()`, and `IntQubit(n)` and
`IntQubitBra(n)` labeled by an integer (`nqubits` fixes the width).
Gates act on them through
[`qapply()`](https://adamrsophie.github.io/symcompR/reference/qapply.md).
`qubit_to_matrix()` and `matrix_to_qubit()` convert to and from column
vectors, `matrix_to_density()` gives a
[`Density()`](https://adamrsophie.github.io/symcompR/reference/E.md) and
`measure_all()`, `measure_partial()`, `measure_all_oneshot()` and
`measure_partial_oneshot()` give measurement outcomes with their
probabilities (the one shot versions draw a random outcome).

## Usage

``` r
QubitState(...)

Qubit(...)

QubitBra(...)

IntQubit(..., nqubits = NULL)

IntQubitBra(..., nqubits = NULL)

matrix_to_qubit(matrix)

matrix_to_density(mat)

qubit_to_matrix(qubit, format = "sympy")

measure_all(qubit, format = "sympy", normalize = TRUE)

measure_partial(qubit, bits, format = "sympy", normalize = TRUE)

measure_partial_oneshot(qubit, bits, format = "sympy")

measure_all_oneshot(qubit, format = "sympy")
```

## Arguments

- ...:

  Qubit values, a string of zeros and ones, or a qubit.

- nqubits:

  The number of qubits of an integer qubit.

- matrix, mat:

  A matrix.

- qubit:

  A qubit expression.

- format:

  The matrix format, only `"sympy"`.

- normalize:

  Normalize the state first.

- bits:

  The qubits to measure.

## Value

A qubit expression, a matrix or a list of outcomes.

## Examples

``` r
q <- Qubit("01")
qubit_to_matrix(q)
#> Matrix(4, 1)
#> [0]
#> [1]
#> [0]
#> [0]
matrix_to_qubit(qubit_to_matrix(q))
#> |01>
measure_all(qapply(HadamardGate(0) * Qubit("0")))
#> [[1]]
#> [[1]][[1]]
#> |0>
#> 
#> [[1]][[2]]
#> 1/2
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> |1>
#> 
#> [[2]][[2]]
#> 1/2
#> 
#> 
```
