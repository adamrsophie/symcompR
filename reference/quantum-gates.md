# Quantum gates

The gates of `sympy.physics.quantum.gate` acting on qubits: the one
qubit gates `IdentityGate()`, `HadamardGate()`, `XGate()`, `YGate()`,
`ZGate()`, `PhaseGate()` (SymPy's `S`) and `TGate()` take a target
qubit; `CNotGate(control, target)` (also `CNOT()`), `SwapGate(i, j)`
(also `SWAP()`), `CPHASE(a, b)`, the controlled gate
`CGate(controls, gate)`, `CGateS()` and the general
`UGate(targets, matrix)` and `Gate()`. Products of gates applied to
[`Qubit()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
states are evaluated by
[`qapply()`](https://adamrsophie.github.io/symcompR/reference/qapply.md);
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
with `nqubits` gives the matrix. `gate_sort()` and `gate_simp()` reorder
and simplify circuits and `random_circuit()` builds a random one.
`represent_zbasis()` is the matrix of a controlled gate.

## Usage

``` r
Gate(...)

CGate(controls, gate)

CGateS(controls, gate)

UGate(targets, mat)

OneQubitGate(...)

TwoQubitGate(...)

IdentityGate(...)

HadamardGate(...)

XGate(...)

YGate(...)

ZGate(...)

PhaseGate(...)

TGate(...)

CNotGate(...)

CNOT(...)

SwapGate(...)

SWAP(...)

CPHASE(a, b)

represent_zbasis(controls, targets, target_matrix, nqubits, format = "sympy")

gate_simp(circuit)

gate_sort(circuit)

random_circuit(
  ngates,
  nqubits,
  gate_space = c("XGate", "YGate", "ZGate", "PhaseGate", "TGate", "HadamardGate",
    "CNotGate", "SwapGate")
)
```

## Arguments

- ...:

  Target qubits, or the arguments of the gate.

- controls, targets:

  Control and target qubits.

- gate:

  A gate.

- mat:

  A matrix.

- a, b:

  Qubits.

- target_matrix:

  The matrix of the gate on its targets.

- format:

  The matrix format, only `"sympy"`.

- circuit:

  A product of gates.

- ngates, nqubits:

  The number of gates and qubits.

- gate_space:

  The gate classes to draw from.

## Value

A gate, a circuit or a matrix.

## Examples

``` r
qapply(HadamardGate(0) * Qubit("0"))
#> sqrt(2)*|0>/2 + sqrt(2)*|1>/2
represent(CNOT(1, 0), nqubits = 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0]
#> [0, 1, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
gate_simp(HadamardGate(0) * HadamardGate(0))
#> 1
```
