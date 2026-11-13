# Measurement gates and gate factories

`Mz(target)` and `Mx(target)` are the mock measurement gates of
`sympy.physics.quantum.circuitplot`. `CreateOneQubitGate(name)` makes a
new one-qubit gate class and returns its constructor, and
`CreateCGate(name)` returns a function `(controls, target)` building a
controlled gate of a new one-qubit gate. `labeller(n)` gives wire
labels. Circuit plotting itself needs matplotlib and is not available.

## Usage

``` r
Mz(...)

Mx(...)

CreateOneQubitGate(name, latexname = NULL)

CreateCGate(name, latexname = NULL)

labeller(n, symbol = "q")
```

## Arguments

- ...:

  The target qubit.

- name, latexname:

  The gate name and its LaTeX form.

- n:

  The number of qubits.

- symbol:

  The label prefix.

## Value

A gate, a gate constructor or a character vector.

## Examples

``` r
Mz(0)
#> Mz(0)
Q <- CreateOneQubitGate("Q")
Q(1)
#> Q(1)
labeller(3)
#> [1] "q_2" "q_1" "q_0"
```
