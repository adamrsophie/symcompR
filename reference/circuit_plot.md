# Plot quantum circuits

`circuit_plot(c, nqubits)` (also `CircuitPlot()`) draws a circuit, a
gate or a product of gates, on `nqubits` wires: one qubit gates as
boxes, controlled gates as a control line with control points (with
`simplify_cgate`, a controlled `X` is drawn as a not point and a
controlled `Z` as a control point), swaps as crosses and measured wires
doubled after the measurement. The gates act from right to left, so the
rightmost factor of a product is drawn first. Options are SymPy's
attributes: `scale`, `fontsize` (here a `cex` multiplier, default 1),
`linewidth`, `control_radius`, `not_radius`, `swap_delta`, `labels`
(wire labels, see
[`labeller()`](https://adamrsophie.github.io/symcompR/reference/circuit-gates.md)),
`inits` (a named list of initial states) and `label_buffer`.

## Usage

``` r
circuit_plot(c, nqubits, ...)

CircuitPlot(c, nqubits, ...)
```

## Arguments

- c:

  A gate or a product of gates.

- nqubits:

  The number of wires.

- ...:

  Options, see Details.

## Value

Invisibly, a `CircuitPlot` object with the grids and gates.

## Examples

``` r
circuit_plot(CNOT(1, 0) * HadamardGate(1), 2, labels = labeller(2))

circuit_plot(SwapGate(0, 1) * Mz(0), 2)
```
