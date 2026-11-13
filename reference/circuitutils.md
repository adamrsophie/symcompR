# Circuit utilities

Tools of `sympy.physics.quantum.circuitutils` for circuits given as
products of gates or lists of gates: `kmp_table()` is the
Knuth-Morris-Pratt table of a list, `find_subcircuit()` finds the
0-based position of a subcircuit (or -1), `replace_subcircuit()`
replaces it, `convert_to_symbolic_indices()` and
`convert_to_real_indices()` switch between numeric and symbolic qubit
labels, and `random_reduce()` and `random_insert()` remove or insert
gate identities. A `seed` given as a list of integers is used as the
sequence of random choices.

## Usage

``` r
kmp_table(word)

find_subcircuit(circuit, subcircuit, start = 0, end = 0)

replace_subcircuit(circuit, subcircuit, replace = NULL, pos = 0)

convert_to_symbolic_indices(seq, start = NULL, gen = NULL, qubit_map = NULL)

convert_to_real_indices(seq, qubit_map)

random_reduce(circuit, gate_ids, seed = NULL)

random_insert(circuit, choices, seed = NULL)
```

## Arguments

- word:

  A list.

- circuit, subcircuit, replace:

  Circuits.

- start, end, pos:

  Positions (0-based).

- seq:

  A circuit.

- gen:

  A list of symbols to draw new indices from.

- qubit_map:

  A dictionary of symbols and qubit indices.

- gate_ids:

  A list of
  [`GateIdentity()`](https://adamrsophie.github.io/symcompR/reference/identitysearch.md)
  objects.

- seed:

  `NULL`, an integer or a list of integers.

- choices:

  A list of circuits.

## Value

A list, a position or a circuit.

## Examples

``` r
x <- XGate(0); y <- YGate(0)
find_subcircuit(x * y * x, y * x)
#> [1] 1
replace_subcircuit(x * y * x, y)
#> [[1]]
#> X(0)
#> 
#> [[2]]
#> X(0)
#> 
```
