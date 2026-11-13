# Gate identities

Tools of `sympy.physics.quantum.identitysearch`.
`is_scalar_nonsparse_matrix()` (also `is_scalar_matrix()`) tests whether
a circuit is a scalar multiple of the identity, `generate_gate_rules()`
and `generate_equivalent_ids()` give the equivalent forms of a gate
identity, `GateIdentity(...)` holds an identity and
`bfs_identity_search()` and `random_identity_search()` look for them.
The scipy sparse version `is_scalar_sparse_matrix()` uses the same exact
computation here.

## Usage

``` r
is_scalar_nonsparse_matrix(circuit, nqubits, identity_only, eps = NULL)

is_scalar_sparse_matrix(circuit, nqubits, identity_only, eps = 1e-11)

is_scalar_matrix(circuit, nqubits, identity_only, eps = NULL)

generate_gate_rules(gate_seq, return_as_muls = FALSE)

generate_equivalent_ids(gate_seq, return_as_muls = FALSE)

GateIdentity(...)

is_degenerate(identity_set, gate_identity)

is_reducible(circuit, nqubits, begin, end)

bfs_identity_search(
  gate_list,
  nqubits,
  max_depth = NULL,
  identity_only = FALSE
)

random_identity_search(gate_list, numgates, nqubits)
```

## Arguments

- circuit, gate_seq:

  A circuit (a product or list of gates).

- nqubits:

  The number of qubits.

- identity_only:

  Accept only the identity, not its multiples.

- eps:

  Unused tolerance of the sparse version.

- return_as_muls:

  Return products instead of lists.

- ...:

  Gates.

- identity_set:

  A list of `GateIdentity` objects.

- gate_identity:

  A list of gates.

- begin, end:

  The range of gates (0-based, end excluded).

- gate_list:

  A list of gates.

- max_depth:

  The maximum circuit length.

- numgates:

  The number of gates.

## Value

A logical, a list (a set) or a `GateIdentity`.

## Examples

``` r
x <- XGate(0)
is_scalar_nonsparse_matrix(list(x, x), 1, TRUE)
#> [1] TRUE
generate_equivalent_ids(list(x, YGate(0), ZGate(0)))
#> [[1]]
#> [[1]][[1]]
#> X(0)
#> 
#> [[1]][[2]]
#> Y(0)
#> 
#> [[1]][[3]]
#> Z(0)
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> Y(0)
#> 
#> [[2]][[2]]
#> Z(0)
#> 
#> [[2]][[3]]
#> X(0)
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> Z(0)
#> 
#> [[3]][[2]]
#> X(0)
#> 
#> [[3]][[3]]
#> Y(0)
#> 
#> 
#> [[4]]
#> [[4]][[1]]
#> Z(0)
#> 
#> [[4]][[2]]
#> Y(0)
#> 
#> [[4]][[3]]
#> X(0)
#> 
#> 
#> [[5]]
#> [[5]][[1]]
#> Y(0)
#> 
#> [[5]][[2]]
#> X(0)
#> 
#> [[5]][[3]]
#> Z(0)
#> 
#> 
#> [[6]]
#> [[6]][[1]]
#> X(0)
#> 
#> [[6]][[2]]
#> Z(0)
#> 
#> [[6]][[3]]
#> Y(0)
#> 
#> 
#> attr(,"qset")
#> [1] TRUE
```
