# State space models

`StateSpace(A, B, C, D)` is the continuous state space model
`x' = A x + B u`, `y = C x + D u` and
`DiscreteStateSpace(A, B, C, D, sampling_time)` the discrete one. The
matrices default to zero matrices of the matching shapes. A model has
the properties `$A`, `$B`, `$C`, `$D` (also `$state_matrix`,
`$input_matrix`, `$output_matrix`, `$feedforward_matrix`),
`$num_states`, `$num_inputs`, `$num_outputs` and `$shape`, and the
methods `$append()`, `$observability_matrix()`,
`$controllability_matrix()`, `$is_observable()`, `$is_controllable()`,
the observable and controllable subspaces and forms,
`$apply_similarity_transform()`, `$dsolve()`,
`$get_asymptotic_stability_conditions()` and
`$rewrite("TransferFunction")`. `+`, `-` and `*` connect models.
[`create_state_space()`](https://adamrsophie.github.io/symcompR/reference/create_transfer_function.md)
picks the class from the sampling time. After
`sympy.physics.control.lti`.

## Usage

``` r
StateSpace(A = NULL, B = NULL, C = NULL, D = NULL)

DiscreteStateSpace(A = NULL, B = NULL, C = NULL, D = NULL, sampling_time = 1)
```

## Arguments

- A, B, C, D:

  The state, input, output and feedforward matrices.

- sampling_time:

  The sampling time.

## Value

A `StateSpace` or `DiscreteStateSpace` object.

## Examples

``` r
A <- Matrix(list(list(-1, 0), list(0, -2)))
ss <- StateSpace(A, Matrix(list(list(1), list(1))), Matrix(list(list(1, 1))), Matrix(list(list(0))))
ss$num_states
#> [1] 2
ss$is_controllable()
#> [1] TRUE
ss$rewrite("TransferFunction")
#> [[1]]
#> [[1]][[1]]
#> TransferFunction(2*s + 3, (s + 1)*(s + 2), s)
#> 
#> 
```
