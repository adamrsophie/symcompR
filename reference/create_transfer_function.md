# Create a transfer function

`create_transfer_function(num, den, var, sampling_time)` is a continuous
[`TransferFunction()`](https://adamrsophie.github.io/symcompR/reference/TransferFunction.md)
when `sampling_time` is 0 and a
[`DiscreteTransferFunction()`](https://adamrsophie.github.io/symcompR/reference/TransferFunction.md)
otherwise; `create_state_space()` does the same for state space models.
After `sympy.physics.control.lti`.

## Usage

``` r
create_transfer_function(num, den, var, sampling_time = 0L)

create_state_space(A, B, C, D, sampling_time = 0L)
```

## Arguments

- num, den:

  The numerator and denominator.

- var:

  The complex variable of the Laplace transform or z-transform.

- sampling_time:

  The sampling time; 0 for a continuous system.

- A, B, C, D:

  The state, input, output and feedforward matrices.

## Value

A transfer function or a state space model.

## Examples

``` r
s <- Symbol("s"); z <- Symbol("z")
create_transfer_function(s + 1, s^2 + 2, s)
#> TransferFunction(s + 1, s^2 + 2, s)
create_transfer_function(z, z + 1, z, 0.1)
#> DiscreteTransferFunction(z, z + 1, z, 0.1)
```
