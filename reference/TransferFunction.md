# Transfer functions

`TransferFunction(num, den, var)` is the continuous transfer function
`num/den` in the complex variable `var` of the Laplace transform, and
`DiscreteTransferFunction(num, den, var, sampling_time)` the discrete
one in the variable of the z-transform.
`PIDController(kp, ki, kd, tf, var)` is the transfer function of a PID
controller with a derivative filter time constant `tf`. The class
methods `TransferFunction$from_rational_expression(expr, var)`,
`$from_coeff_lists(num_list, den_list, var)` and
`$from_zpk(zeros, poles, gain, var)` build transfer functions; the
discrete class also has `$from_gbt()`, `$from_bilinear()`,
`$from_forward_diff()` and `$from_backward_diff()`. A transfer function
has the properties `$num`, `$den`, `$var`, `$sampling_time`,
`$is_proper`, `$is_strictly_proper` and `$is_biproper` and the methods
`$poles()`, `$zeros()`, `$dc_gain()`, `$eval_frequency(other)`,
`$is_stable()`, `$to_expr()`, `$to_standard_form()`, `$simplify()`,
`$expand()` and `$get_asymptotic_stability_conditions()`. The operators
`+`, `-`, `*`, `/` and `^` build
[`Series()`](https://adamrsophie.github.io/symcompR/reference/control-interconnections.md),
[`Parallel()`](https://adamrsophie.github.io/symcompR/reference/control-interconnections.md)
and
[`Feedback()`](https://adamrsophie.github.io/symcompR/reference/control-interconnections.md)
interconnections. After `sympy.physics.control.lti`.

## Usage

``` r
TransferFunction(num, den, var)

DiscreteTransferFunction(num, den, var, sampling_time = 1)

TransferFunctionBase

PIDController(
  kp = Symbol("kp"),
  ki = Symbol("ki"),
  kd = Symbol("kd"),
  tf = 0L,
  var = Symbol("s")
)
```

## Arguments

- num, den:

  The numerator and denominator.

- var:

  The complex variable.

- sampling_time:

  The sampling time.

- kp, ki, kd:

  The proportional, integral and derivative gains.

- tf:

  The derivative filter time constant.

## Value

A `TransferFunction`, `DiscreteTransferFunction` or `PIDController`
object.

## Examples

``` r
s <- Symbol("s")
G <- TransferFunction(s + 1, s^2 + 2, s)
G$poles()
#> [[1]]
#> sqrt(2)*I
#> 
#> [[2]]
#> -sqrt(2)*I
#> 
G * TransferFunction(1, s, s)
#> Series(TransferFunction(s + 1, s^2 + 2, s), TransferFunction(1, s, s))
TransferFunction$from_rational_expression((s + 1) / (s^2 - 1), s)
#> TransferFunction(s + 1, s^2 - 1, s)
```
