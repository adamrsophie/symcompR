# Numerical data and plots of SISO systems

`pole_zero_numerical_data(system)` gives the zeros and the poles of a
SISO system. `impulse_response_numerical_data()`,
`step_response_numerical_data()` and `ramp_response_numerical_data()`
sample the impulse, unit step and ramp responses on `n` equally spaced
times; `bode_magnitude_numerical_data()` and
`bode_phase_numerical_data()` sample the magnitude (dB) and the phase on
`n` logarithmically spaced frequencies. `nyquist_plot_expr()` and
`nichols_plot_expr()` give the expressions of the Nyquist and Nichols
plots. The `*_plot()` functions draw the plots with base R graphics.
After `sympy.physics.control.control_plots`.

## Usage

``` r
pole_zero_numerical_data(system)

impulse_response_numerical_data(
  system,
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  n = 1000L,
  adaptive = FALSE,
  ...
)

step_response_numerical_data(
  system,
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  n = 1000L,
  adaptive = FALSE,
  ...
)

ramp_response_numerical_data(
  system,
  slope = 1L,
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  n = 1000L,
  adaptive = FALSE,
  ...
)

bode_magnitude_numerical_data(
  system,
  initial_exp = -5L,
  final_exp = 5L,
  freq_unit = "rad/sec",
  n = 1000L,
  ...
)

bode_phase_numerical_data(
  system,
  initial_exp = -5L,
  final_exp = 5L,
  freq_unit = "rad/sec",
  phase_unit = "rad",
  phase_unwrap = TRUE,
  n = 1000L,
  ...
)

nyquist_plot_expr(system)

nichols_plot_expr(system)

pole_zero_plot(
  system,
  pole_color = "blue",
  pole_markersize = 10,
  zero_color = "orange",
  zero_markersize = 7,
  grid = TRUE,
  show_axes = TRUE,
  show = TRUE,
  ...
)

step_response_plot(
  system,
  color = "blue",
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  show_axes = FALSE,
  grid = TRUE,
  show = TRUE,
  ...
)

impulse_response_plot(
  system,
  color = "blue",
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  show_axes = FALSE,
  grid = TRUE,
  show = TRUE,
  ...
)

ramp_response_plot(
  system,
  slope = 1L,
  color = "blue",
  prec = 8L,
  lower_limit = 0L,
  upper_limit = 10L,
  show_axes = FALSE,
  grid = TRUE,
  show = TRUE,
  ...
)

bode_magnitude_plot(
  system,
  initial_exp = -5L,
  final_exp = 5L,
  color = "blue",
  show_axes = FALSE,
  grid = TRUE,
  show = TRUE,
  freq_unit = "rad/sec",
  ...
)

bode_phase_plot(
  system,
  initial_exp = -5L,
  final_exp = 5L,
  color = "blue",
  show_axes = FALSE,
  grid = TRUE,
  show = TRUE,
  freq_unit = "rad/sec",
  phase_unit = "rad",
  phase_unwrap = TRUE,
  ...
)

bode_plot(
  system,
  initial_exp = -5L,
  final_exp = 5L,
  grid = TRUE,
  show_axes = FALSE,
  show = TRUE,
  freq_unit = "rad/sec",
  phase_unit = "rad",
  phase_unwrap = TRUE,
  ...
)

nyquist_plot(
  system,
  initial_omega = 0.01,
  final_omega = 100,
  show = TRUE,
  color = "blue",
  n = 1000L,
  ...
)

nichols_plot(
  system,
  initial_omega = 0.01,
  final_omega = 100,
  show = TRUE,
  color = "blue",
  n = 1000L,
  ...
)
```

## Arguments

- system:

  A SISO system.

- prec:

  The precision of the evaluation (kept for compatibility).

- lower_limit, upper_limit:

  The time range.

- n:

  The number of points.

- adaptive:

  Ignored; the points are always equally spaced.

- ...:

  Further arguments.

- slope:

  The slope of the ramp.

- initial_exp, final_exp:

  The frequency range as powers of ten.

- freq_unit:

  `"rad/sec"` or `"Hz"`.

- phase_unit:

  `"rad"` or `"deg"`.

- phase_unwrap:

  Whether to unwrap the phase.

- pole_markersize, zero_markersize:

  The marker sizes.

- grid, show_axes, show:

  Whether to draw a grid, the axes, and the plot.

- color, pole_color, zero_color:

  The colors.

- initial_omega, final_omega:

  The frequency range of the Nyquist and Nichols plots.

## Value

A list of the x and y values (`pole_zero_numerical_data()`: of the zeros
and poles), a list of expressions, or, for the plots, the data
invisibly.

## Examples

``` r
s <- Symbol("s")
G <- TransferFunction(1, s^2 + s + 1, s)
pole_zero_numerical_data(G)
#> ([], [-0.5-0.8660254i, -0.5+0.8660254i]) 
step_response_numerical_data(G, n = 5)
#> [[1]]
#> [1]  0.0  2.5  5.0  7.5 10.0
#> 
#> [[2]]
#> [1] 0.0000000 1.0233596 1.0745906 0.9741517 1.0021701
#> 
```
