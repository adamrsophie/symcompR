# Musculotendon actuators

Musculotendon models after SymPy's `MusculotendonDeGroote2016`: a force
actuator along a pathway whose force follows from the activation
dynamics, the characteristic curves of De Groote et al. (2016) and the
chosen
[`MusculotendonFormulation()`](https://adamrsophie.github.io/symcompR/reference/MusculotendonFormulation.md).
Constants left as `NULL` become symbols named after the model. The rigid
tendon formulation and the explicit fiber length and tendon force
formulations are available; the implicit ones raise
`NotImplementedError` as in SymPy.
`MusculotendonDeGroote2016_with_defaults()` uses the published curve
constants and default fiber velocity, pennation and damping values.
`MusculotendonBase()` is the abstract base class and raises an error.

## Usage

``` r
MusculotendonBase(...)

MusculotendonDeGroote2016(
  name,
  pathway,
  activation_dynamics,
  musculotendon_dynamics = 0L,
  tendon_slack_length = NULL,
  peak_isometric_force = NULL,
  optimal_fiber_length = NULL,
  maximal_fiber_velocity = NULL,
  optimal_pennation_angle = NULL,
  fiber_damping_coefficient = NULL,
  with_defaults = FALSE
)

MusculotendonDeGroote2016_with_defaults(
  name,
  pathway,
  activation_dynamics,
  musculotendon_dynamics = 0L,
  tendon_slack_length = NULL,
  peak_isometric_force = NULL,
  optimal_fiber_length = NULL,
  maximal_fiber_velocity = Float("10.0"),
  optimal_pennation_angle = Float("0.0"),
  fiber_damping_coefficient = Float("0.1")
)
```

## Arguments

- ...:

  Ignored.

- name:

  Name of the musculotendon.

- pathway:

  A pathway, such as
  [`LinearPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md).

- activation_dynamics:

  An activation model, such as
  [`FirstOrderActivationDeGroote2016()`](https://adamrsophie.github.io/symcompR/reference/activation-dynamics.md).

- musculotendon_dynamics:

  A
  [`MusculotendonFormulation()`](https://adamrsophie.github.io/symcompR/reference/MusculotendonFormulation.md)
  value.

- tendon_slack_length, peak_isometric_force, optimal_fiber_length,
  maximal_fiber_velocity, optimal_pennation_angle,
  fiber_damping_coefficient:

  Model constants.

- with_defaults:

  Use the published curve constants.

## Value

A musculotendon object.

## Details

The objects expose SymPy's attributes through `$`, among them `force`,
`pathway`, `activation_dynamics`, `state_vars` (`x`), `input_vars`
(`r`), `constants` (`p`), `M`, `F`,
[`rhs()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) and
`to_loads()`.

## Examples

``` r
N <- ReferenceFrame("N")
O <- PhysPoint("O"); I <- PhysPoint("I")
I$set_pos(O, dynamicsymbols("q") * N$x)
muscle <- MusculotendonDeGroote2016("m", LinearPathway(O, I), FirstOrderActivationDeGroote2016("m"))
muscle$x
#> Matrix(1, 1)
#> [a_m(t)]
muscle$constants
#> Matrix(31, 1)
#> [l_T_slack_m    ]
#> [F_M_max_m      ]
#> [l_M_opt_m      ]
#> [v_M_max_m      ]
#> [alpha_opt_m    ]
#> [beta_m         ]
#> [tau_a_m        ]
#> [tau_d_m        ]
#> [b_m            ]
#> [c_0_fl_T_m     ]
#> [c_1_fl_T_m     ]
#> [c_2_fl_T_m     ]
#> [c_3_fl_T_m     ]
#> [c_0_fl_M_pas_m ]
#> [c_1_fl_M_pas_m ]
#> [c_0_fl_M_act_m ]
#> [c_1_fl_M_act_m ]
#> [c_2_fl_M_act_m ]
#> [c_3_fl_M_act_m ]
#> [c_4_fl_M_act_m ]
#> [c_5_fl_M_act_m ]
#> [c_6_fl_M_act_m ]
#> [c_7_fl_M_act_m ]
#> [c_8_fl_M_act_m ]
#> [c_9_fl_M_act_m ]
#> [c_10_fl_M_act_m]
#> [c_11_fl_M_act_m]
#> [c_0_fv_M_m     ]
#> [c_1_fv_M_m     ]
#> [c_2_fv_M_m     ]
#> [c_3_fv_M_m     ]
```
