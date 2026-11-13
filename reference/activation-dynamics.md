# Activation dynamics

Models of how a muscle's activation follows its neural excitation, as
SymPy's `sympy.physics.biomechanics.activation` classes.
`ZerothOrderActivation()` makes the activation equal to the excitation;
`FirstOrderActivationDeGroote2016()` uses the first order differential
equation of De Groote et al. (2016), with time constants and a smoothing
rate that default to symbols named after the model. The
`*_with_defaults()` variants use the published constants.
`ActivationBase()` is the abstract base class and raises an error, as
SymPy does.

## Usage

``` r
ActivationBase(name)

ZerothOrderActivation(name)

ZerothOrderActivation_with_defaults(name)

FirstOrderActivationDeGroote2016(
  name,
  activation_time_constant = NULL,
  deactivation_time_constant = NULL,
  smoothing_rate = NULL
)

FirstOrderActivationDeGroote2016_with_defaults(name)
```

## Arguments

- name:

  Name of the model, used in the symbols it creates.

- activation_time_constant, deactivation_time_constant, smoothing_rate:

  Constants of the first order model; `NULL` gives symbols.

## Value

An activation dynamics object.

## Details

The objects expose SymPy's attributes through `$`: `name`, `excitation`
(`e`), `activation` (`a`), `order`, `state_vars` (`x`), `input_vars`
(`r`), `constants` (`p`), `M`, `F` and the method
[`rhs()`](https://adamrsophie.github.io/symcompR/reference/Eq.md).

## Examples

``` r
act <- FirstOrderActivationDeGroote2016("bicep")
act$x
#> Matrix(1, 1)
#> [a_bicep(t)]
act$rhs()
#> Matrix(1, 1)
#> [((1/2 - tanh(b_bicep*(-a_bicep(t) + e_bicep(t)))/2)*(3*a_bicep(t)/2 + 1/2)/tau_d_bicep + (tanh(b_bicep*(-a_bicep(t) + e_bicep(t)))/2 + 1/2)/(tau_a_bicep*(3*a_bicep(t)/2 + 1/2)))*(-a_bicep(t) + e_bicep(t))]
ZerothOrderActivation("bicep")$r
#> Matrix(1, 1)
#> [e_bicep(t)]
```
