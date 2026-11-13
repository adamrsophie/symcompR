# Musculotendon characteristic curves

The characteristic curves of De Groote et al. (2016), as SymPy's
`sympy.physics.biomechanics` curve classes: the tendon force length
curve, the passive and active fiber force length curves, the fiber force
velocity curve and the inverses of three of them. Each call builds an
unevaluated function of its argument and constants;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
gives the formula, with `doit(f, evaluate = FALSE)` keeping its
arrangement, and [`diff()`](https://rdrr.io/r/base/diff.html)
differentiates through it. `inverse(f)` returns the constructor of the
inverse curve. The `*_with_defaults()` variants use the published
constants, as SymPy's `with_defaults()` class methods.

## Usage

``` r
CharacteristicCurveFunction(...)

TendonForceLengthDeGroote2016(l_T_tilde, c0, c1, c2, c3)

TendonForceLengthDeGroote2016_with_defaults(l_T_tilde)

TendonForceLengthInverseDeGroote2016(fl_T, c0, c1, c2, c3)

TendonForceLengthInverseDeGroote2016_with_defaults(fl_T)

FiberForceLengthPassiveDeGroote2016(l_M_tilde, c0, c1)

FiberForceLengthPassiveDeGroote2016_with_defaults(l_M_tilde)

FiberForceLengthPassiveInverseDeGroote2016(fl_M_pas, c0, c1)

FiberForceLengthPassiveInverseDeGroote2016_with_defaults(fl_M_pas)

FiberForceLengthActiveDeGroote2016(
  l_M_tilde,
  c0,
  c1,
  c2,
  c3,
  c4,
  c5,
  c6,
  c7,
  c8,
  c9,
  c10,
  c11
)

FiberForceLengthActiveDeGroote2016_with_defaults(l_M_tilde)

FiberForceVelocityDeGroote2016(v_M_tilde, c0, c1, c2, c3)

FiberForceVelocityDeGroote2016_with_defaults(v_M_tilde)

FiberForceVelocityInverseDeGroote2016(fv_M, c0, c1, c2, c3)

FiberForceVelocityInverseDeGroote2016_with_defaults(fv_M)
```

## Arguments

- ...:

  Arguments, none accepted.

- l_T_tilde:

  Normalized tendon length.

- c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11:

  Curve constants.

- fl_T:

  Normalized tendon force.

- l_M_tilde:

  Normalized fiber length.

- fl_M_pas:

  Normalized passive fiber force.

- v_M_tilde:

  Normalized fiber velocity.

- fv_M:

  Fiber force velocity multiplier.

## Value

A `sym` expression.

## Details

`CharacteristicCurveFunction()` is the abstract base class and raises an
error, as SymPy does.

## Examples

``` r
l <- Symbol("l_T_tilde")
f <- TendonForceLengthDeGroote2016_with_defaults(l)
f
#> TendonForceLengthDeGroote2016(l_T_tilde, 0.2, 0.995, 0.25, 33.93669377311689)
doit(f)
#> 4.327167934027389e-16*exp(33.93669377311689*l_T_tilde) - 0.25
diff(f, l)
#> 6.787338754623378*exp(33.93669377311689*(l_T_tilde - 0.995))
inverse(f)
#> function (fl_T, c0, c1, c2, c3) 
#> .bm_curve("TendonForceLengthInverseDeGroote2016", list(fl_T, 
#>     c0, c1, c2, c3))
#> <bytecode: 0x619422045af8>
#> <environment: namespace:symcompR>
```
