# Integrals of fields over regions

`ParametricIntegral(field, region)` integrates a scalar or vector field
over a
[`ParametricRegion()`](https://adamrsophie.github.io/symcompR/reference/ParametricRegion.md):
line integrals for one parameter, surface (flux) integrals for two and
volume integrals for three. When the integral cannot be computed it
stays unevaluated. `vector_integrate()` also accepts implicit regions
and geometric entities, and otherwise falls back to
[`integrate()`](https://rdrr.io/r/stats/integrate.html).

## Usage

``` r
ParametricIntegral(field, parametricregion)

vector_integrate(field, ...)
```

## Arguments

- field:

  A scalar or vector field.

- parametricregion, ...:

  A region, or integration limits.

## Value

An expression.

## Examples

``` r
C <- CoordSys3D("C")
t <- Symbol("t")
ParametricIntegral(C$x * C$i + C$y * C$j, ParametricRegion(list(t, t), list(t, 0, 1)))
#> 1
```
