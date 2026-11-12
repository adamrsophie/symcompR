# Regions of the complex plane

`ComplexRegion(sets, polar = FALSE)` is the set of `x + y*I` with
`(x, y)` in `sets`, a
[`ProductSet()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
of two real sets or a union of them, as SymPy's
`CartesianComplexRegion`. With `polar = TRUE` the pairs are `(r, theta)`
and the points `r*(cos(theta) + I*sin(theta))`, with theta normalized
into `[0, 2*pi)` (`PolarComplexRegion`). `complex_region_from_real(s)`
embeds a real set (SymPy's `from_real`); `a_interval()` and
`b_interval()` are the unions of the first and second factors, and
`psets()` lists the products.
[`contains()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
[`Intersection()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
and
[`Union()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
accept regions.

## Usage

``` r
ComplexRegion(sets, polar = FALSE)

psets(region)

a_interval(region)

b_interval(region)

complex_region_from_real(s)

normalize_theta_set(theta)
```

## Arguments

- sets:

  A product of two real sets, or a union of products.

- polar:

  Polar form.

- region:

  A complex region.

- s:

  A real set.

- theta:

  A set of angles.

## Value

A set.

## Examples

``` r
a <- Interval(2, 3); b <- Interval(4, 6)
C1 <- ComplexRegion(ProductSet(a, b))
contains(C1, Rational(5, 2) + Rational(9, 2)*sym_I)
#> [1] TRUE
ComplexRegion(ProductSet(Interval(0, 1), Interval(0, 2*sym_pi)), polar = TRUE)
#> PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))
Intersection(ComplexRegion(ProductSet(Interval(1, 3), Interval(-1, 1))), Reals)
#> Interval(1, 3)
```
