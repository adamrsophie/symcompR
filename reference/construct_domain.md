# Constructing the minimal domain

`construct_domain(obj)` finds the smallest domain containing the
coefficients `obj` (an expression or a list of them) and converts them,
as SymPy's `construct_domain()`: `ZZ` or `QQ` for rational numbers, `RR`
or `CC` for floats, `ZZ_I` or `QQ_I` for Gaussian rationals, a
polynomial ring or a fraction field over those for expressions in
symbols, and `EX` otherwise; with `extension = TRUE` a single radical
gives an algebraic field.

## Usage

``` r
construct_domain(obj, extension = FALSE, field = FALSE)
```

## Arguments

- obj:

  An expression or a list of expressions.

- extension:

  Use algebraic fields for radicals.

- field:

  Return a field (`QQ` rather than `ZZ`, a fraction field rather than a
  polynomial ring).

## Value

A list of the domain and the converted elements.

## Examples

``` r
construct_domain(list(1, Rational(1, 2)))
#> [[1]]
#> QQ
#> 
#> [[2]]
#> [[2]][[1]]
#> Big Rational ('bigq') :
#> [1] 1
#> 
#> [[2]][[2]]
#> Big Rational ('bigq') :
#> [1] 1/2
#> 
#> 
x <- Symbol("x")
construct_domain(list(x, 1))[[1L]]
#> ZZ[x]
construct_domain(list(1/x))[[1L]]
#> ZZ(x)
```
