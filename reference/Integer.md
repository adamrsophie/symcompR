# Exact integers, rationals and floating point numbers

`Integer()` builds an exact integer of arbitrary size, `Rational()` an
exact fraction reduced to lowest terms and `Float()` an arbitrary
precision binary floating point number, as SymPy's `Float`: 15
significant digits by default, or `dps` digits. Strings are the safe way
to pass integers larger than `2^53`, and decimals with more than 15
digits, since R doubles round them.

## Usage

``` r
Integer(p)

Rational(p, q = 1)

Float(x, dps = NULL, precision = NULL)
```

## Arguments

- p:

  Numerator: a number, string, or `gmp` big integer.

- q:

  Denominator.

- x:

  A number or numeric string.

- dps:

  Decimal digits of precision for `Float()`.

- precision:

  Binary precision in bits for `Float()`; overrides `dps`.

## Value

A `sym` number.

## Examples

``` r
Integer("123456789012345678901234567890")
#> 123456789012345678901234567890
Rational(6, 4)
#> 3/2
Float(0.1)
#> 0.100000000000000
Float("3.14159265358979323846264338327950288", 30)
#> 3.14159265358979323846264338328
```
