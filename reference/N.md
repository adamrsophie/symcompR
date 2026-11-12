# Numerical evaluation

`N()` (alias `evalf()`) evaluates an expression to a floating point
number, optionally substituting values for symbols first. Evaluation is
in arbitrary precision binary arithmetic on gmp integers, as SymPy's
evalf (mpmath): the working precision grows until `digits` significant
digits are certain. Functions without an arbitrary precision evaluator
fall back to R double precision, capped at 15 digits.

## Usage

``` r
N(e, digits = 15, subs = list())

evalf(e, digits = 15, subs = list())
```

## Arguments

- e:

  A `sym` expression.

- digits:

  Significant digits.

- subs:

  Named list of values for symbols.

## Value

A `sym` Float, a complex `a + b*I` with Float parts, or the expression
with numeric subexpressions evaluated when free symbols remain.

## Examples

``` r
N(sym_pi)
#> 3.14159265358979
N(sqrt(S(2)), 5)
#> 1.4142
x <- Symbol("x")
N(sin(x), subs = list(x = 1))
#> 0.841470984807897
```
