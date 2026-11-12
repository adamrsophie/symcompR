# Fu trigonometric simplification

`fu(rv)` simplifies trigonometric expressions with the algorithm of Fu,
Zhong and Zeng, as SymPy's `fu()`: it applies the transformation rules
below in greedy combinations and keeps the result with the fewest
trigonometric functions, then the fewest operations. Each rule is also
exported: `TR0` (normal, factor, expand), `TR1` (sec, csc to cos, sin),
`TR2` (tan, cot to sin, cos), `TR2i` (the reverse), `TR3` (induced
formulas), `TR4` (special angles), `TR5` and `TR6` (sin^2 and cos^2 to
the other), `TR7` (lower the degree of cos^2), `TR8` (products to sums),
`TR9` (sums to products), `TR10` and `TR10i` (expand and contract sums
of angles), `TR11` and `TR11` with a base (double angles), `TR12` and
`TR12i` (tangent of sums), `TR13` (products of tangents), `TRmorrie`
(Morrie's law), `TR14` (factored powers of sin and cos), `TR15`, `TR16`
and `TR22` (powers in terms of cot, tan, sec and csc), `TR111` (negative
powers to the reciprocal functions) and `TRpower` (powers to multiple
angles).

## Usage

``` r
TR0(rv)

TR1(rv)

TR2(rv)

TR2i(rv, half = FALSE)

TR3(rv)

TR4(rv)

TR5(rv, max = 4, pow = FALSE)

TR6(rv, max = 4, pow = FALSE)

TR7(rv)

TR8(rv, first = TRUE)

TR9(rv)

TR10(rv, first = TRUE)

TR10i(rv)

TR11(rv, base = NULL)

TR12(rv, first = TRUE)

TR13(rv)

TRmorrie(rv)

TR12i(rv)

TR14(rv, first = TRUE)

TR15(rv, max = 4, pow = FALSE)

TR16(rv, max = 4, pow = FALSE)

TR111(rv)

TR22(rv, max = 4, pow = FALSE)

TRpower(rv)

fu(rv, measure = .fu_measure)

hyper_as_trig(rv)
```

## Arguments

- rv:

  An expression.

- half:

  For `TR2i()`, also recognize half angle forms.

- max, pow:

  For the power rules: the largest exponent, and whether perfect powers
  are handled.

- first:

  Internal; leave the default.

- base:

  For `TR11()`, the angle to write in terms of.

- measure:

  The function to minimize, returning a numeric vector compared
  lexicographically.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
fu(sin(50)^2 + cos(50)^2 + sin(pi/6))
#> 1.50000000000000
fu(sqrt(6)*cos(x) + sqrt(2)*sin(x))
#> 1.4142135623731*sin(x) + 2.44948974278318*cos(x)
TR10i(cos(1)*cos(3) + sin(1)*sin(3))
#> -0.416146836547142
TR8(cos(x)*sin(y))
#> -sin(x - y)/2 + sin(x + y)/2
```
