# Common factors of terms

`factor_terms(e)` pulls out the factors shared by all terms of every sum
in `e`, without expanding: `3*x*y*(x + 2*y)`. `gcd_terms(e)` does this
for the top level sum and also puts it over a common denominator, as
SymPy's functions of the same names.

## Usage

``` r
gcd_terms(e, isprimitive = FALSE, clear = TRUE, fraction = TRUE)

factor_terms(e, radical = FALSE, clear = FALSE, fraction = FALSE, sign = TRUE)
```

## Arguments

- e:

  An expression.

- isprimitive:

  The terms are already primitive.

- clear:

  Clear a rational coefficient from a sum (`FALSE` keeps a denominator
  inside when some term keeps an integer coefficient).

- fraction:

  Put the result over a common denominator.

- radical:

  Also take out common radicals.

- sign:

  Take out a common negative sign.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
factor_terms(3*x^2*y + 6*x*y^2)
#> 3*x*y*(x + 2*y)
gcd_terms(x / 2 + 1 / x)
#> (x^2 + 2)/(2*x)
```
