# Resultant and discriminant

`resultant(f, g, x)` is the resultant with respect to `x`, computed as
the determinant of the Sylvester matrix. `discriminant(f, x)` is
`(-1)^(n(n-1)/2) * resultant(f, f') / lc(f)`.

## Usage

``` r
resultant(f, g, x = NULL)

discriminant(f, x = NULL)
```

## Arguments

- f, g:

  Expressions.

- x:

  The variable.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
discriminant(a * x^2 + b * x + c, x)
#> -4*a*c + b^2
resultant(x^2 - 1, x - 2, x)
#> 3
```
