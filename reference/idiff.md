# Implicit differentiation

Returns `dy/dx` (or the `n`-th derivative) for `eq = 0` defining `y` as
a function of `x`, like SymPy's `idiff()`.

## Usage

``` r
idiff(eq, y, x, n = 1)
```

## Arguments

- eq:

  An expression equal to zero.

- y:

  The dependent symbol.

- x:

  The independent symbol.

- n:

  Order of the derivative.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
idiff(x^2 + y^2 - 4, y, x)
#> -x/y
```
