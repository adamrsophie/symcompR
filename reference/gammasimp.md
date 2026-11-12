# Simplify gamma functions and combinatorial expressions

`gammasimp(expr)` simplifies products and quotients of gamma functions
using the functional equation, the reflection and duplication formulas
and the multiplication theorem, as SymPy's `gammasimp()`. `combsimp()`
does the same for factorials and binomials, keeping the result in terms
of them when the arguments are integers.

## Usage

``` r
gammasimp(expr)

combsimp(expr)
```

## Arguments

- expr:

  An expression.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); n <- Symbol("n", integer = TRUE)
gammasimp(gamma(x + 1) / gamma(x))
#> x
combsimp(factorial(n) / factorial(n - 3))
#> n*(n - 2)*(n - 1)
```
