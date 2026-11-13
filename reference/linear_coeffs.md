# Coefficients of a linear expression

`linear_coeffs(eq, ...)` returns the coefficients of the symbols `...`
in the linear expression `eq` followed by the constant term, as SymPy's
`linear_coeffs()`; with `dict = TRUE` it returns a dictionary keyed by
the symbols present, with key `1` for the constant. Nonlinear terms
raise a `NonlinearError`.

## Usage

``` r
linear_coeffs(eq, ..., dict = FALSE)
```

## Arguments

- eq:

  An expression.

- ...:

  The symbols.

- dict:

  Return a dictionary.

## Value

A list, or a dictionary.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
linear_coeffs(3*x + 2*y - 1, x, y)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 2
#> 
#> [[3]]
#> -1
#> 
```
