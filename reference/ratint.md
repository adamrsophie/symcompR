# Integration of rational functions

A port of SymPy's `ratint()` (sympy/integrals/rationaltools.py): the
Hermite reduction (`ratint_ratpart()`), the Lazard, Rioboo and Trager
logarithmic part (`ratint_logpart()`) and its conversion to real
functions (`log_to_real()`, `log_to_atan()`).

## Usage

``` r
ratint(f, x, symbol = "t", real = NULL)

ratint_ratpart(f, g, x)

ratint_logpart(f, g, x, t = NULL)

log_to_atan(f, g)

log_to_real(h, q, x, t)
```

## Arguments

- f:

  A rational function, or `list(p, q)`.

- x:

  The variable.

- symbol:

  The name of the variable of the roots in `RootSum` terms.

- real:

  `TRUE` to convert to real functions, `FALSE` for `RootSum` terms,
  `NULL` to decide from the atoms of `f`.

- g:

  The denominator.

- t:

  The variable of the roots.

- h, q:

  Polynomials in `x` and `t`:
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
  objects or dense polynomials.

## Value

The antiderivative.

## Examples

``` r
x <- Symbol("x")
ratint(1/(x^2 + 1), x)
#> atan(x)
ratint(36/(x^5 - 2*x^4 - 2*x^3 + 4*x^2 + x - 2), x)
#> (12*x + 6)/(x^2 - 1) + 4*log(x - 2) - 4*log(x + 1)
```
