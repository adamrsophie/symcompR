# Integrate an expression

`sym_integrate(f, x)` returns an antiderivative;
`sym_integrate(f, list(x, a, b))` a definite integral, evaluated with
limits at the endpoints. Several variables integrate in turn. Rules
cover polynomials, rational functions (partial fractions with logarithms
and arctangents), exponentials, logarithms, trigonometric and hyperbolic
functions and their powers, inverse trigonometric functions,
substitution, integration by parts and a table of special forms (`erf`,
`Si`, `Ei`, `li`, ...). When no rule applies an unevaluated `Integral`
is returned. For an R function as first argument this calls
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

## Usage

``` r
Integral(f, ...)

sym_integrate(
  f,
  ...,
  meijerg = NULL,
  conds = "piecewise",
  risch = NULL,
  heurisch = NULL,
  manual = NULL
)
```

## Arguments

- f:

  An expression, or an R function for
  [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

- ...:

  Variables or `list(var, lower, upper)` specifications, or arguments
  for [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

- meijerg, risch, heurisch, manual:

  `TRUE` to use only that algorithm, `FALSE` to skip it, `NULL` (the
  default) to try it among the others, as SymPy's hints.

- conds:

  How conditions of convergence are returned: `"piecewise"`,
  `"separate"` or `"none"`.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
sym_integrate(x^2 * exp(x), x)
#> (x^2 - 2*x + 2)*exp(x)
sym_integrate(1 / (x^2 + 1), list(x, 0, 1))
#> pi/4
sym_integrate(exp(-x^2), list(x, -oo, oo))
#> sqrt(pi)
```
