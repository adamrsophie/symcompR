# The Risch algorithm

A port of SymPy's `risch_integrate()` (sympy/integrals/risch.py, rde.py
and prde.py) for elementary transcendental functions: exponentials and
logarithms over a tower built by `DifferentialExtension()`. Returns the
antiderivative; a part proven nonelementary is an unevaluated integral.

## Usage

``` r
risch_integrate(
  f,
  x,
  extension = NULL,
  handle_first = "log",
  separate_integral = FALSE,
  rewrite_complex = NULL,
  conds = "piecewise"
)
```

## Arguments

- f:

  The integrand.

- x:

  The variable.

- extension:

  A `DifferentialExtension` to use.

- handle_first:

  `"log"` or `"exp"`: which extensions are built first.

- separate_integral:

  Return `list(elementary, nonelementary)`.

- rewrite_complex:

  Rewrite trigonometric functions with `exp`.

- conds:

  `"piecewise"` for conditional results.

## Value

An expression or a list.

## Examples

``` r
x <- Symbol("x")
risch_integrate(exp(x^2)*x, x)
#> exp(x^2)/2
```
