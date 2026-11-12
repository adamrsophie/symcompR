# Integration by the manual (step-by-step) method

`integral_steps(f, x)` returns the tree of integration rules that
SymPy's `manualintegrate` would apply; `manualintegrate(f, x)` evaluates
it to an antiderivative.

## Usage

``` r
integral_steps(f, var, ...)

manualintegrate(f, var)
```

## Arguments

- f:

  Expression.

- var:

  Integration variable.

- ...:

  Options of the IntegrationSolver (`max_depth`, `branch`).

## Value

A rule object (`integral_steps`) or a `sym` expression.
