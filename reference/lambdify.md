# Turn an expression into an R function

`lambdify(args, e)` returns a fast vectorized R function, like SymPy's
`lambdify()`. The function body is built as an R call from the
expression tree.

## Usage

``` r
lambdify(args, e)
```

## Arguments

- args:

  A symbol or list of symbols.

- e:

  An expression or a list of expressions.

## Value

An R function.

## Examples

``` r
x <- Symbol("x")
f <- lambdify(x, x^2 + sin(x))
f(c(0, 1, 2))
#> [1] 0.000000 1.841471 4.909297
```
