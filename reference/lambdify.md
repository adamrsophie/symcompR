# Turn an expression into an R function

`lambdify(args, e)` returns a fast vectorized R function, like SymPy's
`lambdify()`. The function body is built as an R call from the
expression tree.

## Usage

``` r
lambdify(args, e, cse = FALSE)
```

## Arguments

- args:

  A symbol or list of symbols.

- e:

  An expression or a list of expressions.

- cse:

  If `TRUE`, common subexpressions are computed once, into temporaries,
  before the result. This makes the function slower to build and faster
  to call when the expressions share work. Expressions with
  [`implemented_function()`](https://adamrsophie.github.io/symcompR/reference/implemented_function.md)
  calls or `Piecewise` are built without it, since moving a call out of
  a branch could evaluate it where it is not defined.

## Value

An R function.

## Examples

``` r
x <- Symbol("x")
f <- lambdify(x, x^2 + sin(x))
f(c(0, 1, 2))
#> [1] 0.000000 1.841471 4.909297
y <- Symbol("y")
h <- exp(sin(x + y)) + cos(x * y)
g <- lambdify(list(x, y), list(h^2, h^3), cse = TRUE)
g(1, 2)
#> [[1]]
#> [1] 0.5408367
#> 
#> [[2]]
#> [1] 0.39774
#> 
```
