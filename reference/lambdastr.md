# Python source of a lambda function

`lambdastr(args, expr)` returns the source of a Python lambda function
that evaluates `expr` at `args`, as SymPy's `lambdastr()`.

## Usage

``` r
lambdastr(args, expr, printer = NULL, dummify = NULL)
```

## Arguments

- args:

  A symbol, or a list of symbols (possibly nested).

- expr:

  An expression, or a list of expressions.

- printer:

  A function of an expression giving its source; the default is SymPy's
  `LambdaPrinter`.

- dummify:

  Whether to replace the arguments by dummy symbols. The default
  dummifies when an argument contains a function or a derivative.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
lambdastr(x, x^2)
#> [1] "lambda x: (x**2)"
lambdastr(list(x, list(y)), x + y)
#> [1] "lambda _0,_1: (lambda x,y: (x + y))(_0,_1[0])"
```
