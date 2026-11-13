# Matrix form of a linear system

`linear_eq_to_matrix(equations, ...)` writes the linear equations (or
expressions equated to zero) as `A x = b` for the given symbols and
returns the pair `(A, b)`, as SymPy's `linear_eq_to_matrix()`. The
symbols may be given one by one or as a list; nonlinear terms raise a
`NonlinearError`.

## Usage

``` r
linear_eq_to_matrix(equations, ...)
```

## Arguments

- equations:

  A list of equations or expressions, one of them, or a matrix of
  expressions.

- ...:

  The symbols, or one list of them.

## Value

A tuple of the matrices `A` and `b`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
linear_eq_to_matrix(list(x + 2*y - 3, 3*x - y), list(x, y))
#> (Matrix([[1, 2], [3, -1]]), Matrix([[3], [0]])) 
```
