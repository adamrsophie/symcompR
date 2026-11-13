# Solve a linear system in generated code

`MatrixSolve(A, b)` stands for the solution `x` of `A x = b`, as SymPy's
`sympy.codegen.matrix_nodes.MatrixSolve`; the NumPy printer prints it as
`numpy.linalg.solve(A, b)`. It differentiates as
`MatrixSolve(A, b' - A' MatrixSolve(A, b))`.

## Usage

``` r
MatrixSolve(...)
```

## Arguments

- ...:

  The matrix and the vector.

## Value

A node.

## Examples

``` r
A <- MatrixSymbol("A", 3, 3); b <- MatrixSymbol("b", 3, 1)
printer_doprint(NumPyPrinter(), MatrixSolve(A, b))
#> [1] "numpy.linalg.solve(A, b)"
```
