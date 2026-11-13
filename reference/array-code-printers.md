# Array code printers

Printer objects after SymPy's `CuPyPrinter` and `JaxPrinter`
(`sympy.printing.numpy`), which print Python code for the CuPy and JAX
array libraries as
[`NumPyPrinter()`](https://adamrsophie.github.io/symcompR/reference/code-printer-objects.md)
does for NumPy. The NumPy, SciPy, CuPy, JAX, TensorFlow and PyTorch
printers print array expressions
([`ArrayTensorProduct()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`ArrayContraction()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`ArrayDiagonal()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`PermuteDims()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`ArrayAdd()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`ZeroArray()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md),
[`OneArray()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md))
as einsum, transpose and add calls.

## Usage

``` r
CuPyPrinter(settings = NULL)

JaxPrinter(settings = NULL)
```

## Arguments

- settings:

  A named list of printer settings.

## Value

A printer, for
[`printer_doprint()`](https://adamrsophie.github.io/symcompR/reference/code-printer-objects.md).

## Examples

``` r
x <- Symbol("x")
printer_doprint(CuPyPrinter(), sin(x) + sqrt(x))
#> [1] "cupy.sqrt(x) + cupy.sin(x)"
printer_doprint(JaxPrinter(), And(x > 1, x < 2))
#> [1] "jax.numpy.all(jax.numpy.asarray([jax.numpy.greater(x, 1),jax.numpy.less(x, 2)]), axis=0)"
M <- MatrixSymbol("M", 2, 2)
N <- MatrixSymbol("N", 2, 2)
printer_doprint(NumPyPrinter(), ArrayContraction(ArrayTensorProduct(M, N), c(1, 2)))
#> [1] "numpy.einsum(\"ab,bc->ac\", M,N)"
```
