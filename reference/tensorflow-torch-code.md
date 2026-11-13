# TensorFlow and PyTorch code printing

Code printers after SymPy's `TensorflowPrinter`
(`sympy.printing.tensorflow`) and `TorchPrinter`
(`sympy.printing.pytorch`). `tensorflow_code()` and `torch_code()` print
an expression as Python code calling the TensorFlow or PyTorch
functions: elementary functions, relationals, logic, `Piecewise` (as
[`where()`](https://adamrsophie.github.io/symcompR/reference/given.md)),
derivatives (as gradient calls), matrices and matrix expressions, and
array expressions (as einsum and transpose calls).

## Usage

``` r
TensorflowPrinter(settings = NULL)

tensorflow_code(expr, ...)

TorchPrinter(settings = NULL)

torch_code(
  expr,
  requires_grad = FALSE,
  dtype = "torch.float64",
  assign_to = NULL
)
```

## Arguments

- settings:

  A named list of printer settings: `tensorflow_version` for
  `TensorflowPrinter()`; `requires_grad` and `dtype` for
  `TorchPrinter()`.

- expr:

  An expression.

- ...:

  Printer settings, such as `tensorflow_version = "1.13"`.

- requires_grad:

  Whether printed tensors require gradients.

- dtype:

  The dtype of printed tensors.

- assign_to:

  Optional variable to assign to.

## Value

A printer, or a character scalar of code.

## Examples

``` r
x <- Symbol("x")
y <- Symbol("y")
tensorflow_code(sin(x) + x^2)
#> [1] "tensorflow.math.pow(x, 2) + tensorflow.math.sin(x)"
tensorflow_code(Piecewise(list(x, x < 0), list(y, TRUE)))
#> [1] "tensorflow.where(tensorflow.math.less(x, 0), x, y)"
torch_code(Max(x, y))
#> [1] "torch.maximum(x, y)"
torch_code(Matrix(list(list(x, 1), list(2, y))), requires_grad = TRUE)
#> [1] "torch.tensor([[x, 1], [2, y]], dtype=torch.float64, requires_grad=True)"
printer_doprint(TensorflowPrinter(), exp(x))
#> [1] "tensorflow.math.exp(x)"
```
