# Fast checks of the array, TensorFlow, PyTorch, numexpr and interval
# printers and of the branching strategies; the SymPy battery is in
# test-printers-sympy.R.

test_that("TensorFlow and PyTorch printers print functions and matrices", {
  x <- Symbol("x")
  y <- Symbol("y")
  expect_identical(tensorflow_code(sin(x)), "tensorflow.math.sin(x)")
  expect_identical(tensorflow_code(sqrt(x)), "tensorflow.math.sqrt(x)")
  expect_identical(tensorflow_code(Piecewise(list(x, x < 0), list(y, TRUE)), tensorflow_version = "0.9"),
                   "tensorflow.select(tensorflow.math.less(x, 0), x, y)")
  expect_identical(torch_code(Max(x, y)), "torch.maximum(x, y)")
  expect_identical(torch_code(Heaviside(x)), "torch.heaviside(x, 1/2)")
  expect_identical(torch_code(Matrix(list(list(1, 2), list(3, 4)))),
                   "torch.tensor([[1, 2], [3, 4]], dtype=torch.float64)")
  expect_error(torch_code(polygamma(1L, x)), "digamma")
})

test_that("array expressions print as einsum calls", {
  M <- MatrixSymbol("M", 3L, 3L)
  N <- MatrixSymbol("N", 3L, 3L)
  expect_identical(printer_doprint(NumPyPrinter(), convert_matrix_to_array(M * N)),
                   "numpy.einsum(\"ab,bc->ac\", M,N)")
  expect_identical(printer_doprint(JaxPrinter(), ZeroArray(2L, 3L)), "jax.numpy.zeros((2,3,))")
  expect_identical(printer_doprint(CuPyPrinter(), OneArray(5L)), "cupy.ones((5,))")
  expect_identical(printer_doprint(NumPyPrinter(), np_maximum(Symbol("x"), Symbol("y"))),
                   "functools.reduce(numpy.maximum, [x, y])")
  expect_identical(pycode(List(1L, 2L)), "[1, 2]")
  expect_identical(printer_doprint(NumPyPrinter(), Array(list(list(1L, 2L), list(3L, 5L)))),
                   "numpy.array([[1, 2], [3, 5]])")
})

test_that("lambdarepr, numexpr and interval printers follow SymPy", {
  x <- Symbol("x")
  y <- Symbol("y")
  expect_identical(lambdarepr(sqrt(x)), "sqrt(x)")
  expect_identical(printer_doprint(NumExprPrinter(), sin(x)), "numexpr.evaluate('sin(x)', truediv=True)")
  expect_error(printer_doprint(NumExprPrinter(), Ei(x)), "does not support function")
  expect_identical(printer_doprint(IntervalPrinter(), x / 2), "(mpi('1/2'))*x")
})

test_that("branching strategies give SymPy's results", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L) else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L) else list(x)
  }
  expect_identical(unlist(branch_exhaust(branch5)(5L)), c(10L, 0L))
  expect_identical(unlist(multiplex(posdec, branch5)(5L)), c(4L, 6L))
  expect_identical(unlist(notempty(function(x) list())(3L)), 3L)
  expect_identical(unlist(sfilter(function(x) x %% 2L == 0L, branch5)(5L)), c(6L, 4L))
  out <- capture.output(invisible(branch_debug(posdec)(5L)))
  expect_identical(out, c("Rule: posdec", "In: 5", "Out: 4", ""))
  x <- Symbol("x")
  inc <- function(e) if (is_sym(e) && e$type == "Integer") list(e + 1L) else list()
  expect_identical(vapply(branch_top_down(inc)(x + 1), sstr, ""), "x + 2")
})
