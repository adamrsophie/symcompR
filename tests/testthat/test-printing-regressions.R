# Regression cases for printers and code printers.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("printing and drawing follow SymPy", {
  x <- Symbol("x"); n <- Symbol("n")
  expect_identical(latex(-oo), "-\\infty")
  expect_identical(latex(besselj(n, x)), "J_{n}\\left(x\\right)")
  expect_identical(format(ArrowStringDescription("mm", NULL, NULL, NULL, NULL, "d", "r", "_", "f")), "\\ar[dr]_{f}")
})

test_that("printing_test_numpy: test_issue_17006", {
  skip_unless_full_tests()
  M <- MatrixSymbol("M", 2L, 2L)
  expect_identical(printer_doprint(NumPyPrinter(), M + Identity(2L)), "numpy.eye(2) + M")
})

test_that("printing_test_tensorflow: test_MatrixElement_printing", {
  skip_unless_full_tests()
  A <- MatrixSymbol("A", 1L, 3L)
  B <- MatrixSymbol("B", 1L, 3L)
  F <- new_sym("MatrixElement", list(A - B, S(0L), S(0L)))
  expect_identical(tensorflow_code(F), "(tensorflow.math.add((-1)*B, A))[0, 0]")
})

test_that("printing_test_pycode: Assignment between IndexedBase objects without a shape", {
  skip_unless_full_tests()
  I <- IndexedBase("I")
  J <- IndexedBase("J")
  i <- Idx("i", list(0L, 2L))
  j <- Idx("j", list(0L, 2L))
  k <- Idx("k", list(0L, 2L))
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i, k), I(i, j) * J(j, k))), "I = numpy.einsum(\",->\", I, J)")
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i, j), J(j, i))), "I = J")
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i), J(i, i))), "I = numpy.einsum(\"->\", J)")
  expect_identical(printer_doprint(CuPyPrinter(), Assignment(I(i, k), I(i, j) * J(j, k))), "I = cupy.einsum(\",->\", I, J)")
  expect_identical(printer_doprint(JaxPrinter(), Assignment(I(i, k), I(i, j) * J(j, k))), "I = jax.numpy.einsum(\",->\", I, J)")
})

test_that("negated powers with negative exponents print as fractions", {
  skip_on_cran()
  y <- Symbol("y")
  z <- Symbol("z")
  a <- Symbol("a")
  expect_equal(sstr(-1 / (sqrt(z) + 1)^(2 * a)), "-1/(sqrt(z) + 1)^(2*a)")
  expect_equal(sstr(z - 1 / (sqrt(z) + 1)^(2 * a)), "z - 1/(sqrt(z) + 1)^(2*a)")
  expect_equal(sstr(-y^(-a - 1)), "-y^(-a - 1)")
})

test_that("preview() raises when the output file cannot be written", {
  x <- Symbol("x")
  bad <- file.path(tempdir(), "no-such-directory", "out.tex")
  expect_error(preview(x, output = "tex", viewer = "file", filename = bad), "could not write")
})

test_that("lambdify with cse computes shared work once and gives the same values", {
  x <- Symbol("x"); y <- Symbol("y")
  u <- sin(x + y); v <- cos(x * y); h <- exp(u) + v
  ex <- list(u + h, h^2 + v, h^3 + exp(u))
  f0 <- lambdify(list(x, y), ex); f1 <- lambdify(list(x, y), ex, cse = TRUE)
  a <- c(0.1, 0.5, 2); b <- c(1, -1, 3)
  expect_equal(f1(a, b), f0(a, b))
  expect_true(any(grepl(".cse", deparse(body(f1)), fixed = TRUE)))
  M <- Matrix(list(list(sin(x)^2, sin(x)), list(1, x)))
  expect_equal(lambdify(x, M, cse = TRUE)(1), lambdify(x, M)(1))
  pw <- Piecewise(list(log(x), x > 0), list(log(x) + 1, TRUE))
  expect_false(any(grepl(".cse", deparse(body(lambdify(x, pw, cse = TRUE))), fixed = TRUE)))
})

test_that("implementations with different captured values are different functions", {
  maker <- function(k) function(z) z + k
  x <- Symbol("x")
  f1 <- implemented_function("f", maker(1)); f2 <- implemented_function("f", maker(2))
  expect_false(identical(f1(x)$key, f2(x)$key))
  g <- function(z) z^2
  expect_identical(implemented_function("g", g)(x)$key, implemented_function("g", g)(x)$key)
})

test_that("lambdify calls elementwise numeric implementations once on the vector", {
  x <- Symbol("x")
  v <- c(-0.5, 0.1, 0.7)
  f <- lambdify(x, erfinv(x) + cot(x + 2))
  expect_equal(f(v), stats::qnorm((v + 1) / 2) / sqrt(2) + 1 / tan(v + 2))
  expect_false(any(grepl("mapply", deparse(body(f)), fixed = TRUE)))
})
