# Fast checks of the printer, array expression, matrix equality, set
# membership, free symbol, count_ops and fraction field fixes; the SymPy
# rows are in test-leftovers-sympy.R.

test_that("matrix sums print and sort as SymPy's MatAdd", {
  M <- MatrixSymbol("M", 2L, 2L)
  expect_identical(printer_doprint(NumPyPrinter(), M + Identity(2L)), "numpy.eye(2) + M")
  A <- MatrixSymbol("A", 1L, 3L)
  B <- MatrixSymbol("B", 1L, 3L)
  e <- A - B
  expect_identical(vapply(e$args, sstr, ""), c("-B", "A"))
  expect_identical(sstr(e), "A - B")
  expect_identical(tensorflow_code(new_sym("MatrixElement", list(e, S(0L), S(0L)))),
                   "(tensorflow.math.add((-1)*B, A))[0, 0]")
})

test_that("an IndexedBase without a shape has shape None in array expressions", {
  I <- IndexedBase("I")
  J <- IndexedBase("J")
  i <- Idx("i", list(0L, 2L))
  j <- Idx("j", list(0L, 2L))
  k <- Idx("k", list(0L, 2L))
  r <- convert_indexed_to_array(I(i, j) * J(j, k))
  expect_identical(sstr(r), "ArrayDiagonal(ArrayTensorProduct(I, J), (1, 2))")
  expect_null(get_shape(r))
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i, k), I(i, j) * J(j, k))),
                   "I = numpy.einsum(\",->\", I, J)")
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i), J(i, i))), "I = numpy.einsum(\"->\", J)")
  expect_identical(sstr(ArrayContraction(ArrayTensorProduct(I, J), list(0L, 1L))),
                   "ArrayContraction(ArrayTensorProduct(I, J), (0, 1))")
  A <- IndexedBase(Symbol("A"), Tuple(4L, 4L))
  expect_null(get_shape(ArrayTensorProduct(A, I)))
  expect_length(get_shape(ArrayTensorProduct(A, A)), 4L)
})

test_that("equal explicit matrices are equal and belong to a condition set", {
  w <- ImmutableMatrix(list(list(Rational(11, 39), Rational(16, 39), Rational(4, 13))))
  expect_identical(sstr(Eq(w, w)), "True")
  expect_identical(sstr(w), "Matrix([[11/39, 16/39, 4/13]])")
  T4 <- Matrix(list(list(Rational(1, 5), Rational(2, 5), Rational(2, 5)),
                    list(Rational(1, 10), S_Half, Rational(2, 5)),
                    list(Rational(3, 5), Rational(3, 10), Rational(1, 10))))
  wm <- MatrixSymbol("wm", 1L, 3L)
  cs <- ConditionSet(wm, Eq(wm * T4, wm))
  expect_identical(sstr(doit(contains(cs, w))), "True")
})

test_that("a random indexed symbol is one of its free symbols", {
  B <- BernoulliProcess("B", p = Rational(3, 5), success = 1, failure = 0)
  t <- Symbol("t", positive = TRUE, integer = TRUE)
  expect_identical(vapply(free_symbols(rv_getitem(B, 2 * t)), sstr, ""), c("B[2*t]", "t"))
  expect_identical(vapply(free_symbols(rv_getitem(B, 4)), sstr, ""), "B[4]")
})

test_that("count_ops counts a power of E once and simplify follows", {
  e <- (-18 * exp(S(3)) + 234 + exp(S(9))) * exp(S(-9))
  expect_identical(count_ops(expand_power_exp(expand_mul(e))), 6L)
  expect_identical(sstr(simplify(1 - 18 * (-13 + exp(S(3))) * exp(S(-9)))), "-18*exp(-6) + 234*exp(-9) + 1")
})

test_that("coefficients over a fraction field are a numerator over a denominator", {
  x <- Symbol("x")
  a <- Symbol("a")
  p <- Poly(x^2 * (a + 1) / a + x / (2 * a) - x / a^2, x)
  expect_identical(vapply(all_coeffs(p), sstr, ""), c("(a + 1)/a", "(a - 2)/(2*a^2)", "0"))
})
