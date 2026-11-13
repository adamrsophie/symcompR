# Expected values are SymPy's outputs for the same calls
# (sympy/printing/tests/test_numpy.py, test_tensorflow.py, test_pycode.py and
# sympy/stats/tests/test_stochastic_process.py).

skip_unless_full_tests()

test_that("printing_test_numpy: test_issue_17006", {
  M <- MatrixSymbol("M", 2L, 2L)
  expect_identical(printer_doprint(NumPyPrinter(), M + Identity(2L)), "numpy.eye(2) + M")
})

test_that("printing_test_tensorflow: test_MatrixElement_printing", {
  A <- MatrixSymbol("A", 1L, 3L)
  B <- MatrixSymbol("B", 1L, 3L)
  F <- new_sym("MatrixElement", list(A - B, S(0L), S(0L)))
  expect_identical(tensorflow_code(F), "(tensorflow.math.add((-1)*B, A))[0, 0]")
})

test_that("printing_test_pycode: Assignment between IndexedBase objects without a shape", {
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

test_that("test_DiscreteMarkovChain: limiting and stationary distributions of a symbolic chain", {
  TO4 <- Matrix(list(list(Rational(1, 5), Rational(2, 5), Rational(2, 5)), list(Rational(1, 10), S_Half, Rational(2, 5)),
                     list(Rational(3, 5), Rational(3, 10), Rational(1, 10))))
  w <- ImmutableMatrix(list(list(Rational(11, 39), Rational(16, 39), Rational(4, 13))))
  TS1 <- MatrixSymbol("T", 3, 3)
  Y5 <- DiscreteMarkovChain("Y", trans_probs = TS1)
  expect_identical(sstr(doit(rv_attr(Y5, "limiting_distribution", w, TO4))), "True")
  cs <- subs(rv_attr(Y5, "stationary_distribution", condition_set = TRUE), TS1, TO4)
  expect_identical(sstr(doit(contains(cs, w))), "True")
})

test_that("test_ContinuousMarkovChain: transition probabilities", {
  T2 <- Matrix(list(list(-S_One, S_One, S_Zero), list(S_One, -S_One, S_Zero), list(S_Zero, S_One, -S_One)))
  C2 <- ContinuousMarkovChain("C", list(0, 1, 2), T2)
  A <- rv_attr(C2, "generator_matrix")
  t <- Symbol("t", positive = TRUE)
  expect_identical(sstr(apply_lambda(rv_attr(C2, "transition_probabilities", A), t)),
                   "Matrix([[1/2 + exp(-2*t)/2, 1/2 - exp(-2*t)/2, 0], [1/2 - exp(-2*t)/2, 1/2 + exp(-2*t)/2, 0], [1/2 - exp(-t) + exp(-2*t)/2, 1/2 - exp(-2*t)/2, exp(-t)]])")
})

test_that("test_BernoulliProcess: free symbols", {
  B <- BernoulliProcess("B", p = Float(0.6), success = 1, failure = 0)
  t <- Symbol("t", positive = TRUE, integer = TRUE)
  x <- Symbol("x", positive = TRUE, integer = TRUE)
  fs <- function(e) sort(vapply(free_symbols(e), sstr, ""))
  expect_identical(fs(rv_getitem(B, 2 * t)), sort(c("t", "B[2*t]")))
  expect_identical(fs(rv_getitem(B, 4)), "B[4]")
  expect_identical(fs(rv_getitem(B, x * t)), sort(c("B[t*x]", "t", "x")))
})

test_that("test_PoissonProcess and test_GammaProcess: simplified probabilities", {
  t <- Symbol("t", positive = TRUE)
  d <- Symbol("d", positive = TRUE)
  X <- PoissonProcess("X", 3)
  expect_identical(sstr(simplify(P(Not(Eq(X(t), 2) & (X(d) > 3)), Contains(t, Interval(2, 4, FALSE, TRUE)) & Contains(d, Interval(7, 8, TRUE, FALSE))))),
                   "-18*exp(-6) + 234*exp(-9) + 1")
  X <- GammaProcess("X", 1, 2)
  expect_identical(sstr(simplify(P(Not((X(t) < 5) & (X(d) > 3)), Contains(t, Interval(2, 4, FALSE, TRUE)) & Contains(d, Interval(7, 8, TRUE, FALSE))))),
                   "-4*exp(-3) + 472*exp(-8)/3 + 1")
  expect_identical(sstr(simplify(P((X(t) > 2) | (X(t) < 4), Contains(t, Interval(1, 4, FALSE, TRUE))))),
                   "-643*exp(-4)/15 + 109*exp(-2)/15 + 1")
})

test_that("test_GammaProcess_symbolic: expectation of a sum", {
  t <- Symbol("t", positive = TRUE)
  d <- Symbol("d", positive = TRUE)
  y <- Symbol("y", positive = TRUE)
  g <- Symbol("g", positive = TRUE)
  l <- Symbol("l", positive = TRUE)
  X <- GammaProcess("X", l, g)
  expect_identical(sstr(E(X(t)^2 + X(d) * 2 + X(y)^3, Contains(t, Interval(0, 1, TRUE, FALSE)) & Contains(d, Interval(1, 2, TRUE, FALSE)) &
                            Contains(y, Interval(3, 4, FALSE, TRUE)))),
                   "2*g/l + (g^2 + g)/l^2 + (g^3 + 3*g^2 + 2*g)/l^3")
})

test_that("test_handmade: the set of a finite distribution", {
  FH <- FiniteDistributionHandmade(Dict(list(1, 2), list(S(1) / 2, S(1) / 2)))
  expect_identical(sstr(rv_attr(FH, "set")), "{1, 2}")
})
