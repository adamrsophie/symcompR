# Stochastic processes and Markov chains; expected values are SymPy's.

test_that("discrete Markov chains answer probability queries", {
  T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
  Y <- DiscreteMarkovChain("Y", list(0, 1), T)
  expect_equal(sstr(P(Eq(rv_getitem(Y, 2), 1), Eq(rv_getitem(Y, 1), 0))), "1/2")
  expect_equal(sstr(rv_attr(Y, "stationary_distribution")), "Matrix([[1/3, 2/3]])")
  expect_equal(sstr(rv_attr(Y, "is_regular")), "True")
  expect_equal(sstr(rv_attr(Y, "number_of_states")), "2")
  expect_error(DiscreteMarkovChain(1), "neither a Symbol nor a string")
  expect_error(DiscreteMarkovChain("Z", list(0, 1, 2), T), "must be the same")
})

test_that("communication classes and decomposition follow SymPy", {
  T <- Matrix(list(list(1, 0, 0), list(Rational(1, 2), 0, Rational(1, 2)), list(0, 0, 1)))
  Y <- DiscreteMarkovChain("Y", list(0, 1, 2), T)
  cc <- rv_attr(Y, "communication_classes")
  expect_length(cc, 3L)
  expect_equal(sstr(rv_attr(Y, "is_absorbing_chain")), "True")
  expect_equal(sstr(rv_attr(Y, "absorbing_probabilities")), "Matrix([[1/2, 1/2]])")
})

test_that("continuous Markov chains and assumptions", {
  G <- Matrix(list(list(-1, 1), list(1, -1)))
  C <- ContinuousMarkovChain("C", list(0, 1), G)
  expect_equal(sstr(rv_attr(C, "limiting_distribution")), "Matrix([[1/2, 1/2]])")
  Y <- DiscreteMarkovChain("Y")
  expect_error(GeneratorMatrixOf(Y, G), "only ContinuousMarkovChain")
})

test_that("Poisson processes merge and split", {
  X <- PoissonProcess("X", 2)
  Z <- X + PoissonProcess("Y", 3)
  expect_equal(sstr(rv_attr(Z, "lamda")), "5")
  expect_error(X + 2, "Only instances of Poisson Process")
  expect_error(rv_attr(Z, "split", 3, 3), "Sum of l1 and l2 should be 5")
  expect_equal(sstr(P(Eq(X(1), 2) & Eq(X(2), 3))), "4*exp(-4)")
})

test_that("sampling follows SymPy without scipy", {
  expect_error(sampling_density(Normal("X", 0, 1)), "Failed to import scipy")
  expect_error(sample_stochastic_process(1), "Process must be an instance")
  set.seed(1)
  T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
  g <- sample_stochastic_process(DiscreteMarkovChain("Y", list(0, 1), T))
  expect_true(sstr(g()) %in% c("0", "1"))
})
