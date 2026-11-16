# Regression cases for simplification, powsimp, expand_log and count_ops.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("count_ops counts a power of E once and simplify follows", {
  e <- (-18 * exp(S(3)) + 234 + exp(S(9))) * exp(S(-9))
  expect_identical(count_ops(expand_power_exp(expand_mul(e))), 6L)
  expect_identical(sstr(simplify(1 - 18 * (-13 + exp(S(3))) * exp(S(-9)))), "-18*exp(-6) + 234*exp(-9) + 1")
})

test_that("powsimp joins bases and expand_log takes out perfect powers", {
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  y <- Symbol("y"); p <- Symbol("p", positive = TRUE)
  expect_equal(sstr(powsimp((-1)^x / 2^x)), "(-1/2)^x")
  expect_equal(sstr(powsimp(y^x * 2^x)), "(2*y)^x")
  expect_equal(sstr(powsimp(Integer(-2)^x * 3^x)), "(-6)^x")
  expect_equal(sstr(powsimp(2^x / 3^x)), "(2/3)^x")
  expect_equal(sstr(powsimp(2^(2 * x))), "4^x")
  expect_equal(sstr(powsimp(sqrt(x) * sqrt(y))), "sqrt(x)*sqrt(y)")
  expect_equal(sstr(powsimp(sqrt(p) * sqrt(Symbol("q", positive = TRUE)))), "sqrt(p)*sqrt(q)")
  expect_equal(sstr(expand_log(log(Integer(4)))), "2*log(2)")
  expect_equal(sstr(expand_log(log(Rational(8, 3)))), "-log(3) + 3*log(2)")
  expect_equal(sstr(expand_log(log(Integer(12)))), "log(12)")
})
