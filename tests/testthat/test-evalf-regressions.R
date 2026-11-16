# Regression cases for numeric evaluation, Float precision and sympify of decimals.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("Floats keep their multiprecision value in facts, Integer, Rational and Mod", {
  skip_on_cran()
  expect_false(ask_fact(Float("1e-400", 30), "zero"))
  expect_true(ask_fact(Float("1e400", 30), "finite"))
  expect_equal(sstr(Integer(Float("9007199254740993", 30))), "9007199254740993")
  expect_equal(sstr(Integer(1.9)), "1")
  expect_equal(sstr(Integer(-1.9)), "-1")
  expect_equal(sstr(Integer(Rational(3, 2))), "1")
  expect_equal(sstr(Integer(Float("-2.5"))), "-2")
  expect_equal(sstr(Rational(Float("0.5"))), "1/2")
  expect_equal(sstr(Rational(Float(0.1))), "3602879701896397/36028797018963968")
  expect_equal(sstr(sym_Mod(Float("1.00000000000000000001", 30), 1)), "9.99999999997452451411988619036e-21")
  expect_equal(sstr(sym_Mod(Float(5.5), 2)), "1.50000000000000")
  expect_equal(sstr(sym_Mod(Float("7.25", 30), Rational(1, 2))), "0.250000000000000000000000000000")
})

test_that("sympify keeps decimal precision", {
  skip_on_cran()
  expect_equal(sstr(sympify("1e-400")), "1.00000000000000e-400")
  expect_equal(sstr(sympify("0.123456789123456789123456789")), "0.123456789123456789123456789")
})

test_that("arbitrary precision keeps tiny and cancelling values", {
  skip_on_cran()
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(N(sin(exp(Integer(-1000))), 20)), "5.0759588975494567653e-435")
  expect_equal(sstr(N(sin(Integer(355)), 30)), "-0.0000301443533594884492143302800087")
  expect_equal(sstr(N(log(1 + exp(Integer(-100))), 20)), "3.7200759760208359630e-44")
  expect_equal(sstr(N(S_ComplexInfinity)), "zoo")
  expect_equal(sstr(N(exp(oo * (I + Rational(1, "10000000000000"))))), "zoo")
})

test_that("arbitrary precision keeps tiny values and certifies zeros", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  expect_equal(sstr(N(log(1 + Rational("1e-100")))), "1.00000000000000e-100")
  expect_equal(sstr(N(floor(1 / (1 + exp(Integer(-1000)))))), "0")
  expect_equal(sstr(N(floor(log(Integer(8)) / log(Integer(2))))), "3.00000000000000")
  expect_equal(sstr(N(Rational("1e-400"))), "1.00000000000000e-400")
  expect_equal(sstr(N(x + Rational("1e-400"), 30)), "x + 1.0e-400")
  expect_equal(sstr(Integer(0)^Float("-1e-400", 50)), "zoo")
  expect_equal(sstr(Float(0)^Float("1e-400", 30)), "0.0")
  expect_equal(sstr(N(I * exp(Integer(-1000)), 30)), "5.07595889754945676529180947957e-435*I")
  expect_equal(sstr(Float("1e-200")^Integer(2)), "1.00000000000000e-400")
  expect_equal(sstr(N(log(1 + exp(Integer(-200))), 20)), "1.3838965267367375306e-87")
  expect_equal(sstr(N(Integer(0))), "0")
})

test_that("Float functions at the Floats' precision", {
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  expect_equal(sstr(sinh(Float("1e-100", 30))), "1.00000000000000000000000000000e-100")
  expect_equal(sstr(tanh(Float("1e-100", 30))), "1.00000000000000000000000000000e-100")
  expect_equal(sstr(asinh(Float("1e-100", 30))), "1.00000000000000000000000000000e-100")
  expect_equal(sstr(erfc(Float(20, 30))), "5.39586561160790092893499916791e-176")
  expect_equal(sstr(erfc(Float(200, 30))), "4.68935929598360149803569977699e-17375")
  expect_equal(sstr(log(Float("1.000000000000000000000000000001", 50))), "9.9999999999999999999961694691082510426940898108545e-31")
  expect_equal(sstr(Float("-1", 50)^Float("-0.5", 50)), "-1.0*I")
  expect_equal(sstr(sin(Float("1", 5))), "0.84147")
  expect_equal(sstr(N(sin(Float("0.1", 5)), 20)), "0.099833488464355468750")
  expect_equal(sstr(N(x^y, subs = list(x = 0, y = 0))), "1.00000000000000")
  expect_equal(sstr(N(x + y, 30, subs = list(x = Symbol("z")))), "y + z")
})

test_that("exact values from high precision Floats", {
  skip_on_cran()
  expect_true(sqrt(Integer(10)^800) == Integer(10)^400)
  expect_equal(sstr(nsimplify(Float("100000000000000000001", 50), rational = TRUE)), "100000000000000000001")
  expect_true(nsimplify(Float("1e400", 30)) == Integer(10)^400)
  expect_equal(sstr(FiniteSet(Float("9007199254740993", 30), Integer("9007199254740992"))), "{9007199254740992, 9007199254740993.0}")
  a <- Float("1.00000000000000000001", 30); b <- Float(1, 30)
  expect_true(Dict(list(a, b), list(0, 0)) == Dict(list(b, a), list(0, 0)))
})

test_that("N keeps a symbolic evalf rewrite", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(N(logaddexp(x, y))), "log(exp(x) + exp(y))")
})
