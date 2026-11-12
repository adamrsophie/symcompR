# Expected values are SymPy's outputs for the same calls.

x <- Symbol("x"); y <- Symbol("y"); I <- sym_I; oo <- S_Infinity

test_that("arbitrary precision keeps tiny values and certifies zeros", {
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
  expect_true(sqrt(Integer(10)^800) == Integer(10)^400)
  expect_equal(sstr(nsimplify(Float("100000000000000000001", 50), rational = TRUE)), "100000000000000000001")
  expect_true(nsimplify(Float("1e400", 30)) == Integer(10)^400)
  expect_equal(sstr(FiniteSet(Float("9007199254740993", 30), Integer("9007199254740992"))), "{9007199254740992, 9007199254740993.0}")
  a <- Float("1.00000000000000000001", 30); b <- Float(1, 30)
  expect_true(Dict(list(a, b), list(0, 0)) == Dict(list(b, a), list(0, 0)))
})

test_that("limits on cuts are corrected once", {
  expect_equal(sstr(limit(sqrt(-1 - I * x / Integer(10)^400), x, 0)), "-I")
  expect_equal(sstr(limit(asin(2 + I * x) + floor(x), x, 0)), "pi - asin(2)")
  expect_equal(sstr(limit((-1 - I * x)^Rational(1, 3) + floor(x), x, 0)), "-(-1)^(2/3)")
  expect_equal(sstr(limit((asin(2 + I * x) - asin(Integer(2))) / x, x, 0)), "oo*sign(pi - 2*asin(2))")
})

test_that("series of logs and inverse trigonometric functions on cuts", {
  expect_equal(sstr(series(log(-1 - I * x^6), x, 0, 8)), "-I*pi + I*x^6 + O(x^8)")
  expect_equal(sstr(series(acot(-x), x, 0, 2)), "-pi/2 + x + O(x^2)")
  expect_equal(sstr(series(acot(x), x, 0, 2)), "pi/2 - x + O(x^2)")
  expect_equal(sstr(series(asin(1 + x^3), x, 0, 2)), "pi/2 - sqrt(2)*I*x^(3/2) + O(x^2)")
  expect_equal(sstr(series(asin(2 + I * x), x, 0, 3)), "-asin(2) + pi - sqrt(3)*x/3 + sqrt(3)*I*x^2/9 + O(x^3)")
  expect_equal(sstr(series(acos(2 + I * x), x, 0, 2)), "-acos(2) + sqrt(3)*x/3 + O(x^2)")
  expect_equal(sstr(series(atan(2 * I + x), x, 0, 2)), "I*atanh(2) - x/3 + O(x^2)")
  expect_equal(sstr(series(acsc(Rational(1, 2) - I * x), x, 0, 2)), "-acsc(1/2) + pi - 4*sqrt(3)*x/3 + O(x^2)")
  expect_equal(sstr(series(loggamma(x), x, 0, 0)), "-log(x) + O(1)")
})

test_that("small fixes", {
  expect_equal(sstr(loggamma(I * oo)), "zoo")
  expect_equal(sstr(RootSum(Integer(0), Lambda(x, x^2), x)), "nan")
  r <- symmetrize(x^2 + y^2, list(x, y))
  expect_equal(sstr(r[[1L]]), "-2*x*y + (x + y)^2")
})

test_that("limits that need Limit.doit's leading term step", {
  oo <- S_Infinity
  expect_equal(sstr(limit(x^(-sym_pi), x, 0, "-")), "-oo*(-1)^(1 - pi)")
  expect_equal(sstr(limit((x - 1)^Rational(-1, 3), x, 1, "-")), "-oo*(-1)^(2/3)")
  expect_equal(sstr(limit(x^2 + x * sin(x) + cos(x), x, oo)), "oo")
  expect_equal(sstr(limit(x^2 + x * sin(x) + cos(x), x, -oo)), "oo")
  expect_equal(sstr(limit(x^2 + 2 * x * cos(x) + cos(x)^2, x, oo)), "oo")
  expect_equal(sstr(limit((1 + 1 / (x^2 + cos(x)))^(x^2 + x), x, oo)), "E")
  expect_equal(sstr(limit(x / ((x - 6) * sinh(tanh(Float(0.03) * x)) + tanh(x) - Float(0.5)), x, oo)), "1/sinh(1)")
  expect_equal(sstr(limit(sqrt(4^(x + 1) + 1) / 2^x, x, oo)), "2")
  expect_equal(sstr(limit((Rational(2, 3)^(x + 1) + 1)^(2^x) / 2^(Rational(4, 3)^(x - 1)), x, oo)), "oo")
  expect_equal(sstr(limit((-1)^x / 2^x, x, oo)), "0")
  expect_equal(sstr(limit((x^(-x - 1) * (x^(x + 1) + (x + 1)^x))^x, x, oo)), "exp(E)")
})

test_that("powsimp joins bases and expand_log takes out perfect powers", {
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
