# Regression cases for limits and series, including branch cuts and Limit.doit.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("limits keep tiny exact coefficients and certify signs", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(limit(Float(1e-12) / x, x, 0)), "oo")
  expect_equal(sstr(limit(exp(-Float(1e-12) * x), x, oo)), "0")
  expect_equal(sstr(limit(sin(exp(Integer(-1000))) / x, x, 0)), "oo")
  expect_equal(sstr(limit(I / x, x, 0, "+-")), "zoo")
  expect_true(is_infinite(oo * I))
})

test_that("limits on inverse trigonometric branch cuts follow SymPy", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  expect_equal(sstr(limit((asin(2 + I * x) - asin(Integer(2))) / x, x, 0)), "oo*sign(pi - 2*asin(2))")
  expect_equal(sstr(limit(asin(2 + I * x) / x, x, 0)), "oo*sign(pi - asin(2))")
  expect_equal(sstr(limit(asin(2 + I * Rational(1, "10000000000000") * x), x, 0)), "pi - asin(2)")
  v <- limit(asin(2 + x + I * x^2), x, 0)
  expect_equal(as.complex(N(v)), as.complex(N(sym_pi - asin(Integer(2)))), tolerance = 1e-12)
})

test_that("log at a pole takes the side of the cut", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(limit(log(-1 / x - I) - log(1 / x), x, 0)), "-I*pi")
  expect_equal(sstr(limit(log(-x - I) - log(x), x, oo)), "-I*pi")
  expect_equal(sstr(limit(log(-1 / x + I) - log(1 / x), x, 0)), "I*pi")
  expect_equal(sstr(series(log(-1 / x - I), x, 0, 3)), "-I*pi - log(x) + I*x + x^2/2 + O(x^3)")
})

test_that("gamma and loggamma at poles", {
  skip_on_cran()
  x <- Symbol("x")
  oo <- S_Infinity
  expect_equal(sstr(limit(gamma(x), x, -1)), "-oo")
  expect_equal(sstr(limit(gamma(x - 1), x, 0)), "-oo")
  expect_equal(sstr(limit(x * gamma(x - 1), x, 0)), "-1")
  expect_equal(sstr(limit(gamma(x - 2) * x, x, 0)), "1/2")
  expect_equal(sstr(limit(gamma(x), x, -2, "-")), "-oo")
  expect_equal(sstr(gruntz(loggamma(x) / log(x), x, 0)), "-1")
  expect_equal(sstr(limit(loggamma(-x), x, oo)), "zoo")
  expect_equal(sstr(limit(loggamma(x - 1), x, 0)), "oo")
  expect_equal(sstr(series(loggamma(x), x, 0, 2)), "-log(x) - EulerGamma*x + O(x^2)")
  expect_equal(sstr(series(gamma(x - 2), x, 0, 1)), "1/(2*x) + 3/4 - EulerGamma/2 + O(x)")
  expect_error(series(loggamma(x - 1), x, 0, 2))
  expect_equal(sstr(loggamma(Rational(5, 2))), "log(3*sqrt(pi)/4)")
  expect_equal(sstr(loggamma(Integer(-3))), "oo")
  expect_equal(sstr(loggamma(-oo)), "zoo")
})

test_that("limit forms of SymPy's Limit.doit", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(limit(S_Infinity / (x^2 - 4), x, oo)), "Limit(oo/(x^2 - 4), x, oo, dir='-')")
  expect_equal(sstr(limit(3 * sec(-4 * sym_pi * x + x / 3), x, 3 * sym_pi / (-2 + 24 * sym_pi))), "-oo")
  expect_equal(sstr(limit((I * x - 1)^Rational(1, 3), x, 0, "-")), "-(-1)^(2/3)")
})

test_that("limits on cuts are corrected once", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  expect_equal(sstr(limit(sqrt(-1 - I * x / Integer(10)^400), x, 0)), "-I")
  expect_equal(sstr(limit(asin(2 + I * x) + floor(x), x, 0)), "pi - asin(2)")
  expect_equal(sstr(limit((-1 - I * x)^Rational(1, 3) + floor(x), x, 0)), "-(-1)^(2/3)")
  expect_equal(sstr(limit((asin(2 + I * x) - asin(Integer(2))) / x, x, 0)), "oo*sign(pi - 2*asin(2))")
})

test_that("series of logs and inverse trigonometric functions on cuts", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
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

test_that("limits that need Limit.doit's leading term step", {
  skip_on_cran()
  x <- Symbol("x")
  oo <- S_Infinity
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
