# Reference values from SymPy's gruntz() on cases from its test_gruntz.py.
test_that("gruntz matches SymPy on nested exponential and logarithmic limits", {
  x <- Symbol("x"); a <- Symbol("a", positive = TRUE)
  cases <- list(
    c("exp(x)*(exp(1/x - exp(-x)) - exp(1/x))", "oo", "+", "-1"),
    c("exp(x)*(exp(1/x + exp(-x) + exp(-x**2)) - exp(1/x - exp(-exp(x))))", "oo", "+", "1"),
    c("exp(exp(x - exp(-x))/(1 - 1/x)) - exp(exp(x))", "oo", "+", "oo"),
    c("exp(exp(exp(x + exp(-x)))) / exp(exp(exp(x)))", "oo", "+", "oo"),
    c("exp(exp(exp(x))) / exp(exp(exp(x - exp(-exp(exp(x))))))", "oo", "+", "1"),
    c("exp(exp(x)) / exp(exp(x - exp(-exp(exp(x)))))", "oo", "+", "1"),
    c("(x*log(x)*(log(x*exp(x) - x**2))**2)/(log(log(x**2 + 2*exp(exp(3*x**3*log(x))))))", "oo", "+", "1/3"),
    c("(exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1)))) - exp(x))/x", "oo", "+", "-exp(2)"),
    c("x/log(x**(log(x**(log(2)/log(x)))))", "oo", "+", "oo"),
    c("log(x)*(log(log(x) + log(log(x))) - log(log(x)))/(log(log(x) + log(log(log(x)))))", "oo", "+", "1"),
    c("exp((log(log(x + exp(log(x)*log(log(x))))))/(log(log(log(exp(x) + x + log(x))))))", "oo", "+", "E"),
    c("exp(exp(exp(x)/(1 - 1/x))) - exp(exp(exp(x)/(1 - 1/x - log(x)**(-log(x)))))", "oo", "+", "-oo"),
    c("(exp(exp(-x/(1 + exp(-x))))*exp(-x/(1 + exp(-x/(1 + exp(-x))))) *exp(exp(-x + exp(-x/(1 + exp(-x))))))/ (exp(-x/(1 + exp(-x))))**2 - exp(x) + x", "oo", "+", "2"),
    c("exp(x)*(sin(1/x + exp(-x)) - sin(1/x + exp(-x**2)))", "oo", "+", "1"),
    c("exp(exp(x)) * (exp(sin(1/x + exp(-exp(x)))) - exp(sin(1/x)))", "oo", "+", "1"),
    c("cosh(x)", "oo", "+", "oo"),
    c("sinh(x)", "-oo", "+", "-oo"),
    c("2*cosh(x)*exp(x)", "-oo", "+", "1"),
    c("2*sinh(x)*exp(x)", "-oo", "+", "-1"),
    c("tanh(x)", "-oo", "+", "-1"),
    c("coth(x)", "-oo", "+", "-1"),
    c("x*log(x)", "0", "+", "0"),
    c("x**x", "0", "+", "1"),
    c("x - log(1 + exp(x))", "oo", "+", "0"),
    c("x - log(a + exp(x))", "oo", "+", "0"),
    c("exp(x)/(a + exp(x))", "oo", "+", "1"),
    c("(3**x + 5**x)**(1/x)", "oo", "+", "5"),
    c("(3**(1/x) + 5**(1/x))**x", "0", "+", "5"),
    c("(x + 1)**(1/log(x + 1))", "oo", "+", "E"),
    c("x/log(log(x*exp(x)))", "oo", "+", "oo"),
    c("((x**7 + x + 1)/(2**x + x**2))**(-1/x)", "oo", "+", "2"),
    c("x**2*exp(-x)", "oo", "+", "0"),
    c("(1 + 1/x)**x", "oo", "+", "E"),
    c("1/x", "0", "-", "-oo"),
    c("exp(-1/x)", "0", "+", "0"),
    c("exp(-1/x)", "0", "-", "oo"),
    c("x*sin(1/x)", "oo", "+", "1"),
    c("atan(x)", "oo", "+", "pi/2"))
  for (cs in cases) {
    e <- sympify(cs[1], locals = list(x = x, a = a))
    expect_equal(sstr(gruntz(e, x, sympify(cs[2]), cs[3])), cs[4], label = cs[1])
  }
})

test_that("limit uses the Gruntz algorithm and SymPy's signs at poles", {
  x <- Symbol("x"); y <- Symbol("y"); I <- sym_I
  expect_equal(sstr(limit((3^x + 5^x)^(1/x), x, oo)), "5")
  expect_equal(sstr(limit(((x^7 + x + 1)/(2^x + x^2))^(-1/x), x, oo)), "2")
  expect_equal(sstr(limit(2*cosh(x)*exp(x), x, -oo)), "1")
  expect_equal(sstr(limit(I*x, x, oo)), "oo*I")
  expect_equal(sstr(limit(y*I*x, x, oo)), "oo*I*sign(y)")
  expect_equal(sstr(limit(y/x, x, 0)), "oo*sign(y)")
  expect_equal(sstr(limit(floor(cos(x)), x, 0)), "0")
  expect_equal(sstr(limit(atan(x), x, oo)), "pi/2")
})

test_that("special functions are rewritten or refused, never mis-limited", {
  x <- Symbol("x")
  # SymPy: gamma -> exp(loggamma) with Stirling's series at infinity
  expect_equal(sstr(gruntz(gamma(x)*exp(-x), x, oo)), "oo")
  expect_equal(sstr(gruntz(exp(-2*x)*gamma(x)^2, x, oo)), "oo")
  expect_equal(sstr(gruntz(loggamma(x)/x, x, oo)), "oo")
  expect_equal(sstr(gruntz(1/gamma(x), x, oo)), "0")
  expect_equal(sstr(gruntz(gamma(Rational(1, 7) + 1/x), x, oo)), "gamma(1/7)")
  expect_equal(sstr(limit(factorial(x)*exp(-x), x, oo)), "oo")
  expect_equal(sstr(limit(binomial(2*x, x)*exp(-x), x, oo)), "oo")
  # Ei is made tractable through _eis, as in SymPy
  expect_equal(sstr(gruntz(x*exp(-x)*Ei(x), x, oo)), "1")
  expect_equal(sstr(limit(x*exp(-x)*Ei(x), x, oo)), "1")
  expect_equal(sstr(limit(erf(x), x, oo)), "1")
  expect_equal(sstr(limit(sqrt(x)*acosh(1 + 1/x), x, oo)), "sqrt(2)")
})
