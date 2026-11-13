# Core fixes behind solveset and ODE results (expected values from SymPy)

test_that("arg and atan2 follow SymPy's eval", {
  expect_equal(sstr(arg(1 + sym_I)), "pi/4")
  expect_equal(sstr(arg(-1 - sym_I)), "-3*pi/4")
  expect_equal(sstr(atan2(as_sym(1), as_sym(-1))), "3*pi/4")
  expect_equal(sstr(atan2(as_sym(0), as_sym(-1))), "pi")
  expect_equal(sstr(atan2(as_sym(-1), as_sym(0))), "-pi/2")
})

test_that("Eq decides infinite sides as is_eq does", {
  z <- Symbol("z", real = TRUE)
  w <- Symbol("w", extended_real = TRUE)
  expect_identical(Eq(S_NegativeInfinity, z), S_false)
  expect_identical(Eq(S_NegativeInfinity, z + 1), S_false)
  expect_identical(Eq(S_Infinity, S_NegativeInfinity), S_false)
  expect_equal(sstr(Eq(S_Infinity, w)), "Eq(oo, w)")
})

test_that("Pow is imaginary and conjugates as SymPy's Pow", {
  expect_true(ask_fact(sqrt(as_sym(-2)), "imaginary"))
  r <- sqrt(-sqrt(as_sym(2)) - Rational(1, 2))
  expect_equal(sstr(conjugate(r)), "-I*sqrt(1/2 + sqrt(2))")
  expect_equal(sstr(Abs(sqrt(as_sym(2)) / 2 + r)), "sqrt(1 + sqrt(2))")
})

test_that("log of a finite nonzero argument is finite", {
  d <- Dummy("d", finite = TRUE)
  expect_true(ask_fact(log(as_sym(3)^(-d)), "finite"))
})

test_that("Float coefficients factor and integrate over RR", {
  x <- Symbol("x"); y <- Symbol("y"); u <- Symbol("u")
  e <- Float(1.0) * x - Float(1.0) * log(y - Float(2.5)) + Float(1.0) * log(y + Float(2.5))
  expect_equal(sstr(factor_sym(e)), "1.0*(x - log(y - 2.5) + log(y + 2.5))")
  expect_equal(sstr(gcd_terms(Float(1.0) * x - Float(1.0) * y)), "1.0*(x - y)")
  expect_equal(sstr(factor_terms(Float(2.0) * x + Float(4.0) * y)), "2.0*x + 4.0*y")
  expect_equal(sstr(integrate(1 / (Float(51.05) * u * y - Float(1.0)), u)),
               "0.019588638589618*log(51.05*u*y - 1.0)/y")
})

test_that("ODE terms_gcd keeps Float and Gaussian coefficients as SymPy", {
  a <- Symbol("a"); y <- Symbol("y"); C1 <- Symbol("C1"); x <- Symbol("x")
  expect_equal(sstr(.cs_terms_gcd(C1 * y - Float(51.05) * a * y, clear = FALSE, deep = TRUE)),
               "y*(1.0*C1 - 51.05*a)")
  expect_equal(sstr(.cs_terms_gcd(2 * sym_I * x + 4 * sym_I * x^2, clear = FALSE, deep = TRUE)),
               "2*x*(2*I*x + I)")
  f <- Function("f")(x); C3 <- Symbol("C3")
  expect_equal(sstr(constantsimp(Eq(f, (-C3 * sym_I + sym_I * x^4) / (C3 + x^4)), list(C3))),
               "Eq(f(x), (-I*C3 + I*x^4)/(C3 + x^4))")
})

test_that("keep_coeff and Add printing follow SymPy", {
  expect_equal(sstr(4 * (Rational(1, 4) - sqrt(as_sym(17)) / 4)), "1 - sqrt(17)")
  s <- sin(S_Pi / 9); c <- cos(S_Pi / 9); b <- s^2 + c^2; r3 <- sqrt(as_sym(3))
  expect_equal(sstr(-r3 * s / 2 - r3 * s / (2 * b) - c / 2 - c / (2 * b)),
               "-cos(pi/9)/2 - cos(pi/9)/(2*sin(pi/9)^2 + 2*cos(pi/9)^2) - sqrt(3)*sin(pi/9)/2 - sqrt(3)*sin(pi/9)/(2*sin(pi/9)^2 + 2*cos(pi/9)^2)")
})

test_that("solve skips a non-polynomial generator", {
  C1 <- Symbol("C1"); x <- Symbol("x")
  expect_equal(vapply(solve(sqrt(C1 + 1) - x, C1), sstr, ""), "x^2 - 1")
})

test_that("solveset of exponentials and Float image sets", {
  y <- Symbol("y")
  expect_equal(sstr(solveset(4 * as_sym(2)^y * as_sym(2)^(log(as_sym(3)^(-y)) / log(as_sym(3))) - as_sym(4)^y, y)), "{1}")
})
