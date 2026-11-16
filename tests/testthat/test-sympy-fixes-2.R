# Fixes found by comparing with SymPy; reference values are SymPy's.

test_that("a non elementary Risch result prints", {
  x <- Symbol("x")
  r <- risch_integrate((1 + 2*x^2 + x^4 + 2*x^3*exp(2*x^2))/(x^4*exp(x^2) + 2*x^2*exp(x^2) + exp(x^2)), x)
  expect_equal(sstr(r), "Integral(exp(-x^2), x) + exp(x^2)/(x^2 + 1)")
})

test_that("a sum evaluates an unevaluated numeric power", {
  x <- Symbol("x")
  s <- Pow(2, Rational(-1, 2), evaluate = FALSE)
  expect_equal(sstr(s - sqrt(Integer(2)) / 2), "0")
  expect_equal(sstr(s + x), "x + sqrt(2)/2")
  expect_equal(sstr(Pow(2, -3, evaluate = FALSE) + x), "x + 1/8")
  expect_equal(sstr(s + s), "sqrt(2)")
})

test_that("N keeps a symbolic evalf rewrite", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(N(logaddexp(x, y))), "log(exp(x) + exp(y))")
})

test_that("equations between noncommutative objects stay undecided", {
  B <- Symbol("B", commutative = FALSE)
  expect_equal(sstr(Eq(B, 1)), "Eq(B, 1)")
  expect_equal(sstr(Ne(B, 1)), "Ne(B, 1)")
  expect_equal(sstr(Eq(B^2, 1)), "Eq(B^2, 1)")
})

test_that("a number as integration variable is rejected", {
  x <- Symbol("x")
  expect_error(sym_integrate(x^2, x, 0, 1), "Invalid limits")
  expect_equal(sstr(sym_integrate(x^2, list(x, 0, 1))), "1/3")
})

test_that("Poly floor division, remainder and equality", {
  x <- Symbol("x"); y <- Symbol("y")
  p <- Poly(x^2 + 2, x); q <- Poly(x, x)
  expect_equal(sstr(p %/% q), "Poly(x, x, domain='ZZ')")
  expect_equal(sstr(p %% q), "Poly(2, x, domain='ZZ')")
  expect_false(Poly(x + y, x) == Poly(x + y, x, y))
  expect_true(Poly(x + y, x) == Poly(x + y, x))
})

test_that("factorials are exact", {
  expect_equal(sstr(sym_factorial(25)), "15511210043330985984000000")
  expect_equal(sstr(factorial2(Integer(-3))), "-1")
  expect_equal(sstr(factorial2(Integer(-5))), "1/3")
  expect_error(factorial2(Integer(-4)), "negative odd")
})

test_that("invalid input is refused", {
  expect_error(Integer(NA), "ValueError")
  expect_error(Integer(Inf), "ValueError")
  expect_error(normalized(Matrix(list(3, 4))), "TRUE or FALSE")
  x <- Symbol("x")
  expect_warning(m <- as.matrix(Matrix(list(list(1, 2), list(3, x)))), "not numbers")
  expect_equal(m[1, 2], 2)
})

test_that("Poly() of a Poly applies requested generators and options", {
  x <- Symbol("x"); y <- Symbol("y")
  p <- Poly(x + 2, x)
  expect_identical(Poly(p, x), p)
  expect_equal(sstr(Poly(p, x, modulus = 2)), "Poly(x, x, modulus=2)")
  expect_equal(sstr(Poly(Poly(x * y + 1, x, y), y, x)), "Poly(y*x + 1, y, x, domain='ZZ')")
})

test_that("arithmetic between Polys unifies their generators", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(Poly(x, x) + Poly(y, y)), "Poly(x + y, x, y, domain='ZZ')")
  expect_equal(sstr(Poly(y, y) * Poly(x, x)), "Poly(y*x, y, x, domain='ZZ')")
})

test_that("a non integer power of a higher order holonomic function raises", {
  x <- Symbol("x")
  h <- expr_to_holonomic(sin(x), x)
  expect_error(h^Rational(1, 2), "TypeError")
  expect_equal(sstr(to_expr_holonomic(h^2)), "sin(x)^2")
})

test_that("preview() raises when the output file cannot be written", {
  x <- Symbol("x")
  bad <- file.path(tempdir(), "no-such-directory", "out.tex")
  expect_error(preview(x, output = "tex", viewer = "file", filename = bad), "could not write")
})
