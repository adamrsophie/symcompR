# Expected values are SymPy's outputs for the same calls.

x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); I <- sym_I

test_that("cancel works over the Gaussian integers and rationals", {
  expect_equal(sstr(cancel((x + I) / (x^2 + 1))), "1/(x - I)")
  expect_equal(sstr(cancel((1 + I) / (2 * x))), "I*(1 - I)/(2*x)")
  expect_equal(sstr(cancel((x^2 + 1) / (x - I))), "x + I")
  expect_equal(sstr(cancel((x^2 + 1) / (2 * x - 2 * I))), "x/2 + I/2")
  expect_equal(sstr(cancel((x / 2 + I / 3) / (x^2 + Rational(4, 9)))), "3/(6*x - 4*I)")
  expect_equal(sstr(cancel(1 / (x - I / 2))), "2/(2*x - I)")
})

test_that("expand takes the fraction of a product as Mul._eval_expand_mul", {
  expect_equal(sstr(expand((sqrt(z) / 2)^(2 - 2 * a))), "2^(2*a)*z^(1 - a)/4")
  expect_equal(sstr(expand(exp(-x) * (y + 1))), "y*exp(-x) + exp(-x)")
  fr <- fraction(exp(-x) * y)
  expect_equal(sstr(fr[[2L]]), "exp(x)")
  fr <- fraction(y / x^2, exact = TRUE)
  expect_equal(sstr(fr[[2L]]), "x^2")
})

test_that("products rebuild factor exponents as Mul.flatten does", {
  expect_equal(sstr(powsimp(Rational(1, 2)^(2 - 2 * a) * y, combine = "base")), "(1/2)^(2 - 2*a)*y")
  expect_equal(sstr(powdenest((sqrt(z) / 2)^(2 - 2 * a), polar = TRUE)), "(1/2)^(2 - 2*a)*z^(1 - a)")
  expect_equal(sstr(simplify((sqrt(z) / 2)^(2 - 2 * a))), "2^(2*a - 2)*z^(1 - a)")
})

test_that("negated powers with negative exponents print as fractions", {
  expect_equal(sstr(-1 / (sqrt(z) + 1)^(2 * a)), "-1/(sqrt(z) + 1)^(2*a)")
  expect_equal(sstr(z - 1 / (sqrt(z) + 1)^(2 * a)), "z - 1/(sqrt(z) + 1)^(2*a)")
  expect_equal(sstr(-y^(-a - 1)), "-y^(-a - 1)")
})

test_that("ordered() counts a Float as half a node", {
  expect_equal(sstr(hyper(list(), list(Float(2) / 3, Rational(1, 3)), z)), "hyper((), (0.666666666666667, 1/3), z)")
})

test_that("sums and integrals that go through apart", {
  i <- Symbol("i", integer = TRUE, positive = TRUE)
  expect_equal(sstr(summation(1 / (i * (i + 1)), list(i, 1, oo))), "1")
  expect_equal(sstr(sym_integrate(x / (1 + x), x)), "x - log(x + 1)")
})

test_that("binomial, harmonic, bernoulli and nsimplify evaluate as in SymPy", {
  n <- Symbol("n", integer = TRUE, positive = TRUE); m <- Symbol("m")
  expect_equal(sstr(binomial(n, n + 1)), "0")
  expect_equal(sstr(binomial(Rational(1, 2), Integer(3))), "1/16")
  expect_equal(sstr(binomial(m, Integer(-1))), "0")
  expect_equal(sstr(harmonic(oo, Integer(2))), "pi^2/6")
  expect_equal(sstr(harmonic(m, Integer(-2))), "m/6 + (m + 1)^3/3 - (m + 1)^2/2 + 1/6")
  expect_equal(sstr(bernoulli(Integer(3), m + 1)), "m/2 + (m + 1)^3 - 3*(m + 1)^2/2 + 1/2")
  expect_equal(sstr(nsimplify(Float(0.99999))), "99999/100000")
})

test_that("powers of powers, signs of sums and commutativity follow SymPy", {
  k <- Symbol("k", integer = TRUE); p <- Symbol("p", positive = TRUE)
  expect_equal(sstr((Integer(-2)^y)^(k + 1)), "(-2)^(y*(k + 1))")
  expect_equal(sstr((p^2)^y), "p^(2*y)")
  d <- Dummy("d", integer = TRUE, positive = TRUE)
  expect_true(isTRUE(ask_fact(d - 1, "nonnegative")))
  expect_equal(sstr(Ge(d, 1)), "True")
  B <- Symbol("B", commutative = FALSE)
  expect_equal(sstr(conjugate(B) + conjugate(B)^2 + conjugate(B)^3), "conjugate(B) + conjugate(B)^2 + conjugate(B)^3")
  f <- Function("f")
  expect_equal(sstr(f(Integer(10)) + f(Integer(2)) + f(Integer(3))), "f(2) + f(3) + f(10)")
  m <- Symbol("m")
  expect_equal(sstr(apart(-2 * a^2 * m + 2 * a * m + 4 * m, m)), "-m*(2*a^2 - 2*a - 4)")
})
