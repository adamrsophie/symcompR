# Reference values from SymPy.
x <- Symbol("x"); y <- Symbol("y"); h <- Symbol("h"); a <- Symbol("a"); s <- Symbol("s")
f <- Function("f"); g <- Function("g")

test_that("differentiate_finite", {
  expect_equal(sstr(differentiate_finite(f(x)*g(x), x, points = list(x - h, x + h))),
               "-f(-h + x)*g(-h + x)/(2*h) + f(h + x)*g(h + x)/(2*h)")
  expect_equal(sstr(differentiate_finite(f(x)*g(x), x)), "-f(x - 1/2)*g(x - 1/2) + f(x + 1/2)*g(x + 1/2)")
  expect_equal(sstr(differentiate_finite(f(x, y), x, 2)), "-2*f(x, y) + f(x - 1, y) + f(x + 1, y)")
})

test_that("collect_const and kroneckersimp", {
  r3 <- sqrt(Integer(3)); r7 <- sqrt(Integer(7))
  expect_equal(sstr(collect_const(r3 + r3*(1 + sqrt(Integer(2))))), "sqrt(3)*(sqrt(2) + 2)")
  expect_equal(sstr(collect_const(r3*s + r7*s + r3 + r7)), "(sqrt(3) + sqrt(7))*(s + 1)")
  expect_equal(sstr(collect_const(2*r3 + 4*a*sqrt(Integer(5)))), "4*sqrt(5)*a + 2*sqrt(3)")
  expect_equal(sstr(collect_const(2*x - 2*y - 2*a, 2)), "2*(-a + x - y)")
  expect_equal(sstr(collect_const(2*x - 2*y - 2*a, -2)), "2*x - 2*(a + y)")
  expect_equal(latex(collect_const(2*x - 2*y - 2*a, -2)), "2 x - 2 \\left(a + y\\right)")
  i <- Symbol("i"); j <- Symbol("j")
  expect_equal(sstr(kroneckersimp(1 + KroneckerDelta(0, j)*KroneckerDelta(1, j))), "1")
  expect_equal(sstr(kroneckersimp(KroneckerDelta(i, j)*KroneckerDelta(i, 0))), "KroneckerDelta(0, i)*KroneckerDelta(i, j)")
})

test_that("p12 review fixes for simplify extras and coeff", {
  w <- Symbol("w"); z <- Symbol("z")
  expect_equal(sstr(collect_const(2*x + 2*y + z) + w), "w + z + 2*(x + y)")
  j <- Symbol("j")
  expect_identical(kroneckersimp(S_Infinity*KroneckerDelta(0, j)*KroneckerDelta(1, j)), S_NaN)
  expect_equal(sstr(differentiate_finite(x^2)), "-(x - 1/2)^2 + (x + 1/2)^2")
  expect_equal(sstr(differentiate_finite(f(x)*diff(g(x), x), x)), "(-g(x) + g(x + 1))*f(x + 1/2) - (g(x) - g(x - 1))*f(x - 1/2)")
  expect_equal(sstr(differentiate_finite(f(x), x, points = h)), "-f(-h/2 + x)/h + f(h/2 + x)/h")
  expect_equal(sstr(coeff(x - 2*y, -1)), "2*y")
  p <- Symbol("p", commutative = FALSE); q <- Symbol("q", commutative = FALSE); r <- Symbol("r", commutative = FALSE)
  expect_equal(sstr(coeff(p*q + r*q, q)), "p + r")
  expect_equal(sstr(coeff(q*p + p*q*p, q)), "1 + p")
  ar <- Symbol("a", real = TRUE)
  expect_equal(solveset(sin(x) - ar, x, Reals)$type, "ConditionSet")
  expect_identical(solveset(sin(x) - Rational("10000000000001", "10000000000000"), x, Reals), S_EmptySet)
})

test_that("p13 h2 fix: embedded derivatives follow SymPy's Subs stencils", {
  x <- Symbol("x"); h <- Symbol("h"); f <- Function("f")
  expect_equal(sstr(differentiate_finite(diff(f(x), x), x, points = list(x - 1, x, x + 1))), "-2*f(x) + f(x - 1) + f(x + 1)")
  expect_equal(sstr(differentiate_finite(f(x)*diff(f(x), x), x, points = list(x - 1, x, x + 1))),
               "(-2*f(x) + f(x - 1)/2 + 3*f(x + 1)/2)*f(x + 1)/2 - (2*f(x) - 3*f(x - 1)/2 - f(x + 1)/2)*f(x - 1)/2")
  expect_equal(sstr(differentiate_finite(f(x)*diff(f(x), x), x, points = list(x, x + 1))),
               "(-f(x) + f(x + 1))*f(x + 1) - f(x)*Derivative(f(x), x)")
  expect_equal(sstr(differentiate_finite(f(x)*diff(f(x), x), x, points = h)),
               "(-f(x)/h + f(h + x)/h)*f(h/2 + x)/h - (f(x)/h - f(-h + x)/h)*f(-h/2 + x)/h")
  expect_equal(sstr(differentiate_finite(diff(f(x), x), x, points = h)), "-2*f(x)/h^2 + f(-h + x)/h^2 + f(h + x)/h^2")
})

test_that("p14 h2 fixes: mixed embedded derivatives and a step evaluated at x0", {
  x <- Symbol("x"); y <- Symbol("y"); f <- Function("f"); g <- Function("g"); h <- Function("h")
  expect_equal(sstr(differentiate_finite(g(x)*Derivative(f(x, y), x, y), x)),
               "(-Derivative(f(x, y), y) + Derivative(f(x + 1, y), y))*g(x + 1/2) - (Derivative(f(x, y), y) - Derivative(f(x - 1, y), y))*g(x - 1/2)")
  expect_equal(sstr(differentiate_finite(f(x), x, points = h(x), x0 = 0)), "-f(-h(0)/2)/h(0) + f(h(0)/2)/h(0)")
})

test_that("p15 review fix: repeated variables merge in differentiate_finite", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(simplify(differentiate_finite(x^3*y, x, y, x, points = x))), "6*x")
})
