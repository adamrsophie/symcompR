x <- Symbol("x"); f <- Function("f"); y <- Function("y")
n <- Symbol("n", integer = TRUE)
D1 <- Derivative(f(x), x); D2 <- Derivative(f(x), x, 2)
s <- function(e) sstr(e)

check_ode <- function(eq, sol) {
  rhs <- sol$args[[2L]]
  res <- doit(subs(eq, f(x), rhs))
  v <- subs(res, list(x = Rational(7, 10), C1 = Rational(3, 10), C2 = Rational(11, 10)))
  expect_lt(abs(as.numeric(N(v))), 1e-9)
}

test_that("dsolve matches SymPy", {
  expect_equal(s(dsolve(D1 - f(x), f(x))), "Eq(f(x), C1*exp(x))")
  expect_equal(s(dsolve(D2 + f(x), f(x))), "Eq(f(x), C1*sin(x) + C2*cos(x))")
  expect_equal(s(dsolve(D2 + f(x) - x, f(x))), "Eq(f(x), C1*sin(x) + C2*cos(x) + x)")
  expect_equal(s(dsolve(D1 + f(x) - x, f(x))), "Eq(f(x), C1*exp(-x) + x - 1)")
  expect_equal(s(dsolve(D1 - f(x)^2, f(x))), "Eq(f(x), -1/(C1 + x))")
  expect_equal(s(dsolve(x^2 * D2 - 2 * f(x), f(x))), "Eq(f(x), C1/x + C2*x^2)")
  expect_equal(s(dsolve(D2 + f(x), f(x), ics = list(list(0, 0, 0), list(0, 1, 1)))), "Eq(f(x), sin(x))")
  for (eq in list(D2 - 3*D1 + 2*f(x), D2 + 2*D1 + f(x), D2 - f(x) - exp(x), D2 + f(x) - sin(x), D1 + f(x) - f(x)^2)) {
    check_ode(eq, dsolve(eq, f(x)))
  }
})

test_that("rsolve matches SymPy", {
  fib <- rsolve(y(n + 2) - y(n + 1) - y(n), y(n), n, init = list(list(0, 0), list(1, 1)))
  expect_equal(as.numeric(N(subs(fib, n, 10))), 55, tolerance = 1e-9)
  expect_equal(s(rsolve(y(n + 1) - y(n) - 2^n, y(n), n)), "2^n + C0")
  expect_equal(s(rsolve(y(n + 2) - 2*y(n + 1) + y(n), y(n), n)), "C0 + C1*n")
})

test_that("final review regressions stay fixed", {
  expect_equal(sstr(dsolve(Derivative(f(x), x) - sym_I * f(x), f(x))), "Eq(f(x), C1*exp(I*x))")
  sol <- dsolve((x^2 + 2*f(x)) * Derivative(f(x), x) + 2*x*f(x) + 1, f(x), ics = list(list(0, 0, 1)))
  expect_equal(sstr(sol$args[[1L]]), "f(x)")
  res <- doit(subs((x^2 + 2*f(x)) * Derivative(f(x), x) + 2*x*f(x) + 1, f(x), sol$args[[2L]]))
  expect_lt(abs(as.numeric(N(subs(res, x, Rational(1, 3))))), 1e-9)
  expect_equal(sstr(rsolve(y(n + 2) - 2*y(n + 1) + y(n), y(n), n, init = list(list(0, 0)))), "C1*n")
  x0 <- Symbol("x0"); yy <- Symbol("y")
  r <- cse(list(Sum(sin(yy), list(x0, 0, 1)), sin(yy)))
  expect_false(any(vapply(r$replacements, function(p) p[[1L]]$name, "") == "x0"))
  expect_equal(sstr(doit(subs(r$reduced[[1L]], lapply(r$replacements, identity)))), "2*sin(y)")
  k <- Symbol("k", integer = TRUE)
  expect_equal(sstr(gosper_sum(k^2 - 1, list(k, 1, 3))), "11")
  expect_equal(sstr(summation(x, list(x, 3, 1))), "-2")
  t <- Symbol("t", real = TRUE)
  expect_equal(sstr(subs(P(Die("D", 6) > t), t, 1)), "5/6")
  X <- Pareto("X", 1, Rational(1, 2))
  expect_equal(sstr(E(X)), "oo")
  expect_equal(sstr(variance(X)), "oo")
  expect_equal(lambdify(x, LambertW(x, -1))(-0.1), -3.5771520639573, tolerance = 1e-10)
  u <- Wild("u")
  expect_equal(sstr(sym_match(sin(x + yy) + x + yy, sin(u) + u)$u), "x + y")
  aw <- Wild("a", properties = list(function(z) z$type == "Symbol"))
  bw <- Wild("b", properties = list(function(z) z$type == "Symbol"))
  m <- sym_match(x + yy, aw + bw)
  expect_setequal(vapply(m, sstr, ""), c("x", "y"))
})
