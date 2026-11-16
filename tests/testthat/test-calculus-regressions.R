# Regression cases for derivatives, integrals, sums, transforms and holonomic functions.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("calculus helpers follow SymPy", {
  x <- Symbol("x"); t <- Symbol("t")
  expect_identical(sstr(pade(exp(x), x, 2, 2)), "(x^2/4 + 3*x/2 + 3)/(x^2/4 - 3*x/2 + 3)")
  expect_false(is_monotonic(x^3))
  expect_true(is_monotonic(exp(x)))
  expect_identical(sstr(line_integrate(x*Symbol("y"), Curve(list(t, 2*t), list(t, 0, 1)), list(x, Symbol("y")))), "2*sqrt(5)/3")
  out <- utils::capture.output(textplot(x^2, 0, 2, W = 30, H = 8))
  expect_identical(out[4L], "2.28571 |---------------------...------")
  expect_identical(out[9L], "         0              1.03448       2")
})

test_that("Integral differentiates under the integral sign", {
  x <- Symbol("x"); y <- Symbol("y"); t <- Symbol("t"); f <- Function("f")
  expect_equal(sstr(diff(Integral(x + y, y, list(y, 1, x)), x)),
               "Integral(x + y, (y, x)) + Integral(1, y, (y, 1, x))")
  expect_equal(sstr(diff(Integral(x + y, y, list(y, 1, x)), y)), "0")
  expect_equal(sstr(diff(Integral(sin(t * x), list(x, 0, t)), t)),
               "sin(t^2) + Integral(x*cos(t*x), (x, 0, t))")
  expect_equal(sstr(diff(Integral(f(x) * y, x), y)), "Integral(f(x), x)")
  expect_equal(sstr(diff(Integral(x, x), x)), "x")
})

test_that("an undone integral keeps its limit and factors its integrand", {
  skip_on_cran()
  x <- Symbol("x"); n <- Symbol("n"); y <- Dummy("y"); f <- Function("f")
  expect_equal(sstr(doit(Integral(exp(Integral(y^2 + n * y, y)), list(y, f(x))))),
               "Integral(exp(_y^3/3)*exp(_y^2*n/2), (_y, f(x)))")
  expect_equal(sstr(sym_integrate(exp(x^3 / 3 + x^2 * n / 2), x)), "Integral(exp(x^3/3)*exp(n*x^2/2), x)")
})

test_that("Integral, tryit, particular solve and noncommutative subs follow SymPy", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_equal(sstr(Integral(sympify("nan"), x)), "nan")
  expect_equal(sstr(Integral(Integral(x, x), x)), "Integral(x, x, x)")
  expect_error(Integral(1), "specify dummy variables")
  h <- tryit(function(e) stop("TypeError: bad"), "ValueError")
  expect_error(h(x), "TypeError")
  s <- solve(list(x + y - 1), list(x, y), particular = TRUE)
  expect_equal(vapply(s, sstr, ""), c(x = "0", y = "1"))
  A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE); C <- Symbol("C", commutative = FALSE)
  expect_equal(sstr(subs(A * B * C, A * B, z)), "z*C")
  expect_equal(sstr(subs(A * B * A * B, A * B, z)), "z^2")
  expect_equal(sstr(subs(A * B^2 * C, A * B, z)), "z*B*C")
})

test_that("holonomic powers, diffgeom coordinates and nsolve follow SymPy", {
  skip_on_cran()
  x <- Symbol("x")
  h <- expr_to_holonomic(exp(x), x)
  expect_equal(sstr(to_expr_holonomic(h^Rational(1, 2))), "exp(x/2)")
  expect_equal(sstr(evalf_holonomic(h, 0.02)), "1.05127109375000")
  C <- CoordSystem("C", Patch("P", Manifold("M", 1)), symbols = list(x))
  expect_equal(sstr(diff(coord_symbols(C)[[1]], x)), "0")
  expect_equal(sstr(nsolve(list(x - 1, 2 * x - 2), list(x), list(0))[[1]]), "1.00000000000000")
  expect_equal(sstr(nsolve(list(x / Integer(10)^400 - 1), list(x), list(0))[[1]]), "1.00000000000000e+400")
})

test_that("Integral as_sum, principal_value and transform follow SymPy", {
  skip_on_cran()
  x <- Symbol("x"); u <- Symbol("u"); a <- Symbol("a")
  I2 <- Integral(x^2, list(x, 0, 2))
  expect_equal(sstr(integral_as_sum(I2, 4)), "21/8")
  expect_equal(sstr(integral_as_sum(I2, 4, method = "trapezoid")), "11/4")
  expect_equal(sstr(integral_as_sum(I2, 4, method = "left")), "7/4")
  expect_equal(sstr(principal_value(Integral(1/x, list(x, -1, 1)))), "0")
  expect_equal(sstr(principal_value(Integral(1/(x - 1), list(x, 0, 3)))), "log(2)")
  expect_equal(sstr(integral_transform(Integral(x, list(x, 0, 1)), x, 2*u)), "Integral(4*u, (u, 0, 1/2))")
  expect_equal(sstr(integral_transform(Integral(x*exp(x^2), list(x, 0, 1)), x^2, u)), "Integral(exp(u)/2, (u, 0, 1))")
  expect_equal(sstr(integral_transform(Integral(x, list(x, a, 1)), x, 1 - u)), "Integral(u - 1, (u, 1 - a, 0))")
  expect_equal(sstr(integral_transform(Integral(1/x, list(x, 1, 2)), log(x), u)), "Integral(1, (u, 0, log(2)))")
})

test_that("AccumBounds, hankel_transform and ordinals follow SymPy", {
  skip_on_cran()
  r <- Symbol("r", positive = TRUE)
  expect_equal(sstr(Float(0) / AccumBounds(-1, 1)), "0")
  expect_equal(sstr(1 / AccumBounds(Float(0), 2)), "AccumBounds(1/2, oo)")
  expect_equal(sstr(log(AccumBounds(Float(0), 2))), "AccumBounds(-oo, log(2))")
  expect_equal(sstr(hankel_transform(r^(-4), r, -1, 1)), "-oo")
  expect_equal(sstr(hankel_transform(r^(-4), r, 0, 1)), "0")
  w <- ord_omega
  expect_equal(format(Ordinal()), "ord0")
  expect_equal(format(w + Integer(1)), "w + 1")
  expect_equal(format(Integer(2) + w), "w")
  expect_true(OmegaPower(1, 1) < OmegaPower(2, 1))
  expect_true(OmegaPower(1, 1) == OmegaPower(1, 1))
})

test_that("sums and integrals that go through apart", {
  skip_on_cran()
  x <- Symbol("x")
  i <- Symbol("i", integer = TRUE, positive = TRUE)
  expect_equal(sstr(summation(1 / (i * (i + 1)), list(i, 1, oo))), "1")
  expect_equal(sstr(sym_integrate(x / (1 + x), x)), "x - log(x + 1)")
})

test_that("a non elementary Risch result prints", {
  skip_on_cran()
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

test_that("a number as integration variable is rejected", {
  x <- Symbol("x")
  expect_error(sym_integrate(x^2, x, 0, 1), "Invalid limits")
  expect_equal(sstr(sym_integrate(x^2, list(x, 0, 1))), "1/3")
})

test_that("a non integer power of a higher order holonomic function raises", {
  skip_on_cran()
  x <- Symbol("x")
  h <- expr_to_holonomic(sin(x), x)
  expect_error(h^Rational(1, 2), "TypeError")
  expect_equal(sstr(to_expr_holonomic(h^2)), "sin(x)^2")
})
