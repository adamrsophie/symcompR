# Fixes found by comparing dsolve and solveset with SymPy; the reference
# values are SymPy's outputs.

test_that("exp and power factors read as powers of the generators", {
  x <- Symbol("x"); n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(Poly(exp(Integer(6)) + S_E, S_E)), "Poly(E^6 + E, E, domain='ZZ')")
  expect_equal(sstr(Poly(exp(Integer(6)) * x + 1, x)), "Poly(exp(6)*x + 1, x, domain='ZZ[E]')")
  expect_equal(sstr(as_expr(Poly(exp(2 * x) + exp(x), exp(x)))), "exp(2*x) + exp(x)")
  expect_equal(degree(Poly(2^(2 * n) + 2^n, 2^n)), 2)
})

test_that("exp(y) is never zero for finite y and an integral of zero is zero", {
  y <- Symbol("y"); x <- Symbol("x")
  expect_false(ask_fact(exp(y), "zero"))
  expect_true(ask_fact(Integral(0, x), "zero"))
  expect_equal(sstr(exp(-Integral(0, x))), "1")
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
  x <- Symbol("x"); n <- Symbol("n"); y <- Dummy("y"); f <- Function("f")
  expect_equal(sstr(doit(Integral(exp(Integral(y^2 + n * y, y)), list(y, f(x))))),
               "Integral(exp(_y^3/3)*exp(_y^2*n/2), (_y, f(x)))")
  expect_equal(sstr(sym_integrate(exp(x^3 / 3 + x^2 * n / 2), x)), "Integral(exp(x^3/3)*exp(n*x^2/2), x)")
})

test_that("subs replaces a function that is an integration variable", {
  x <- Symbol("x"); yy <- Symbol("y"); f <- Function("f")
  expect_equal(sstr(subs(Integral(x * f(x), f(x)), f(x), yy)), "Integral(x*y, y)")
})

test_that("dsolve agrees with SymPy on the fixed cases", {
  x <- Symbol("x"); f <- Function("f")
  expect_equal(sstr(dsolve(x * diff(f(x), x) - f(x) - x^2 * f(x)^2, f(x))),
               "Eq(f(x), -3*x/(C1 + x^3))")
  eq <- diff(f(x), x) + (3 + 2 * f(x)) / (x + 3)
  expect_equal(sstr(dsolve(eq, f(x), hint = "linear_coefficients")), "Eq(f(x), C1/(x^2 + 6*x + 9) - 3/2)")
  expect_equal(sstr(dsolve(eq, f(x), hint = "linear_coefficients_Integral", simplify = FALSE)),
               "Eq(log(f(x) + 3/2), log(C1) + Integral(-2/(3*_u1), (_u1, (x + 3)/(f(x) + 3/2))))")
  expect_equal(sstr(dsolve(diff(f(x), x), f(x), hint = "1st_linear_Integral")), "Eq(f(x), C1 + Integral(0, x))")
  expect_equal(sstr(dsolve(x * diff(f(x), x) - f(x) - x^2 * f(x)^2, f(x), hint = "1st_exact_Integral")),
               paste0("Eq(Subs(Integral(_y*(-_y*x^2 - 1)*exp(Integral(-2/_y, _y)), x) + ",
                      "Integral(x*exp(Integral(-2/_y, _y)) - Integral(exp(Integral(-2/_y, _y)), x), _y), _y, f(x)), C1)"))
})

test_that("Poly reads I as a coefficient in ZZ_I and QQ_I", {
  x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a")
  expect_equal(sstr(Poly(S_I * x + 1, x)), "Poly(I*x + 1, x, domain='ZZ_I')")
  expect_equal(sstr(Poly(S_I * x + y)), "Poly(I*x + y, x, y, domain='ZZ_I')")
  expect_equal(sstr(Poly(S_I * x / 2 + 1)), "Poly(I/2*x + 1, x, domain='QQ_I')")
  expect_equal(sstr(Poly(S_I * a * x + 1, x)), "Poly(I*a*x + 1, x, domain='ZZ_I[a]')")
  expect_equal(sstr(Poly(S_I * sqrt(Integer(3)) * x + 1, x)), "Poly(sqrt(3)*I*x + 1, x, domain='EX')")
})

test_that("terms_gcd follows SymPy's options", {
  x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a"); I <- S_I
  expect_equal(sstr(terms_gcd(2 * x + 2 * I * x)), "x*(2 + 2*I)")
  expect_equal(sstr(terms_gcd(I * x + I * y)), "I*x + I*y")
  expect_equal(sstr(terms_gcd(x^2 * I + x)), "x*(I*x + 1)")
  expect_equal(sstr(terms_gcd((1 + I) * x + (1 + I) * y)), "(1 + I)*(x + y)")
  expect_equal(sstr(terms_gcd(a * x + a * y, x)), "a*(x + y)")
  expect_equal(sstr(terms_gcd(x / 2 + x^2 / 4)), "x*(x + 2)/4")
  expect_equal(sstr(terms_gcd(x * (2 * y + 2 * x), deep = TRUE)), "2*x*(x + y)")
  expect_equal(sstr(terms_gcd(Eq(2 * x + 2 * y, 4 * x * y))), "Eq(2*(x + y), 4*x*y)")
  expect_equal(sstr(terms_gcd(x / 2 + y / 3, clear = FALSE)), "(3*x + 2*y)/6")
  expect_equal(sstr(terms_gcd(x / 2 + 1, clear = FALSE)), "x/2 + 1")
  expect_error(terms_gcd(Lt(x, y)), "Inequalities")
})

test_that("constantsimp factors a Gaussian rational content", {
  x <- Symbol("x"); f <- Function("f")
  s <- dsolve(1 + f(x)^3 + x^2 * f(x)^2 * diff(f(x), x), f(x), hint = "almost_linear")
  expect_equal(vapply(s, sstr, ""),
               c("Eq(f(x), (C1*exp(3/x) - 1)^(1/3))", "Eq(f(x), (-1 - sqrt(3)*I)*(C1*exp(3/x) - 1)^(1/3)/2)",
                 "Eq(f(x), (-1 + sqrt(3)*I)*(C1*exp(3/x) - 1)^(1/3)/2)"))
})
