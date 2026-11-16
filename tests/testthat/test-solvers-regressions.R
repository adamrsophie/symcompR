# Regression cases for equation and ODE solvers.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("dsolve agrees with SymPy on the fixed cases", {
  skip_on_cran()
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

test_that("constantsimp factors a Gaussian rational content", {
  skip_on_cran()
  x <- Symbol("x"); f <- Function("f")
  s <- dsolve(1 + f(x)^3 + x^2 * f(x)^2 * diff(f(x), x), f(x), hint = "almost_linear")
  expect_equal(vapply(s, sstr, ""),
               c("Eq(f(x), (C1*exp(3/x) - 1)^(1/3))", "Eq(f(x), (-1 - sqrt(3)*I)*(C1*exp(3/x) - 1)^(1/3)/2)",
                 "Eq(f(x), (-1 + sqrt(3)*I)*(C1*exp(3/x) - 1)^(1/3)/2)"))
})

test_that("an equation keeps lhs before rhs when the functions are extracted", {
  x <- Symbol("x"); f <- Function("f"); g <- Function("g")
  eqs <- list(Eq(Derivative(f(x), list(x, 2)), -2 * g(x) + 2 * x * diff(g(x), x)),
              Eq(Derivative(g(x), list(x, 2)), 2 * f(x) - 2 * x * diff(f(x), x)))
  expect_equal(vapply(symcompR:::.ode_extract_funcs(eqs), sstr, ""), c("f(x)", "g(x)"))
})
