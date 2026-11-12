# Reference values from SymPy.
x <- Symbol("x"); y <- Symbol("y"); t <- Symbol("t")
f <- Function("f"); u <- Function("u"); F <- Function("F")
X <- Function("X"); T <- Function("T"); Y <- Function("Y")
fx <- diff(f(x, y), x); fy <- diff(f(x, y), y)
s <- function(e) if (is.null(e)) "NULL" else if (is_sym(e)) sstr(e) else vapply(e, sstr, "")

test_that("classify_pde finds first order linear methods", {
  expect_equal(classify_pde(fx + fy + f(x, y)), "1st_linear_constant_coeff_homogeneous")
  expect_equal(classify_pde(fx + fy + f(x, y) - x),
               c("1st_linear_constant_coeff", "1st_linear_constant_coeff_Integral"))
  expect_equal(classify_pde(x*fx + fy), "1st_linear_variable_coeff")
  expect_length(classify_pde(fx^2 + fy), 0L)
  expect_length(classify_pde(diff(f(x, y), x, 2) + fy), 0L)
})

test_that("pdsolve solves first order linear PDEs", {
  expect_equal(s(pdsolve(fx + fy + f(x, y))), "Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))")
  expect_equal(s(pdsolve(2*fx + 3*fy + 5*f(x, y))), "Eq(f(x, y), F(3*x - 2*y)*exp(-10*x/13 - 15*y/13))")
  expect_equal(s(pdsolve(fx - fy)), "Eq(f(x, y), F(-x - y))")
  expect_equal(s(pdsolve(fx + fy + f(x, y) - x)), "Eq(f(x, y), x + F(x - y)*exp(-x/2 - y/2) - 1)")
  expect_equal(s(pdsolve(Eq(fx + fy, 1))), "Eq(f(x, y), x/2 + y/2 + F(x - y))")
  expect_equal(s(pdsolve(x*fx + y*fy + f(x, y))), "Eq(f(x, y), F(y/x)/x)")
  expect_equal(s(pdsolve(fx + x*fy)), "Eq(f(x, y), F(-x^2/2 + y))")
  expect_equal(s(pdsolve(x*fx - 2)), "Eq(f(x, y), F(y) + 2*log(x))")
  expect_equal(s(pdsolve(fy - x*y)), "Eq(f(x, y), x*y^2/2 + F(x))")
  expect_equal(s(pdsolve(x*fx + f(x, y) - x^2)), "Eq(f(x, y), (x^3/3 + F(y))/x)")
  expect_error(pdsolve(f(x, y)*fx + fy), "psolve: Cannot solve")
  sol <- pdsolve(fx + 2*fy - 3*f(x, y) - exp(x + y))
  expect_true(checkpdesol(fx + 2*fy - 3*f(x, y) - exp(x + y), sol)[[1L]])
})

test_that("checkpdesol and the chain rule for undefined functions", {
  expect_identical(checkpdesol(fx + fy + f(x, y), Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))), list(TRUE, S_Zero))
  expect_equal(sstr(checkpdesol(fx + fy + f(x, y), Eq(f(x, y), x))[[2L]]), "x + 1")
  expect_equal(sstr(diff(F(x^2), x)), "2*x*Subs(Derivative(F(_xi_1), _xi_1), _xi_1, x^2)")
  expect_equal(sstr(diff(F(x, x*y), x)),
               "y*Subs(Derivative(F(x, _xi_2), _xi_2), _xi_2, x*y) + Subs(Derivative(F(_xi_1, x*y), _xi_1), _xi_1, x)")
  expect_identical(diff(F(x - y), x) + diff(F(x - y), y), S_Zero)
})

test_that("pde_separate splits variables", {
  expect_equal(s(pde_separate(Eq(diff(u(x, t), x), diff(u(x, t), t)), u(x, t), list(X(x), T(t)))),
               c("Derivative(X(x), x)/X(x)", "Derivative(T(t), t)/T(t)"))
  expect_equal(s(pde_separate_add(Eq(diff(u(x, t), x), exp(diff(u(x, t), t))), u(x, t), list(X(x), T(t)))),
               c("Derivative(X(x), x)", "exp(Derivative(T(t), t))"))
  expect_equal(s(pde_separate_mul(Eq(diff(u(x, y), x, 2) + diff(u(x, y), y, 2), 0), u(x, y), list(X(x), Y(y)))),
               c("Derivative(X(x), (x, 2))/X(x)", "-Derivative(Y(y), (y, 2))/Y(y)"))
  expect_equal(s(pde_separate_mul(Eq(diff(u(x, t), t), x*diff(u(x, t), x, 2)), u(x, t), list(X(x), T(t)))),
               c("x*Derivative(X(x), (x, 2))/X(x)", "Derivative(T(t), t)/T(t)"))
})

test_that("p10 review fixes for PDEs", {
  u <- f(x, y)
  expect_equal(sstr(doit(pdsolve(diff(u, x) - 1, u, hint = "1st_linear_constant_coeff_Integral"))), "Eq(f(x, y), x + F(-y))")
  r <- checkpdesol(diff(u, x), list(Eq(u, y), Eq(u, x)))
  expect_true(r[[1L]][[1L]])
  expect_false(r[[2L]][[1L]])
  expect_equal(s(pdsolve((diff(u, x) + diff(u, y))^2)), "Eq(f(x, y), F(x - y))")
  sol <- pdsolve(y*x^2*u + y*diff(u, x) + diff(u, y))
  expect_true(checkpdesol(y*x^2*u + y*diff(u, x) + diff(u, y), sol)[[1L]])
  z <- Symbol("C1")
  v <- f(z, y)
  expect_equal(s(pdsolve(z*diff(v, z) + z*diff(v, y))), "Eq(f(C1, y), F(-C1 + y))")
})
