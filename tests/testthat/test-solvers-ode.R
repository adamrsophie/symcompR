x <- Symbol("x")
t <- Symbol("t")
f <- Function("f")
g <- Function("g")
s <- function(e) sstr(e)

test_that("dsolve solves linear constant coefficient ODEs", {
  eq <- Derivative(f(x), x, 2) - 3 * Derivative(f(x), x) + 2 * f(x)
  expect_equal(s(dsolve(eq)), "Eq(f(x), (C1 + C2*exp(x))*exp(x))")
  expect_equal(s(dsolve(Derivative(f(x), x) - f(x))), "Eq(f(x), C1*exp(x))")
})

test_that("dsolve applies initial conditions", {
  ics <- .pydict(list(f(0)), list(1))
  expect_equal(s(dsolve(Derivative(f(x), x) - f(x), f(x), ics = ics)), "Eq(f(x), exp(x))")
})

test_that("classify_ode lists the matching hints", {
  h <- classify_ode(Derivative(f(x), x) - f(x), f(x))
  expect_true("separable" %in% h)
  expect_true("1st_linear" %in% h)
})

test_that("ode_order and checkodesol work", {
  expect_equal(ode_order(Derivative(f(x), x, 3) + f(x), f(x)), 3L)
  chk <- checkodesol(Derivative(f(x), x) - f(x), Eq(f(x), Symbol("C1") * exp(x)))
  expect_true(isTRUE(chk[[1L]]))
})

test_that("dsolve solves a linear system", {
  eqs <- list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), -f(t)))
  expect_equal(vapply(dsolve(eqs), s, ""),
               c("Eq(f(t), -C1*sin(t) + C2*cos(t))", "Eq(g(t), -C1*cos(t) - C2*sin(t))"))
})

test_that("linear_ode_to_matrix and canonical_odes work", {
  eqs <- list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), f(t)))
  r <- linear_ode_to_matrix(eqs, list(f(t), g(t)), t, 1)
  expect_equal(s(r[[2L]]), "Matrix([[0], [0]])")
  ce <- canonical_odes(eqs, list(f(t), g(t)), t)
  expect_length(ce, 1L)
})

test_that("homogeneous_order finds the order", {
  y <- Symbol("y")
  expect_equal(s(homogeneous_order(x^2 + x * y, x, y)), "2")
})
