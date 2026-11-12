# Expected strings are SymPy outputs for the same calls.

x <- Symbol("x"); y <- Symbol("y")
pr <- function(v) paste(capture.output(print(v)), collapse = "\n")

test_that("domains print and compute as SymPy's", {
  expect_true(ZZ(3) == 3)
  expect_true(QQ(1, 2) == gmp::as.bigq(1, 2))
  expect_equal(pr(GF(7)(10)), "3 mod 7")
  expect_equal(pr(GF(7)(3) * GF(7)(5)), "1 mod 7")
  expect_equal(pr(ZZ[x]), "ZZ[x]")
  expect_equal(pr(frac_field(QQ, x)), "QQ(x)")
  expect_equal(pr(algebraic_field(QQ, sqrt(Integer(2)))), "QQ<sqrt(2)>")
  expect_equal(pr(construct_domain(list(1, Rational(1, 2)))[[1L]]), "QQ")
  expect_equal(pr(construct_domain(list(x, 1))[[1L]]), "ZZ[x]")
  expect_equal(pr(construct_domain(list(1 / x))[[1L]]), "ZZ(x)")
})

test_that("rings, fields and monomials", {
  r <- ring("a,b", ZZ)
  expect_equal(pr(r[[1L]]), "Polynomial ring in a, b over ZZ with lex order")
  expect_equal(pr((r[[2L]] + r[[3L]])^2), "a^2 + 2*a*b + b^2")
  f <- field("u,v", QQ); u <- f[[2L]]; v <- f[[3L]]
  expect_equal(pr((u^2 - v^2) / (u - v)), "u + v")
  m <- Monomial(c(1, 2))
  expect_equal(pr(m), "Monomial((1, 2))")
  expect_equal(pr(m * Monomial(c(2, 1))), "Monomial((3, 3))")
})

test_that("apart_list and assemble_partfrac_list", {
  pl <- apart_list((2 * x^3 - 2 * x) / (x^2 - 2 * x + 1), x)
  expect_equal(sstr(pl[[1L]]), "1")
  expect_equal(pr(pl[[2L]]), "Poly(2*x + 4, x, domain='ZZ')")
  expect_length(pl[[3L]], 1L)
  expect_equal(sstr(assemble_partfrac_list(pl)), "2*x + 4 + 4/(x - 1)")
  pl <- apart_list((x^2 + 1) / (x^3 - x), x)
  expect_length(pl[[3L]], 1L)
  expect_equal(pr(pl[[3L]][[1L]][[1L]]), "Poly(_w^3 - _w, _w, domain='ZZ')")
  expect_equal(sstr(assemble_partfrac_list(pl)), "1/(x + 1) + 1/(x - 1) - 1/x")
})

test_that("sqf_list gives one factor per multiplicity", {
  fmt <- function(f) { s <- sqf_list(expand(f))
    paste(sstr(s$coeff), paste(vapply(s$factors, function(p) paste0("(", sstr(p[[1L]]), ", ", p[[2L]], ")"), ""), collapse = " ")) }
  expect_equal(fmt(x^3 - x), "1 (x^3 - x, 1)")
  expect_equal(fmt(x^5 - x^3), "1 (x^2 - 1, 1) (x, 3)")
  expect_equal(fmt(x * y * (x + y)^2), "1 (x*y, 1) (x + y, 2)")
  expect_equal(fmt(2 * x^3 * (x + 1)), "2 (x + 1, 1) (x, 3)")
  expect_equal(fmt((x^2 + y) * (x + 1)^2 * y^3), "1 (x^2 + y, 1) (x + 1, 2) (y, 3)")
  expect_equal(fmt(3 * x^4 - 3 * x^2), "3 (x^2 - 1, 1) (x, 2)")
})

test_that("SymPy exceptions are classed conditions", {
  e <- tryCatch(sympify("x +* 1"), error = function(err) err)
  expect_s3_class(e, "SympifyError")
  expect_s3_class(e, "sympy_error")
  e <- tryCatch(stop(MultivariatePolynomialError("two variables")), PolynomialError = function(err) err)
  expect_s3_class(e, "BasePolynomialError")
  expect_equal(conditionMessage(e), "two variables")
  e <- tryCatch(series(loggamma(x - 1), x, 0, 2), PoleError = function(err) "pole")
  expect_equal(e, "pole")
})

test_that("utility helpers", {
  expect_equal(filldedent("
    A long docstring line
    continues here."), "\nA long docstring line continues here.")
  expect_equal(capture(function() print(Symbol("x")^2)), "x^2\n")
  d <- threaded(function(e) diff(e, x))
  expect_equal(vapply(d(list(x^2, sin(x))), sstr, ""), c("2*x", "cos(x)"))
  expect_equal(sstr(d(x^3 + x)), "3*x^2 + 1")
  calls <- 0
  sq <- memoize_property(function(v) { calls <<- calls + 1; v^2 })
  sq(x); sq(x)
  expect_equal(calls, 1)
  tm <- timed(function() NULL, limit = 10)
  expect_true(tm$number >= 1 && tm$unit %in% c("s", "ms", "us", "ns"))
})
