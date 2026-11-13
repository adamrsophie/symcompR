# Expected values are SymPy's outputs for the same calls (sympy/series/tests/test_kauers.py, test_approximants.py and acceleration);
# the rows come from the battery that compares this port with SymPy.

test_that("kauers: finite_diff", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(z^2 - 2*z + 3, z) })), "2*z - 1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(w^2 + 3*w - 2, w) })), "2*w + 4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x^2 - 2*x + 3, x, 2) })), "4*x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(n^2 - 2*n + 3, n, 3) })), "6*n + 3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x^3, x) })), "3*x^2 + 3*x + 1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x*y^2, x, y) })), "y^3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(2^x, x) })), "2^x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x^4 - x, x, -1) })), "-4*x^3 + 6*x^2 - 4*x + 2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff((x + 1)^3, x) })), "3*x^2 + 9*x + 7")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(1/x, x) })), "1/(x + 1) - 1/x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x^2, x, Rational(1, 2)) })), "x + 1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(5, x) })), "ERR AttributeError: 'int' object has no attribute 'expand'")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff(x*sin(x), x) })), "-x*sin(x) + x*sin(x + 1) + sin(x + 1)")
})

test_that("kauers: finite_diff_kauers", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(x^2, list(x, 1, n))) })), "(n + 1)^2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(y, list(y, 1, m))) })), "m + 1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum((x*y), list(x, 1, m), list(y, 1, n))) })), "(m + 1)*(n + 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum((x*y^2), list(x, 1, m), list(y, 1, n))) })), "(m + 1)*(n + 1)^2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(1/k, list(k, 1, n))) })), "1/(n + 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(sin(k), list(k, 0, n))) })), "sin(n + 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(2^k*x, list(k, 0, m))) })), "2^(m + 1)*x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(k^3, list(k, 1, 2*n))) })), "(2*n + 1)^3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(x^k, list(k, 0, n))) })), "x^(n + 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); finite_diff_kauers(Sum(k*z, list(k, 1, n), list(z, 0, m))) })), "(m + 1)*(n + 1)")
})

test_that("acceleration: richardson", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson((1 + 1/n)^n, n, 10, 3) })), "64760917612327084418953891207076443/23824634392675154069619793920000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(1/n, n, 1, 2) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(1 - 1/n^2, n, 5, 4) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(n/(n + 1), n, 3, 0) })), "3/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(n/(n + 1), n, 3, 1) })), "19/20")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(x + 1/n, n, 2, 2) })), "x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); richardson(Sum(1/factorial(k), list(k, 0, n)), n, 2, 2) })), "8/3")
})

test_that("acceleration: shanks", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); shanks(1 - 1/n, n, 4) })), "7/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); m <- Symbol("m"); n <- Symbol("n"); w <- Symbol("w"); k <- Symbol("k"); shanks(n^2, n, 2) })), "-7/2")
})

test_that("approximants: pade_approximants_gcdex", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); lapply(pade_approximants_gcdex(Poly(x, x), 2), function(pq) .pyT(lapply(pq, as_expr))) })), "[(x, 1)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); lapply(pade_approximants_gcdex(Poly(1 - x, x), 1), function(pq) .pyT(lapply(pq, as_expr))) })), "[(1 - x, 1), (1, x + 1)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); lapply(pade_approximants_gcdex(Poly(1 + x + x^2 + x^3 + x^4, x), 2), function(pq) .pyT(lapply(pq, as_expr))) })), "[(x^2 + x + 1, 1), (1, 1 - x)]")
})

test_that("approximants: pade_approximant_gcdex", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); .pyT(lapply(pade_approximant_gcdex(Poly(1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120, x), -1, 2), as_expr)) })), "ERR ValueError: m=-1 and n=2 must be non-negative integers.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); .pyT(lapply(pade_approximant_gcdex(Poly(x^2 + x^3, x), 1, 1), as_expr)) })), "ERR ValueError: polynomial f=Poly(x^3 + x^2, x, domain='QQ') has zero of order 2, which is greater than 1, the requested order of the numerator.")
})

test_that("approximants: pade_approximant", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); pade_approximant(1/x, x, 0, 1, 0) })), "ERR ValueError: Expression f=1/x has pole of order 1, which is greater than 0, the requested order of the denominator.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); pade_approximant(sqrt(x), x, 0, 1, 1) })), "ERR ValueError: f=sqrt(x) does not have a laurent expansion in x around x0.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); pade_approximant(exp(x), x, 0, -1, 1) })), "ERR ValueError: m=-1 and n=1 must be non-negative integers")
})

