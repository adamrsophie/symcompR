# Reference values from SymPy.
x <- Symbol("x", integer = TRUE); y <- Symbol("y", integer = TRUE)
z <- Symbol("z", integer = TRUE); w <- Symbol("w", integer = TRUE)
s <- function(r) paste(sort(vapply(if (r$type == "FiniteSet") r$args else list(), sstr, ""), method = "radix"),
                       collapse = " ; ")

test_that("linear diophantine equations", {
  expect_equal(s(diophantine(2*x + 3*y - 5)), "(3*t_0 - 5, 5 - 2*t_0)")
  expect_equal(s(diophantine(3*x - 6*y + 7)), "")
  expect_equal(s(diophantine(x + 3*y - 4*z + 5)), "(t_0, 5*t_0 + 4*t_1 + 5, 4*t_0 + 3*t_1 + 5)")
  expect_equal(s(diophantine(4*x + 6*y + 3*z - 2)), "(3*t_0 + 2, t_0 + t_1, -6*t_0 - 2*t_1 - 2)")
  expect_equal(s(diophantine(2*x + 4*y + 6*z + 8*w - 10)),
               "(t_0, t_0 + t_1, 8*t_0 + 4*t_1 + 3*t_2 - 5, -7*t_0 - 3*t_1 - 2*t_2 + 5)")
  expect_equal(s(diophantine(x^2 - 5*x + 6)), "(2,) ; (3,)")
})

test_that("binary quadratic and Pell equations", {
  expect_equal(s(diophantine(x*y + 3*x - 2*y - 7)), "(1, -4) ; (3, -2)")
  expect_equal(s(diophantine(x^2 + 2*x*y + y^2 + x + y - 6)), "(t_0, -t_0 - 3) ; (t_0, 2 - t_0)")
  expect_equal(s(diophantine(2*x^2 + 3*x*y + y^2 - 5)), "(-4, 3) ; (-4, 9) ; (4, -3) ; (4, -9)")
  expect_equal(s(diophantine(x^2 + 2*y^2 - 9)), "(-1, -2) ; (-1, 2) ; (-3, 0) ; (1, -2) ; (1, 2) ; (3, 0)")
  expect_equal(s(diophantine(3*x^2 + 2*x*y - y^2 + 7*x - 4*y + 1)), "(4, -7) ; (4, 11)")
  expect_equal(s(diophantine(x^2 - 2*y^2 + 2)),
               paste("(-sqrt(2)*(3 - 2*sqrt(2))^t/2 + sqrt(2)*(2*sqrt(2) + 3)^t/2, (3 - 2*sqrt(2))^t/2 + (2*sqrt(2) + 3)^t/2) ;",
                     "(sqrt(2)*(3 - 2*sqrt(2))^t/2 - sqrt(2)*(2*sqrt(2) + 3)^t/2, -(3 - 2*sqrt(2))^t/2 - (2*sqrt(2) + 3)^t/2)"))
  expect_equal(s(diophantine(x^2 + x*y - 3*y^2 + 2*x - 5)), "")
})

test_that("ternary, Pythagorean and power sums", {
  expect_equal(s(diophantine(x^2 + y^2 - z^2)), "(2*p*q, p^2 - q^2, p^2 + q^2)")
  expect_equal(s(diophantine(x^2 + 2*y^2 - 3*z^2)), "(2*p^2 - 3*q^2, -2*p^2 + 6*p*q - 3*q^2, 2*p^2 - 4*p*q + 3*q^2)")
  expect_equal(s(diophantine(x*y - 2*z^2)), "(-2*q^2, -p^2 + 4*p*q - 4*q^2, -p*q + 2*q^2)")
  expect_equal(s(diophantine(x^2 - 4*x*y + 2*y^2 - z^2)), "(2*p^2 - q^2, 4*p^2 + 2*p*q, 2*p^2 + 4*p*q + q^2)")
  expect_equal(s(diophantine(4*x^2 + 9*y^2 - 36*z^2)), "(18*p*q, 2*p^2 - 18*q^2, p^2 + 9*q^2)")
  expect_equal(s(diophantine(x^2 + y^2 + z^2 - 14)), "(1, 2, 3)")
  expect_equal(s(diophantine(x^4 + y^4 + z^4 - 98)), "(1, 2, 3)")
})

test_that("products, denominators, permutations and assumptions", {
  expect_equal(s(diophantine((x - 1)*(y - 2))), "(1, n1) ; (n1, 2)")
  expect_equal(s(diophantine((x^2 + y^2 - 1)*(x - y))), "(-1, 0) ; (0, -1) ; (0, 1) ; (1, 0) ; (t_0, t_0)")
  expect_equal(s(diophantine(x/y - 3)), "(3*t_0, t_0)")
  expect_equal(s(diophantine(x^2 + y^2 - 5, permute = TRUE)),
               "(-1, -2) ; (-1, 2) ; (-2, -1) ; (-2, 1) ; (1, -2) ; (1, 2) ; (2, -1) ; (2, 1)")
  expect_equal(s(diophantine(2*x + 3*y - 5, syms = list(y, x))), "(5 - 2*t_0, 3*t_0 - 5)")
  p <- Symbol("p", positive = TRUE, integer = TRUE); q <- Symbol("q", positive = TRUE, integer = TRUE)
  expect_equal(s(diophantine(p*q - 6)), "(1, 6) ; (2, 3) ; (3, 2) ; (6, 1)")
  expect_error(diophantine(x^2 + y^2 - z^2 - w^2), "not yet recognized")
  cl <- classify_diop(x^2 + y^2 - x*y + x + 5)
  expect_equal(cl[[3L]], "binary_quadratic")
  expect_equal(sstr(diop_solve(2*x + 3*y - 5)), "(3*t_0 - 5, 5 - 2*t_0)")
})

test_that("p11 review fixes", {
  expect_equal(s(diophantine(9*x^2 + 6*x*y + y^2 + x + y)), "(2*t^2 - 3*t + 1, -6*t^2 + 7*t - 2) ; (2*t^2 - t, -6*t^2 + t)")
  expect_equal(s(diophantine(9*x^2 + 6*x*y + y^2 + x)), "(-t^2, 3*t^2 + t)")
  b <- gmp::as.bigz(2)^100 + gmp::as.bigz(2)^47 + 1
  expect_equal(as.character(.zroot(b^4, 4)[[1L]]), as.character(b))
  expect_equal(sstr(diop_solve(2*x + 3*y - 5, param = NULL)), "(-5, 5)")
  expect_identical(diop_solve(x^2 + y^2 + 1), S_EmptySet)
  expect_equal(.cf_length(0, -1, 2), 4L)
  expect_equal(.cf_length(1, 1, 4), 1L)
  pr <- vapply(.power_representation(10, 1, 3), function(t) paste(vapply(t, as.character, ""), collapse = ","), "")
  expect_equal(pr, c("1,1,8", "1,2,7", "1,3,6", "1,4,5", "2,2,6", "2,3,5", "2,4,4", "3,3,4"))
  expect_equal(s(diophantine(x/2 + y/2 - 1)), "(t_0, 2 - t_0)")
  expect_equal(sstr(diop_solve(x^1024 + y^1024 + z^1024 - as_sym(gmp::as.bigz(2)^1024 + 2))), "{(1, 1, 2)}")
})
