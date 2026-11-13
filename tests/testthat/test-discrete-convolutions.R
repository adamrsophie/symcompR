# Expected values are SymPy's outputs for the same calls (sympy/discrete/tests/test_convolutions.py);
# the rows come from the battery that compares this port with SymPy.

test_that("convolutions: conv", {
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(b, d, dps = 5, prime = q) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(b, d, dps = 2, dyadic = TRUE) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(b, cc, prime = p, dyadic = TRUE) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(a, b, subset = TRUE, dyadic = TRUE) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(cc, d, subset = TRUE, dps = 6) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(a, cc, subset = TRUE, prime = q) })), "ERR TypeError: Ambiguity in determining the type of convolution")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(Rational(1,2)), list(Rational(1,2))) })), "[1/4]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(b, e) })), "[-9, 10, 239/15, 34/3, 32/3, 43/5, 113/15, 14/5]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 0, 1), list(1, 1)) })), "[1, 1, 1, 1]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(Integer(0), Integer(1), Integer(0), Rational(-1, 6), Integer(0), Rational(1, 120)), list(Integer(1), Integer(0), Rational(-1, 2), Integer(0), Rational(1, 24), Integer(0))) })), "[0, 1, 0, -2/3, 0, 2/15, 0, -1/90, 0, 1/2880, 0]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 0), list(1, 1)) })), "[1, 1, 0]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(0, 0), list(1, 2)) })), "[0, 0, 0]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 0), list(1, 1), cycle = 3) })), "[1, 1, 0]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 2, 3), list(4, 5, 6), cycle = 5) })), "[4, 13, 28, 27, 18]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 2, 3), list(4, 5, 6), cycle = 3) })), "[31, 31, 28]")
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5)); b <- c(9, 5, 5, 4, 3, 2); cc <- c(3, 5, 3, 7, 8); d <- c(1422, 6572, 3213, 5552); e <- list(-1, Rational(5, 3), Rational(7, 5)); p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; x <- Symbol("x"); convolution(list(1, 2, 3), list(4, 5, 6), cycle = -1) })), "ERR ValueError: The length for cyclic convolution must be non-negative")
})

test_that("convolutions: cyc", {
  expect_identical(.pyrun(local({ a <- list(Rational(1, 3), Rational(7, 3), Rational(5, 9), Rational(2, 7), Rational(5, 8)); b <- list(Rational(3, 5), Rational(4, 7), Rational(7, 8), Rational(8, 9)); convolution(a, b, cycle = 0) })), "[1/5, 167/105, 47/24, 3053/1080, 16397/5292, 2497/2268, 3229/4032, 5/9]")
})

test_that("convolutions: cyc2", {
  expect_identical(.pyrun(local({ a <- list(2313, 5323532, Integer(3232), 42142, 42242421); b <- list(Integer(33456), 56757, 45754, 432423); u <- Symbol("u"); v <- Symbol("v"); w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); p <- Symbol("p"); q <- Symbol("q"); r <- Symbol("r"); s <- Symbol("s"); t <- Symbol("t"); cc <- list(u, v, w, x, y); d <- list(p, q, r, s, t); convolution(cc, d) })), "ERR ValueError: Expected non-symbolic coefficients")
})

test_that("convolutions: fft", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); convolution_fft(list(), list(), dps = NULL) })), "[]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); convolution_fft(list(), list(1), dps = 3) })), "[]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); convolution_fft(x, y) })), "ERR TypeError: 'Symbol' object is not subscriptable")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); convolution_fft(list(x, y), list(y, x)) })), "ERR ValueError: Expected non-symbolic coefficients")
})

test_that("convolutions: ntt", {
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(), list(1), prime = p) })), "[]")
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(2), list(3), r) })), "[6]")
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(2, 3), list(4), r) })), "[8, 12]")
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(2, 3), list(4, 5), r) })), "ERR ValueError: Expected prime modulus of the form (m*2^k + 1)")
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(x, y), list(y, x), q) })), "ERR ValueError: x is not an integer")
  expect_identical(.pyrun(local({ p <- 7*17*2^23 + 1; q <- 19*2^10 + 1; r <- 2*500000003 + 1; x <- Symbol("x"); y <- Symbol("y"); convolution_ntt(list(1, 2), list(3, 4), 15) })), "ERR ValueError: Expected prime modulus for Number Theoretic Transform")
})

test_that("convolutions: fwht", {
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), Rational(7, 5), 4 + 5*sym_I); b <- c(94, 51, 53, 45, 31, 27, 13); cc <- list(3 + 4*sym_I, 5 + 7*sym_I, 3, Rational(7, 6), 8); u <- Symbol("u"); v <- Symbol("v"); w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); convolution_fwht(list(), list(1)) })), "[]")
})

test_that("convolutions: subset", {
  expect_identical(.pyrun(local({ a <- list(1, Rational(5, 3), sqrt(Integer(3)), 4 + 5*sym_I); b <- c(64, 71, 55, 47, 33, 29, 15); cc <- list(3 + sym_I*2/3, 5 + 7*sym_I, 7, Rational(7, 5), 9); u <- Symbol("u"); v <- Symbol("v"); w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); convolution_subset(list(), list()) })), "[]")
})

test_that("convolutions: int", {
  expect_identical(.pyrun(local({ convolution_int(c(1), c(1)) })), "[1]")
  expect_identical(.pyrun(local({ convolution_int(c(1, 1), c(0)) })), "[0]")
  expect_identical(.pyrun(local({ convolution_int(c(1, 2, 3), c(4, 5, 6)) })), "[4, 13, 28, 27, 18]")
  expect_identical(.pyrun(local({ convolution_int(c(1), c(5, 6, 7)) })), "[5, 6, 7]")
  expect_identical(.pyrun(local({ convolution_int(c(1, 3), c(5, 6, 7)) })), "[5, 21, 25, 21]")
  expect_identical(.pyrun(local({ convolution_int(c(10, -5, 1, 3), c(-5, 6, 7)) })), "[-50, 85, 35, -44, 25, 21]")
  expect_identical(.pyrun(local({ convolution_int(c(0, 1, 0, -1), c(1, 0, -1, 0)) })), "[0, 1, 0, -2, 0, 1]")
  expect_identical(.pyrun(local({ convolution_int(c(0, 0, 1), c(1, 0)) })), "[0, 0, 1]")
  expect_identical(.pyrun(local({ convolution_int(c(-1, 2), c(3, 0, 0)) })), "[-3, 6]")
  expect_identical(.pyrun(local({ convolution_int(c(123456789012345, 3), c(987654321, -2)) })), "[121932631124827861592745, -246910615061727, -6]")
})

