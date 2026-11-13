# Expected values are SymPy's outputs for the same calls (sympy/concrete/tests/test_guess.py);
# the rows come from the battery that compares this port with SymPy.

test_that("guess: fsrv", {
  expect_identical(.pyrun(local({ find_simple_recurrence_vector(c(0, 0, 0, 0, 0, 0, 0, 0)) })), "[1]")
  expect_identical(.pyrun(local({ find_simple_recurrence_vector(c(1, 2, 3)) })), "[0]")
  expect_identical(.pyrun(local({ find_simple_recurrence_vector(c(1, 1)) })), "[0]")
})

test_that("guess: fsr", {
  expect_identical(.pyrun(local({ find_simple_recurrence(c(0, 0, 0, 0, 0, 0, 0, 0)) })), "0")
  expect_identical(.pyrun(local({ find_simple_recurrence(c(1, 1)) })), "0")
})

test_that("guess: guess", {
  expect_identical(.pyrun(local({ guess(c(0, 0, 0, 0, 0, 0, 0, 0)) })), "[]")
  expect_identical(.pyrun(local({ guess(c(0, 0, 1, 2, 3, 4, 5, 6, 7, 8)) })), "[]")
  expect_identical(.pyrun(local({ guess(list(Rational(3, 2), Rational(11, 2), 0, Rational(-121, 2), Rational(-363, 2), 121, Rational(4719, 2), Rational(11495, 2), -8712, Rational(-178717, 2))) })), "[]")
  expect_identical(.pyrun(local({ guess(c(1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496)) })), "[]")
  expect_identical(.pyrun(local({ guess(c(0, 2, 8, 24, 64, 160, 384, 896, 2048, 4608)) })), "[]")
  expect_identical(.pyrun(local({ guess(c(1, 1)) })), "[]")
  expect_identical(.pyrun(local({ guess(c(1, 0, 2)) })), "[]")
  expect_identical(.pyrun(local({ guess(5) })), "[]")
  expect_identical(.pyrun(local({ guess(c(1, 0, 0)) })), "[]")
})

test_that("guess: ggfr", {
  expect_identical(.pyrun(local({ guess_generating_function_rational(c(1, 1)) })), "None")
})

test_that("guess: rationalize", {
  expect_identical(.pyrun(local({ rationalize(cos(pi/3)) })), "1/2")
  expect_identical(.pyrun(local({ rationalize(-0.333333333333333) })), "-1/3")
  expect_identical(.pyrun(local({ rationalize(pi, maxcoeff = 250) })), "355/113")
  expect_identical(.pyrun(local({ rationalize(0.1) })), "1/10")
  expect_identical(.pyrun(local({ rationalize(0.5) })), "1/2")
  expect_identical(.pyrun(local({ rationalize(2) })), "2")
  expect_identical(.pyrun(local({ rationalize(-2.5) })), "-5/2")
  expect_identical(.pyrun(local({ rationalize(0) })), "0")
  expect_identical(.pyrun(local({ rationalize(1/7) })), "1/7")
  expect_identical(.pyrun(local({ rationalize(22/7, maxcoeff = 5) })), "3")
  expect_identical(.pyrun(local({ rationalize(355/113) })), "355/113")
  expect_identical(.pyrun(local({ rationalize(0.142857142857) })), "1/7")
  expect_identical(.pyrun(local({ rationalize(1234.5678) })), "6172839/5000")
  expect_identical(.pyrun(local({ rationalize(2/3) })), "2/3")
  expect_identical(.pyrun(local({ rationalize(-5/12) })), "-5/12")
  expect_identical(.pyrun(local({ rationalize(sqrt(2), maxcoeff = 3) })), "209064253/147830751")
  expect_identical(.pyrun(local({ rationalize(pi, maxcoeff = 10) })), "22/7")
  expect_identical(.pyrun(local({ rationalize(pi, maxcoeff = 20) })), "355/113")
  expect_identical(.pyrun(local({ rationalize(1e-5) })), "0")
  expect_identical(.pyrun(local({ rationalize(123456) })), "123456")
  expect_identical(.pyrun(local({ rationalize(0.75, maxcoeff = 1) })), "0")
  expect_identical(.pyrun(local({ rationalize(1/3, maxcoeff = 2) })), "0")
  expect_identical(.pyrun(local({ rationalize(Rational(3, 7)) })), "3/7")
  expect_identical(.pyrun(local({ rationalize(Rational(-22, 7)) })), "-22/7")
  expect_identical(.pyrun(local({ rationalize(Float(0.5)) })), "ERR TypeError: Invalid comparison of non-real zoo")
})

