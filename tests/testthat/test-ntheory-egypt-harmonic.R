# Expected values are SymPy's outputs for the same calls (sympy/ntheory/tests/test_egyptian_fraction.py (egypt_harmonic));
# the rows come from the battery that compares this port with SymPy.

test_that("egypt_harmonic: egypt", {
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(1, 3))) })), "([], 1/3)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(5, 2))) })), "([1, 2, 3, 4, 5, 6], 1/20)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Integer(2))) })), "([1, 2, 3], 1/6)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Integer(1))) })), "([1], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(11, 6))) })), "([1, 2, 3], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(25, 12))) })), "([1, 2, 3, 4], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Integer(3))) })), "([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 179/2520)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(7, 3))) })), "([1, 2, 3, 4, 5], 1/20)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(3, 7))) })), "([], 3/7)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(1, 1000))) })), "([], 1/1000)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(99, 100))) })), "([], 99/100)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(3, 2))) })), "([1, 2], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Integer(0))) })), "([], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(-1, 2))) })), "([], -1/2)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(137, 60))) })), "([1, 2, 3, 4, 5], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(49, 20))) })), "([1, 2, 3, 4, 5, 6], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(363, 140))) })), "([1, 2, 3, 4, 5, 6, 7], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(761, 280))) })), "([1, 2, 3, 4, 5, 6, 7, 8], 0)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Integer(4))) })), "([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30], 11675421053/2329089562800)")
  expect_identical(.pyrun(local({ .pyT(egypt_harmonic(Rational(101, 37))) })), "([1, 2, 3, 4, 5, 6, 7, 8], 123/10360)")
})

