# Expected values are SymPy's outputs for the same calls (sympy/polys/tests/test_appellseqs.py and test_polytools.py (gcdex_steps));
# the rows come from the battery that compares this port with SymPy.

test_that("appellseqs: bernoulli_poly", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_poly(-1, x) })), "ERR ValueError: Cannot generate Bernoulli polynomial of index -1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_poly(4, x) })), "x^4 - 2*x^3 + x^2 - 1/30")
  expect_identical(.pyrun(local({ bernoulli_poly(3) })), "_0^3 - 3*_0^2/2 + _0/2")
  expect_identical(.pyrun(local({ bernoulli_poly(4, Symbol("t")) })), "t^4 - 2*t^3 + t^2 - 1/30")
})

test_that("appellseqs: bernoulli_c_poly", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_c_poly(0, x) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_c_poly(2, x) })), "x^2 - 1/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_c_poly(3, x) })), "x^3 - x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_c_poly(4, x) })), "x^4 - 2*x^2 + 7/15")
  expect_identical(.pyrun(local({ x <- Symbol("x"); bernoulli_c_poly(5, x) })), "x^5 - 10*x^3/3 + 7*x/3")
  expect_identical(.pyrun(local({ bernoulli_c_poly(3) })), "_0^3 - _0")
  expect_identical(.pyrun(local({ bernoulli_c_poly(4, Symbol("t")) })), "t^4 - 2*t^2 + 7/15")
  expect_identical(.pyrun(local({ bernoulli_c_poly(4, Rational(1, 2)) })), "7/240")
})

test_that("appellseqs: genocchi_poly", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); genocchi_poly(-1, x) })), "ERR ValueError: Cannot generate Genocchi polynomial of index -1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); genocchi_poly(0, x) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); genocchi_poly(2, x) })), "1 - 2*x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); genocchi_poly(3, x) })), "-3*x^2 + 3*x")
  expect_identical(.pyrun(local({ genocchi_poly(3) })), "-3*_0^2 + 3*_0")
  expect_identical(.pyrun(local({ genocchi_poly(4, Symbol("t")) })), "-4*t^3 + 6*t^2 - 1")
  expect_identical(.pyrun(local({ genocchi_poly(4, Rational(1, 2)) })), "0")
})

test_that("appellseqs: euler_poly", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(-1, x) })), "ERR ValueError: Cannot generate Euler polynomial of index -1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(0, x) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(1, x) })), "x - 1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(2, x) })), "x^2 - x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(3, x) })), "x^3 - 3*x^2/2 + 1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); euler_poly(4, x) })), "x^4 - 2*x^3 + x")
  expect_identical(.pyrun(local({ euler_poly(3) })), "_0^3 - 3*_0^2/2 + 1/4")
  expect_identical(.pyrun(local({ euler_poly(4, Symbol("t")) })), "t^4 - 2*t^3 + t")
  expect_identical(.pyrun(local({ euler_poly(4, Rational(1, 2)) })), "5/16")
})

test_that("appellseqs: andre_poly", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); andre_poly(-1, x) })), "ERR ValueError: Cannot generate Andre polynomial of index -1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); andre_poly(0, x) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); andre_poly(2, x) })), "x^2 - 1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); andre_poly(3, x) })), "x^3 - 3*x")
  expect_identical(.pyrun(local({ andre_poly(3) })), "_0^3 - 3*_0")
  expect_identical(.pyrun(local({ andre_poly(4, Symbol("t")) })), "t^4 - 6*t^2 + 5")
  expect_identical(.pyrun(local({ andre_poly(4, Rational(1, 2)) })), "57/16")
})

test_that("gcdex_steps: steps", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lapply(gcdex_steps(x, x), .pyT) })), "[(0, 1, x)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lapply(gcdex_steps(12, 18), .pyT) })), "[(0, 1, 18), (1, 0, 12), (-1, 1, 6)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lapply(gcdex_steps(240, 46), .pyT) })), "[(0, 1, 46), (1, -5, 10), (-4, 21, 6), (5, -26, 4), (-9, 47, 2)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lapply(gcdex_steps(Rational(1, 2), Rational(3, 4)), .pyT) })), "[(0, 1, 3/4)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lapply(gcdex_steps(17, 5), .pyT) })), "[(0, 1, 5), (1, -3, 2), (-2, 7, 1)]")
})

