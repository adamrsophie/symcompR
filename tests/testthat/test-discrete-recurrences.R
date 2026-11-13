# Expected values are SymPy's outputs for the same calls (sympy/discrete/tests/test_recurrences.py);
# the rows come from the battery that compares this port with SymPy.

test_that("linrec: num", {
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 1), init = list(1, 1), n = 20) })), "10946")
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 2, 3, 4, 5), init = list(1, 1, 0, 2), n = 10) })), "1040")
  expect_identical(.pyrun(local({ linrec(coeffs = list(0, 0, 11, 13), init = list(23, 27), n = 25) })), "59628567384")
  expect_identical(.pyrun(local({ linrec(coeffs = list(0, 0, 1, 1, 2), init = list(1, 5, 3), n = 15) })), "165")
  expect_identical(.pyrun(local({ linrec(coeffs = list(11, 13, 15, 17), init = list(1, 2, 3, 4), n = 70) })), "56889923441670659718376223533331214868804815612050381493741233489928913241")
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 1), init = list(0, 1), n = 100) })), "354224848179261915075")
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 1), init = list(1, 1), n = 600) })), "178684461669052552311410692812805706249615844217278044703496837914086683543763273909969771627106004287604844670397177991379601")
  expect_identical(.pyrun(local({ linrec(coeffs = list(), init = list(1, 1), n = 20) })), "0")
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 1), init = list(1, 1), n = 1) })), "1")
  expect_identical(.pyrun(local({ linrec(coeffs = list(1, 1), init = list(3), n = 1) })), "0")
  expect_identical(.pyrun(local({ linrec(coeffs = list(11, 13, 15, 17), init = list(1, 2, 3, 4, 5), n = 1) })), "ERR TypeError: Count of initial values should not exceed the order of the recurrence")
  expect_identical(.pyrun(local({ linrec(coeffs = list(11, 13, 15, 17), init = list(1, 2, 3, 4), n = -1) })), "ERR ValueError: Point of evaluation of recurrence must be a non-negative integer")
  expect_identical(.pyrun(local({ linrec(coeffs = list(11, 13, 15, 17), init = list(1, 2, 3, 4), n = 2.5) })), "ERR ValueError: 2.5 is not an integer")
})

test_that("linrec: big", {
  expect_identical(.pyrun(local({ linrec(c(2, 3), c(1, 1), 10^4) %% 1000003 })), "109376")
  expect_identical(.pyrun(local({ linrec(c(rep(0, 55), 1, 1, 2, 3), c(rep(0, 50), 1, 2, 3), 4000) })), "702633573874937994980598979769135096432444135301118916539")
})

test_that("linrec: rat", {
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(a[1:5], b[1:4], 80) })), "1726244235456268979436592226626304376013002142588105090705187189/1960143456748895967474334873705475211264")
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(a[1:4], b[1:4], 50) })), "368949940033050147080268092104304441/504857282956046106624")
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(a[4:6], b[1:3], 35) })), "97409272177295731943657945116791049305244422833125109/814315512679031689453125")
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(c(rep(list(0), 60), list(Rational(2, 3), Rational(4, 5))), b, 3000) })), "26777668739896791448594650497024/48084516708184142230517578125")
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(a[1:4], b[1:5], 10000) })), "ERR TypeError: Count of initial values should not exceed the order of the recurrence")
  expect_identical(.pyrun(local({ a <- list(Rational(1, 2), Rational(3, 4), Rational(5, 6), 7, Rational(8, 9), Rational(3, 5)); b <- list(1, 2, 8, Rational(5, 7), Rational(3, 7), Rational(2, 9), 6); linrec(a[1:4], b[1:4], -10000) })), "ERR ValueError: Point of evaluation of recurrence must be a non-negative integer")
})

test_that("linrec: sym", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(x, c(1, 2), 10000) })), "ERR TypeError: Expected a sequence of coefficients for the recurrence")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(c(1, 2), y, 10000) })), "ERR TypeError: Expected a sequence of values for the initialization of the recurrence")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(c(1, 2, 1), list(x, y, z), 20) })), "269542*x + 664575*y + 578949*z")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(c(0, 3, 1, 2), list(x, y), 30) })), "58516436*x + 56372788*y")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(list(x, y, z), c(1, 2, 3), 2) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(list(x, y), c(1, 1), 5) })), "x^2*y + x*(x^3 + 2*x*y) + y^2 + y*(x^3 + 2*x*y)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(list(x, 1), c(0, 1), 6) })), "x*(x^2 + 1)^2 + 2*x*(x^2 + 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linrec(list(2, x), list(y, 1), 4) })), "4*x + y*(x^2 + 4*x) + 8")
})

