# The caches take keys of any length, and Wang's algorithm reuses its
# univariate diophantine solutions within one factorization.

test_that("expressions with keys over 10000 bytes are cached", {
  x <- Symbol("x"); y <- Symbol("y")
  e <- powsimp(expand((x + y + 1)^25))
  expect_length(e$args, 351L)
  expect_equal(substr(sstr(e), 1, 51), "x^25 + 25*x^24*y + 25*x^24 + 300*x^23*y^2 + 600*x^2")
})

test_that("factor of x^20 - y^20 matches SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(sym_factor(x^20 - y^20)),
               paste0("(x - y)*(x + y)*(x^2 + y^2)*(x^4 - x^3*y + x^2*y^2 - x*y^3 + y^4)*",
                      "(x^4 + x^3*y + x^2*y^2 + x*y^3 + y^4)*(x^8 - x^6*y^2 + x^4*y^4 - x^2*y^6 + y^8)"))
})
