# Reference layouts from SymPy's pretty().
x <- Symbol("x"); y <- Symbol("y")
pa <- function(e) pretty(e, use_unicode = FALSE)
pu <- function(e) pretty(e, use_unicode = TRUE)

test_that("ASCII pretty printing", {
  expect_equal(pa(x^2 + 1), " 2    \nx  + 1")
  expect_equal(pa((x + 1)/(y - 2)), "x + 1\n-----\ny - 2")
  expect_equal(pa(sqrt(x)), "  ___\n\\/ x ")
  expect_equal(pa(sqrt(x^2 + 1)), "   ________\n  /  2     \n\\/  x  + 1 ")
  expect_equal(pa(sin(x)^2), "   2   \nsin (x)")
  expect_equal(pa(x^(y + 1)), " y + 1\nx     ")
  expect_equal(pa(1/x^2), "1 \n--\n 2\nx ")
  expect_equal(pa(x^Rational(1, 3)), "3 ___\n\\/ x ")
  expect_equal(pa((x + y)^2), "       2\n(x + y) ")
  expect_equal(pa(sin(x/2)), "   /x\\\nsin|-|\n   \\2/")
  expect_equal(pa(-x/2 + y), "  x    \n- - + y\n  2    ")
  expect_equal(pa(2*x*y), "2*x*y")
})

test_that("Unicode pretty printing", {
  expect_equal(pu(x/y), "x\n\u2500\ny")
  expect_equal(pu(sqrt(x)), "\u221ax")
  expect_equal(pu(exp(x)), " x\n\u212f ")
  expect_equal(pu(2*x*y), "2\u22c5x\u22c5y")
  expect_equal(pu(sin(x/2)), "   \u239bx\u239e\nsin\u239c\u2500\u239f\n   \u239d2\u23a0")
  expect_output(pprint(x^2, use_unicode = FALSE), " 2\nx ")
})
