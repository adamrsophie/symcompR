# N() through the port of sympy/core/evalf.py: the evalf_table dispatch,
# complex accuracy tracking and the options. Expected strings are SymPy's.

P <- function(s) parse_expr(s)
x <- Symbol("x")

test_that("constants, rationals and cancellation carry SymPy's digits", {
  expect_identical(sstr(N(P("pi"), 50)), "3.1415926535897932384626433832795028841971693993751")
  expect_identical(sstr(N(P("355/113-pi"), 6)), "2.66764e-7")
  expect_identical(sstr(N(P("2**(1/10**50)-1"))), "6.93147180559945e-51")
  expect_identical(sstr(N(P("Add(pi, 1/10**1000, -pi, evaluate=False)"), maxn = 1200)), "1.00000000000000e-1000")
  expect_identical(sstr(N(P("Add(pi, 1/10**1000, -pi, evaluate=False)"))), "-0.e-124")
  expect_identical(sstr(N(P("exp(pi*sqrt(163)) - 640320**3 - 744"), 10)), "-7.499274028e-13")
})

test_that("complex parts keep their own accuracy", {
  expect_identical(sstr(N(P("(pi + 1/10**8 + pi*I)**2"))), "6.2831853e-8 + 19.7392088650106*I")
  expect_identical(sstr(N(P("(pi + pi*I)**2"))), "-0.e-21 + 19.7392088021787*I")
  expect_identical(sstr(N(P("(pi + pi*I)**2"), chop = TRUE)), "19.7392088021787*I")
  expect_identical(sstr(N(P("(-pi)**pi"))), "-32.9138577418939 - 15.6897116534332*I")
  expect_identical(sstr(N(P("log(-1 + 0.00001)"), 2)), "-1.0e-5 + 3.1*I")
})

test_that("functions, subs and symbols", {
  expect_identical(sstr(N(P("sin(exp(pi*sqrt(163))*pi)"))), "-2.35596641936785e-12")
  expect_identical(sstr(N(P("x + pi + E"))), "x + 5.85987448204884")
  expect_identical(sstr(N(P("Max(y, 4)"))), "Max(4.0, y)")
  expect_identical(sstr(N(P("ceiling(x)"), subs = list(x = P("2 + 3*I")))), "2.0 + 3.0*I")
  expect_identical(sstr(N(P("log(100, 10, evaluate=False)"))), "2.00000000000000")
  expect_identical(sstr(N(P("sin(160*pi, evaluate=False)"))), "-0.e-177")
  expect_identical(sstr(N(P("1/x"), subs = list(x = 0))), "zoo")
})

test_that("sums, products and CRootOf at full precision", {
  expect_identical(sstr(N(P("Sum(1/factorial(n), (n, 0, oo))"))), "2.71828182845905")
  expect_identical(sstr(N(P("CRootOf(x**5 - x + 1, 0)"), 40)), "-1.167303978261418684256045899854842180721")
  expect_identical(sstr(N(P("CRootOf(x**5 - x + 1, 1)"), 40)),
                   "-0.1812324444698753839018002377811206399687 - 1.083954101317710668430344492980766574274*I")
  expect_identical(sstr(N(CRootOf(x^2 + 1, x, 0), 4)), "-1.0*I")
})

test_that("strict raises when precision runs out", {
  expect_error(N(P("sin(160*pi, evaluate=False)"), strict = TRUE), class = "PrecisionExhausted")
})
