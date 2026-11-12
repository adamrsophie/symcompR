# Expected values are SymPy 1.15 outputs (sympy.simplify.epathtools).

skip_unless_full_tests()

test_that("epath selects and modifies nested lists", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expr <- list(list(list(x, 1), 2), list(list(3, y), z))
  f <- function(a) if (is_sym(a)) sstr(a) else if (is.list(a)) paste0("(", paste(vapply(a, f, ""), collapse = ", "), ")") else format(a)
  expect_equal(vapply(epath("/*/[0]/Symbol", expr), sstr, ""), c("x", "y"))
  expect_equal(f(epath("/*/[0]/Symbol", expr, function(e) e^2)), "(((x^2, 1), 2), ((3, y^2), z))")
  expect_equal(f(epath("/*/[-1]", expr)), "(2, z)")
  expect_equal(f(epath("/[::-1]", expr)), "(((3, y), z), ((x, 1), 2))")
  expect_equal(f(epath("/[1:]", expr)), "(((3, y), z))")
  expect_equal(f(epath("/*/*", expr)), "((x, 1), 2, (3, y), z)")
})

test_that("epath walks expressions in printing order", {
  x <- Symbol("x"); y <- Symbol("y"); t <- Symbol("t")
  e <- t + sin(x + 1) + cos(x + y + sym_E)
  expect_equal(vapply(epath("/*/*/Symbol", e), sstr, ""), c("x", "x", "y"))
  expect_equal(sstr(epath("/*/*/Symbol", e, function(a) 2 * a)), "t + sin(2*x + 1) + cos(2*x + 2*y + E)")
  e <- 3*x*sin(x + 1)^2 + exp(y)/2 + cos(x + y + sym_E) + Rational(1, 2)*sym_pi + 7
  sel <- function(p) paste(vapply(epath(p, e), sstr, ""), collapse = " ; ")
  expect_equal(sel("/TrigonometricFunction"), "cos(x + y + E)")
  expect_equal(sel("/*/Number"), "3 ; 1/2 ; 1/2")
  expect_equal(sel("/Mul/*"), "3 ; x ; sin(x + 1)^2 ; 1/2 ; exp(y) ; 1/2 ; pi")
  expect_equal(sel("/AtomicExpr"), "7")
  expect_equal(sstr(epath("/*/*", e, function(a) a + 1)),
               "4*(x + 1)*(sin(x + 1)^2 + 1) + 3*exp(y)/2 + cos(x + y + 1 + E) + 3*pi/2 + 10")
  expect_equal(sstr(epath("/*/Add/Symbol", e, function(a) a + 1)),
               "3*x*sin(x + 1)^2 + exp(y)/2 + cos(x + y + 2 + E) + pi/2 + 7")
})

test_that("EPath raises SymPy's errors", {
  msg <- function(p) tryCatch(epath(p, list()), error = function(e) paste(class(e)[1], conditionMessage(e)))
  expect_equal(msg(""), "ValueError empty EPath")
  expect_equal(msg("x"), "NotImplementedError non-root EPath")
  expect_equal(msg("//"), "ValueError empty selector")
  expect_equal(msg("/[0"), "ValueError expected ']', got EOL")
  expect_equal(msg("/[a]"), "ValueError invalid literal for int() with base 10: 'a'")
  expect_equal(msg("/Symbol[0]x"), "ValueError trailing characters in selector")
  expect_equal(msg("/[1:2:3:4]"), "TypeError slice expected at most 3 arguments, got 4")
  expect_equal(msg("/a||b"), "ValueError empty element")
  expect_output(print(EPath("/*/[0]")), "EPath('/*/[0]')", fixed = TRUE)
})
