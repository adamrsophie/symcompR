# Expected values are SymPy's outputs for the same calls (sympy/calculus/tests/test_util.py and test_singularities.py);
# the rows come from the battery that compares this port with SymPy.

test_that("function_range: fr", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); function_range(cos(x), x, EmptySet) })), "EmptySet")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); function_range(x, x, FiniteSet(1, 2, 3)) })), "{1, 2, 3}")
})

test_that("lcim: lcim", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(sym_pi/2, sym_pi/4, sym_pi)) })), "pi")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(2*sym_pi, sym_pi/2)) })), "2*pi")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(Integer(1), 2*sym_pi)) })), "None")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(Integer(4), Integer(6))) })), "12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(Rational(2, 3), Rational(3, 4))) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(sqrt(Integer(2)), 3*sqrt(Integer(2))/2)) })), "3*sqrt(2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(sym_pi, sym_E)) })), "None")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(sym_pi/3, 2*sym_pi/5)) })), "2*pi")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(x, Integer(2))) })), "None")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(Integer(5))) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); lcim(list(sym_pi)) })), "pi")
})

test_that("is_convex: convex", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(1/x, x, domain = Interval(0, oo, TRUE, TRUE)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(x^2, x, domain = Interval(0, oo)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(-1/x^3, x, domain = Interval(-oo, 0, FALSE, TRUE)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(log(x) ,x) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(exp(x), x) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(-x^2, x) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); is_convex(x*log(x), x, domain = Interval(0, oo, TRUE, TRUE)) })), "True")
})

test_that("stationary_points: stationary", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); stationary_points((x+3)*(x-2), x) })), "{-1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); stationary_points(exp(x), x) })), "EmptySet")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); stationary_points(y, x, Reals) })), "Reals")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); stationary_points(y, x, EmptySet) })), "EmptySet")
})

test_that("monotonicity_helper: mono", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(4*x + 5, function(d) Ge(d, 0), Reals) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(4*x + 5, function(d) Gt(d, 0), Reals) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(4*x + 5, function(d) Le(d, 0), Reals) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(4*x + 5, function(d) Lt(d, 0), Reals) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(1/x, function(d) Ge(d, 0), Interval(-2, 2)) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(1/x, function(d) Le(d, 0), Interval(-2, 2)) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(1/x, function(d) Lt(d, 0), Interval(-2, 2)) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(exp(x), function(d) Ge(d, 0), Reals) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(exp(x), function(d) Gt(d, 0), Reals) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(exp(x), function(d) Le(d, 0), Reals) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(exp(x), function(d) Lt(d, 0), Reals) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(y*x, function(d) Ge(d, 0), Reals) })), "ERR NotImplementedError: The function has not yet been implemented for all multivariate expressions.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(y*x, function(d) Gt(d, 0), Reals) })), "ERR NotImplementedError: The function has not yet been implemented for all multivariate expressions.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(y*x, function(d) Le(d, 0), Reals) })), "ERR NotImplementedError: The function has not yet been implemented for all multivariate expressions.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(y*x, function(d) Lt(d, 0), Reals) })), "ERR NotImplementedError: The function has not yet been implemented for all multivariate expressions.")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); monotonicity_helper(Integer(5), function(d) Ge(d, 0)) })), "True")
})

