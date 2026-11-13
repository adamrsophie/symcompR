# Expected values are SymPy's outputs for the same calls (sympy/integrals/tests/test_intpoly.py);
# the rows come from the battery that compares this port with SymPy.

test_that("intpoly: decompose separate", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x + y, TRUE)) })), "{x, y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^2 + y, TRUE)) })), "{x^2, y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(3*x^2*y, TRUE)) })), "{3*x^2*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^5, TRUE)) })), "{x^5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x + y + z + 1, TRUE)) })), "{1, x, y, z}")
})

test_that("intpoly: decompose", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(Integer(5))) })), "{0: 5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(3*x^2*y)) })), "{3: 3*x^2*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^5)) })), "{5: x^5}")
})

test_that("intpoly: best_origin", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 2), 3, Segment(Point(0, Integer(3)/2), Point(3, Integer(3)/2)), x ^ 7)) })), "(0, 1.5)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 2, Segment(Point(2, 0), Point(1, 1)), x ^ 9 * y ^ 2)) })), "(2.0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), 5)) })), "(-1, 2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 1), 1, Segment(Point(-2, 1), Point(2, 1)), x)) })), "(0, 1.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 1), 1, Segment(Point(Integer(3)/2, 0), Point(Integer(3)/2, 3)), x)) })), "(3/2, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 0), 2, Segment(Point(0, Integer(3)/2), Point(3, Integer(3)/2)), y)) })), "(0, 3/2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 0), 2, Segment(Point(-2, 1), Point(2, 1)), y^2)) })), "(-2, 1)")
})

test_that("intpoly: point_sort list", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(0, 1))) })), "[Point2D(0, 1)]")
})

test_that("intpoly: hrep 3d", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(1, 0, 0), 1), list(c(0, 1, 0), 1)), 1) })), "ERR NotImplementedError: Integration for H-representation 3Dcase not implemented yet.")
})

test_that("intpoly: lineseg_integrate", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), 0, 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), 0, 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), z, 1) })), "0")
})

test_that("intpoly: gradient_terms", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(0), .pyL) })), "[[1, 0, 0, 0]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(1), .pyL) })), "[[1, 0, 0, 0], [y, 0, 1, 0], [x, 1, 0, 0]]")
})

test_that("intpoly: gradient_terms 3d", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(0, 3), function(a) lapply(a, function(b) lapply(b, .pyL))) })), "[[[[1, 0, 0, 0, 0, 0, 0, 0]]]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(1, 3), function(a) lapply(a, function(b) lapply(b, .pyL))) })), "[[[[1, 0, 0, 0, 0, 0, 0, 0]]], [[[y, 0, 1, 0, 1, 0, 0, 0], [z, 0, 0, 1, 1, 0, 1, 0]], [[x, 1, 0, 0, 1, 1, 0, 0]]]]")
})

test_that("intpoly: misc", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(c(1, 2)) })), "sqrt(5)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(c(3, 4, 12)) })), "13")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(list(1, 1)) })), "sqrt(2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(Point(3, 4)) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(Point(1, 2, 3)) })), "2*sqrt(2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(list(Rational(1, 2), sqrt(Integer(3))/2)) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_norm(list(x, y)) })), "sqrt(x^2 + y^2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); is_vertex(2) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); is_vertex(c(2, 3)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); is_vertex(Point(2, 3)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); is_vertex(c(2, 3, 4)) })), "True")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); is_vertex(c(2, 3, 4, 5)) })), "False")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(7)) })), "(7, 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(0)) })), "(0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(x)) })), "(1, x)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_cross_product(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0)) })), "[0, 0, 1]")
})

