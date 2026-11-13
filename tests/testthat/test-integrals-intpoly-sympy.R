skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/integrals/tests/test_intpoly.py);
# the rows come from the battery that compares this port with SymPy.

test_that("intpoly: decompose matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x)) })), "{1: x}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^2)) })), "{2: x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x*y)) })), "{2: x*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x + y)) })), "{1: x + y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^2 + y)) })), "{1: y, 2: x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(8*x^2 + 4*y + 7)) })), "{0: 7, 1: 4*y, 2: 8*x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^2 + 3*y*x)) })), "{2: x^2 + 3*x*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(9*x^2 + y + 4*x + x^3 + y^2*x + 3)) })), "{0: 3, 1: 4*x + y, 2: 9*x^2, 3: x^3 + x*y^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(Integer(5))) })), "{0: 5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^3*y^2 + z)) })), "{1: z, 5: x^3*y^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x*y*z + x^3 - 2)) })), "{0: -2, 3: x^3 + x*y*z}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(3*x^2*y)) })), "{3: 3*x^2*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^5)) })), "{5: x^5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x + y + z + 1)) })), "{0: 1, 1: x + y + z}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(x^2*y^2 + x*y + x + 1)) })), "{0: 1, 1: x, 2: x*y, 4: x^2*y^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(intpoly_decompose(y^4 - x^4 + x^2*y^2)) })), "{4: -x^4 + x^2*y^2 + y^4}")
})

test_that("intpoly: decompose separate matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x, TRUE)) })), "{x}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^2, TRUE)) })), "{x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x*y, TRUE)) })), "{x*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x + y, TRUE)) })), "{x, y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^2 + y, TRUE)) })), "{x^2, y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(8*x^2 + 4*y + 7, TRUE)) })), "{4*y, 7, 8*x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^2 + 3*y*x, TRUE)) })), "{3*x*y, x^2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(9*x^2 + y + 4*x + x^3 + y^2*x + 3, TRUE)) })), "{3, 4*x, 9*x^2, x^3, x*y^2, y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(Integer(5), TRUE)) })), "{5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^3*y^2 + z, TRUE)) })), "{x^3*y^2, z}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x*y*z + x^3 - 2, TRUE)) })), "{-2, x^3, x*y*z}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(3*x^2*y, TRUE)) })), "{3*x^2*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^5, TRUE)) })), "{x^5}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x + y + z + 1, TRUE)) })), "{1, x, y, z}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(x^2*y^2 + x*y + x + 1, TRUE)) })), "{1, x, x^2*y^2, x*y}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyset(intpoly_decompose(y^4 - x^4 + x^2*y^2, TRUE)) })), "{-x^4, x^2*y^2, y^4}")
})

test_that("intpoly: best_origin matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(2, 1), 3, Segment(Point(0, 3), Point(1, 1)), y ^ 2 * x ^ 5 + y ^ 5 * x ^ 7 + 7 * x + x ^ 12 + y ^ 7 * x)) })), "(0, 3.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(2, 0), 3, Segment(Point(Integer(3)/2, 0), Point(Integer(3)/2, 3)), x ^ 7)) })), "(1.5, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 2), 3, Segment(Point(0, Integer(3)/2), Point(3, Integer(3)/2)), x ^ 7)) })), "(0, 1.5)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 2, Segment(Point(0, 2), Point(2, 0)), x ^ 7 * y ^ 3)) })), "(0, 2.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 2, Segment(Point(0, 2), Point(2, 0)), x ^ 3 * y ^ 7)) })), "(2.0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 2, Segment(Point(0, 2), Point(1, 1)), x ^ 2 * y ^ 9)) })), "(0, 2.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 2, Segment(Point(2, 0), Point(1, 1)), x ^ 9 * y ^ 2)) })), "(2.0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), x ^ 2 * y)) })), "(0, 1.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), x * y ^ 2)) })), "(1.0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), 5)) })), "(-1, 2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), x ^ 3)) })), "(0, 1.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), y)) })), "(1.0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 1), 1, Segment(Point(-1, 2), Point(3, -2)), 3*x*y)) })), "(0, 1.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, -1), 0, Segment(Point(1, 1), Point(3, 3)), x^2 + y)) })), "(1, 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 1), 1, Segment(Point(-2, 1), Point(2, 1)), x)) })), "(0, 1.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(0, 1), 1, Segment(Point(Integer(3)/2, 0), Point(Integer(3)/2, 3)), x)) })), "(3/2, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 0), 2, Segment(Point(0, Integer(3)/2), Point(3, Integer(3)/2)), y)) })), "(0, 3/2)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(1, 0), 2, Segment(Point(-2, 1), Point(2, 1)), y^2)) })), "(-2, 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(2, 1), 3, Segment(Point(0, 3), Point(1, 1)), x^2*y^2)) })), "(0, 3.0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(best_origin(c(2, 1), 3, Segment(Point(0, 3), Point(1, 1)), x^3*y + y^4)) })), "(0, 3.0)")
})

test_that("intpoly: polytope tri1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), 1) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x) })), "16/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), y) })), "8/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x*y) })), "8/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x^2 + x*y + y^2) })), "16")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), 6*x^2 - 40*y) })), "-128/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x^3) })), "128/5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x^2*y^3 + 2*x) })), "1376/105")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x^4 - y^4 + 7) })), "92")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(4, 0)), x + y + 1) })), "12")
})

test_that("intpoly: polytope sq matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), 1) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), y) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x*y) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x^2 + x*y + y^2) })), "11/12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), 6*x^2 - 40*y) })), "-18")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x^3) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x^2*y^3 + 2*x) })), "13/12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x^4 - y^4 + 7) })), "7")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x + y + 1) })), "2")
})

test_that("intpoly: polytope tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), 1) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x) })), "10")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), y) })), "35/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x*y) })), "145/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x^2 + x*y + y^2) })), "235/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), 6*x^2 - 40*y) })), "-935/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x^3) })), "78")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x^2*y^3 + 2*x) })), "19195/42")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x^4 - y^4 + 7) })), "113")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x + y + 1) })), "80/3")
})

test_that("intpoly: polytope sq3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), 1) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), y) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x*y) })), "9/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x^2 + x*y + y^2) })), "33/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), 6*x^2 - 40*y) })), "18 - 60*sqrt(3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x^3) })), "9*sqrt(3)/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x^2*y^3 + 2*x) })), "21*sqrt(3)/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x^4 - y^4 + 7) })), "21")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0)), x + y + 1) })), "3 + 3*sqrt(3)")
})

test_that("intpoly: polytope hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), 1) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), y) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x^2 + x*y + y^2) })), "sqrt(3)*(-sqrt(3)/3 - sqrt(3)*(3/2 - sqrt(3))/12 - sqrt(3)*(1/2 - sqrt(3))/12 + 4/3)/8 + sqrt(3)*(sqrt(3)*(1/2 + sqrt(3))/12 + sqrt(3)*(3/2 + sqrt(3))/12 + sqrt(3)/3 + 4/3)/8 + sqrt(3)*(-7*sqrt(3)/24 - sqrt(3)*(3/2 - sqrt(3))/12 + 37/12)/4 + sqrt(3)*(5*sqrt(3)/12 + 10/3)/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), 6*x^2 - 40*y) })), "-465*sqrt(3)/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x^3) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x^2*y^3 + 2*x) })), "71*sqrt(3)/160")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x^4 - y^4 + 7) })), "57*sqrt(3)/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x + y + 1) })), "sqrt(3)*(1 - sqrt(3)/2)/6 + sqrt(3)*(sqrt(3)/2 + 1)/6 + sqrt(3)*(7/4 - sqrt(3)/4)/3 + sqrt(3)*(sqrt(3)/4 + 7/4)/3 + 3*sqrt(3)/2")
})

test_that("intpoly: polytope nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), 1) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x) })), "-2/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x^2 + x*y + y^2) })), "2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), 6*x^2 - 40*y) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x^3) })), "-2/5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x^2*y^3 + 2*x) })), "-4/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x^4 - y^4 + 7) })), "311/15")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x + y + 1) })), "7/3")
})

test_that("intpoly: polytope nc2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), 1) })), "2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x^2 + x*y + y^2) })), "4/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), 6*x^2 - 40*y) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x^3) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x^2*y^3 + 2*x) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x^4 - y^4 + 7) })), "218/15")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(0, 0), Point(1, 1), Point(1, -1), Point(0, 0)), x + y + 1) })), "2")
})

test_that("intpoly: polytope rect matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), 1) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x) })), "9")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), y) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x*y) })), "9")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x^2 + x*y + y^2) })), "35")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), 6*x^2 - 40*y) })), "-132")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x^3) })), "81/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x^2*y^3 + 2*x) })), "54")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x^4 - y^4 + 7) })), "120")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 2), Point(3, 2), Point(3, 0)), x + y + 1) })), "21")
})

test_that("intpoly: polytope col matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), 1) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), y) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x*y) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x^2 + x*y + y^2) })), "44/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), 6*x^2 - 40*y) })), "-128")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x^3) })), "8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x^2*y^3 + 2*x) })), "56/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x^4 - y^4 + 7) })), "28")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0)), x + y + 1) })), "12")
})

test_that("intpoly: polytope dup matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), 1) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), y) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x*y) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x^2 + x*y + y^2) })), "11/12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), 6*x^2 - 40*y) })), "-18")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x^3) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x^2*y^3 + 2*x) })), "13/12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x^4 - y^4 + 7) })), "7")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0)), x + y + 1) })), "2")
})

test_that("intpoly: polytope trap matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), 1) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x) })), "12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), y) })), "16/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x*y) })), "32/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x^2 + x*y + y^2) })), "139/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), 6*x^2 - 40*y) })), "-118/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x^3) })), "78")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x^2*y^3 + 2*x) })), "7064/105")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x^4 - y^4 + 7) })), "3772/15")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0)), x + y + 1) })), "70/3")
})

test_that("intpoly: polytope pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), 1) })), "13/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x) })), "25/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), y) })), "25/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x*y) })), "39/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x^2 + x*y + y^2) })), "605/24")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), 6*x^2 - 40*y) })), "-1775/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x^3) })), "173/20")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x^2*y^3 + 2*x) })), "5629/210")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x^4 - y^4 + 7) })), "91/30")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x + y + 1) })), "19")
})

test_that("intpoly: polytope half matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), 1) })), "1/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x) })), "1/36")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), y) })), "1/24")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x*y) })), "1/144")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x^2 + x*y + y^2) })), "35/1296")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), 6*x^2 - 40*y) })), "-44/27")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x^3) })), "1/648")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x^2*y^3 + 2*x) })), "289/5184")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x^4 - y^4 + 7) })), "9059/7776")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), x + y + 1) })), "17/72")
})

test_that("intpoly: polytope fig1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(1.220, -0.827), Point(-1.490, -4.503), Point(-3.766, -1.622), Point(-4.240, -0.091), Point(-3.160, 4), Point(-0.981, 4.447), Point(0.132, 4.027)), x^2 + x*y + y^2) })), "2031627344735367/8000000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(1.220, -0.827), Point(-1.490, -4.503), Point(-3.766, -1.622), Point(-4.240, -0.091), Point(-3.160, 4), Point(-0.981, 4.447), Point(0.132, 4.027)), 1) })), "64737657/2000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(1.220, -0.827), Point(-1.490, -4.503), Point(-3.766, -1.622), Point(-4.240, -0.091), Point(-3.160, 4), Point(-0.981, 4.447), Point(0.132, 4.027)), x) })), "-99066199641/2000000000")
})

test_that("intpoly: polytope fig2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(4.561, 2.317), Point(1.491, -1.315), Point(-3.310, -3.164), Point(-4.845, -3.110), Point(-4.569, 1.867)), x^2 + x*y + y^2) })), "517091313866043/1600000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(4.561, 2.317), Point(1.491, -1.315), Point(-3.310, -3.164), Point(-4.845, -3.110), Point(-4.569, 1.867)), 1) })), "65914981/2000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(4.561, 2.317), Point(1.491, -1.315), Point(-3.310, -3.164), Point(-4.845, -3.110), Point(-4.569, 1.867)), x) })), "-9143108543/250000000")
})

test_that("intpoly: polytope fig3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-2.740, -1.888), Point(-3.292, 4.233), Point(-2.723, -0.697), Point(-0.643, -3.151)), x^2 + x*y + y^2) })), "147449361647041/8000000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-2.740, -1.888), Point(-3.292, 4.233), Point(-2.723, -0.697), Point(-0.643, -3.151)), 1) })), "3280487/2000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-2.740, -1.888), Point(-3.292, 4.233), Point(-2.723, -0.697), Point(-0.643, -3.151)), x) })), "-22047837983/6000000000")
})

test_that("intpoly: polytope fig4 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0.211, -4.622), Point(-2.684, 3.851), Point(0.468, 4.879), Point(4.630, -1.325), Point(-0.411, -1.044)), x^2 + x*y + y^2) })), "180742845225803/1000000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0.211, -4.622), Point(-2.684, 3.851), Point(0.468, 4.879), Point(4.630, -1.325), Point(-0.411, -1.044)), 1) })), "5295863/200000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0.211, -4.622), Point(-2.684, 3.851), Point(0.468, 4.879), Point(4.630, -1.325), Point(-0.411, -1.044)), x) })), "1595137013/120000000")
})

test_that("intpoly: polytope fig5 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-4.165, -0.832), Point(-3.668, 1.568), Point(-3.266, 1.279), Point(-1.090, -2.080), Point(3.313, -0.683), Point(3.033, -4.845), Point(-4.395, 4.840), Point(-1.007, -3.328)), x^2 + x*y + y^2) })), "1633405224899363/24000000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-4.165, -0.832), Point(-3.668, 1.568), Point(-3.266, 1.279), Point(-1.090, -2.080), Point(3.313, -0.683), Point(3.033, -4.845), Point(-4.395, 4.840), Point(-1.007, -3.328)), 1) })), "7446737/2000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-4.165, -0.832), Point(-3.668, 1.568), Point(-3.266, 1.279), Point(-1.090, -2.080), Point(3.313, -0.683), Point(3.033, -4.845), Point(-4.395, 4.840), Point(-1.007, -3.328)), x) })), "1071729601/100000000")
})

test_that("intpoly: polytope fig6 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-3.018, -4.473), Point(-0.103, 2.378), Point(-1.605, -2.308), Point(4.516, -0.771), Point(4.203, 0.478)), x^2 + x*y + y^2) })), "88161333955921/3000000000000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-3.018, -4.473), Point(-0.103, 2.378), Point(-1.605, -2.308), Point(4.516, -0.771), Point(4.203, 0.478)), 1) })), "194049/100000")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-3.018, -4.473), Point(-0.103, 2.378), Point(-1.605, -2.308), Point(4.516, -0.771), Point(4.203, 0.478)), x) })), "-19340473267/1500000000")
})

test_that("intpoly: clockwise ccw matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 0), Point(1, 1)), x*y, clockwise = TRUE) })), "1/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 0), Point(1, 1)), 1, clockwise = TRUE) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(1, 0), Point(1, 1)), x^2, clockwise = TRUE) })), "1/4")
})

test_that("intpoly: point_sort ccw matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(1, 0), Point(1, 1))) })), "[Point2D(1, 1), Point2D(1, 0), Point2D(0, 0)]")
})

test_that("intpoly: point_sort ccw ccw matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(1, 0), Point(1, 1)), clockwise = FALSE) })), "[Point2D(0, 0), Point2D(1, 0), Point2D(1, 1)]")
})

test_that("intpoly: clockwise sq matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x*y, clockwise = TRUE) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), 1, clockwise = TRUE) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x^2, clockwise = TRUE) })), "1/3")
})

test_that("intpoly: point_sort sq matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0))) })), "[Point2D(1, 1), Point2D(1, 0), Point2D(0, 0), Point2D(0, 1)]")
})

test_that("intpoly: point_sort ccw sq matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), clockwise = FALSE) })), "[Point2D(0, 1), Point2D(0, 0), Point2D(1, 0), Point2D(1, 1)]")
})

test_that("intpoly: clockwise tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x*y, clockwise = TRUE) })), "145/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), 1, clockwise = TRUE) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), x^2, clockwise = TRUE) })), "155/6")
})

test_that("intpoly: point_sort tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 3), Point(5, 3), Point(1, 1))) })), "[Point2D(5, 3), Point2D(1, 1), Point2D(0, 3)]")
})

test_that("intpoly: point_sort ccw tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), clockwise = FALSE) })), "[Point2D(0, 3), Point2D(1, 1), Point2D(5, 3)]")
})

test_that("intpoly: clockwise pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x*y, clockwise = TRUE) })), "39/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), 1, clockwise = TRUE) })), "13/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), x^2, clockwise = TRUE) })), "25/4")
})

test_that("intpoly: point_sort pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0))) })), "[Point2D(2, 2), Point2D(2, 0), Point2D(0, 0), Point2D(-1, 1), Point2D(0, 3)]")
})

test_that("intpoly: point_sort ccw pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), clockwise = FALSE) })), "[Point2D(0, 3), Point2D(-1, 1), Point2D(0, 0), Point2D(2, 0), Point2D(2, 2)]")
})

test_that("intpoly: clockwise hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x*y, clockwise = TRUE) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), 1, clockwise = TRUE) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), x^2, clockwise = TRUE) })), "5*sqrt(3)/16")
})

test_that("intpoly: point_sort hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2)))) })), "[Point2D(0, 2), Point2D(sqrt(3)/2, 3/2), Point2D(sqrt(3)/2, 1/2), Point2D(0, 0), Point2D(-sqrt(3)/2, 1/2), Point2D(-sqrt(3)/2, 3/2)]")
})

test_that("intpoly: point_sort ccw hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), clockwise = FALSE) })), "[Point2D(-sqrt(3)/2, 3/2), Point2D(-sqrt(3)/2, 1/2), Point2D(0, 0), Point2D(sqrt(3)/2, 1/2), Point2D(sqrt(3)/2, 3/2), Point2D(0, 2)]")
})

test_that("intpoly: clockwise nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x*y, clockwise = TRUE) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), 1, clockwise = TRUE) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), x^2, clockwise = TRUE) })), "7/6")
})

test_that("intpoly: point_sort nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1))) })), "[Point2D(1, 1), Point2D(1, -1), Point2D(0, 0), Point2D(-1, -1), Point2D(-1, 1)]")
})

test_that("intpoly: point_sort ccw nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1)), clockwise = FALSE) })), "[Point2D(-1, 1), Point2D(-1, -1), Point2D(0, 0), Point2D(1, -1), Point2D(1, 1)]")
})

test_that("intpoly: point_sort list matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(0, 0), Point(1, 0), Point(1, 1))) })), "[Point2D(1, 1), Point2D(1, 0), Point2D(0, 0)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(2, 0), Point(0, 0), Point(1, 1), Point(1, -1))) })), "[Point2D(1, 1), Point2D(2, 0), Point2D(1, -1), Point2D(0, 0)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(0, 1))) })), "[Point2D(0, 1)]")
})

test_that("intpoly: point_sort 3d matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(1, 1, 0)), normal = c(0, 0, 1)) })), "ERR TypeError: '<' not supported between instances of 'NoneType' and 'int'")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); point_sort(list(Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)), normal = c(1, 1, 1)) })), "[Point3D(0, 0, 1), Point3D(0, 1, 0), Point3D(1, 0, 0)]")
})

test_that("intpoly: clockwise error matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), 1, clockwise = TRUE) })), "ERR TypeError: clockwise=True works for only 2-PolytopeV-representation input")
})

test_that("intpoly: no expr matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0))) })), "ERR TypeError: Input expression must be a valid SymPy expression")
})

test_that("intpoly: hrep h1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), 1) })), "4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), x*y) })), "8/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), 6*x^2 - 40*y) })), "-128/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), x^2 + y) })), "40/3")
})

test_that("intpoly: hrep h2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(0, 1), 1), list(c(1, 0), 1), list(c(0, -1), 0)), 1) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(0, 1), 1), list(c(1, 0), 1), list(c(0, -1), 0)), x*y) })), "1/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(0, 1), 1), list(c(1, 0), 1), list(c(0, -1), 0)), 6*x^2 - 40*y) })), "-18")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(c(0, 1), 1), list(c(1, 0), 1), list(c(0, -1), 0)), x^2 + y) })), "5/6")
})

test_that("intpoly: hrep h3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 1), 3), list(c(1, -2), -1), list(c(-2, -1), -3)), 1) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 1), 3), list(c(1, -2), -1), list(c(-2, -1), -3)), x*y) })), "145/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 1), 3), list(c(1, -2), -1), list(c(-2, -1), -3)), 6*x^2 - 40*y) })), "-935/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 1), 3), list(c(1, -2), -1), list(c(-2, -1), -3)), x^2 + y) })), "75/2")
})

test_that("intpoly: hrep h4 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(list(0, sqrt(Integer(3))), 3), list(list(sqrt(Integer(3)), 0), 3), list(c(0, -1), 0)), 1) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(list(0, sqrt(Integer(3))), 3), list(list(sqrt(Integer(3)), 0), 3), list(c(0, -1), 0)), x*y) })), "9/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(list(0, sqrt(Integer(3))), 3), list(list(sqrt(Integer(3)), 0), 3), list(c(0, -1), 0)), 6*x^2 - 40*y) })), "18 - 60*sqrt(3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 0), list(list(0, sqrt(Integer(3))), 3), list(list(sqrt(Integer(3)), 0), 3), list(c(0, -1), 0)), x^2 + y) })), "3*sqrt(3)/2 + 3")
})

test_that("intpoly: hrep hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Rational(-1, 2), -sqrt(Integer(3)) / 2), 0), list(c(-1, 0), sqrt(Integer(3)) / 2), list(list(Rational(-1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(list(Rational(1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(c(1, 0), sqrt(Integer(3)) / 2), list(list(Rational(1, 2), -sqrt(Integer(3)) / 2), 0)), 1) })), "3*sqrt(3)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Rational(-1, 2), -sqrt(Integer(3)) / 2), 0), list(c(-1, 0), sqrt(Integer(3)) / 2), list(list(Rational(-1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(list(Rational(1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(c(1, 0), sqrt(Integer(3)) / 2), list(list(Rational(1, 2), -sqrt(Integer(3)) / 2), 0)), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Rational(-1, 2), -sqrt(Integer(3)) / 2), 0), list(c(-1, 0), sqrt(Integer(3)) / 2), list(list(Rational(-1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(list(Rational(1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(c(1, 0), sqrt(Integer(3)) / 2), list(list(Rational(1, 2), -sqrt(Integer(3)) / 2), 0)), 6*x^2 - 40*y) })), "-465*sqrt(3)/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Rational(-1, 2), -sqrt(Integer(3)) / 2), 0), list(c(-1, 0), sqrt(Integer(3)) / 2), list(list(Rational(-1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(list(Rational(1, 2), sqrt(Integer(3)) / 2), sqrt(Integer(3))), list(c(1, 0), sqrt(Integer(3)) / 2), list(list(Rational(1, 2), -sqrt(Integer(3)) / 2), 0)), x^2 + y) })), "29*sqrt(3)/16")
})

test_that("intpoly: hrep nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(0, 1), 1), list(c(1, -1), 0), list(c(1, 1), 0), list(c(0, -1), 1)), 1) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(0, 1), 1), list(c(1, -1), 0), list(c(1, 1), 0), list(c(0, -1), 1)), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(0, 1), 1), list(c(1, -1), 0), list(c(1, 1), 0), list(c(0, -1), 1)), 6*x^2 - 40*y) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(0, 1), 1), list(c(1, -1), 0), list(c(1, 1), 0), list(c(0, -1), 1)), x^2 + y) })), "5/6")
})

test_that("intpoly: hrep nc2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(1, 1), 0), list(c(-1, 1), 0), list(c(1, 0), 1), list(c(-1, -1), 0), list(c(1, -1), 0)), 1) })), "2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(1, 1), 0), list(c(-1, 1), 0), list(c(1, 0), 1), list(c(-1, -1), 0), list(c(1, -1), 0)), x*y) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(1, 1), 0), list(c(-1, 1), 0), list(c(1, 0), 1), list(c(-1, -1), 0), list(c(1, -1), 0)), 6*x^2 - 40*y) })), "6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(-1, 0), 1), list(c(1, 1), 0), list(c(-1, 1), 0), list(c(1, 0), 1), list(c(-1, -1), 0), list(c(1, -1), 0)), x^2 + y) })), "1")
})

test_that("intpoly: hrep 3d matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(1, 0, 0), 1), list(c(0, 1, 0), 1)), 1) })), "ERR NotImplementedError: Integration for H-representation 3Dcase not implemented yet.")
})

test_that("intpoly: max_degree matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), list(x^9*y + x^7*y^3 + 2*x^2*y^8, x^6*y^4 + x^5*y^5 + 2*y^10, x^10 + x^9*y + x^8*y^2 + x^5*y^5), max_degree = 10)) })), "{x^10 + x^9*y + x^8*y^2 + x^5*y^5: 1946257153/924, x^6*y^4 + x^5*y^5 + 2*y^10: 13062161/27, x^9*y + x^7*y^3 + 2*x^2*y^8: 615780107/594}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), list(x^7*y^1 + 2*x^2*y^6, x^6*y^4 + x^5*y^5 + 2*y^10), max_degree = 9)) })), "{x^7*y + 2*x^2*y^6: 489262/9}")
})

test_that("intpoly: max_degree all matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), max_degree = 1)) })), "{0: 0, 1: 1, x: 1/2, y: 1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), max_degree = 2)) })), "{0: 0, 1: 1, x: 1/2, x^2: 1/3, x*y: 1/4, y: 1/2, y^2: 1/3}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), max_degree = 3)) })), "{0: 0, 1: 1, x: 1/2, x^2: 1/3, x^2*y: 1/6, x^3: 1/4, x*y: 1/4, x*y^2: 1/6, y: 1/2, y^2: 1/3, y^3: 1/4}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), max_degree = 4)) })), "{0: 0, 1: 1, x: 1/2, x^2: 1/3, x^2*y: 1/6, x^2*y^2: 1/9, x^3: 1/4, x^3*y: 1/8, x^4: 1/5, x*y: 1/4, x*y^2: 1/6, x*y^3: 1/8, y: 1/2, y^2: 1/3, y^3: 1/4, y^4: 1/5}")
})

test_that("intpoly: max_degree polys matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(c(0, 0), c(0, 1), c(1, 1), c(1, 0)), list(1, x, y, x*y, x^2*y, x*y^2), max_degree = 1)) })), "{1: 1, x: 1/2, y: 1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(c(0, 0), c(0, 1), c(1, 1), c(1, 0)), list(1, x, y, x*y, x^2*y, x*y^2), max_degree = 2)) })), "{1: 1, x: 1/2, x*y: 1/4, y: 1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(c(0, 0), c(0, 1), c(1, 1), c(1, 0)), list(1, x, y, x*y, x^2*y, x*y^2), max_degree = 3)) })), "{1: 1, x: 1/2, x^2*y: 1/6, x*y: 1/4, x*y^2: 1/6, y: 1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(c(0, 0), c(0, 1), c(1, 1), c(1, 0)), list(1, x, y, x*y, x^2*y, x*y^2), max_degree = 4)) })), "{1: 1, x: 1/2, x^2*y: 1/6, x*y: 1/4, x*y^2: 1/6, y: 1/2}")
})

test_that("intpoly: max_degree all tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), max_degree = 3)) })), "{0: 0, 1: 5, x: 10, x^2: 155/6, x^2*y: 389/6, x^3: 78, x*y: 145/6, x*y^2: 181/3, y: 35/3, y^2: 85/3, y^3: 71}")
})

test_that("intpoly: max_degree list tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)), list(x^2 + 3, 0, x*y + y, 2*x^3), max_degree = 3)) })), "{0: 0, 2*x^3: 156, x^2 + 3: 245/6, x*y + y: 215/6}")
})

test_that("intpoly: max_degree all pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), max_degree = 3)) })), "{0: 0, 1: 13/2, x: 25/6, x^2: 25/4, x^2*y: 71/10, x^3: 173/20, x*y: 39/8, x*y^2: 47/6, y: 25/3, y^2: 169/12, y^3: 271/10}")
})

test_that("intpoly: max_degree list pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0)), list(x^2 + 3, 0, x*y + y, 2*x^3), max_degree = 3)) })), "{0: 0, 2*x^3: 173/10, x^2 + 3: 103/4, x*y + y: 317/24}")
})

test_that("intpoly: max_degree all hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), max_degree = 3)) })), "{0: 0, 1: 3*sqrt(3)/2, x: 0, x^2: 5*sqrt(3)/16, x^2*y: 5*sqrt(3)/16, x^3: 0, x*y: 0, x*y^2: 0, y: 3*sqrt(3)/2, y^2: 29*sqrt(3)/16, y^3: 39*sqrt(3)/16}")
})

test_that("intpoly: max_degree list hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2))), list(x^2 + 3, 0, x*y + y, 2*x^3), max_degree = 3)) })), "{0: 0, 2*x^3: 0, x^2 + 3: 77*sqrt(3)/16, x*y + y: 3*sqrt(3)/2}")
})

test_that("intpoly: max_degree all half matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), max_degree = 3)) })), "{0: 0, 1: 1/6, x: 1/36, x^2: 1/162, x^2*y: 1/648, x^3: 1/648, x*y: 1/144, x*y^2: 1/432, y: 1/24, y^2: 1/72, y^3: 1/192}")
})

test_that("intpoly: max_degree list half matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0)), list(x^2 + 3, 0, x*y + y, 2*x^3), max_degree = 3)) })), "{0: 0, 2*x^3: 1/324, x^2 + 3: 41/81, x*y + y: 7/144}")
})

test_that("intpoly: max_degree hrep matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(-1, 0), Integer(0)), list(c(0, 1), Integer(1)), list(c(1, 0), Integer(1)), list(c(0, -1), Integer(0))), max_degree = 2)) })), "{0: 0, 1: 1, x: 1/2, x^2: 1/3, x*y: 1/4, y: 1/2, y^2: 1/3}")
})

test_that("intpoly: issue 19234 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), list(1, x, y, x*y, x^2*y, x*y^2))) })), "{1: 1, x: 1/2, x^2*y: 1/6, x*y: 1/4, x*y^2: 1/6, y: 1/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), list(1, x, y, x*y, 3 + x^2*y, x + x*y^2))) })), "{1: 1, x: 1/2, x^2*y + 3: 19/6, x*y: 1/4, x*y^2 + x: 2/3, y: 1/2}")
})

test_that("intpoly: max_degree not list matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x, max_degree = 2) })), "ERR TypeError: 'Symbol' object is not iterable")
})

test_that("intpoly: polytope3d cube1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), 1) })), "162")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "6561")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), x) })), "378")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), x*y*z) })), "4131")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6)), c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d cube1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(1, 2, 3, 4), c(3, 2, 5, 6), c(1, 7, 5, 2), c(0, 6, 5, 7), c(1, 4, 0, 7), c(0, 4, 3, 6)), list(c(0, 0, 0), c(0, 6, 6), c(6, 6, 6), c(3, 6, 0), c(0, 6, 0), c(6, 0, 6), c(3, 0, 0), c(0, 0, 6))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([0, 1, 0], 6), ([2, 0, -1], 6), ([0, 0, 1], 6), ([0, -1, 0], 0), ([-1, 0, 0], 0), ([0, 0, -1], 0)]")
})

test_that("intpoly: polytope3d cube2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), 1) })), "125")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "15625/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), x) })), "625/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), x*y*z) })), "15625/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d cube2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([0, 1, 0], 5), ([0, 0, 1], 5), ([1, 0, 0], 5), ([0, -1, 0], 0), ([-1, 0, 0], 0), ([0, 0, -1], 0)]")
})

test_that("intpoly: polytope3d cube3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), 1) })), "105")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "33835/12")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), x) })), "1415/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), x*y*z) })), "28675/24")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5)), c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d cube3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(6, 11, 5, 0), c(1, 7, 6, 0), c(5, 4, 3, 2, 1, 0), c(11, 10, 4, 5), c(10, 9, 3, 4), c(9, 8, 2, 3), c(8, 7, 1, 2), c(7, 8, 9, 10, 11, 6)), list(c(0, 0, 0), c(5, 0, 0), c(5, 4, 0), c(3, 2, 0), c(3, 5, 0), c(0, 5, 0), c(0, 0, 5), c(5, 0, 5), c(5, 4, 5), c(3, 2, 5), c(3, 5, 5), c(0, 5, 5))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([-1, 0, 0], 0), ([0, -1, 0], 0), ([0, 0, -1], 0), ([0, 1, 0], 5), ([1, 0, 0], 3), ([-1, 1, 0], -1), ([1, 0, 0], 5), ([0, 0, 1], 5)]")
})

test_that("intpoly: polytope3d cube4 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), 1) })), "1/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "37/960")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), x) })), "11/384")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), x*y*z) })), "19/30720")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d cube4 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([0, 0, -1], 0), ([0, -1, 0], 0), ([2, 1, 1], 1), ([1, 2, 1], 1), ([0, 0, 1], 0), ([1, 1, 2], 1), ([-1, 0, 0], 0)]")
})

test_that("intpoly: polytope3d octa matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), 1) })), "sqrt(2)/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "sqrt(2)/20")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), x) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), x*y*z) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0)), c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d octa matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(3, 4, 5), c(3, 5, 1), c(3, 1, 0), c(3, 0, 4), c(4, 0, 2), c(4, 2, 5), c(2, 0, 1), c(5, 2, 1)), list(list(Integer(-1) / sqrt(Integer(2)), 0, 0), list(0, Integer(1) / sqrt(Integer(2)), 0), list(0, 0, Integer(-1) / sqrt(Integer(2))), list(0, 0, Integer(1) / sqrt(Integer(2))), list(0, Integer(-1) / sqrt(Integer(2)), 0), list(Integer(1) / sqrt(Integer(2)), 0, 0))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([1, -1, 1], sqrt(2)/2), ([1, 1, 1], sqrt(2)/2), ([-1, 1, 1], sqrt(2)/2), ([-1, -1, 1], sqrt(2)/2), ([-1, -1, -1], sqrt(2)/2), ([1, -1, -1], sqrt(2)/2), ([-1, 1, -1], sqrt(2)/2), ([1, 1, -1], sqrt(2)/2)]")
})

test_that("intpoly: polytope3d cubec matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), 1) })), "-125")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "-15625/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), x) })), "-625/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), x*y*z) })), "-15625/8")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d cubec matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)), list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([0, -1, 0], -5), ([0, 0, -1], -5), ([-1, 0, 0], -5), ([0, 1, 0], 0), ([1, 0, 0], 0), ([0, 0, 1], 0)]")
})

test_that("intpoly: polytope3d tet matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), 1) })), "1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), x ^ 2 + y ^ 2 + x * y + z ^ 2) })), "17/10")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), x) })), "1/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), x*y*z) })), "1/20")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), z^2 + 1) })), "nan")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), x^3 + y) })), "nan")
})

test_that("intpoly: hyperplanes3d tet matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(list(c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1))), function(h) .pytuple(h[[1]], h[[2]])) })), "[([0, 0, -1], 0), ([0, -1, 0], 0), ([-1, 0, 0], 0), ([3, 2, 6], 6)]")
})

test_that("intpoly: polytope3d max_degree cube2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), max_degree = 1)) })), "{1: 125, x: 625/2, y: 625/2, z: 625/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), max_degree = 2)) })), "{1: 125, x: 625/2, x^2: 3125/3, x*y: 3125/4, x*z: 3125/4, y: 625/2, y^2: 3125/3, y*z: 3125/4, z: 625/2, z^2: 3125/3}")
})

test_that("intpoly: polytope3d max_degree tet matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(0, 0, 0), c(2, 0, 0), c(0, 3, 0), c(0, 0, 1)), c(0, 2, 1), c(0, 1, 3), c(0, 3, 2), c(1, 2, 3)), max_degree = 2)) })), "{1: 1, x: 1/2, x^2: 2/5, x*y: 3/10, x*z: 1/10, y: 3/4, y^2: 9/10, y*z: 3/20, z: 1/4, z^2: 1/10}")
})

test_that("intpoly: polytope3d max_degree cube4 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), list(Integer(1) / 4, Integer(1) / 4, Integer(1) / 4)), c(0, 2, 1), c(1, 3, 0), c(4, 2, 3), c(4, 3, 1), c(0, 1, 2), c(2, 4, 1), c(0, 3, 2)), max_degree = 2)) })), "{1: 1/8, x: 11/384, x^2: 43/3840, x*y: 19/3840, x*z: 19/3840, y: 11/384, y^2: 43/3840, y*z: 19/3840, z: 11/384, z^2: 43/3840}")
})

test_that("intpoly: polytope3d list matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pySD(polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0)), list(x^2, y*z), max_degree = 2)) })), "{x^2: 3125/3, y*z: 3125/4}")
})

test_that("intpoly: polytope3d noexpr matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); polytope_integrate(list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0))) })), "ERR TypeError: Input expression must be a valid SymPy expression")
})

test_that("intpoly: helpers3d matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); main_integrate3d(1, faces, vertices, hp) })), "-125")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); .pySD(main_integrate3d(1, faces, vertices, hp, max_degree = 1)) })), "{1: -125, x: -625/2, y: -625/2, z: -625/2}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); main_integrate3d(x*y, faces, vertices, hp) })), "-3125/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); polygon_integrate(faces[[1]], list(c(0, 1, 0), 5), 0, faces, vertices, 1, 0) })), "-25")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); polygon_integrate(faces[[1]], hp[[1]], 0, faces, vertices, x, 1) })), "125/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); polygon_integrate(faces[[3]], hp[[3]], 2, faces, vertices, y*z, 2) })), "625/4")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); polygon_integrate(faces[[2]], hp[[2]], 1, faces, vertices, 0, 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); left_integral3D(faces, 0, x, vertices, hp[[1]], 1) })), "375/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)), c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0), c(3, 1, 0, 2), c(0, 4, 6, 2)); vertices <- cube[[1]]; faces <- cube[-1]; hp <- hyperplane_parameters(faces, vertices); left_integral3D(faces, 3, 1, vertices, hp[[4]], 0) })), "50")
})

test_that("intpoly: distance_to_side matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); distance_to_side(c(0, 0, 0), list(c(0, 0, 1), c(0, 1, 0)), c(1, 0, 0)) })), "-sqrt(2)/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); distance_to_side(c(1, 2, 3), list(c(0, 0, 1), c(4, 1, 0)), c(1, 2, 2)) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); distance_to_side(c(0, 0, 0), list(c(1, 0, 0), c(0, 1, 0)), c(0, 0, 1)) })), "sqrt(2)/2")
})

test_that("intpoly: lineseg_integrate matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), 1, 0) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), 1, 0) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), 0, 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), 0, 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), x, 1) })), "25/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), x, 1) })), "25/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), x^2 + y*z, 2) })), "125/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), x^2 + y*z, 2) })), "95")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0, list(c(0, 5, 0), c(5, 5, 0)), z, 1) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lineseg_integrate(NULL, 0, list(c(1, 2, 3), c(4, 6, 3)), z, 1) })), "15")
})

test_that("intpoly: helpers2d matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; main_integrate(x^2 + y^2, facets, hp) })), "325/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pySD(main_integrate(x^2 + y^2, facets, hp, max_degree = 1)) })), "{0: 0, 1: 5, x: 10, y: 35/3}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pySD(main_integrate(x^2 + y^2, facets, hp, max_degree = 2)) })), "{0: 0, 1: 5, x: 10, x^2: 155/6, x*y: 145/6, y: 35/3, y^2: 85/3}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pySD(main_integrate(list(x, y^2, 1), facets, hp)) })), "{1: 5, x: 10, y^2: 85/3}")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction(facets, 0, a, b, 1, list(x, y), 0) })), "5")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction(facets, 0, a, b, 0, list(x, y), 0) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction(facets, 1, a, b, x, list(x, y), 1) })), "6*sqrt(5)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction(facets, 2, a, b, x*y, list(x, y), 2) })), "5*sqrt(5)/6")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction(facets, 2, a, b, x^2 + y^2, list(x, y), 2) })), "14*sqrt(5)/3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction_dynamic(facets, 0, a, b, x, 1, list(x, y), 1, 0, 1, Point(0, 3), list(list(0, 0, 0, 0), list(1, 0, 0, 5), list(y, 0, 1, 15), list(x, 1, 0, NULL)), 3) })), "25/2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction_dynamic(facets, 0, a, b, 0L, 1, list(x, y), 1, 0, 1, Point(0, 3), list(list(0, 0, 0, 0), list(1, 0, 0, 5), list(y, 0, 1, 15), list(x, 1, 0, NULL)), 3) })), "0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; integration_reduction_dynamic(facets, 0, a, b, y, 1, list(x, y), 0, 1, 1, Point(0, 3), list(list(0, 0, 0, 0), list(1, 0, 0, 5), list(y, 0, 1, NULL), list(x, 1, 0, NULL)), 2) })), "15")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; left_integral2D(3, 0, facets, Point(0, 3), x*y, list(x, y)) })), "75")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; left_integral2D(3, 1, facets, Point(5, 3), 1, list(x, y)) })), "2*sqrt(5)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; lapply(hp, function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((0, 1), 3), ((1, -2), -1), ((-2, -1), -3)]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pyT(intpoly_intersection(facets[[1]], facets[[2]], "segment2D")) })), "(5, 3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pyT(intpoly_intersection(facets[[3]], facets[[1]], "segment2D")) })), "(0, 3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; intpoly_intersection(facets[[1]], facets[[1]], "segment2D") })), "None")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pyT(intpoly_intersection(hp[[1]], hp[[2]], "plane2D")) })), "(5, 3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; .pyT(intpoly_intersection(hp[[2]], hp[[3]], "plane2D")) })), "(1, 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1)); facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)), Segment(Point(1, 1), Point(0, 3))); hp <- hyperplane_parameters(tri); a <- hp[[1]][[1]]; b <- hp[[1]][[2]]; intpoly_intersection(hp[[1]], hp[[1]], "plane2D") })), "None")
})

test_that("intpoly: hyperplanes tri1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, 2), Point(4, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((1, 2), 4), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes sq matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((0, 1), 1), ((1, 0), 1), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes tri2 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 3), Point(5, 3), Point(1, 1))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((0, 1), 3), ((1, -2), -1), ((-2, -1), -3)]")
})

test_that("intpoly: hyperplanes sq3 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, sqrt(Integer(3))), Point(sqrt(Integer(3)), sqrt(Integer(3))), Point(sqrt(Integer(3)), 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((0, sqrt(3)), 3), ((sqrt(3), 0), 3), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes hex matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(-sqrt(Integer(3)) / 2, Rational(1, 2)), Point(-sqrt(Integer(3)) / 2, Integer(3) / 2), Point(0, 2), Point(sqrt(Integer(3)) / 2, Integer(3) / 2), Point(sqrt(Integer(3)) / 2, Rational(1, 2)))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1/2, -sqrt(3)/2), 0), ((-1, 0), sqrt(3)/2), ((-1/2, sqrt(3)/2), sqrt(3)), ((1/2, sqrt(3)/2), sqrt(3)), ((1, 0), sqrt(3)/2), ((1/2, -sqrt(3)/2), 0)]")
})

test_that("intpoly: hyperplanes nc1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(-1, -1), Point(-1, 1), Point(1, 1), Point(0, 0), Point(1, -1))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 1), ((0, 1), 1), ((1, -1), 0), ((1, 1), 0), ((0, -1), 1)]")
})

test_that("intpoly: hyperplanes col matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, 1), Point(0, 2), Point(2, 2), Point(2, 0), Point(1, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((0, 1), 2), ((1, 0), 2), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes dup matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((0, 1), 1), ((1, 0), 1), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes trap matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(1, 2), Point(3, 2), Point(4, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-2, 1), 0), ((0, 1), 2), ((2, 1), 8), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes pent matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(-1, 1), Point(0, 3), Point(2, 2), Point(2, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, -1), 0), ((-2, 1), 3), ((1, 2), 6), ((1, 0), 2), ((0, -1), 0)]")
})

test_that("intpoly: hyperplanes fig1 matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(1.220, -0.827), Point(-1.490, -4.503), Point(-3.766, -1.622), Point(-4.240, -0.091), Point(-3.160, 4), Point(-0.981, 4.447), Point(0.132, 4.027))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((367600, -271000), 672589), ((-1440500, -1138000), 7270759), ((-765500, -237000), 3267287), ((-102275, 27000), 431189), ((-11175, 54475), 253213), ((20000, 53000), 216071), ((606750, 136000), 627763)]")
})

test_that("intpoly: hyperplanes half matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(hyperplane_parameters(Polygon(Point(0, 0), Point(0, Rational(1, 2)), Point(Integer(1)/3, Rational(1, 2)), Point(Integer(1)/3, 0))), function(h) .pytuple(.pyT(h[[1]]), h[[2]])) })), "[((-1, 0), 0), ((0, 2), 1), ((3, 0), 1), ((0, -1), 0)]")
})

test_that("intpoly: gradient_terms matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(0), .pyL) })), "[[1, 0, 0, 0]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(1), .pyL) })), "[[1, 0, 0, 0], [y, 0, 1, 0], [x, 1, 0, 0]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(2), .pyL) })), "[[1, 0, 0, 0], [y, 0, 1, 0], [y^2, 0, 2, 0], [x, 1, 0, 0], [x*y, 1, 1, 0], [x^2, 2, 0, 0]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(3), .pyL) })), "[[1, 0, 0, 0], [y, 0, 1, 0], [y^2, 0, 2, 0], [y^3, 0, 3, 0], [x, 1, 0, 0], [x*y, 1, 1, 0], [x*y^2, 1, 2, 0], [x^2, 2, 0, 0], [x^2*y, 2, 1, 0], [x^3, 3, 0, 0]]")
})

test_that("intpoly: gradient_terms 3d matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(0, 3), function(a) lapply(a, function(b) lapply(b, .pyL))) })), "[[[[1, 0, 0, 0, 0, 0, 0, 0]]]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(1, 3), function(a) lapply(a, function(b) lapply(b, .pyL))) })), "[[[[1, 0, 0, 0, 0, 0, 0, 0]]], [[[y, 0, 1, 0, 1, 0, 0, 0], [z, 0, 0, 1, 1, 0, 1, 0]], [[x, 1, 0, 0, 1, 1, 0, 0]]]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); lapply(gradient_terms(2, 3), function(a) lapply(a, function(b) lapply(b, .pyL))) })), "[[[[1, 0, 0, 0, 0, 0, 0, 0]]], [[[y, 0, 1, 0, 1, 0, 0, 0], [z, 0, 0, 1, 1, 0, 1, 0]], [[x, 1, 0, 0, 1, 1, 0, 0]]], [[[y^2, 0, 2, 0, 2, 0, 0, 0], [y*z, 0, 1, 1, 2, 0, 1, 0], [z^2, 0, 0, 2, 2, 0, 2, 0]], [[x*y, 1, 1, 0, 2, 1, 0, 0], [x*z, 1, 0, 1, 2, 1, 1, 0]], [[x^2, 2, 0, 0, 2, 2, 0, 0]]]]")
})

test_that("intpoly: misc matches SymPy", {
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
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(3*x^2*y)) })), "(3, x^2*y)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(7)) })), "(7, 1)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(0)) })), "(0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(x)) })), "(1, x)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); .pyT(intpoly_strip(-x*y^3/2)) })), "(-1/2, x*y^3)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_cross_product(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0)) })), "[0, 0, 1]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); intpoly_cross_product(c(1, 2, 3), c(4, 5, 6), c(7, 8, 10)) })), "[3, -3, 0]")
})

