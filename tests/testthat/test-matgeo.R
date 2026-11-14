# Matrix powers, decompositions and geometry entities; reference values
# from SymPy.
x <- Symbol("x"); n <- Symbol("n", integer = TRUE, nonnegative = TRUE)
G <- function(v) if (inherits(v, "GeoEntity")) format(v) else if (is_sym(v) || inherits(v, "symmatrix")) sstr(v) else
  paste0("[", paste(vapply(v, G, ""), collapse = ", "), "]")

test_that("matrix powers follow SymPy", {
  J <- Matrix(list(list(1, 1), list(0, 1)))
  expect_equal(sstr(J %^% n), "Matrix([[1, n], [0, 1]])")
  F <- Matrix(list(list(1, 1), list(1, 0)))
  expect_equal(sstr(F %^% 10), "Matrix([[89, 55], [55, 34]])")
  expect_equal(sstr(F %^% -2), "Matrix([[1, -1], [-1, 2]])")
  D <- Matrix(list(list(4, 0), list(0, 9)))
  expect_equal(sstr(D %^% Rational(1, 2)), "Matrix([[2, 0], [0, 3]])")
})

test_that("Gram-Schmidt, QR and lambdify of matrices", {
  expect_error(GramSchmidt(list(Matrix(list(1, 2, 3)), Matrix(list(2, 4, 6)))),
               "GramSchmidt: vector set not linearly independent")
  A <- Matrix(list(list(1, 2), list(3, 4)))
  qr <- QRdecomposition(A)
  expect_equal(sstr(simplify(qr[[1L]] %*% qr[[2L]])), sstr(A))
  f <- lambdify(x, Matrix(list(list(x, 1), list(0, x))))
  expect_equal(f(2), matrix(c(2, 1, 0, 2), 2, 2, byrow = TRUE))
})

test_that("points and lines", {
  expect_equal(G(Point(0.5, 0.25)), "Point2D(1/2, 1/4)")
  expect_equal(G(Point(1, 2) + Point(3, 4)), "Point2D(4, 6)")
  expect_equal(G(Point(1, 2) * 3), "Point2D(3, 6)")
  expect_equal(G(Segment(Point(1, 1), Point(1, 1))), "Point2D(1, 1)")
  expect_true(Line(Point(0, 0), Point(1, 1)) == Line(Point(0, 0), Point(1, 1)))
  expect_true(equals(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 0), Point(2, 2))))
  expect_equal(G(intersection(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 2), Point(2, 0)))), "[Point2D(1, 1)]")
})

test_that("circles, ellipses, polygons and planes", {
  expect_equal(G(arbitrary_point(Circle(Point(0, 0), 2))), "Point2D(2*cos(t), 2*sin(t))")
  expect_equal(G(arbitrary_point(Ellipse(Point(0, 0), 3, 2))), "Point2D(3*cos(t), 2*sin(t))")
  expect_equal(G(equation(Circle(Point(1, 2), 3))), "(x - 1)^2 + (y - 2)^2 - 9")
  expect_equal(G(Circle(Point(0, 0), Point(1, 0), Point(0, 1))), "Circle(Point2D(1/2, 1/2), sqrt(2)/2)")
  C <- Circle(Point(0, 0), 1)
  expect_equal(G(intersection(C, C)), "[Circle(Point2D(0, 0), 1)]")
  expect_equal(G(Polygon(Point(0, 0), Point(1, 0), Point(2, 0), Point(1, 1))),
               "Triangle(Point2D(0, 0), Point2D(2, 0), Point2D(1, 1))")
  expect_equal(G(Triangle(Point(0, 0), Point(1, 1), Point(2, 2))), "Segment2D(Point2D(0, 0), Point2D(2, 2))")
  p1 <- Plane(Point3D(0, 0, 0), normal_vector = c(1, 0, 0))
  p2 <- Plane(Point3D(0, 0, 0), normal_vector = c(0, 1, 0))
  expect_equal(G(intersection(p1, p2)), "[Line3D(Point3D(0, 0, 0), Point3D(0, 0, 1))]")
  expect_equal(G(area(Triangle(Point(0, 0), Point(4, 0), Point(0, 3)))), "6")
})
