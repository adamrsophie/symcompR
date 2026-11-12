# Reference values from SymPy's ComplexRegion and imageset.
I <- sym_I; pi <- sym_pi
a <- Interval(2, 3); b <- Interval(4, 6); c <- Interval(1, 8)

test_that("ComplexRegion construction, membership and set operations follow SymPy", {
  C1 <- ComplexRegion(ProductSet(a, b))
  C2 <- ComplexRegion(Union(ProductSet(a, b), ProductSet(b, c)))
  P <- ComplexRegion(ProductSet(Interval(0, 1), Interval(0, 2*pi)), polar = TRUE)
  cases <- list(
    list(C1, "CartesianComplexRegion(ProductSet(Interval(2, 3), Interval(4, 6)))"),
    list(C2, "CartesianComplexRegion(Union(ProductSet(Interval(2, 3), Interval(4, 6)), ProductSet(Interval(4, 6), Interval(1, 8))))"),
    list(P, "PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))"),
    list(ComplexRegion(ProductSet(Reals, Reals)), "Complexes"),
    list(ComplexRegion(ProductSet(FiniteSet(1, 2), FiniteSet(0, 1))), "{1, 2, 1 + I, 2 + I}"),
    list(a_interval(C1), "Interval(2, 3)"), list(b_interval(C2), "Interval(1, 8)"),
    list(Intersection(C1, ComplexRegion(ProductSet(Interval(2, 5), Interval(5, 7)))), "CartesianComplexRegion(ProductSet(Interval(2, 3), Interval(5, 6)))"),
    list(Intersection(C1, Reals), "EmptySet"),
    list(Intersection(ComplexRegion(ProductSet(Interval(1, 3), Interval(-1, 1))), Reals), "Interval(1, 3)"),
    list(Intersection(ComplexRegion(ProductSet(Interval(1, 2), Interval(0, pi)), polar = TRUE), Reals), "Union(Interval(-2, -1), Interval(1, 2))"),
    list(Union(C1, ComplexRegion(ProductSet(c, a))), "CartesianComplexRegion(Union(ProductSet(Interval(1, 8), Interval(2, 3)), ProductSet(Interval(2, 3), Interval(4, 6))))"),
    list(Union(C1, Interval(0, 1)), "CartesianComplexRegion(Union(ProductSet(Interval(0, 1), {0}), ProductSet(Interval(2, 3), Interval(4, 6))))"),
    list(complex_region_from_real(Interval(1, 2)), "CartesianComplexRegion(ProductSet(Interval(1, 2), {0}))"),
    list(ComplexRegion(ProductSet(Interval(0, 1), Interval(-pi/2, pi/2)), polar = TRUE),
         "PolarComplexRegion(ProductSet(Interval(0, 1), Union(Interval(0, pi/2), Interval.Ropen(3*pi/2, 2*pi))))"),
    list(ComplexRegion(ProductSet(Interval(0, 1), Interval(0, 3*pi)), polar = TRUE), "PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]), cs[[2L]])
  expect_true(contains(C1, Rational(5, 2) + Rational(9, 2)*I))
  expect_false(contains(C1, Rational(5, 2) + Rational(13, 2)*I))
  expect_true(contains(P, Rational(1, 2) + Rational(1, 2)*I))
  expect_false(contains(P, 1 + 2*I))
  expect_true(contains(P, 0))
  expect_equal(latex(C1), "\\left\\{x + y i\\; \\middle|\\; x, y \\in \\left[2, 3\\right] \\times \\left[4, 6\\right] \\right\\}")
})

test_that("imageset evaluates images of intervals under real polynomials", {
  x <- Symbol("x", real = TRUE)
  expect_equal(sstr(imageset(Lambda(x, -x), Interval(1, 2))), "Interval(-2, -1)")
  expect_equal(sstr(imageset(Lambda(x, 2*x + 1), Interval(0, 1, TRUE, FALSE))), "Interval.Lopen(1, 3)")
  expect_equal(sstr(imageset(Lambda(x, -3*x + 1), Interval(0, Inf))), "Interval(-oo, 1)")
  expect_equal(sstr(imageset(Lambda(x, x^2), Interval(-1, 1, TRUE, TRUE))), "Interval.Ropen(0, 1)")
  expect_equal(sstr(Intersection(ImageSet(Lambda(x, -x), Interval(1, 2)), Reals)), "Interval(-2, -1)")
})
