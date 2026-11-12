# Expected values are SymPy outputs for the cases of sympy/sets/tests/test_fancysets.py
# (ComplexRegion, normalize_theta_set).

test_that("ComplexRegion matches SymPy", {
  x <- Symbol("x"); r <- Symbol("r", real = TRUE); I <- sym_I; pi <- sym_pi
  Iv <- function(a, b, lo = FALSE, ro = FALSE) Interval(a, b, lo, ro)
  PS <- ProductSet
  a <- Iv(2, 3); b <- Iv(4, 6); c <- Iv(7, 9)
  c1 <- ComplexRegion(PS(a, b)); c2 <- ComplexRegion(Union(PS(a, b), PS(c, a)))
  c3 <- ComplexRegion(PS(Iv(0, 1), Iv(0, 2 * pi)), polar = TRUE)
  c4 <- ComplexRegion(PS(Iv(0, 3), Iv(pi, 2 * pi, TRUE)), polar = TRUE)
  out <- character()
  P <- function(e) out <<- c(out, tryCatch(if (is.list(e) && !is_sym(e)) paste0("(", paste(vapply(e, sstr, ""), collapse = ", "), ")") else sstr(e), error = function(err) paste("ERR", conditionMessage(err))))
  for (v in list(Float(2.5) + Float(4.5) * I, 2 + 4 * I, 3 + 4 * I)) P(contains(c1, v))
  P(contains(c2, 8 + Float(2.5) * I)); P(contains(c1, Float(2.5) + Float(6.1) * I)); P(contains(c1, Float(4.5) + Float(3.2) * I))
  P(contains(c1, x)); P(contains(c1, r)); P(contains(c2, r))
  for (v in list(Float(0.5) + I * Rational(6, 10), Rational(1, 2) + Float(0.6) * I, I, Integer(1), Integer(0), 1 + I, 1 - I)) P(contains(c3, v))
  P(contains(c3, r + 2 * I)); P(contains(c3, 1 / (1 + r^2)))
  for (v in list(Integer(0), 2 + I, -2 + I, -2 - I, 2 - I, Integer(-2), Integer(2))) P(contains(c4, v))
  P(contains(c4, 3 / (1 + r^2)))
  P(c3); P(c4)
  X <- ComplexRegion(PS(Iv(0, oo), FiniteSet(0, pi)), polar = TRUE)
  ud <- ComplexRegion(PS(Iv(0, 1), Iv(0, 2 * pi)), polar = TRUE)
  uhd <- ComplexRegion(PS(Iv(0, oo), Iv(0, pi)), polar = TRUE)
  lhd <- ComplexRegion(PS(Iv(0, oo), Iv(pi, 2 * pi)), polar = TRUE)
  rhd <- ComplexRegion(PS(Iv(0, oo), Iv(-pi / 2, pi / 2)), polar = TRUE)
  fqd <- ComplexRegion(PS(Iv(0, oo), Iv(0, pi / 2)), polar = TRUE)
  P(Intersection(uhd, ud)); P(Intersection(rhd, fqd)); P(Intersection(uhd, rhd)); P(Intersection(uhd, lhd)); P(X); P(rhd)
  k1 <- ComplexRegion(PS(Iv(0, 4), Iv(0, 2 * pi)), polar = TRUE)
  P(Intersection(k1, Iv(1, 5))); P(Intersection(k1, Iv(4, 9))); P(Intersection(k1, Iv(5, 12)))
  us <- ComplexRegion(PS(Iv(-1, 1), Iv(-1, 1))); uhp <- ComplexRegion(PS(Iv(-oo, oo), Iv(0, oo)))
  lhp <- ComplexRegion(PS(Iv(-oo, oo), Iv(-oo, 0))); rhp <- ComplexRegion(PS(Iv(0, oo), Iv(-oo, oo)))
  fqp <- ComplexRegion(PS(Iv(0, oo), Iv(0, oo)))
  P(Intersection(uhp, us)); P(Intersection(rhp, fqp)); P(Intersection(uhp, rhp)); P(Intersection(uhp, lhp))
  k2 <- ComplexRegion(PS(Iv(-5, 5), Iv(-10, 10)))
  P(Intersection(k2, Iv(2, 7))); P(Intersection(k2, Iv(5, 7))); P(Intersection(k2, Iv(6, 9)))
  P(Intersection(ud, us))
  cc2 <- ComplexRegion(PS(Iv(0, 1), Iv(0, pi)), polar = TRUE)
  P(Union(ud, cc2)); P(Union(uhd, lhd))
  c5 <- ComplexRegion(PS(Iv(2, 5), Iv(6, 9))); c6 <- ComplexRegion(PS(Iv(4, 6), Iv(10, 12)))
  c7 <- ComplexRegion(PS(Iv(0, 10), Iv(-10, 0))); c8 <- ComplexRegion(PS(Iv(12, 16), Iv(14, 20)))
  P(Union(c5, c6)); P(Union(c7, c8)); P(Union(ud, Iv(2, 4))); P(Union(c5, Iv(2, 4)))
  P(complex_region_from_real(Iv(-1, 1)))
  cm1 <- ComplexRegion(PS(Iv(2, 5), Iv(4, 8)))
  cm2 <- ComplexRegion(Union(PS(Iv(2, 5), Iv(0, 2 * pi)), PS(Iv(4, 8), Iv(0, pi))), polar = TRUE)
  P(measure(cm1)); P(measure(cm2))
  R <- Rational
  for (t in list(Iv(pi, 2 * pi), Iv(pi * R(9, 2), 5 * pi), Iv(pi * R(-3, 2), pi / 2),
    Iv(pi * R(-3, 2), pi / 2, TRUE, TRUE), Iv(pi * R(-7, 2), pi * R(-3, 2), TRUE, TRUE),
    Iv(-pi / 2, pi / 2), Iv(-pi / 2, pi / 2, TRUE, TRUE), Iv(-4 * pi, 3 * pi), Iv(pi * R(-3, 2), -pi / 2),
    Iv(0, 2 * pi, TRUE, TRUE), Iv(-pi / 2, pi / 2, FALSE, TRUE), Iv(-pi / 2, pi / 2, TRUE, FALSE), Iv(4 * pi, pi * R(9, 2), TRUE, TRUE),
    Iv(4 * pi, pi * R(9, 2), TRUE, FALSE), Iv(4 * pi, pi * R(9, 2), FALSE, TRUE), Iv(3 * pi, 5 * pi, TRUE, TRUE),
    FiniteSet(0, pi, 3 * pi), FiniteSet(0, pi / 2, pi, 2 * pi), FiniteSet(0, -pi / 2, -pi, -2 * pi), FiniteSet(pi * R(-3, 2), pi / 2),
    FiniteSet(2 * pi), Union(Iv(0, pi / 3), Iv(pi / 2, pi)), Union(Iv(0, pi), Iv(2 * pi, pi * R(7, 3)))))
    P(normalize_theta_set(t))
  y <- Symbol("y"); z <- Symbol("z"); aa <- Symbol("a"); bb <- Symbol("b"); cc <- Symbol("c")
  P(ComplexRegion(PS(FiniteSet(aa, bb, cc), FiniteSet(x, y, z)))); P(ComplexRegion(PS(FiniteSet(2), FiniteSet(3))))
  P(a_interval(c2)); P(b_interval(c2)); P(psets(c2)); P(c3$args[[1]]); P(c2$args[[1]])
  P(Union(Complexes, Iv(1, 2)))
  P(Intersection(c1, Reals)); P(Intersection(c3, Reals))
  expected <- c(
    "True",
    "True",
    "True",
    "True",
    "False",
    "False",
    "Contains(x, CartesianComplexRegion(ProductSet(Interval(2, 3), Interval(4, 6))))",
    "False",
    "False",
    "True",
    "True",
    "True",
    "True",
    "True",
    "False",
    "False",
    "Contains(r + 2*I, PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi))))",
    "Contains(1/(r^2 + 1), PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi))))",
    "True",
    "False",
    "False",
    "True",
    "True",
    "False",
    "True",
    "Contains(3/(r^2 + 1), PolarComplexRegion(ProductSet(Interval(0, 3), Union({0}, Interval.open(pi, 2*pi)))))",
    "PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))",
    "PolarComplexRegion(ProductSet(Interval(0, 3), Union({0}, Interval.open(pi, 2*pi))))",
    "PolarComplexRegion(ProductSet(Interval(0, 1), Interval(0, pi)))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), Interval(0, pi/2)))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), Interval(0, pi/2)))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), {0, pi}))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), {0, pi}))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), Union(Interval(0, pi/2), Interval.Ropen(3*pi/2, 2*pi))))",
    "Interval(1, 4)",
    "{4}",
    "EmptySet",
    "CartesianComplexRegion(ProductSet(Interval(-1, 1), Interval(0, 1)))",
    "CartesianComplexRegion(ProductSet(Interval(0, oo), Interval(0, oo)))",
    "CartesianComplexRegion(ProductSet(Interval(0, oo), Interval(0, oo)))",
    "CartesianComplexRegion(ProductSet(Interval(-oo, oo), {0}))",
    "Interval(2, 5)",
    "{5}",
    "EmptySet",
    "Intersection(CartesianComplexRegion(ProductSet(Interval(-1, 1), Interval(-1, 1))), PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi))))",
    "PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))",
    "PolarComplexRegion(ProductSet(Interval(0, oo), Interval.Ropen(0, 2*pi)))",
    "CartesianComplexRegion(Union(ProductSet(Interval(2, 5), Interval(6, 9)), ProductSet(Interval(4, 6), Interval(10, 12))))",
    "CartesianComplexRegion(Union(ProductSet(Interval(0, 10), Interval(-10, 0)), ProductSet(Interval(12, 16), Interval(14, 20))))",
    "Union(Interval(2, 4), PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi))))",
    "CartesianComplexRegion(Union(ProductSet(Interval(2, 4), {0}), ProductSet(Interval(2, 5), Interval(6, 9))))",
    "CartesianComplexRegion(ProductSet(Interval(-1, 1), {0}))",
    "12",
    "9*pi",
    "Union({0}, Interval.Ropen(pi, 2*pi))",
    "Interval(pi/2, pi)",
    "Interval.Ropen(0, 2*pi)",
    "Union(Interval.Ropen(0, pi/2), Interval.open(pi/2, 2*pi))",
    "Union(Interval.Ropen(0, pi/2), Interval.open(pi/2, 2*pi))",
    "Union(Interval(0, pi/2), Interval.Ropen(3*pi/2, 2*pi))",
    "Union(Interval.Ropen(0, pi/2), Interval.open(3*pi/2, 2*pi))",
    "Interval.Ropen(0, 2*pi)",
    "Interval(pi/2, 3*pi/2)",
    "Interval.open(0, 2*pi)",
    "Union(Interval.Ropen(0, pi/2), Interval.Ropen(3*pi/2, 2*pi))",
    "Union(Interval(0, pi/2), Interval.open(3*pi/2, 2*pi))",
    "Interval.open(0, pi/2)",
    "Interval.Lopen(0, pi/2)",
    "Interval.Ropen(0, pi/2)",
    "Union(Interval.Ropen(0, pi), Interval.open(pi, 2*pi))",
    "{0, pi}",
    "{0, pi/2, pi}",
    "{0, pi, 3*pi/2}",
    "{pi/2}",
    "{0}",
    "Union(Interval(0, pi/3), Interval(pi/2, pi))",
    "Interval(0, pi)",
    "{a + I*x, a + I*y, a + I*z, b + I*x, b + I*y, b + I*z, c + I*x, c + I*y, c + I*z}",
    "{2 + 3*I}",
    "Union(Interval(2, 3), Interval(7, 9))",
    "Union(Interval(2, 3), Interval(4, 6))",
    "(ProductSet(Interval(2, 3), Interval(4, 6)), ProductSet(Interval(7, 9), Interval(2, 3)))",
    "ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi))",
    "Union(ProductSet(Interval(2, 3), Interval(4, 6)), ProductSet(Interval(7, 9), Interval(2, 3)))",
    "Complexes",
    "EmptySet",
    "Interval(-1, 1)"
  )
  expect_equal(out, expected)
})
