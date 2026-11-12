# Reference values from SymPy.
x <- Symbol("x"); y <- Symbol("y")

test_that("B-splines and interpolating splines", {
  expect_equal(sstr(bspline_basis(2, c(0, 1, 2, 3), 0, x)),
               "Piecewise((x^2/2, (x >= 0) & (x <= 1)), (-x^2 + 3*x - 3/2, (x >= 1) & (x <= 2)), (x^2/2 - 3*x + 9/2, (x >= 2) & (x <= 3)), (0, True))")
  bs <- bspline_basis_set(1, c(0, 1, 2, 3), x)
  expect_equal(sstr(bs[[2L]]), "Piecewise((x - 1, (x >= 1) & (x <= 2)), (3 - x, (x >= 2) & (x <= 3)), (0, True))")
  expect_equal(sstr(interpolating_spline(3, x, c(-2, 0, 1, 3, 4), c(4, 2, 1, 1, 3))),
               "Piecewise((7*x^3/117 + 7*x^2/117 - 131*x/117 + 2, (x >= -2) & (x <= 1)), (10*x^3/117 - 2*x^2/117 - 122*x/117 + 77/39, (x >= 1) & (x <= 4)))")
})

test_that("small helpers", {
  p <- polarify(exp(x)*y)
  expect_equal(sstr(p[[1L]]), "_y*exp(_x)")
  expect_identical(polarify(x, subs = FALSE), x)
  expect_equal(sstr(sym_prod(list(1, 2, 3))), "6")
  expect_equal(sym_reshape(1:6, 2), list(list(1L, 2L), list(3L, 4L), list(5L, 6L)))
  expect_error(sym_reshape(1:6, 4), "divides the sequence length")
  r <- sqf_norm(x^2 - 2, sqrt(Integer(3)))
  expect_equal(sstr(r[[2L]]), "x^2 - 2*sqrt(3)*x + 1")
  expect_equal(sstr(r[[3L]]), "x^4 - 10*x^2 + 1")
  expect_error(gff(x^3 - x), "symbolic falling factorial")
  expect_equal(igrevlex(c(1, 2)), c(-3, 2, 1))
  expect_equal(carmichael$find_first_n_carmichaels(3), c(561, 1105, 1729))
})

test_that("p14 review fixes: sqf_norm, polarify, splines, reshape", {
  r <- sqf_norm(x - sqrt(Integer(2)), sqrt(Integer(2)))
  expect_equal(sstr(r[[3L]]), "x^2 - 2")
  expect_equal(sstr(polarify(x, lift = TRUE)), "polar_lift(x)")
  expect_equal(sstr(polarify(-1, subs = FALSE)), "exp_polar(I*pi)")
  expect_equal(sstr(polarify(x + 1, lift = TRUE)), "polar_lift(x + 1)")
  expect_equal(sstr(polarify(2*x + 3)[[1L]]), "2*_x + 3")
  expect_equal(sstr(polarify(exp(x)*y, lift = TRUE)), "exp(polar_lift(x))*polar_lift(y)")
  expect_equal(sstr(polarify(sin(x) + x^2)[[1L]]), "_x^2 + sin(_x)")
  expect_equal(sym_reshape(1:2, list(1, list(0))), list(list(1L, list()), list(2L, list())))
  expect_equal(sym_reshape(1:2, list(1, list(0), 1)), list(list(1L, list(), 2L)))
  expect_equal(sstr(bspline_basis(1, c(0, 1, 2), 0, 1)), "1")
  xp <- Symbol("x", positive = TRUE)
  expect_equal(sstr(subs(bspline_basis(1, c(0, 1, 2), 0, xp), xp, Rational(1, 2))), "1/2")
})

test_that("p15 review fixes: sqf_norm over composite extensions, polarify of integrals", {
  s2 <- sqrt(Integer(2)); s3 <- sqrt(Integer(3))
  r <- sqf_norm(x - s2, 2*s2)
  expect_equal(sstr(r[[3L]]), "x^2 - 2")
  r <- sqf_norm(x + s2, s2 + s3)
  expect_equal(sstr(r[[1L]][[1L]]), "2")
  expect_equal(sstr(r[[2L]]), "x - 2*sqrt(3) - sqrt(2)")
  expect_equal(sstr(r[[3L]]), "x^4 - 28*x^2 + 100")
  f <- Function("f")
  expect_equal(sstr(polarify(Integral(x, list(x, 1, y)), lift = TRUE)), "Integral(polar_lift(x), (x, exp_polar(0), polar_lift(y)))")
  expect_equal(sstr(polarify(Integral(f(x), x), lift = TRUE)), "Integral(f(polar_lift(x)), x)")
})
