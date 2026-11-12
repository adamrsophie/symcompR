# Reference values from SymPy's AccumBounds.
AB <- AccumBounds
x <- Symbol("x"); y <- Symbol("y", real = TRUE)
oo <- S_Infinity

test_that("construction and printing", {
  expect_equal(sstr(AB(1, 2)), "AccumBounds(1, 2)")
  expect_equal(sstr(AB(1, 1)), "1")
  expect_error(AB(2, 1), "Lower limit")
  expect_error(AB(sym_I, 2), "Only real")
  expect_equal(latex(AB(1, 2)), "\\left\\langle 1, 2\\right\\rangle")
})

test_that("interval arithmetic", {
  cases <- list(
    list(AB(1, 2) + 3, "AccumBounds(4, 5)"), list(AB(1, 2) + AB(2, 3), "AccumBounds(3, 5)"),
    list(AB(1, 2) + x, "x + AccumBounds(1, 2)"), list(x + AB(1, 2) + 1, "x + AccumBounds(2, 3)"),
    list(AB(1, 2) - AB(0, 5), "AccumBounds(-4, 2)"), list(-AB(1, 2), "AccumBounds(-2, -1)"),
    list(AB(1, 2) + oo, "oo"), list(AB(-Inf, 2) + oo, "AccumBounds(-oo, oo)"), list(AB(1, Inf) - 5, "AccumBounds(-4, oo)"),
    list(AB(0, 1) + sym_pi/2, "AccumBounds(pi/2, 1 + pi/2)"), list(AB(1, 2) + y, "AccumBounds(y + 1, y + 2)"),
    list(2*AB(1, 2), "AccumBounds(2, 4)"), list(-3*AB(1, 2), "AccumBounds(-6, -3)"), list(x*AB(1, 2), "x*AccumBounds(1, 2)"),
    list(AB(1, 2)*sym_pi, "AccumBounds(pi, 2*pi)"), list(AB(-1, 2)*AB(3, 4), "AccumBounds(-4, 8)"),
    list(AB(0, 2)*oo, "AccumBounds(0, oo)"), list(AB(-2, 0)*oo, "AccumBounds(-oo, 0)"),
    list(0*AB(1, 2), "0"), list(0*AB(1, Inf), "AccumBounds(0, oo)"),
    list(AB(1, 2)/2, "AccumBounds(1/2, 1)"), list(AB(1, 2)/AB(3, 4), "AccumBounds(1/4, 2/3)"),
    list(1/AB(1, 2), "AccumBounds(1/2, 1)"), list(1/AB(-1, 2), "AccumBounds(-oo, oo)"),
    list(1/AB(0, 2), "AccumBounds(1/2, oo)"), list(2/AB(-2, 0), "AccumBounds(-oo, -1)"),
    list(AB(1, 2)^2, "AccumBounds(1, 4)"), list(AB(-1, 2)^2, "AccumBounds(0, 4)"), list(AB(-2, -1)^3, "AccumBounds(-8, -1)"),
    list(AB(-1, 2)^3, "AccumBounds(-1, 8)"), list(AB(0, 2)^(-2), "AccumBounds(1/4, oo)"),
    list(AB(1, 4)^Rational(1, 2), "AccumBounds(1, 2)"), list(AB(1, 2)^oo, "AccumBounds(0, oo)"),
    list(AB(0, Rational(1, 2))^oo, "0"), list(2^AB(1, 3), "AccumBounds(2, 8)"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]), cs[[2L]])
})

test_that("functions, comparisons, sets and limits", {
  cases <- list(
    list(exp(AB(0, 1)), "AccumBounds(1, E)"), list(sym_E^AB(0, 1), "AccumBounds(1, E)"), list(log(AB(1, sym_E)), "AccumBounds(0, 1)"),
    list(sin(oo), "AccumBounds(-1, 1)"), list(cos(-oo), "AccumBounds(-1, 1)"), list(tan(oo), "AccumBounds(-oo, oo)"),
    list(sin(AB(0, sym_pi/2)), "AccumBounds(0, 1)"), list(sin(AB(0, 2*sym_pi)), "AccumBounds(-1, 1)"),
    list(sin(AB(-Inf, 2)), "AccumBounds(-1, 1)"), list(cos(AB(0, sym_pi)), "AccumBounds(-1, 1)"),
    list(cos(AB(0, sym_pi/2)), "AccumBounds(0, 1)"), list(tan(AB(0, sym_pi/4)), "AccumBounds(0, 1)"),
    list(cot(AB(sym_pi/4, sym_pi/2)), "AccumBounds(0, 1)"), list(Abs(AB(-3, 2)), "Abs(AccumBounds(-3, 2))"),
    list(floor(AB(1, 2)), "floor(AccumBounds(1, 2))"),
    list(Lt(AB(1, 2), 3), "True"), list(Lt(AB(1, 2), 2), "AccumBounds(1, 2) < 2"), list(Gt(AB(1, 2), AB(3, 4)), "False"),
    list(intersection(AB(1, 3), AB(2, 4)), "AccumBounds(2, 3)"), list(intersection(AB(1, 3), FiniteSet(2, 5)), "{2}"),
    list(ab_union(AB(1, 2), AB(3, 4)), "AccumBounds(1, 4)"),
    list(limit(sin(x), x, Inf), "AccumBounds(-1, 1)"), list(limit(sin(1/x), x, 0), "AccumBounds(-1, 1)"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]), cs[[2L]])
  expect_true(contains(AB(1, 3), 2))
  expect_false(contains(AB(1, 3), 4))
})

test_that("oo absorbs factors of known sign, as Mul.flatten", {
  xp <- Symbol("x", positive = TRUE)
  expect_equal(sstr(S_Infinity*sym_pi), "oo")
  expect_equal(sstr(-S_Infinity/(2*sym_pi)), "-oo")
  expect_equal(sstr(S_Infinity*xp), "oo")
  expect_equal(sstr(-S_Infinity*xp*x), "-oo*x")
  expect_equal(sstr(floor(S_NegativeInfinity)), "-oo")
})

test_that("p15 review fixes for AccumBounds", {
  y <- Symbol("y", real = TRUE)
  s1 <- sin(Integer(1))
  cases <- list(
    list(Float(0)*AB(1, Inf), "AccumBounds(0, oo)"), list(0/AB(-1, 1), "0"), list(AB(1, Inf)/(-S_Infinity), "AccumBounds(-oo, 0)"),
    list(AB(-1, 0)/AB(0, 1), "AccumBounds(-oo, 0)"), list(AB(-Inf, 0)^Integer("9007199254740993"), "AccumBounds(-oo, 0)"),
    list(exp(log(AB(-1, 1))), "nan"), list(1/subs(AB(-y^2, y^2), y, 0), "zoo"), list(subs(AB(-y^2, y^2), y, 0), "0"),
    list(limit(AB(1, 2)*x, x, 0), "0"), list(limit(sin(x)*log(x)/x, x, Inf), "0"), list(limit(sin(x)/sqrt(x), x, Inf), "0"),
    list(AB(s1, 2*s1)*AB(1, 2), "AccumBounds(sin(1), 4*sin(1))"), list(AB(Symbol("q", extended_real = TRUE, finite = FALSE), Inf), "AccumBounds(q, oo)"),
    list(Min(s1, 2*s1), "sin(1)"), list(Gt(s1, 0), "True"), list(sin(S_NaN), "nan"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]), cs[[2L]])
})
