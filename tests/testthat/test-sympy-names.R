# Names of SymPy's public API added for completeness; reference values are
# SymPy's outputs.

test_that("kinds follow SymPy", {
  x <- Symbol("x"); A <- Symbol("A", commutative = FALSE)
  expect_equal(format(expr_kind(x + 1)), "NumberKind")
  expect_equal(format(expr_kind(sin(x))), "UndefinedKind")
  expect_equal(format(expr_kind(A * x)), "UndefinedKind")
  expect_equal(format(expr_kind(x > 1)), "BooleanKind")
  expect_equal(format(expr_kind(Matrix(list(1, x)))), "MatrixKind(NumberKind)")
  expect_equal(format(expr_kind(FiniteSet(1, x))), "SetKind(NumberKind)")
  expect_equal(format(expr_kind(S_EmptySet)), "SetKind()")
  expect_equal(format(expr_kind(Tuple(1, 2))), "TupleKind(NumberKind, NumberKind)")
  expect_equal(format(expr_kind(Ket("k") * x)), "KetKind")
  expect_equal(format(expr_kind(Mod(x, 2))), "NumberKind")
  expect_true(MatrixKind() == MatrixKind(NumberKind))
})

test_that("exception classes inherit as in SymPy", {
  r <- tryCatch(stop(DMNotAField("no")), DMDomainError = function(e) "domain", error = function(e) "other")
  expect_equal(r, "domain")
  r <- tryCatch(stop(MissingUnityError("u")), StructureError = function(e) "structure")
  expect_equal(r, "structure")
  r <- tryCatch(stop(HilbertSpaceError("h")), QuantumError = function(e) "quantum")
  expect_equal(r, "quantum")
})

test_that("aliases and constants", {
  expect_equal(sstr(Px), "Px")
  expect_equal(length(allhints), 46L)
  expect_equal(allhints[1:3], c("factorable", "nth_algebraic", "separable"))
  expect_equal(length(all_assumptions), 29L)
  x <- Symbol("x")
  expect_equal(sstr(Relational(x, 1, "<")), "x < 1")
  expect_equal(sstr(RealNumber(1.5)), "1.50000000000000")
})

test_that("binary relation predicates", {
  x <- Symbol("x")
  expect_equal(sstr(AppliedBinaryRelation(StrictLessThanPredicate(), x, 2)), "Q.lt(x, 2)")
  p <- binrel_properties(Q$gt)
  expect_false(p$is_reflexive)
  expect_equal(format(p$reversed), "Q.lt")
  expect_equal(format(p$negated), "Q.le")
  expect_true(test_closed_group(x + 1, Q$positive(x), Q$positive))
})

test_that("stochastic process base classes and sample_iter", {
  s <- StochasticProcess("X", c(1, 2))
  expect_equal(sstr(attr(s, "node")), "StochasticProcess(X, {1, 2})")
  d <- DiscreteTimeStochasticProcess("Y", c(1, 2))
  expect_equal(sstr(rv_getitem(d, 1)), "Y[1]")
  v <- sample_iter(Die("X", 6), numsamples = 4, seed = 1)
  expect_length(v, 4L)
  expect_true(all(vapply(v, function(e) as.numeric(e) %in% 1:6, NA)))
})

test_that("predefined manifolds of diffgeom.rn", {
  expect_equal(sstr(R2), "R^2")
  expect_equal(sstr(R2_r), "rectangular")
  expect_equal(sstr(R3_s), "spherical")
  expect_equal(sstr(BaseScalarField(R2_r, 0)), "x")
  expect_equal(sstr(BaseVectorField(R2_r, 1)), "e_y")
})

test_that("interval arithmetic follows intervalmath", {
  i <- interval(1, 2)
  expect_equal(format(i + 3), "interval(4.000000, 5.000000)")
  expect_equal(format(3 - i), "interval(1.000000, 2.000000)")
  expect_equal(format(i * interval(-1, 2)), "interval(-2.000000, 4.000000)")
  expect_equal(format(i / interval(-2, -1)), "interval(-2.000000, -0.500000)")
  s <- sqrt(interval(-1, 4))
  expect_true(is.na(s$is_valid))
  expect_equal(format(cos(interval(0, 3))), "interval(-0.989992, 1.000000)")
  expect_equal(format(imax(interval(0, 1), interval(-1, 3))), "interval(0.000000, 3.000000)")
  expect_equal(format(i > 0), "intervalMembership(True, True)")
  expect_equal(format(i < interval(1.5, 3)), "intervalMembership(None, True)")
  expect_equal(format(interval(-1, 2)^2), "interval(0.000000, 4.000000)")
  expect_equal(format(log(interval(-1, 1))), "interval(-inf, inf)")
  expect_equal(format(ceil(interval(0.5, 1.5))), "interval(1.000000, 2.000000)")
})

test_that("number field utilities", {
  x <- Symbol("x")
  zp <- AlgIntPowers(Poly(cyclotomic_poly(5, x)))
  expect_equal(as.integer(zp[[4]]), c(-1L, -1L, -1L, -1L))
  expect_equal(as.integer(zp[[7]]), c(0L, 0L, 1L, 0L))
  nxt <- coeff_search(2, 1)
  got <- t(sapply(1:13, function(i) nxt()))
  expect_equal(got[, 1], c(1, 1, 1, 0, 2, 2, 2, 2, 2, 1, 1, 0, 3))
  expect_equal(got[, 2], c(1, 0, -1, 1, 2, 1, 0, -1, -2, 2, -2, 2, 3))
})

test_that("multiple dispatch", {
  ns <- new.env()
  f <- dispatch("Integer", fname = "f", fun = function(x) x + 1, namespace = ns)
  f <- dispatch("Rational", fname = "f", fun = function(x) x * 2, namespace = ns)
  expect_equal(sstr(f(Integer(3))), "4")
  expect_equal(sstr(f(Rational(1, 3))), "2/3")
  expect_error(f("a"), "Could not find signature")
})

test_that("evaluate(FALSE) builds unevaluated nodes", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(srepr(evaluate(FALSE, x + x)), "Add(Symbol('x'), Symbol('x'))")
  expect_equal(srepr(evaluate(FALSE, x - y)), "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))")
  expect_equal(srepr(evaluate(FALSE, x / y)), "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))")
  expect_equal(sstr(evaluate(FALSE, sin(Integer(0)))), "sin(0)")
  expect_equal(sstr(evaluate(FALSE, x + 1 + 2)), "(x + 1) + 2")
  expect_equal(sstr(evaluate(FALSE, x - (x + 1))), "x - (x + 1)")
  expect_equal(sstr(x + x), "2*x")
})

test_that("pretty symbology helpers", {
  expect_equal(pretty_symbol("alpha_1"), "α₁")
  expect_equal(hobj("-", 5), strrep("─", 5))
  expect_equal(center("ab", 6), "  ab  ")
  expect_equal(center_pad(2, 6), c("  ", "  "))
  expect_equal(annotated("F")[[5]], "┌─\n├─\n╵")
})

test_that("small utilities", {
  x <- Symbol("x"); y <- Symbol("y")
  X <- DeferredVector("X")
  expect_equal(lambdify(X, X[0] + 2 * X[1])(c(3, 4)), 11)
  vdiff <- vectorize(0, 1)(function(f, v) diff(f, v))
  r <- vdiff(list(x^2, x * y), list(x, y))
  expect_equal(sstr(r[[1]][[1]]), "2*x")
  expect_equal(sstr(r[[2]][[2]]), "x")
  seterr(TRUE)
  expect_error(Rational(0, 0), "Indeterminate 0/0")
  seterr(FALSE)
  expect_equal(sstr(Rational(0, 0)), "nan")
  M <- Matrix(list(list(1, sym_I), list(0, 2)))
  expect_equal(sympy_to_numpy(M), matrix(c(1, 0, 1i, 2), 2))
  expect_true(isinstance(x, "Atom"))
  expect_false(isinstance(x + 1, "Number"))
  F <- vfree_group("a, b")
  expect_equal(fp_str(a * b), "a*b")
  expect_error(Morphism(Object("A"), Object("B")), "Cannot instantiate Morphism")
  tex <- preview(x^2, output = "tex", viewer = "file", filename = tempfile(fileext = ".tex"))
  expect_true(any(grepl("displaystyle x^{2}", readLines(tex), fixed = TRUE)))
})

test_that("power_series_ring matches sympy.polys.series", {
  r <- power_series_ring("x", QQ, 6); R <- r[[1]]; x <- r[[2]]
  expect_equal(format(R), "Power Series Ring in x over QQ of size 6")
  expect_equal(format(1 + x + x^2), "1 + x + x^2")
  expect_equal(format(R$inverse(1 + x + x^2)), "1 - x + x^3 - x^4 + O(x^6)")
  expect_equal(format(R$reversion(R$sin(x))), "x + 1/6*x^3 + 3/40*x^5 + O(x^6)")
  expect_equal(format(R$tan(x + x^2)), "x + x^2 + 1/3*x^3 + x^4 + 17/15*x^5 + O(x^6)")
  expect_equal(format(R$asin(x - x^3)), "x - 5/6*x^3 - 17/40*x^5 + O(x^6)")
  expect_equal(format(R$sqrt(4 + x)), "2 + 1/4*x - 1/64*x^2 + 1/512*x^3 - 5/16384*x^4 + 7/131072*x^5 + O(x^6)")
  expect_equal(format(R$sqrt((1 + x)^2)), "1 + x")
  expect_equal(format(R$compose(R$exp(x) - 1, x + x^2)), "x + 3/2*x^2 + 7/6*x^3 + 25/24*x^4 + 27/40*x^5 + O(x^6)")
  expect_equal(format(x / (x + x^2)), "1 - x + x^2 - x^3 + x^4 - x^5 + O(x^6)")
  expect_equal(format(R$truncate(R$exp(x), 3)), "1 + x + 1/2*x^2 + O(x^3)")
  expect_null(x$prec); expect_equal(R$exp(x)$prec, 6L)
  X <- Symbol("x")
  expect_equal(format(R$from_expr(X^2 + 3 * X + O(X^3))), "3*x + x^2 + O(x^3)")
  expect_equal(sstr(R$to_expr(R$tan(x))), "x + x^3/3 + 2*x^5/15 + O(x^6)")
  expect_equal(format(x^7), "O(x^6)")
  s <- power_series_ring("y", ZZ, 4); y <- s[[2]]
  expect_equal(format(s[[1]]$inverse(1 - y)), "1 + y + y^2 + y^3 + O(y^4)")
  expect_error(s[[1]]$inverse(2 + y), class = "NotReversible")
})

test_that("set closure, interior, powerset and graph_certificate match SymPy", {
  expect_equal(sstr(closure(Interval(0, 1, TRUE, TRUE))), "Interval(0, 1)")
  expect_equal(sstr(interior(Interval(0, 1))), "Interval.open(0, 1)")
  expect_equal(sstr(powerset(FiniteSet(1, 2))), "FiniteSet(EmptySet, {1}, {2}, {1, 2})")
  expect_equal(sstr(symmetric_difference(Interval(0, 2), Interval(1, 3))),
               "Union(Interval.Ropen(0, 1), Interval.Lopen(2, 3))")
  expect_identical(OrdinalZero(), ord0)
  gr1 <- list(`0` = c(1, 2, 3, 5), `1` = c(0, 2, 4), `2` = c(0, 1, 3, 4),
              `3` = c(0, 2, 4), `4` = c(1, 2, 3, 5), `5` = c(0, 4))
  gr2 <- list(`0` = c(1, 5), `1` = c(0, 2, 3, 4), `2` = c(1, 3, 5),
              `3` = c(1, 2, 4, 5), `4` = c(1, 3, 5), `5` = c(0, 2, 3, 4))
  expect_equal(graph_certificate(gr1), c(0L, 2L, 4L, 6L, 1L, 8L, 10L, 12L, 3L, 14L, 16L, 18L, 5L, 9L, 15L, 7L, 11L, 17L, 13L, 19L, 20L, 21L))
  expect_identical(graph_certificate(gr1), graph_certificate(gr2))
})

test_that("power series elements combine with sym and gmp scalars", {
  r <- power_series_ring("x", QQ, 6); x <- r[[2]]
  expect_equal(format(Rational(1, 2) * x), "1/2*x")
  expect_equal(format(x + Integer(2)), "2 + x")
  expect_equal(format(gmp::as.bigq(1, 2) * x), "1/2*x")
  expect_equal(format(Integer(1) / (1 - x)), "1 + x + x^2 + x^3 + x^4 + x^5 + O(x^6)")
  expect_error(PowerSeriesRingRing(QQ), "Unsupported ground domain")
})
