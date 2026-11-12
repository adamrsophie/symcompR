# Regression tests for the p16 review findings; values from SymPy.
test_that("trig reduction is exact for huge multiples of pi", {
  x <- Symbol("x")
  expect_equal(sstr(tan(Integer(10)^100 * sym_pi)), "0")
  expect_equal(sstr(tan(x + Integer("18014398509481985") * sym_pi)), "tan(x)")
  expect_equal(sstr(cot(x + Integer("18014398509481985") * sym_pi)), "cot(x)")
})

test_that("numeric sign certification refuses unsafe cases", {
  u <- sin(Integer(1))^2 + cos(Integer(1))^2 - exp(Integer(-100))
  # SymPy decides these with mpmath at adaptive precision
  expect_equal(sstr(Gt(floor(u), 0)), "False")
  expect_equal(sstr(cos(Integer("36028797019032683")) > 0), "True")
  expect_equal(sstr(Gt(exp(sym_pi) - sym_pi^sym_E, 0)), "True")
})

test_that("noncommutative products keep their order", {
  x <- Symbol("x"); y <- Symbol("y")
  A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE)
  X <- Symbol("X", commutative = FALSE); Y <- Symbol("Y", commutative = FALSE)
  expect_equal(sstr(expand((1 + X)^(-1) * (1 + Y)^(-1))), "(1 + X)^(-1)*(1 + Y)^(-1)")
  expect_equal(sstr(expand(X^(-1) * (Y + 1)^(-1))), "X^(-1)*(1 + Y)^(-1)")
  expect_equal(sstr(TR0(A*B - B*A)), "A*B - B*A")
  expect_equal(sstr(TR9(cos(x)*A*B + cos(y)*B*A)), "cos(x)*A*B + cos(y)*B*A")
  expect_equal(sstr(x*A + B), "x*A + B")
})

test_that("as_numer_denom, count_ops and TableForm follow SymPy", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_equal(sstr(TR0(x^oo)), "x^oo")
  expect_equal(count_ops(Max(x, y, z)), 2)
  expect_equal(sstr(count_ops(Min(x, y, z), visual = TRUE)), "2*MIN")
  expect_equal(count_ops(Piecewise(list(1, Lt(x, 0)), list(2, TRUE))), 1)
  v <- gmp::as.bigz("12345678901234567890")
  expect_equal(.py_format_int("%.3d", v), "12345678901234567890")
  expect_equal(.py_format_int("%.5d", gmp::as.bigz(12)), "00012")
  expect_equal(.py_format_int("id=%d", v), "id=12345678901234567890")
  expect_equal(.py_format_int("%05d", gmp::as.bigz(-42)), "-0042")
})

test_that("AccumBounds, hankel_transform and ordinals follow SymPy", {
  r <- Symbol("r", positive = TRUE)
  expect_equal(sstr(Float(0) / AccumBounds(-1, 1)), "0")
  expect_equal(sstr(1 / AccumBounds(Float(0), 2)), "AccumBounds(1/2, oo)")
  expect_equal(sstr(log(AccumBounds(Float(0), 2))), "AccumBounds(-oo, log(2))")
  expect_equal(sstr(hankel_transform(r^(-4), r, -1, 1)), "-oo")
  expect_equal(sstr(hankel_transform(r^(-4), r, 0, 1)), "0")
  w <- ord_omega
  expect_equal(format(Ordinal()), "ord0")
  expect_equal(format(w + Integer(1)), "w + 1")
  expect_equal(format(Integer(2) + w), "w")
  expect_true(OmegaPower(1, 1) < OmegaPower(2, 1))
  expect_true(OmegaPower(1, 1) == OmegaPower(1, 1))
})

test_that("blockcut, slices, ask, python() and sqf_norm follow SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  n <- Symbol("n", integer = TRUE, positive = TRUE)
  B <- blockcut(MatrixSymbol("A", n + 1, 2), list(n, 1), c(1, 1))
  expect_equal(vapply(mat_shape(B), sstr, ""), c("n + 1", "2"))
  N <- Integer("9007199254740993")
  expect_equal(sstr(mat_shape(blockcut(MatrixSymbol("A", N, 1), list(N), list(1)))[[1L]]), "9007199254740993")
  m <- as_explicit(MatrixSlice(ZeroMatrix(2, n), c(0, 1), c(0, 1)))
  expect_equal(c(m$nrow, m$ncol, sstr(m$data[[1L]])), c("1", "1", "0"))
  expect_false(ask(Q$commutative(x), !Q$commutative(x)))
  expect_false(ask(Q$negative(x + y), Q$positive(x + y)))
  expect_error(ask(Q$positive(x + y), Q$positive(x + y), !Q$positive(x + y)), "inconsistent")
  expect_equal(python(Symbol("x_{1}") + 1), "x_1 = Symbol('x_{1}')\ne = x_1 + 1")
  expect_equal(python(Symbol("a^{2}") * Symbol("lambda")), "a^2 = Symbol('a^{2}')\nlambda_ = Symbol('lambda')\ne = a**2*lambda_")
  a <- sqrt(Integer(2)) + sqrt(Integer(3))
  s <- sqf_norm(x + 1/a, a)
  expect_equal(c(sstr(s[[1L]][[1L]]), sstr(s[[2L]]), sstr(s[[3L]])), c("0", "x - sqrt(2) + sqrt(3)", "x^4 - 10*x^2 + 1"))
  s <- sqf_norm(x*y + 1, sqrt(Integer(2)))
  expect_equal(c(vapply(s[[1L]], sstr, ""), sstr(s[[2L]]), sstr(s[[3L]])),
               c("1", "0", "x*y - sqrt(2)*y + 1", "x^2*y^2 + 2*x*y - 2*y^2 + 1"))
})
