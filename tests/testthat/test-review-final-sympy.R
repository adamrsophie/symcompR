# Reference values from SymPy 1.14 (count_roots, evaluate = False, Set operators).
test_that("count_roots takes SymPy's numeric bounds", {
  x <- Symbol("x"); f <- x^4 - 5 * x^2 + 4
  expect_equal(count_roots(f, 0, 3), 2L)
  expect_equal(count_roots(f), 4L)
  expect_equal(count_roots(f, -3, 0), 2L)
  expect_equal(count_roots(f, 0), 2L)
  expect_equal(count_roots(x^4 - 1, x, -2, 0), 1L)
})

test_that("Add, Mul and Pow honor evaluate = FALSE", {
  x <- Symbol("x")
  expect_equal(latex(Mul(2, x, evaluate = FALSE)), "2 x")
  expect_equal(sstr(Add(x, 1, x, evaluate = FALSE)), "x + x + 1")
  expect_equal(sstr(Pow(2, 2, evaluate = FALSE)), "2^2")
  expect_equal(sstr(Mul(x, 2, evaluate = FALSE)), "x*2")
})

test_that("Set operators follow Set.__and__, __or__, __sub__, __xor__ and __pow__", {
  A <- Interval(0, 2); B <- Interval(1, 3)
  expect_equal(sstr(A & B), "Interval(1, 2)")
  expect_equal(sstr(A | B), "Interval(0, 3)")
  expect_equal(sstr(A + B), "Interval(0, 3)")
  expect_equal(sstr(Interval(0, 3) - Interval(1, 2)), "Union(Interval.Ropen(0, 1), Interval.Lopen(2, 3))")
  expect_equal(sstr(A ^ B), "Union(Interval.Ropen(0, 1), Interval.Lopen(2, 3))")
  expect_equal(sstr(FiniteSet(1, 2)^2), "ProductSet({1, 2}, {1, 2})")
  expect_equal(sstr(FiniteSet(1, 2)^0), "{()}")
})

test_that("cholesky, LDLdecomposition and dot check their input as SymPy", {
  m <- function(r) paste(capture.output(print(r)), collapse = " ")
  expect_error(cholesky(Matrix(list(list(1, 2), list(3, 4)))), "Matrix must be Hermitian.", fixed = TRUE)
  expect_error(cholesky(Matrix(list(list(1, 2), list(2, 1)))), "Matrix must be positive-definite", fixed = TRUE)
  expect_equal(sstr(cholesky(Matrix(list(list(1, 2), list(2, 1))), hermitian = FALSE)), "Matrix([[1, 0], [2, sqrt(3)*I]])")
  A <- Matrix(list(list(4, 2 * S_I), list(-2 * S_I, 3)))
  expect_equal(sstr(cholesky(A)), "Matrix([[2, 0], [-I, sqrt(2)]])")
  ld <- LDLdecomposition(A)
  expect_equal(sstr(ld$L), "Matrix([[1, 0], [-I/2, 1]])")
  expect_equal(sstr(ld$D), "Matrix([[4, 0], [0, 2]])")
  expect_error(dot(Matrix(list(1, 2, 3)), Matrix(list(1, 2))), "Dimensions incorrect for dot product: (3, 1), (2, 1)", fixed = TRUE)
  expect_equal(sstr(dot(Matrix(list(1, 2, 3)), Matrix(list(list(1, 2, 3))))), "14")
})

test_that("Floats keep their multiprecision value in facts, Integer, Rational and Mod", {
  expect_false(ask_fact(Float("1e-400", 30), "zero"))
  expect_true(ask_fact(Float("1e400", 30), "finite"))
  expect_equal(sstr(Integer(Float("9007199254740993", 30))), "9007199254740993")
  expect_equal(sstr(Integer(1.9)), "1")
  expect_equal(sstr(Integer(-1.9)), "-1")
  expect_equal(sstr(Integer(Rational(3, 2))), "1")
  expect_equal(sstr(Integer(Float("-2.5"))), "-2")
  expect_equal(sstr(Rational(Float("0.5"))), "1/2")
  expect_equal(sstr(Rational(Float(0.1))), "3602879701896397/36028797018963968")
  expect_equal(sstr(sym_Mod(Float("1.00000000000000000001", 30), 1)), "9.99999999997452451411988619036e-21")
  expect_equal(sstr(sym_Mod(Float(5.5), 2)), "1.50000000000000")
  expect_equal(sstr(sym_Mod(Float("7.25", 30), Rational(1, 2))), "0.250000000000000000000000000000")
})

test_that("real_root takes the real branch for a symbolic argument", {
  x <- Symbol("x", real = TRUE); z <- Symbol("z"); n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(real_root(x, 3)), "Abs(x)^(1/3)*sign(x)")
  expect_equal(sstr(subs(real_root(x, 3), x, -8)), "-2")
  expect_equal(sstr(real_root(z, 3)), "Piecewise((Abs(z)^(1/3)*sign(z), Eq(im(z), 0)), (z^(1/3), True))")
  expect_equal(sstr(real_root(z, n)), "Piecewise((z^(1/n), Eq(n, -1) | Eq(n, 1)), (Abs(z)^(1/n)*sign(z), Eq(im(z), 0) & Eq(Mod(n, 2), 1)), (z^(1/n), True))")
  expect_equal(sstr(real_root(x, 2)), "sqrt(x)")
  expect_equal(sstr(real_root(Integer(-8), 3)), "-2")
})
