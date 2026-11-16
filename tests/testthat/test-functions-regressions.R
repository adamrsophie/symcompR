# Regression cases for elementary and special functions, factorials and exact powers.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("special functions follow SymPy", {
  x <- Symbol("x")
  expect_identical(sstr(subfactorial(5)), "44")
  expect_identical(sstr(subfactorial(10)), "1334961")
  expect_identical(sstr(Heaviside(oo)), "1")
  expect_identical(sstr(elliptic_f(x, 0)), "x")
  expect_identical(sstr(elliptic_pi(0, x)), "elliptic_k(x)")
  expect_identical(sstr(diff(asech(x), x)), "-sqrt(1/(x + 1))/(x*sqrt(1 - x))")
  expect_identical(sstr(Li(x)), "Li(x)")
})

test_that("real_root takes the real branch for a symbolic argument", {
  skip_on_cran()
  x <- Symbol("x", real = TRUE); z <- Symbol("z"); n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(real_root(x, 3)), "Abs(x)^(1/3)*sign(x)")
  expect_equal(sstr(subs(real_root(x, 3), x, -8)), "-2")
  expect_equal(sstr(real_root(z, 3)), "Piecewise((Abs(z)^(1/3)*sign(z), Eq(im(z), 0)), (z^(1/3), True))")
  expect_equal(sstr(real_root(z, n)), "Piecewise((z^(1/n), Eq(n, -1) | Eq(n, 1)), (Abs(z)^(1/n)*sign(z), Eq(im(z), 0) & Eq(Mod(n, 2), 1)), (z^(1/n), True))")
  expect_equal(sstr(real_root(x, 2)), "sqrt(x)")
  expect_equal(sstr(real_root(Integer(-8), 3)), "-2")
})

test_that("factorials, Matrix size, is_zero_matrix and mod_inverse follow SymPy", {
  skip_on_cran()
  expect_equal(sstr(RisingFactorial(-2, 5001)), "0")
  expect_equal(sstr(FallingFactorial(2, 101)), "0")
  expect_equal(sstr(FallingFactorial(5, -2)), "1/42")
  expect_equal(sstr(RisingFactorial(2, -2)), "zoo")
  expect_error(Matrix(2, 2, list(1)), "List length should be equal to rows\\*columns")
  x <- Symbol("x"); p <- Symbol("p", positive = TRUE)
  expect_identical(is_zero_matrix(Matrix(list(x))), NA)
  expect_false(is_zero_matrix(Matrix(list(p))))
  expect_true(is_zero_matrix(Matrix(list(0, 0))))
  expect_equal(sstr(mod_inverse(2, Rational(5, 2))), "1/2")
  expect_error(mod_inverse(Rational(2, 7), 0), "does not exist")
  expect_false(RealField(200) == RealField(53))
  expect_true(RealField(dps = 30) == RealField(103))
})

test_that("trig reduction is exact for huge multiples of pi", {
  skip_on_cran()
  x <- Symbol("x")
  expect_equal(sstr(tan(Integer(10)^100 * sym_pi)), "0")
  expect_equal(sstr(tan(x + Integer("18014398509481985") * sym_pi)), "tan(x)")
  expect_equal(sstr(cot(x + Integer("18014398509481985") * sym_pi)), "cot(x)")
})

test_that("exp of an infinite multiple decides by the exact real part", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(simplify(exp(oo * (-Rational(1, 10^13) + I)))), "0")
  expect_equal(sstr(subs(acot(-x), x, 0)), "pi/2")
  expect_equal(sstr(acot(Integer(-2))), "-acot(2)")
})

test_that("loggamma at imaginary infinity, an empty RootSum and symmetrize", {
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  I <- sym_I
  oo <- S_Infinity
  expect_equal(sstr(loggamma(I * oo)), "zoo")
  expect_equal(sstr(RootSum(Integer(0), Lambda(x, x^2), x)), "nan")
  r <- symmetrize(x^2 + y^2, list(x, y))
  expect_equal(sstr(r[[1L]]), "-2*x*y + (x + y)^2")
})

test_that("binomial, harmonic, bernoulli and nsimplify evaluate as in SymPy", {
  skip_on_cran()
  n <- Symbol("n", integer = TRUE, positive = TRUE); m <- Symbol("m")
  expect_equal(sstr(binomial(n, n + 1)), "0")
  expect_equal(sstr(binomial(Rational(1, 2), Integer(3))), "1/16")
  expect_equal(sstr(binomial(m, Integer(-1))), "0")
  expect_equal(sstr(harmonic(oo, Integer(2))), "pi^2/6")
  expect_equal(sstr(harmonic(m, Integer(-2))), "m/6 + (m + 1)^3/3 - (m + 1)^2/2 + 1/6")
  expect_equal(sstr(bernoulli(Integer(3), m + 1)), "m/2 + (m + 1)^3 - 3*(m + 1)^2/2 + 1/2")
  expect_equal(sstr(nsimplify(Float(0.99999))), "99999/100000")
})

test_that("Mod of partition and binomial follows SymPy's hooks", {
  n <- Symbol("n", integer = TRUE, nonnegative = TRUE); k <- Symbol("k", integer = TRUE)
  got <- vapply(list(Mod(partition(5 * n + 4), 5), Mod(partition(7 * n + 5), 7),
                     Mod(partition(11 * n + 6), 11), Mod(partition(5 * n + 4), 7),
                     Mod(partition(n), 5)), sstr, "")
  expect_equal(got, c("0", "0", "0", "Mod(partition(5*n + 4), 7)", "Mod(partition(n), 5)"))
  expect_error(Mod(binomial(Rational(1, 2), k), 3), "Integers expected for binomial Mod")
})

test_that("factorials are exact", {
  expect_equal(sstr(sym_factorial(25)), "15511210043330985984000000")
  expect_equal(sstr(factorial2(Integer(-3))), "-1")
  expect_equal(sstr(factorial2(Integer(-5))), "1/3")
  expect_error(factorial2(Integer(-4)), "negative odd")
})
