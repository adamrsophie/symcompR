# Rows checked against SymPy 1.15: assumptions of re, im and Abs, the
# unevaluated MatAdd and MatMul constructors, and the print order of sums
# with matrix scalars.

test_that("re, im and Abs are extended real, finite only for a finite argument", {
  x <- Symbol("x"); r <- Symbol("r", real = TRUE); f <- Symbol("f", finite = TRUE)
  facts <- function(e) vapply(c("real", "extended_real", "finite", "nonnegative", "extended_nonnegative", "zero"),
                              function(k) ask_fact(e, k), NA)
  na <- NA
  expect_equal(unname(facts(re(x))), c(na, TRUE, na, na, na, na))
  expect_equal(unname(facts(im(x))), c(na, TRUE, na, na, na, na))
  expect_equal(unname(facts(re(f))), c(TRUE, TRUE, TRUE, na, na, na))
  expect_equal(unname(facts(Abs(x))), c(na, TRUE, na, na, TRUE, na))
  expect_equal(unname(facts(Abs(r))), c(TRUE, TRUE, TRUE, TRUE, TRUE, na))
  expect_equal(sstr(Contains(x, Complexes)), "Contains(x, Complexes)")
})

test_that("Min and Max compare by the extended sign", {
  x <- Symbol("x"); w <- Symbol("w")
  expect_equal(sstr(Max(0, 2 * Abs(im(w)))), "2*Abs(im(w))")
  expect_equal(sstr(Min(0, -Abs(x))), "-Abs(x)")
  expect_equal(sstr(Max(0, Abs(x))), "Abs(x)")
})

test_that("MatAdd and MatMul keep their arguments until doit", {
  A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
  got <- vapply(list(MatAdd(A, A), MatMul(2, A, 3), MatMul(A, A), MatMul(2, 3),
                     doit(MatMul(2, A, 3)), doit(MatAdd(A, A)), A + A, -A, A - B, A * 2, A / 2), sstr, "")
  expect_equal(got, c("A + A", "2*A*3", "A*A", "6", "6*A", "2*A", "2*A", "-A", "A - B", "2*A", "(1/2)*A"))
})

test_that("matrix scalars sort before numbers in a sum", {
  A <- MatrixSymbol("A", 2, 2); x <- Symbol("x")
  got <- vapply(list(Trace(A) + 5, MatrixElement(A, 0, 0) + 1, x + MatrixElement(A, 0, 0) + 1,
                     Trace(A) - 1, -Trace(A) + x), sstr, "")
  expect_equal(got, c("Trace(A) + 5", "A[0, 0] + 1", "x + A[0, 0] + 1", "Trace(A) - 1", "x - Trace(A)"))
})

test_that("relationals of Abs evaluate by extended sign facts", {
  x <- Symbol("x")
  got <- vapply(list(Gt(-Abs(x), 0), Lt(-Abs(x), 0), Ge(Abs(x), 0), Le(-Abs(x), 0)), sstr, "")
  expect_equal(got, c("False", "-Abs(x) < 0", "True", "True"))
})

test_that("Mod of partition and binomial follows SymPy's hooks", {
  n <- Symbol("n", integer = TRUE, nonnegative = TRUE); k <- Symbol("k", integer = TRUE)
  got <- vapply(list(Mod(partition(5 * n + 4), 5), Mod(partition(7 * n + 5), 7),
                     Mod(partition(11 * n + 6), 11), Mod(partition(5 * n + 4), 7),
                     Mod(partition(n), 5)), sstr, "")
  expect_equal(got, c("0", "0", "0", "Mod(partition(5*n + 4), 7)", "Mod(partition(n), 5)"))
  expect_error(Mod(binomial(Rational(1, 2), k), 3), "Integers expected for binomial Mod")
})

test_that("Subs placeholders keep symbol facts once noncommutative symbols exist", {
  x <- Symbol("x")
  Symbol("A_nc", commutative = FALSE)
  expect_equal(sstr(doit(Subs(sin(x) * x, x, 0))), "0")
  expect_equal(sstr(doit(Subs(x^2 + 1, x, 2))), "5")
})
