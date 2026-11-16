# Regression cases for assumption facts: re, im, Abs, Min, Max, relationals and numeric sign certification.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("exp(y) is never zero for finite y and an integral of zero is zero", {
  y <- Symbol("y"); x <- Symbol("x")
  expect_false(ask_fact(exp(y), "zero"))
  expect_true(ask_fact(Integral(0, x), "zero"))
  expect_equal(sstr(exp(-Integral(0, x))), "1")
})

test_that("numeric sign certification refuses unsafe cases", {
  skip_on_cran()
  u <- sin(Integer(1))^2 + cos(Integer(1))^2 - exp(Integer(-100))
  # SymPy decides these with mpmath at adaptive precision
  expect_equal(sstr(Gt(floor(u), 0)), "False")
  expect_equal(sstr(cos(Integer("36028797019032683")) > 0), "True")
  expect_equal(sstr(Gt(exp(sym_pi) - sym_pi^sym_E, 0)), "True")
})

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

test_that("relationals of Abs evaluate by extended sign facts", {
  x <- Symbol("x")
  got <- vapply(list(Gt(-Abs(x), 0), Lt(-Abs(x), 0), Ge(Abs(x), 0), Le(-Abs(x), 0)), sstr, "")
  expect_equal(got, c("False", "-Abs(x) < 0", "True", "True"))
})

test_that("Subs placeholders keep symbol facts once noncommutative symbols exist", {
  x <- Symbol("x")
  Symbol("A_nc", commutative = FALSE)
  expect_equal(sstr(doit(Subs(sin(x) * x, x, 0))), "0")
  expect_equal(sstr(doit(Subs(x^2 + 1, x, 2))), "5")
})
