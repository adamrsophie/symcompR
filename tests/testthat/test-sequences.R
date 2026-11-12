# Reference values from SymPy.
n <- Symbol("n"); m <- Symbol("m"); x <- Symbol("x")
L <- function(v) paste0("[", paste(vapply(v, sstr, ""), collapse = ", "), "]")

test_that("SeqFormula and SeqPer terms and printing", {
  s <- sym_sequence(n^2)
  expect_equal(sstr(s), "SeqFormula(n^2, (n, 0, oo))")
  expect_equal(latex(s), "\\left[0, 1, 4, 9, \\ldots\\right]")
  expect_equal(L(sequence_terms(s, 5)), "[0, 1, 4, 9, 16]")
  expect_equal(sstr(coeff(s, 3)), "9")
  p <- sym_sequence(list(1, 2, 3))
  expect_equal(sstr(p), "SeqPer((1, 2, 3), (_k, 0, oo))")
  expect_equal(L(sequence_terms(p, 7)), "[1, 2, 3, 1, 2, 3, 1]")
  expect_equal(L(sequence_terms(SeqPer(list(n, n^2), list(n, 0, Inf)), 5)), "[0, 1, 2, 9, 4]")
  expect_equal(L(sequence_terms(sym_sequence(list(1, 2, 3), list(n, -Inf, 0)), 5)), "[1, 2, 3, 1, 2]")
  expect_equal(latex(sym_sequence(list(1, 2, 3), list(n, -Inf, 0))), "\\left[\\ldots, 1, 3, 2, 1\\right]")
  expect_equal(latex(SeqFormula(n^2, list(n, -Inf, 0))), "\\left[\\ldots, 9, 4, 1, 0\\right]")
  expect_equal(latex(SeqFormula(n^2, list(n, 0, 3))), "\\left[0, 1, 4, 9\\right]")
  expect_identical(sym_sequence(n^2, list(n, 3, 2)), EmptySequence)
  expect_length(free_symbols(sym_sequence(n^2)), 0L)
  expect_equal(vapply(free_symbols(sym_sequence(m*n^2, list(n, 0, Inf))), sstr, ""), "m")
  expect_error(coeff(SeqFormula(n^2, list(n, 0, 5)), 7), "out of bounds")
})

test_that("sequence arithmetic reduces like SymPy", {
  expect_equal(sstr(sym_sequence(n^2) + sym_sequence(n^3)), "SeqFormula(n^3 + n^2, (n, 0, oo))")
  expect_equal(sstr(sym_sequence(n^2) * sym_sequence(n)), "SeqFormula(n^3, (n, 0, oo))")
  expect_equal(sstr(-sym_sequence(n)), "SeqFormula(-n, (n, 0, oo))")
  expect_equal(sstr(sym_sequence(n) - sym_sequence(n^2)), "SeqFormula(-n^2 + n, (n, 0, oo))")
  expect_equal(sstr(SeqPer(list(1, 2), list(n, 0, Inf)) + SeqPer(list(3, 4, 5), list(n, 0, Inf))),
               "SeqPer((4, 6, 6, 5, 5, 7), (n, 0, oo))")
  expect_equal(sstr(SeqPer(list(1, 2), list(n, 0, Inf)) * SeqPer(list(3, 4, 5), list(n, 0, Inf))),
               "SeqPer((3, 8, 5, 6, 4, 10), (n, 0, oo))")
  s2 <- SeqFormula(n^2, list(n, 0, 5)) + SeqPer(list(1, 2), list(n, 0, Inf))
  expect_equal(sstr(s2), "SeqAdd(SeqFormula(n^2, (n, 0, 5)), SeqPer((1, 2), (n, 0, oo)))")
  expect_equal(L(sequence_terms(s2, 4)), "[1, 3, 5, 11]")
  expect_error(2 * sym_sequence(n), "cannot multiply sequence")
})

test_that("find_linear_recurrence", {
  expect_equal(L(find_linear_recurrence(sym_sequence(fibonacci(n)), 10)), "[1, 1]")
  r <- find_linear_recurrence(sym_sequence(fibonacci(n)), 10, gfvar = x)
  expect_equal(sstr(r[[2L]]), "-x/(x^2 + x - 1)")
  expect_equal(sstr(Min(oo, 5)), "5")
  expect_equal(sstr(Max(-oo, x)), "x")
})

test_that("approximants match SymPy", {
  t <- Symbol("t")
  expect_equal(L(approximants(list(1, 1, 2, 3, 5, 8, 13, 21))), "[1, -1/(x - 1), (x - 1)/(2*x - 1), -1/(x^2 + x - 1)]")
  expect_equal(L(approximants(list(1, t, t^2, t^3), X = x)), "[1, -1/(t*x - 1)]")
  expect_equal(L(approximants(list(1, 2, 3, 4, 5, 6))), "[1, -1/(2*x - 1), (x + 2)/(2 - 3*x), 1/(x^2 - 2*x + 1)]")
  expect_equal(sstr(approximants(list(1, 1, 2, 5, 14, 42, 132))[[4L]]), "(1 - 2*x)/(x^2 - 3*x + 1)")
})
