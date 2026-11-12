# Reference values from sympy.sets.ordinals (powers print as ^ here).
w <- ord_omega

test_that("ordinal arithmetic and printing follow SymPy", {
  cases <- list(
    list(w, "w"), list(ord0, "ord0"), list(Ordinal(OmegaPower(3, 2), OmegaPower(1, 4), OmegaPower(0, 7)), "w^3*2 + w*4 + 7"),
    list(w + 1, "w + 1"), list(1 + w, "w"), list(w + w, "w*2"), list(w*2, "w*2"), list(2*w, "w"), list(w*w, "w^2"),
    list(w^w, "w^(w)"), list((w + 3)*(w + 2), "w^2 + w*2 + 3"), list((w*2 + 1)*3, "w*6 + 1"),
    list(w^2 + w*3 + 5, "w^2 + w*3 + 5"), list(Ordinal(OmegaPower(w + 1, 1)), "w^(w + 1)"),
    list(Ordinal(OmegaPower(w, 2)), "w^(w)*2"), list(ord0 + w, "w"), list(w*ord0, "ord0"),
    list((w^2*3 + w + 5) + (w*2 + 1), "w^2*3 + w*3 + 1"), list((w + 2)*w, "w^2"), list(ordinal_degree(w^2 + 3), "2"))
  for (cs in cases) expect_equal(format(cs[[1L]]), cs[[2L]])
  expect_true(w < w + 1)
  expect_false(w + 1 < w)
  expect_true(3 < w)
  expect_true(w == 1 + w)
  expect_true(is_successor_ordinal(w + 1))
  expect_true(is_limit_ordinal(w))
  expect_error(Ordinal(OmegaPower(1, 1), OmegaPower(2, 1)), "decreasing order")
  expect_error(OmegaPower(1, 0), "positive integer")
})

test_that("tensordiagonal and sparse array aliases follow SymPy", {
  flat <- function(a) vapply(a$data, sstr, "")
  A <- Array(as.list(0:26), c(3, 3, 3)); B <- Array(as.list(0:8), c(3, 3)); C <- Array(as.list(0:15), c(2, 2, 2, 2))
  expect_equal(flat(tensordiagonal(B, c(1, 2))), c("0", "4", "8"))
  expect_equal(flat(tensordiagonal(A, c(1, 2))), as.character(c(0, 12, 24, 1, 13, 25, 2, 14, 26)))
  expect_equal(flat(tensordiagonal(A, c(2, 3))), as.character(c(0, 4, 8, 9, 13, 17, 18, 22, 26)))
  expect_equal(flat(tensordiagonal(A, c(1, 2, 3))), c("0", "13", "26"))
  expect_equal(flat(tensordiagonal(A, c(1, 3))), as.character(c(0, 10, 20, 3, 13, 23, 6, 16, 26)))
  r <- tensordiagonal(C, c(1, 2), c(3, 4))
  expect_equal(r$shape, c(2L, 2L)); expect_equal(flat(r), c("0", "3", "12", "15"))
  expect_error(tensordiagonal(B, 1), "at least two axes")
  expect_equal(ImmutableSparseNDimArray(list(1, 0, 0, 2), c(2, 2))$shape, c(2L, 2L))
})
