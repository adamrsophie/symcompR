# Reference values from SymPy's Idx, IndexedBase, get_indices and
# get_contraction_structure.
test_that("Indexed objects print as in SymPy", {
  i <- Idx("i"); j <- Idx("j"); n <- Symbol("n", integer = TRUE)
  A <- IndexedBase("A"); x <- IndexedBase("x")
  expect_equal(srepr(A(i)), "Indexed(IndexedBase(Symbol('A')), Idx(Symbol('i', integer=True)))")
  expect_equal(srepr(Idx("i", n)), "Idx(Symbol('i', integer=True), Tuple(Integer(0), Add(Integer(-1), Symbol('n', integer=True))))")
  expect_equal(latex(A(i, j)), "{A}_{i,j}")
  expect_equal(sstr(Indexed("A", i, j)), "A[i, j]")
  expect_equal(sstr(A(i + 1)), "A[i + 1]")
  expect_equal(sstr(x(i) * A(i, i)), "A[i, i]*x[i]")
  expect_equal(srepr(Symbol("x", positive = TRUE, integer = TRUE)), "Symbol('x', integer=True, positive=True)")
})

test_that("get_indices and get_contraction_structure follow SymPy", {
  i <- Idx("i"); j <- Idx("j")
  A <- IndexedBase("A"); x <- IndexedBase("x"); y <- IndexedBase("y")
  gi <- function(e) vapply(get_indices(e)$indices, sstr, "")
  expect_length(gi(A(i, i)), 0L)
  expect_equal(gi(A(i, j) * x(j)), "i")
  expect_equal(gi(x(i) + y(i)), "i")
  expect_equal(gi(A(i, j)), c("i", "j"))
  expect_length(gi(exp(x(i) * y(i))), 0L)
  expect_error(get_indices(x(i) + y(j)), class = "IndexConformanceException")
  f <- function(e) format_contraction_structure(get_contraction_structure(e))
  expect_equal(f(x(i)*y(i) + A(j, j)), "{(i,): {x[i]*y[i]}, (j,): {A[j, j]}}")
  expect_equal(f(x(i)*y(j)), "{None: {x[i]*y[j]}}")
  expect_equal(f(1 + A(i, i)*x(i)), "{None: {1, A[i, i]*x[i]}, A[i, i]*x[i]: [{(i,): {A[i, i]}}]}")
})
