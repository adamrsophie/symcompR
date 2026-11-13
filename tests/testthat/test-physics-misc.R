test_that("hydrogen energies match SymPy", {
  n <- Symbol("n"); Z <- Symbol("Z")
  expect_identical(sstr(hydrogen_E_nl(n, Z)), "-Z^2/(2*n^2)")
  expect_identical(sstr(hydrogen_E_nl(2)), "-1/8")
  expect_error(hydrogen_E_nl(-1), "'n' must be positive integer")
})

test_that("Wigner symbols match SymPy", {
  expect_identical(sstr(wigner_3j(2, 6, 4, 0, 0, 0)), "sqrt(715)/143")
  expect_identical(sstr(clebsch_gordan(1, 1, 1, 1, 0, 1)), "sqrt(2)/2")
  expect_identical(sstr(wigner_6j(3, 3, 3, 3, 3, 3)), "-1/14")
})

test_that("Pauli and Dirac matrices match SymPy", {
  expect_identical(sstr(msigma(3)), "Matrix([[1, 0], [0, -1]])")
  expect_identical(sstr(minkowski_tensor), "Matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])")
})
