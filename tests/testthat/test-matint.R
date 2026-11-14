test_that("inverse dispatches to the SymPy methods", {
  A <- Matrix(list(list(1, 2), list(3, 4)))
  for (m in list(NULL, "GE", "LU", "ADJ", "BLOCK"))
    expect_identical(sstr(inverse(A, method = m)), "Matrix([[-2, 1], [3/2, -1/2]])")
  x <- Symbol("x")
  expect_identical(sstr(inverse(Matrix(list(list(x, 1), list(1, x))))),
                   "Matrix([[x/(x^2 - 1), -1/(x^2 - 1)], [-1/(x^2 - 1), x/(x^2 - 1)]])")
  expect_error(inverse(A, method = "XYZ"), "Inversion method unrecognized")
  expect_error(inverse(Matrix(list(list(1, 2), list(2, 4)))), class = "NonInvertibleMatrixError")
})

test_that("rref returns the reduced matrix and pivots", {
  r <- rref(Matrix(list(list(1, 2, 3), list(2, 4, 7))))
  expect_identical(sstr(r$matrix), "Matrix([[1, 2, 0], [0, 0, 1]])")
  expect_identical(r$pivots, c(1L, 3L))
  x <- Symbol("x")
  expect_identical(sstr(nullspace(Matrix(list(list(1, x), list(x, x^2))))[[1]]), "Matrix([[-x], [1]])")
})

test_that("eigenvals and eigenvects follow SymPy", {
  A <- Matrix(list(list(1, 2), list(3, 4)))
  ev <- eigenvals(A)
  expect_identical(vapply(ev, sstr, ""), c("5/2 - sqrt(33)/2", "5/2 + sqrt(33)/2"))
  expect_identical(attr(ev, "multiplicity"), c(1L, 1L))
  expect_length(eigenvals(Matrix(list(list(2, 0), list(0, 2))), multiple = TRUE), 2L)
  vs <- eigenvects(A)
  expect_identical(sstr(vs[[1]]$vectors[[1]]), "Matrix([[-sqrt(33)/6 - 1/2], [1]])")
  expect_error(eigenvals(Matrix(list(list(1, 2)))), class = "NonSquareMatrixError")
  d <- diagonalize(Matrix(list(list(1, 2), list(2, 1))))
  expect_identical(sstr(d$P), "Matrix([[-1, 1], [1, 1]])")
  expect_identical(sstr(d$D), "Matrix([[-1, 0], [0, 3]])")
  expect_false(is_diagonalizable(Matrix(list(list(1, 1), list(0, 1)))))
  expect_identical(vapply(singular_values(Matrix(list(list(1, 0), list(0, 2)))), sstr, ""), c("2", "1"))
  expect_identical(sstr(left_eigenvects(Matrix(list(list(0, 1), list(1, 0))))[[1]]$vectors[[1]]), "Matrix([[-1, 1]])")
})

test_that("Poly with Float coefficients keeps RR coefficients", {
  n <- Symbol("n"); x <- Symbol("x")
  expect_identical(vapply(all_coeffs(Poly(n * S_Pi - Float(0.5), n)), sstr, ""), c("1.0*pi", "-0.500000000000000"))
  expect_identical(vapply(all_coeffs(Poly(x^2 + Float(1.0), x)), sstr, ""), c("1.00000000000000", "0.0", "1.00000000000000"))
})

test_that("Poly over RR prints its domain and gives Float coefficients", {
  x <- Symbol("x")
  expect_identical(capture.output(print(Poly(Float(0.25) * x^2 + 2 * x + 1, x))),
                   "Poly(0.25*x^2 + 2.0*x + 1.0, x, domain='RR')")
  expect_identical(sstr(as_expr(Poly(x^2 / 3 + 1, x, domain = RR))), "0.333333333333333*x^2 + 1.0")
  expect_identical(.pyfmt(all_coeffs(Poly(x^2 / 4 + 1, x, domain = "RR"))),
                   "[0.250000000000000, 0.0, 1.00000000000000]")
})

test_that("CRootOf follows ComplexRootOf.__new__", {
  x <- Symbol("x")
  got <- vapply(list(CRootOf(x^2 - 2, 0), CRootOf(x^3 - 8, 0), CRootOf(x^3 - 8, 1),
                     CRootOf((x^3 - x - 1) * (x^2 - 2), 1), CRootOf((x^3 - x - 1) * (x^2 - 2), 3),
                     CRootOf(2 * x - 1, 0), CRootOf(x^5 - 9 * x^3 + 243, 4), CRootOf(-2 * x^3 + 2 * x + 2, 0),
                     CRootOf(x^2 - 2, 1, radicals = TRUE), CRootOf(x^3 - x - 1, -1),
                     CRootOf(x^3 - 2, 0, radicals = TRUE)), sstr, "")
  expect_identical(got, c("CRootOf(x^2 - 2, 0)", "2", "2*CRootOf(x^2 + x + 1, 0)", "CRootOf(x^3 - x - 1, 0)",
                          "CRootOf(x^3 - x - 1, 1)", "1/2", "3*CRootOf(x^5 - x^3 + 1, 4)", "CRootOf(x^3 - x - 1, 0)",
                          "sqrt(2)", "CRootOf(x^3 - x - 1, 2)", "2^(1/3)"))
  expect_error(CRootOf(x^3 - x - 1, 3), "root index out of \\[-3, 2\\] range, got 3")
})
