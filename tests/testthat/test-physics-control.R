test_that("transfer functions print and combine as in SymPy", {
  s <- Symbol("s"); p <- Symbol("p")
  tf1 <- TransferFunction(s + 6, s - 5, s)
  tf2 <- TransferFunction(s + 3, s + 1, s)
  expect_identical(sstr(tf1), "TransferFunction(s + 6, s - 5, s)")
  expect_identical(sstr(tf1 + tf2),
                   "Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s))")
  poles <- tf1$poles()
  expect_length(poles, 1L)
  expect_identical(sstr(poles[[1L]]), "5")
  expect_error(1 - p^3 + tf1, "cannot be added with")
})

test_that("Routh Hurwitz tables match SymPy", {
  s <- Symbol("s"); b1 <- Symbol("b1")
  r <- RouthHurwitz(Poly(s^3 + b1 * s^2 + 2 * s + 1, s), s)
  expect_identical(sstr(r), "Matrix([[1, 2], [b1, 1], [(2*b1 - 1)/b1, 0], [1, 0]])")
})

test_that("state space systems connect in series as in SymPy", {
  s <- Symbol("s")
  ss1 <- StateSpace(Matrix(list(list(0L, 1L), list(1L, 0L))), Matrix(list(list(0L), list(1L))),
                    Matrix(list(list(0L, 1L))), Matrix(list(list(0L))))
  ser <- Series(TransferFunction(s, s + 1, s), ss1)
  expect_identical(ss1$num_inputs, 1L)
  expect_identical(sstr(ser$rewrite(TransferFunction)), "TransferFunction(s^2, s^3 + s^2 - s - 1, s)")
})
