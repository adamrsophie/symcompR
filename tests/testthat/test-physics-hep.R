# Expected values are SymPy outputs (sympy.physics.hep.gamma_matrices).

test_that("Kahane's algorithm and traces of gamma matrices", {
  ii <- tensor_indices("i0:4", LorentzIndex)
  i0 <- ii[[1]]; i1 <- ii[[2]]; i2 <- ii[[3]]; i3 <- ii[[4]]
  G <- GammaMatrix
  expect_equal(sstr(kahane_simplify(G(i0) * G(i1) * G(-i0))), "-2*GammaMatrix(i1)")
  expect_equal(sstr(gamma_trace(G(i0) * G(i1))), "4*metric(i0, i1)")
  expect_equal(sstr(kahane_simplify(G(i0) * G(-i0))),
               "Matrix([[4, 0, 0, 0], [0, 4, 0, 0], [0, 0, 4, 0], [0, 0, 0, 4]])")
  expect_equal(sstr(gamma_trace(G(i0) * G(i1) * G(i2) * G(i3))),
               "-4*metric(i0, i2)*metric(i1, i3) + 4*metric(i0, i1)*metric(i2, i3) + 4*metric(i0, i3)*metric(i1, i2)")
})
