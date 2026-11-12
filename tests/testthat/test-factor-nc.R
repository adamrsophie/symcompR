# Reference values from SymPy's factor_nc() and factor().
test_that("factor_nc follows SymPy", {
  x <- Symbol("x")
  A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE); C <- Symbol("C", commutative = FALSE)
  cases <- list(
    list(x^2 + 2*A*x + A^2, "(x + A)^2"), list(expand((x + A)*(x + B)), "(x + A)*(x + B)"),
    list(A*B + A*C, "A*(B + C)"), list(B*A + C*A, "(B + C)*A"), list(2*A*B + 4*A*C, "2*A*(B + 2*C)"),
    list(x*A*B + x*A*C, "x*A*(B + C)"), list(A^2*B + A^3*C, "A^2*(A*C + B)"), list(A*B - B*A, "A*B - B*A"),
    list((A + B)^2, "(A + B)^2"), list(A*B*C + A*C*C, "A*(B + C)*C"), list(sin(x)*(A + B), "sin(x)*(A + B)"),
    list(A*x + x*B, "x*(A + B)"), list(expand((A + B)*(A - B)), "(A + B)*(A - B)"))
  for (cs in cases) expect_equal(sstr(factor_nc(cs[[1L]])), cs[[2L]])
  expect_equal(sstr(sym_factor(A*B + A*x)), "A*(x + B)")
})
