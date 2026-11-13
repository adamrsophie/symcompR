test_that("ray transfer matrices print as in SymPy", {
  f <- Symbol("f"); d <- Symbol("d"); R <- Symbol("R")
  expect_identical(sstr(ThinLens(f)), "Matrix([[1, 0], [-1/f, 1]])")
  expect_identical(sstr(FreeSpace(d)), "Matrix([[1, d], [0, 1]])")
  expect_identical(sstr(CurvedMirror(R)), "Matrix([[1, 0], [-2/R, 1]])")
})

test_that("Gaussian beam helpers match SymPy", {
  w <- Symbol("w"); wavelen <- Symbol("wavelen"); a <- Symbol("a"); b <- Symbol("b")
  expect_identical(sstr(waist2rayleigh(w, wavelen)), "pi*w^2/wavelen")
  expect_identical(sstr(geometric_conj_ab(a, b)), "a*b/(a + b)")
})

test_that("mirror formula and polarization match SymPy", {
  u <- Symbol("u"); v <- Symbol("v"); f <- Symbol("f")
  expect_identical(sstr(mirror_formula(focal_length = f, u = u)), "f*u/(-f + u)")
  expect_identical(sstr(mirror_formula(u = u, v = v)), "u*v/(u + v)")
  expect_identical(sstr(jones_vector(0, 0)), "Matrix([[1], [0]])")
  expect_identical(sstr(stokes_vector(0, 0)), "Matrix([[1], [1], [0], [0]])")
})
