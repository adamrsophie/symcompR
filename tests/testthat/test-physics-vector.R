test_that("frames orient and give direction cosine matrices", {
  q <- dynamicsymbols("q")
  N <- ReferenceFrame("N")
  A <- ReferenceFrame("A")
  A$orient_axis(N, N$z, q)
  expect_equal(sstr(A$dcm(N)), "Matrix([[cos(q(t)), sin(q(t)), 0], [-sin(q(t)), cos(q(t)), 0], [0, 0, 1]])")
  expect_equal(format(A$ang_vel_in(N)), "Derivative(q(t), t)*N.z")
  expect_equal(format(N$x + 2 * N$y), "N.x + 2*N.y")
})

test_that("vectors, dyadics and points follow SymPy", {
  N <- ReferenceFrame("N")
  a <- Symbol("a")
  v <- a * N$x + N$y
  expect_equal(sstr(dot(v, N$x)), "a")
  expect_equal(format(cross(N$x, N$y)), "N.z")
  expect_equal(format(outer_sym(N$x, N$y)), "(N.x|N.y)")
  O <- PhysPoint("O")
  P <- O$locatenew("P", a * N$x)
  expect_equal(format(P$pos_from(O)), "a*N.x")
  O$set_vel(N, 0)
  expect_equal(format(P$vel(N)), "0")
})

test_that("dynamicsymbols expands ranges and prints derivatives with primes", {
  qs <- dynamicsymbols("q1:3")
  expect_length(qs, 2L)
  expect_equal(sstr(qs[[2L]]), "q2(t)")
  expect_equal(length(dynamicsymbols("q0:2_J")), 2L)
  expect_equal(sstr(dynamicsymbols("q0:2_J")[[2L]]), "q1_J(t)")
  expect_equal(vsprint(dynamicsymbols("q", 1)), "q'")
})

test_that("sympy.vector functions keep working beside the physics ones", {
  C <- CoordSys3D("C")
  expect_equal(sstr(dot(C$i, C$i)), "1")
})
