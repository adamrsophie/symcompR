# Fast checks of the solveset port (sympy/solvers/solveset.py); the expected
# strings are SymPy's outputs for the same calls.

x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")

test_that("solveset solves single equations as SymPy", {
  expect_identical(.pyrun(solveset(cos(x), x, Reals)),
                   "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers))")
  expect_identical(.pyrun(solveset(x^2 - 1, x)), "{-1, 1}")
  expect_identical(.pyrun(solveset(exp(x) - 1, x, Reals)), "{0}")
  expect_identical(.pyrun(solveset(exp(x) - 1, x)), "ImageSet(Lambda(_n, 2*_n*I*pi), Integers)")
  expect_identical(.pyrun(solveset(Abs(x) - 2, x, Reals)), "{-2, 2}")
  expect_identical(.pyrun(solveset(sqrt(x) - 2, x, Reals)), "{4}")
  expect_identical(.pyrun(solveset(log(x) - 1, x, Reals)), "{E}")
  expect_identical(.pyrun(solveset(x^2 > 4, x, Reals)), "Union(Interval.open(-oo, -2), Interval.open(2, oo))")
  expect_identical(.pyrun(solveset(Mod(x, 4) - 1, x, Integers)), "ImageSet(Lambda(_n, 4*_n + 1), Integers)")
  expect_identical(.pyrun(solveset_real(x^2 + 1, x)), "EmptySet")
  expect_identical(.pyrun(solveset_complex(x^2 + 1, x)), "{-I, I}")
})

test_that("the inverters and helpers match SymPy", {
  expect_identical(.pyrun(invert_real(exp(x), y, x)), "(x, Intersection({log(y)}, Reals))")
  expect_identical(.pyrun(invert_complex(exp(x), y, x)), "(x, ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers))")
  expect_identical(.pyrun(linear_eq_to_matrix(list(x + 2*y - 3, 3*x - y - 1), list(x, y))),
                   "(Matrix([\\n[1,  2],\\n[3, -1]]), Matrix([\\n[3],\\n[1]]))")
  expect_identical(.pyrun(linear_coeffs(2*x + 3*y - 1, x, y)), "[2, 3, -1]")
  expect_identical(.pyrun(domain_check(1/x, x, 0)), "False")
  expect_identical(.pyrun(solvify(sin(x), x, Interval(0, 2*S_Pi))), "[0, pi, 2*pi]")
})

test_that("linsolve, nonlinsolve and substitution return SymPy's sets", {
  expect_identical(.pyrun(linsolve(list(x + y + z - 1, x + y + 2*z - 3), list(x, y, z))), "{(-y - 1, y, 2)}")
  expect_identical(.pyrun(linsolve(list(x + y - 2, x - y), x, y)), "{(1, 1)}")
  expect_identical(.pyrun(nonlinsolve(list(x^2 + y^2 - 1, x - y), list(x, y))),
                   "{(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)}")
  expect_identical(.pyrun(nonlinsolve(list(x*y - 1, x - 2), list(x, y))), "{(2, 1/2)}")
  expect_identical(.pyrun(substitution(list(x + y), list(x))), "{(-y,)}")
  expect_identical(.pyrun(singular_points(ImplicitRegion(list(x, y, z), x*y*z + y^4 - x^2*z^2))),
                   "{(0, 0, z), (x, 0, 0)}")
})
