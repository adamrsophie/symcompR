# Reference values from SymPy.
s <- function(e) if (is.null(e)) "NULL" else vapply(e, function(v) if (is_sym(v)) sstr(v) else
  paste0("(", paste(vapply(v, sstr, ""), collapse = ", "), ")"), "")

test_that("decompogen splits compositions", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(s(decompogen(sin(cos(x)), x)), c("sin(x)", "cos(x)"))
  expect_equal(s(decompogen(sin(x)^2 + sin(x) + 1, x)), c("x^2 + x + 1", "sin(x)"))
  expect_equal(s(decompogen(sin(sqrt(cos(x^2 + 1))), x)), c("sin(x)", "sqrt(x)", "cos(x)", "x^2 + 1"))
  expect_equal(s(decompogen(x^4 + 2*x^3 - x - 1, x)), c("x^2 - x - 1", "x^2 + x"))
  expect_equal(s(decompogen(Max(sin(x)^2, sin(x)), x)), c("Max(x, x^2)", "sin(x)"))
  expect_equal(s(decompogen(exp(x^2), x)), c("exp(x)", "x^2"))
  expect_equal(s(decompogen(1/x, x)), "1/x")
  expect_equal(s(decompogen(y, x)), "y")
})

test_that("homogeneous_order finds the degree of homogeneity", {
  x <- Symbol("x"); y <- Symbol("y"); f <- Function("f")
  expect_equal(sstr(homogeneous_order(x^2 + x*y, x, y)), "2")
  expect_null(homogeneous_order(x + y^2, x, y))
  expect_equal(sstr(homogeneous_order(sqrt(x^2 + y^2), x, y)), "1")
  expect_equal(sstr(homogeneous_order(sin(x/y), x, y)), "0")
  expect_equal(sstr(homogeneous_order(x^3*f(x)^2, x, f(x))), "5")
  expect_null(homogeneous_order(log(x/y) + x, x, y))
  expect_equal(sstr(homogeneous_order(2, x)), "0")
})

test_that("polynomial system helpers match SymPy", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b")
  r <- solve_linear_system_LU(Matrix(list(list(1, 2, 0, 1), list(3, 2, 2, 1), list(2, 0, 0, 1))), list(x, y, z))
  expect_equal(vapply(r, sstr, ""), c(x = "1/2", y = "1/4", z = "-1/2"))
  expect_equal(s(solve_triangulated(list(x^2 + y^2 + z^2 - 1, x + y + z - 1, x - y), list(x, y, z))),
               c("(0, 0, 1)", "(2/3, 2/3, -1/3)"))
  expect_length(solve_triangulated(list(x^2 - 2, y - x), list(x, y)), 0L)
  expect_equal(s(factor_system(list(x^2 - 1, y - 1))), c("(x + 1, y - 1)", "(x - 1, y - 1)"))
  expect_equal(s(factor_system(list((x - 1)*(y - 2), (x - 1)*(x + y)))), c("(x + y, y - 2)", "(x - 1)"))
  expect_equal(s(factor_system(list(x*y, x*z, y*z))), c("(x, y)", "(x, z)", "(y, z)"))
  expect_equal(s(factor_system(list(x*(x - 1)*(x + 1), x*y))), c("(x + 1, y)", "(x - 1, y)", "(x)"))
  expect_length(factor_system(list(x^2 - 4, a*y, b), list(x, y)), 0L)
  expect_length(factor_system(list(1)), 0L)
  expect_identical(factor_system(list(0)), list(list()))
})

test_that("rational and absolute value inequalities match SymPy", {
  x <- Symbol("x")
  q <- function(n, d, r) list(list(n, d), r)
  expect_equal(sstr(solve_rational_inequalities(list(list(q(-x + 1, 1, ">="), q(-x + 1, 1, "<="))))), "{1}")
  expect_equal(sstr(solve_rational_inequalities(list(list(q(x, 1, "!="), q(-x + 1, 1, ">="))))),
               "Union(Interval.open(-oo, 0), Interval.Lopen(0, 1))")
  expect_equal(sstr(solve_rational_inequalities(list(list(q(x - 1, x + 1, "<"))))), "Interval.open(-1, 1)")
  expect_equal(sstr(solve_rational_inequalities(list(list(q(x - 1, x, ">=")), list(q(x - 5, 1, ">"))))),
               "Union(Interval.open(-oo, 0), Interval(1, oo))")
  expect_equal(sstr(reduce_abs_inequality(Abs(x - 5) - 3, "<", x)), "(2 < x) & (x < 8)")
  expect_equal(sstr(reduce_abs_inequalities(list(list(Abs(3*x - 5) - 7, "<="), list(Abs(x + 25) - 13, ">")), x)),
               "(-2/3 <= x) & (x <= 4) & (((-oo < x) & (x < -38)) | ((-12 < x) & (x < oo)))")
  expect_equal(sstr(reduce_abs_inequalities(list(list(Abs(x - 4) + Abs(3*x - 5) - 7, "<")), x)),
               "(1/2 < x) & (x < 4)")
})

test_that("substitution returns solution tuples", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(substitution(list(x + y), list(x), list(list(y = 1)), list(y), list(), list(x, y))), "{(-1, 1)}")
  expect_equal(sstr(substitution(list(x + y), list(x))), "{(-y,)}")
  expect_equal(sstr(substitution(list(x^2 - 1, y - x), list(x, y))), "{(-1, -1), (1, 1)}")
  expect_equal(sstr(substitution(list(x*y - 1, x - 2), list(x, y))), "{(2, 1/2)}")
})

test_that("rsolve_poly, rsolve_ratio and rsolve_hyper match SymPy", {
  n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(rsolve_poly(list(-1, 1), 4*n^3, n)), "C0 + n^4 - 2*n^3 + n^2")
  expect_equal(sstr(rsolve_poly(list(1, -2, 1), n^2, n)), "C0 + C1*n + n^4/12 - n^3/3 + 5*n^2/12 - n/6")
  expect_equal(sstr(rsolve_poly(list(-n - 1, n), 1, n)), "C0*n - 1")
  expect_equal(sstr(rsolve_poly(list(-(n + 1), n + 2, -1), 0, n)), "C0")
  expect_identical(rsolve_poly(list(n^2 - 1, -2*n, 1), 0, n), S_Zero)
  expect_equal(sstr(rsolve_ratio(list(n, -n - 1), 0, n)), "C0/n")
  r <- rsolve_ratio(list(-2*n^3 + n^2 + 2*n - 1, 2*n^3 + n^2 - 6*n, -2*n^3 - 11*n^2 - 18*n - 9,
                         2*n^3 + 13*n^2 + 22*n + 8), 0, n)
  expect_identical(simplify(r - Symbol("C0")*(2*n - 3)/(2*(n^2 - 1))), S_Zero)
  expect_equal(sstr(rsolve_hyper(list(-1, -1, 1), 0, n)), "C0*(1/2 - sqrt(5)/2)^n + C1*(1/2 + sqrt(5)/2)^n")
  expect_equal(sstr(rsolve_hyper(list(-n - 1, 1), 0, n)), "C0*factorial(n)")
  expect_equal(sstr(rsolve_hyper(list(-1, 1), n, n)), "C0 + n*(n - 1)/2")
  expect_equal(sstr(rsolve_hyper(list(-2, 1), n, n)), "2^n*C0 - n - 1")
  expect_identical(rsolve_hyper(list(-(n + 1)^2, 1), 0, n), S_Zero)
})

test_that("rsolve handles polynomial coefficients", {
  n <- Symbol("n", integer = TRUE)
  y <- Function("y")
  expect_equal(sstr(rsolve((n - 1)*y(n + 2) - (n^2 + 3*n - 2)*y(n + 1) + 2*n*(n + 1)*y(n), y(n), n,
                           init = list(list(0, 0), list(1, 3)))), "3*2^n - 3*factorial(n)")
  expect_equal(sstr(rsolve(y(n + 1) - (n + 1)*y(n), y(n), n)), "C0*factorial(n)")
  expect_null(rsolve(n*y(n + 1) - (n + 1)*y(n) - n, y(n), n))
})

test_that("p10 review fixes for the new solvers", {
  x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a")
  n <- Symbol("n", integer = TRUE)
  M <- Integer("9007199254740992")
  expect_equal(sstr(solve_rational_inequalities(list(list(list(list((x - M)*(x - M - 1), 1), "<"))))),
               "Interval.open(9007199254740992, 9007199254740993)")
  expect_identical(substitution(list(x + y), list(x), list(list(y = 1)), list(y), list(x + 1), list(x, y)), S_EmptySet)
  expect_equal(sstr(substitution(list(exp(x) - 1), list(x))), "{(0,), (ImageSet(Lambda(_n, 2*_n*I*pi), Integers),)}")
  expect_equal(sstr(substitution(list(sin(x)), list(x))),
               "{(ImageSet(Lambda(_n, 2*_n*pi), Integers),), (ImageSet(Lambda(_n, 2*_n*pi + pi), Integers),)}")
  expect_equal(s(solve_triangulated(list(x - a), list(x))), "(a)")
  expect_identical(rsolve_poly(list(-(n - Integer("2147483648") + 1)^2, (n - Integer("2147483648"))^2), 0, n), S_Zero)
  Y <- Function("Y")
  expect_null(rsolve(Y(n + 1) - (n + 1)*Y(n), Y(n), n, init = list(list(0, 1), list(1, 2))))
  expect_equal(sstr(solveset(exp(x) - 2, x)), "ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers)")
  expect_equal(sstr(solveset(sin(x), x)), "Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))")
})
