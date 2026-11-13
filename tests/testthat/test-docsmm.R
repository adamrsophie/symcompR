test_that("ArrayComprehensionMap builds and expands", {
  i <- Symbol("i")
  expect_identical(sstr(ArrayComprehensionMap(function() 1, list(i, 1, 2))), "ArrayComprehensionMap((i, 1, 2))")
  expect_identical(sstr(doit(ArrayComprehensionMap(function(i) i + 1, list(i, 1, 3)))), "[2, 3, 4]")
})

test_that("ArrayDerivative rejects an array symbol as variable", {
  A <- ArraySymbol("A", list(2, 2))
  expect_error(ArrayDerivative(A, A), "Can't calculate derivative wrt A.", fixed = TRUE, class = "ValueError")
})

test_that("PIDController prints its gains", {
  s <- Symbol("s")
  expect_identical(capture.output(print(PIDController(1, 2, 3, 0, s))), "PIDController(1, 2, 3, 0, s)")
})

test_that("special sets as relations and Contains rebuilding", {
  x <- Symbol("x")
  expect_identical(sstr(as_relational(Naturals, x)), "(x >= 1) & (x < oo) & Eq(floor(x), x)")
  expect_identical(sstr(as_relational(Integers, x)), "(-oo < x) & (x < oo) & Eq(floor(x), x)")
  expect_identical(sstr(xreplace(Contains(x, Integers), list(list(x, 2L)))), "True")
})

test_that("subfactorial rewrites as uppergamma", {
  n <- Symbol("n")
  expect_identical(sstr(rewrite(subfactorial(n), "uppergamma")), "exp(-1)*uppergamma(n + 1, -1)")
})

test_that("factorization predicates and shapes", {
  M <- MatrixSymbol("M", 3, 3)
  expect_true(ask(Q$orthogonal(QofQR(M))))
  expect_true(ask(Q$upper_triangular(RofQR(M))))
  expect_true(ask(Q$diagonal(SofSVD(M))))
  expect_identical(format(shape(EigenVectors(M))), "(3, 3)")
})

test_that("dictionaries are read with $", {
  t <- Symbol("t"); f <- Function("f"); g <- Function("g")
  cs <- classify_sysode(list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), -f(t))))
  expect_true(cs$is_linear)
  expect_identical(sstr(cs$no_of_equation), "2")
})

test_that("intpoly_decompose keeps SymPy's term order", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_identical(.pySD(intpoly_decompose(x^2 + x * y + x + y + x^3 * y^2 + y^5)),
                   "{1: x + y, 2: x^2 + x*y, 5: x^3*y^2 + y^5}")
})

test_that("log_to_atan accepts the EX domain", {
  x <- Symbol("x")
  expect_identical(sstr(log_to_atan(Poly(x, x, domain = "EX"), Poly(1, x, domain = "EX"))), "2*atan(x)")
})

test_that("codegen Result datatypes", {
  x <- Symbol("x")
  expect_identical(ast_attr(Result(x^2), "get_datatype")("C"), "double")
  expect_error(ast_attr(Result(x^2), "datatype"), class = "AttributeError")
})
