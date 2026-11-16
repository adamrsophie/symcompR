# Regression cases for matrices and matrix expressions.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("matrix sums print and sort as SymPy's MatAdd", {
  M <- MatrixSymbol("M", 2L, 2L)
  expect_identical(printer_doprint(NumPyPrinter(), M + Identity(2L)), "numpy.eye(2) + M")
  A <- MatrixSymbol("A", 1L, 3L)
  B <- MatrixSymbol("B", 1L, 3L)
  e <- A - B
  expect_identical(vapply(e$args, sstr, ""), c("-B", "A"))
  expect_identical(sstr(e), "A - B")
  expect_identical(tensorflow_code(new_sym("MatrixElement", list(e, S(0L), S(0L)))),
                   "(tensorflow.math.add((-1)*B, A))[0, 0]")
})

test_that("an IndexedBase without a shape has shape None in array expressions", {
  I <- IndexedBase("I")
  J <- IndexedBase("J")
  i <- Idx("i", list(0L, 2L))
  j <- Idx("j", list(0L, 2L))
  k <- Idx("k", list(0L, 2L))
  r <- convert_indexed_to_array(I(i, j) * J(j, k))
  expect_identical(sstr(r), "ArrayDiagonal(ArrayTensorProduct(I, J), (1, 2))")
  expect_null(get_shape(r))
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i, k), I(i, j) * J(j, k))),
                   "I = numpy.einsum(\",->\", I, J)")
  expect_identical(printer_doprint(NumPyPrinter(), Assignment(I(i), J(i, i))), "I = numpy.einsum(\"->\", J)")
  expect_identical(sstr(ArrayContraction(ArrayTensorProduct(I, J), list(0L, 1L))),
                   "ArrayContraction(ArrayTensorProduct(I, J), (0, 1))")
  A <- IndexedBase(Symbol("A"), Tuple(4L, 4L))
  expect_null(get_shape(ArrayTensorProduct(A, I)))
  expect_length(get_shape(ArrayTensorProduct(A, A)), 4L)
})

test_that("equal explicit matrices are equal and belong to a condition set", {
  w <- ImmutableMatrix(list(list(Rational(11, 39), Rational(16, 39), Rational(4, 13))))
  expect_identical(sstr(Eq(w, w)), "True")
  expect_identical(sstr(w), "Matrix([[11/39, 16/39, 4/13]])")
  T4 <- Matrix(list(list(Rational(1, 5), Rational(2, 5), Rational(2, 5)),
                    list(Rational(1, 10), S_Half, Rational(2, 5)),
                    list(Rational(3, 5), Rational(3, 10), Rational(1, 10))))
  wm <- MatrixSymbol("wm", 1L, 3L)
  cs <- ConditionSet(wm, Eq(wm * T4, wm))
  expect_identical(sstr(doit(contains(cs, w))), "True")
})

test_that("cholesky, LDLdecomposition and dot check their input as SymPy", {
  skip_on_cran()
  m <- function(r) paste(capture.output(print(r)), collapse = " ")
  expect_error(cholesky(Matrix(list(list(1, 2), list(3, 4)))), "Matrix must be Hermitian.", fixed = TRUE)
  expect_error(cholesky(Matrix(list(list(1, 2), list(2, 1)))), "Matrix must be positive-definite", fixed = TRUE)
  expect_equal(sstr(cholesky(Matrix(list(list(1, 2), list(2, 1))), hermitian = FALSE)), "Matrix([[1, 0], [2, sqrt(3)*I]])")
  A <- Matrix(list(list(4, 2 * S_I), list(-2 * S_I, 3)))
  expect_equal(sstr(cholesky(A)), "Matrix([[2, 0], [-I, sqrt(2)]])")
  ld <- LDLdecomposition(A)
  expect_equal(sstr(ld$L), "Matrix([[1, 0], [-I/2, 1]])")
  expect_equal(sstr(ld$D), "Matrix([[4, 0], [0, 2]])")
  expect_error(dot(Matrix(list(1, 2, 3)), Matrix(list(1, 2))), "Dimensions incorrect for dot product: (3, 1), (2, 1)", fixed = TRUE)
  expect_equal(sstr(dot(Matrix(list(1, 2, 3)), Matrix(list(list(1, 2, 3))))), "14")
})

test_that("mutable, immutable, dense and sparse matrices follow SymPy's classes", {
  skip_on_cran()
  x <- Symbol("x")
  I <- ImmutableMatrix(list(list(1, 2), list(3, 4))); M <- Matrix(list(list(1, 0), list(0, 1)))
  S <- SparseMatrix(2, 2, list(list(0, 1, x))); IS <- ImmutableSparseMatrix(S)
  k <- function(A) symcompR:::.mclass(A)
  expect_equal(c(k(I + M), k(M %*% I), k(t(I)), k(inverse(I)), k(x * I), k(I %^% 2)), rep("ImmutableDenseMatrix", 6))
  expect_equal(c(k(S + M), k(M + S), k(as_mutable(I)), k(Matrix(I)), k(LUsolve(I, M)), k(row_join(M, I))), rep("MutableDenseMatrix", 6))
  expect_equal(c(k(S %*% S), k(t(S)), k(x * S), k(row_join(S, M))), rep("MutableSparseMatrix", 4))
  expect_equal(c(k(IS + S), k(S + IS), k(subs(IS, x, 1)), k(diff(IS, x)), k(as_immutable(S))), rep("ImmutableSparseMatrix", 5))
  expect_equal(k(I + S), "ImmutableDenseMatrix")
  expect_equal(sstr(applyfunc(S, function(z) z + 1)), "Matrix([[0, x + 1], [0, 0]])")
  expect_equal(sstr(applyfunc(I, function(z) z + 1)), "Matrix([[2, 3], [4, 5]])")
  expect_equal(srepr(I), "ImmutableDenseMatrix([[Integer(1), Integer(2)], [Integer(3), Integer(4)]])")
  expect_equal(srepr(S), "MutableSparseMatrix([[Integer(0), Symbol('x')], [Integer(0), Integer(0)]])")
  expect_error(I[1, 1] <- 5, "Cannot set values of <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>")
  expect_error(IS[1, 1] <- 5, "ImmutableSparseMatrix")
  S2 <- S; S2[2, 2] <- 7
  expect_equal(c(sstr(S2), k(S2), nnz(S2)), c("Matrix([[0, x], [0, 7]])", "MutableSparseMatrix", "2"))
  expect_equal(names(todok(S)), "(0, 1)")
  expect_equal(vapply(row_list(S2), function(t) paste(t[[1L]], t[[2L]], sstr(t[[3L]])), ""), c("0 1 x", "1 1 7"))
  expect_equal(sstr(SparseMatrix(2, 3, function(i, j) i + j)), "Matrix([[0, 1, 2], [1, 2, 3]])")
  expect_equal(sstr(SparseMatrix(2, 2, list(1, 0, 0, 1))), "Matrix([[1, 0], [0, 1]])")
  expect_true(I == Matrix(list(list(1, 2), list(3, 4))))
})

test_that("blockcut, slices, ask, python() and sqf_norm follow SymPy", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y")
  n <- Symbol("n", integer = TRUE, positive = TRUE)
  B <- blockcut(MatrixSymbol("A", n + 1, 2), list(n, 1), c(1, 1))
  expect_equal(vapply(mat_shape(B), sstr, ""), c("n + 1", "2"))
  N <- Integer("9007199254740993")
  expect_equal(sstr(mat_shape(blockcut(MatrixSymbol("A", N, 1), list(N), list(1)))[[1L]]), "9007199254740993")
  m <- as_explicit(MatrixSlice(ZeroMatrix(2, n), c(0, 1), c(0, 1)))
  expect_equal(c(m$nrow, m$ncol, sstr(m$data[[1L]])), c("1", "1", "0"))
  expect_false(ask(Q$commutative(x), !Q$commutative(x)))
  expect_false(ask(Q$negative(x + y), Q$positive(x + y)))
  expect_error(ask(Q$positive(x + y), Q$positive(x + y), !Q$positive(x + y)), "inconsistent")
  expect_equal(python(Symbol("x_{1}") + 1), "x_1 = Symbol('x_{1}')\ne = x_1 + 1")
  expect_equal(python(Symbol("a^{2}") * Symbol("lambda")), "a^2 = Symbol('a^{2}')\nlambda_ = Symbol('lambda')\ne = a**2*lambda_")
  a <- sqrt(Integer(2)) + sqrt(Integer(3))
  s <- sqf_norm(x + 1/a, a)
  expect_equal(c(sstr(s[[1L]][[1L]]), sstr(s[[2L]]), sstr(s[[3L]])), c("0", "x - sqrt(2) + sqrt(3)", "x^4 - 10*x^2 + 1"))
  s <- sqf_norm(x*y + 1, sqrt(Integer(2)))
  expect_equal(c(vapply(s[[1L]], sstr, ""), sstr(s[[2L]]), sstr(s[[3L]])),
               c("1", "0", "x*y - sqrt(2)*y + 1", "x^2*y^2 + 2*x*y - 2*y^2 + 1"))
})

test_that("MatAdd and MatMul keep their arguments until doit", {
  A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
  got <- vapply(list(MatAdd(A, A), MatMul(2, A, 3), MatMul(A, A), MatMul(2, 3),
                     doit(MatMul(2, A, 3)), doit(MatAdd(A, A)), A + A, -A, A - B, A * 2, A / 2), sstr, "")
  expect_equal(got, c("A + A", "2*A*3", "A*A", "6", "6*A", "2*A", "2*A", "-A", "A - B", "2*A", "(1/2)*A"))
})

test_that("matrix scalars sort before numbers in a sum", {
  A <- MatrixSymbol("A", 2, 2); x <- Symbol("x")
  got <- vapply(list(Trace(A) + 5, MatrixElement(A, 0, 0) + 1, x + MatrixElement(A, 0, 0) + 1,
                     Trace(A) - 1, -Trace(A) + x), sstr, "")
  expect_equal(got, c("Trace(A) + 5", "A[0, 0] + 1", "x + A[0, 0] + 1", "Trace(A) - 1", "x - Trace(A)"))
})

test_that("exact number and triangular determinants take the direct routes", {
  m <- matrix(c(2, -1, 0, 3, 4, 1, -2, 5, 1, 0, 3, 2, 7, -3, 1, 1), 4)
  M <- Matrix(lapply(1:4, function(i) as.list(m[i, ])))
  expect_equal(sstr(det(M)), as.character(round(base::det(m))))
  Q <- Matrix(list(list(Rational(1, 2), 1, 2), list(Rational(1, 3), 4, 1), list(5, Rational(-1, 4), 3)))
  expect_equal(sstr(det(Q)), sstr(cancel(symcompR:::.det_berkowitz(Q))))
  expect_equal(sstr(det(Matrix(list(list(1, 2, 3), list(2, 4, 6), list(1, 1, 1))))), "0")
  x <- Symbol("x"); y <- Symbol("y")
  T <- Matrix(list(list(x + 1, y, 3), list(0, x * y, x), list(0, 0, y - 2)))
  expect_equal(sstr(det(T)), sstr(expand((x + 1) * x * y * (y - 2))))
})

test_that("rational functions with a shared denominator add over it", {
  x <- Symbol("x"); y <- Symbol("y")
  d <- (x + y + 1)^2 + 1
  M <- Matrix(list(lapply(1:3, function(i) (x^i + y) / d)))
  N <- Matrix(lapply(1:3, function(i) list((y^i - x) / d)))
  P <- DomainMatrix_from_Matrix(M)$matmul(DomainMatrix_from_Matrix(N))
  expect_equal(sstr(P$to_Matrix()), paste0("Matrix([[", sstr(cancel((M %*% N)[1, 1])), "]])"))
})
