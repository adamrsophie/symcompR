# Matrix expression additions: generic identity and zero, matrix calculus,
# Hadamard and Kronecker products, structured matrices and dotprodsimp.

test_that("empty MatMul and MatAdd give the generic identity and zero", {
  expect_equal(sstr(MatMul()), "I")
  expect_equal(sstr(MatAdd()), "0")
  expect_error(mat_shape(GenericIdentity()), class = "TypeError")
})

test_that("matrix expressions are differentiated through array expressions", {
  k <- Symbol("k"); i <- Symbol("i"); x <- Symbol("x")
  A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
  expect_equal(sstr(diff(x * B, x)), "B")
  expect_equal(sstr(diff(i^2 * A * B, i)), "(2*i)*A*B")
  expect_equal(sstr(diff(A * B, x)), "0")
  expect_equal(sstr(diff(Determinant(A), A)), "Determinant(A)*A.T^(-1)")
  expect_equal(sstr(diff(Trace(A * B), A)), "B.T")
  expect_equal(sstr(diff(MatrixElement(A, 0, 1), A)), "MatrixUnit(k, k, 0, 1)")
})

test_that("Hadamard and Kronecker products canonicalize like SymPy", {
  k <- Symbol("k")
  A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
  expect_equal(sstr(hadamard_product(A, HadamardProduct(A, B))), "HadamardProduct(HadamardPower(A, 2), B)")
  expect_equal(sstr(t(HadamardProduct(A, B))), "HadamardProduct(A.T, B.T)")
  expect_equal(sstr(kronecker_product(A, KroneckerProduct(A, B))), "KroneckerProduct(A, A, B)")
  expect_error(hadamard_product(), "Empty Hadamard product is undefined", class = "TypeError")
})

test_that("structured matrices have SymPy entries and printing", {
  k <- Symbol("k")
  expect_equal(sstr(t(MatrixUnit(k, k, 0, 1))), "MatrixUnit(k, k, 1, 0)")
  expect_equal(sstr(MatrixElement(DFT(4), 1, 1)), "-I/2")
  A <- MatrixSymbol("A", k, k)
  expect_equal(sstr(applyfunc(A, sin)), "Lambda(_d, sin(_d)).(A)")
})

test_that("dotprodsimp simplifies sums of products", {
  x <- Symbol("x")
  expect_equal(sstr(dotprodsimp(x / (x - 1) - 1 / (x - 1))), "1")
})

test_that("block_collapse, DotProduct and Derivative routing follow SymPy", {
  A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
  Z <- ZeroMatrix(2, 2)
  expect_equal(sstr(block_collapse(BlockMatrix(list(list(A, Z), list(Z, B)))^2)),
               "Matrix([\n[A^2,   0],\n[  0, B^2]])")
  expect_equal(sstr(inverse(BlockDiagMatrix(A, B))), "Matrix([\n[A^(-1),      0],\n[     0, B^(-1)]])")
  d <- DotProduct(Matrix(list(1, 2, 3)), Matrix(list(4, 5, 6)))
  expect_equal(sstr(doit(d)), "32")
  expect_error(DotProduct(A, Matrix(list(1, 2))), "Argument 1 of DotProduct is not a vector", class = "TypeError")
  X <- MatrixSymbol("X", 3, 3); C <- MatrixSymbol("C", 3, 3)
  expect_equal(sstr(Derivative(Trace(C * X), X, evaluate = TRUE)), "C.T")
  expect_equal(sstr(diff(Trace(C), as_explicit(C))), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
})

test_that("Eq, simplify and explicit factors of matrix expressions follow SymPy", {
  expect_equal(sstr(Eq(MatrixSymbol("A", 3, 1), MatrixSymbol("B", 2, 1))), "False")
  expect_equal(sstr(Eq(OneMatrix(1, 1), Identity(1))), "True")
  expect_error(doit(MatPow(ZeroMatrix(3, 3), -1)), "Matrix det == 0; not invertible.", fixed = TRUE,
               class = "NonInvertibleMatrixError")
  n <- Symbol("n", integer = TRUE); m <- Symbol("m", integer = TRUE)
  expect_error(as_explicit(MatrixSymbol("A", n, m)), "Matrix with symbolic shape cannot be represented explicitly.",
               fixed = TRUE, class = "ValueError")
  d <- Dummy("x"); X <- MatrixSymbol("X", 2, 2)
  expect_equal(sstr(simplify(ElementwiseApplyFunction(Lambda(d, d^2), as_explicit(X)))),
               "Lambda(_x, _x^2).(Matrix([\n[X[0, 0], X[0, 1]],\n[X[1, 0], X[1, 1]]]))")
  x <- MatrixSymbol("x", 2, 1)
  m1 <- Matrix(list(list(1, 2), list(4, 7))); m3 <- Matrix(list(list(3), list(5)))
  fx <- MatAdd(MatMul(6, Identity(1)), MatMul(Transpose(x), m3), MatMul(Transpose(x), m1, x))
  expect_equal(sstr(diff(fx, x)),
               "Matrix([\n[3],\n[5]]) + Matrix([\n[1, 2],\n[4, 7]])*x + Matrix([\n[1, 4],\n[2, 7]])*x")
})
