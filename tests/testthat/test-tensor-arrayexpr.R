# Array expressions: derivatives, conversions to and from indexed form,
# collecting tensor products and array comprehensions.

test_that("array_derive gives SymPy's array expressions", {
  k <- Symbol("k")
  A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
  expect_equal(sstr(array_derive(A * B, A)), "PermuteDims(ArrayContraction(ArrayTensorProduct(I, I, B), (3, 4)), (3)(1 2))")
  expect_equal(sstr(diff(A * B, A)), "PermuteDims(ArrayTensorProduct(I, B), (3)(1 2))")
})

test_that("indexed and array expressions convert both ways", {
  k <- Symbol("k"); i <- Symbol("i"); j <- Symbol("j"); d <- Symbol("d")
  M <- MatrixSymbol("M", k, k); N <- MatrixSymbol("N", k, k)
  e <- Sum(MatrixElement(M, i, d) * MatrixElement(N, d, j), list(d, 0, k - 1))
  expect_equal(sstr(convert_indexed_to_array(e, list(i, j))), "ArrayContraction(ArrayTensorProduct(M, N), (1, 2))")
  A <- ArraySymbol("A", list(3, 3))
  expect_equal(sstr(convert_array_to_indexed(A, list(i, j))), "A[i, j]")
})

test_that("collect_tensor_products merges sums of tensor products", {
  k <- Symbol("k")
  A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k); C <- MatrixSymbol("C", k, k)
  expect_equal(sstr(collect_tensor_products(ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C)))),
               "ArrayTensorProduct(A, B + C)")
})

test_that("array comprehensions expand with numeric bounds", {
  i <- Symbol("i"); j <- Symbol("j")
  a <- ArrayComprehension(i * j, list(i, 1, 2), list(j, 1, 2))
  expect_equal(sstr(doit(a)), "[[1, 2], [2, 4]]")
  expect_error(ArrayComprehension(i, list(i, 3, 1)), "Lower bound should be inferior to upper bound", class = "ValueError")
})

test_that("array expression doit, reshape and indexing follow SymPy", {
  i <- Symbol("i"); j <- Symbol("j"); k <- Symbol("k")
  M <- ArraySymbol("M", list(k, k))
  expect_equal(sstr(doit(ArrayAdd(ArrayTensorProduct(i, j), ArrayContraction(M, Tuple(Integer(0), Integer(1)))))),
               "i*j + ArrayContraction(M, (0, 1))")
  expect_equal(sstr(doit(Reshape(Array(list(Integer(1), Integer(2), Integer(3), Integer(4)), 4), Tuple(Integer(2), Integer(2))))), "[[1, 2], [3, 4]]")
  A <- ArraySymbol("A", list(3, 3, 3))
  expect_equal(sstr(convert_array_to_indexed(ArrayContraction(A, Tuple(Integer(0), Integer(2))), list(i))), "Sum(A[_d0, i, _d0], (_d0, 0, 2))")
  K <- ArraySymbol("K", list(k, k))
  expect_equal(sstr(ArrayElement(K, list(k, k))), "K[k, k]")
  expect_error(ArrayElement(ArraySymbol("B", list(2, 2)), list(2, 0)), "shape is out of bounds", class = "ValueError")
})
