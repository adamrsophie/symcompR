# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_graycode.py);
# the rows come from the battery that compares this port with SymPy.

test_that("graycode: a", {
  expect_identical(.pyrun(local({ a <- GrayCode(6); a$current })), "000000")
  expect_identical(.pyrun(local({ a <- GrayCode(6); a$rank })), "0")
  expect_identical(.pyrun(local({ a <- GrayCode(3); .pyL(generate_gray(a)) })), "['000', '001', '011', '010', '110', '111', '101', '100']")
  expect_identical(.pyrun(local({ a <- GrayCode(6); a$rank; graycode_next(a, 2)$current })), "000011")
  expect_identical(.pyrun(local({ a <- GrayCode(6); graycode_next(a, -1)$current })), "100000")
  expect_identical(.pyrun(local({ a <- GrayCode(6); graycode_next(a, -1) })), "GrayCode(6)")
})

test_that("graycode: rank", {
  expect_identical(.pyrun(local({ GrayCode(5, start = "10010")$rank })), "28")
  expect_identical(.pyrun(local({ GrayCode(6, start = "101000")$rank })), "48")
  expect_identical(.pyrun(local({ GrayCode(6, rank = 4)$current })), "000110")
  expect_identical(.pyrun(local({ GrayCode(6, rank = 4)$rank })), "4")
  expect_identical(.pyrun(local({ lapply(generate_gray(GrayCode(4)), function(s) GrayCode(4, start = s)$rank) })), "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]")
  expect_identical(.pyrun(local({ GrayCode(15, rank = 15)$current })), "000000000001000")
  expect_identical(.pyrun(local({ GrayCode(3, rank = 11)$current })), "010")
  expect_identical(.pyrun(local({ GrayCode(3, rank = -1)$current })), "100")
  expect_identical(.pyrun(local({ GrayCode(4, start = "11")$current })), "0011")
  expect_identical(.pyrun(local({ GrayCode(4, start = "11")$rank })), "2")
})

test_that("graycode: str", {
  expect_identical(.pyrun(local({ GrayCode(4) })), "GrayCode(4)")
  expect_identical(.pyrun(local({ GrayCode(4)$n })), "4")
  expect_identical(.pyrun(local({ GrayCode(4)$selections })), "16")
  expect_identical(.pyrun(local({ GrayCode(4) == GrayCode(4, start = "0110") })), "True")
})

test_that("graycode: unrank", {
  expect_identical(.pyrun(local({ graycode_unrank(4, 10) })), "1111")
  expect_identical(.pyrun(local({ lapply(0:7, function(i) graycode_unrank(3, i)) })), "['000', '001', '011', '010', '110', '111', '101', '100']")
})

test_that("graycode: b2g", {
  expect_identical(.pyrun(local({ bin_to_gray("111") })), "100")
  expect_identical(.pyrun(local({ bin_to_gray("0100") })), "0110")
  expect_identical(.pyrun(local({ bin_to_gray("0101") })), "0111")
  expect_identical(.pyrun(local({ gray_to_bin(bin_to_gray("0100")) })), "0100")
  expect_identical(.pyrun(local({ gray_to_bin(bin_to_gray("0101")) })), "0101")
  expect_identical(.pyrun(local({ gray_to_bin("111") })), "101")
  expect_identical(.pyrun(local({ gray_to_bin("1") })), "1")
})

test_that("graycode: rbs", {
  expect_identical(.pyrun(local({ a <- random_bitstring(5); .pytuple(is.character(a), nchar(a), all(strsplit(a, "")[[1]] %in% c("0", "1"))) })), "(True, 5, True)")
})

test_that("graycode: sub", {
  expect_identical(.pyrun(local({ .pyL(get_subset_from_bitstring(c("a", "b", "c", "d"), "0011")) })), "['c', 'd']")
  expect_identical(.pyrun(local({ .pyL(get_subset_from_bitstring("abcd", "1001")) })), "['a', 'd']")
  expect_identical(.pyrun(local({ lapply(graycode_subsets(c("a", "b", "c")), .pyL) })), "[[], ['c'], ['b', 'c'], ['b'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'c'], ['a']]")
  expect_identical(.pyrun(local({ lapply(graycode_subsets(1:2), .pyL) })), "[[], [2], [1, 2], [1]]")
})

test_that("graycode: err", {
  expect_identical(.pyrun(local({ GrayCode(0) })), "ERR ValueError: Gray code dimension must be a positive integer, not 0")
  expect_identical(.pyrun(local({ GrayCode(2.2) })), "ERR ValueError: Gray code dimension must be a positive integer, not 2")
  expect_identical(.pyrun(local({ GrayCode(2, start = c(1, 1, 0)) })), "ERR ValueError: Gray code start has length 3 but should not be greater than 2")
  expect_identical(.pyrun(local({ GrayCode(2, rank = 2.5) })), "ERR ValueError: Gray code rank must be a positive integer, not 2")
  expect_identical(.pyrun(local({ get_subset_from_bitstring(c("c", "a", "c"), "1100") })), "ERR ValueError: The sizes of the lists are not equal")
  expect_identical(.pyrun(local({ generate_gray(GrayCode(3), start = "1111") })), "ERR ValueError: Gray code start has length 4 but should not be greater than 3")
})

