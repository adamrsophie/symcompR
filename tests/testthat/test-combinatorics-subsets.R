# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_subsets.py);
# the rows come from the battery that compares this port with SymPy.

test_that("subsets: a", {
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_lexicographic(a)$subset) })), "['d']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_prev_lexicographic(a)$subset) })), "['c']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_gray(a)$subset) })), "['c']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_binary })), "3")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_lexicographic })), "14")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_gray })), "2")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$cardinality })), "16")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$size })), "2")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$superset_size })), "4")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(a$superset) })), "['a', 'b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); bitlist_from_subset(a, c("a", "b", "c", "d")) })), "0011")
  expect_identical(.pyrun(local({ bitlist_from_subset(c("a", "c"), c("a", "b", "c", "d")) })), "1010")
  expect_identical(.pyrun(local({ a <- Subset(c("a", "b", "c"), c("a", "b", "c", "d")); .pyL(subset_iterate_binary(a, 2)$subset) })), "[]")
  expect_identical(.pyrun(local({ a <- Subset(c("a", "b", "c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_binary(a)$subset) })), "[]")
  expect_identical(.pyrun(local({ a <- Subset(character(), c("a", "b", "c", "d")); .pyL(subset_prev_binary(a)$subset) })), "['a', 'b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_iterate_graycode(a, 3)$subset) })), "['b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_iterate_graycode(a, -2)$subset) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_next_lexicographic(Subset("d", c("a", "b", "c", "d")))$subset) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_prev_lexicographic(Subset("a", c("a", "b", "c", "d")))$subset) })), "[]")
  expect_identical(.pyrun(local({ Subset(character(), c("a", "b", "c", "d"))$rank_lexicographic })), "0")
})

test_that("subsets: b", {
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$rank_binary })), "37")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$rank_lexicographic })), "93")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$cardinality })), "128")
  expect_identical(.pyrun(local({ a <- Subset(c(7L, 2L), 1:7); .pyL(subset_next_lexicographic(a)$subset) })), "[3]")
})

test_that("subsets: unrank", {
  expect_identical(.pyrun(local({ subset_unrank_binary(4, c("a", "b", "c", "d"))$rank_binary })), "4")
  expect_identical(.pyrun(local({ subset_unrank_gray(10, c("a", "b", "c", "d"))$rank_gray })), "10")
  expect_identical(.pyrun(local({ subset_unrank_binary(33, 1:9)$rank_binary })), "33")
  expect_identical(.pyrun(local({ .pyL(subset_unrank_binary(4, c("a", "b", "c", "d"))$subset) })), "['b']")
  expect_identical(.pyrun(local({ .pyL(subset_unrank_gray(10, c("a", "b", "c", "d"))$subset) })), "['a', 'b', 'c', 'd']")
})

test_that("subsets: err", {
  expect_identical(.pyrun(local({ Subset(c("a", "b"), "a") })), "ERR ValueError: Invalid arguments have been provided. The superset must be larger than the subset.")
  expect_identical(.pyrun(local({ subset_from_bitlist(c("a", "b"), "010") })), "ERR ValueError: The sizes of the lists are not equal")
})

test_that("subsets: eq", {
  expect_identical(.pyrun(local({ Subset("a", c("a", "b")) != Subset("b", c("a", "b")) })), "True")
  expect_identical(.pyrun(local({ Subset("a", c("a", "b")) != Subset("a", c("a", "c")) })), "True")
})

test_that("subsets: fbl", {
  expect_identical(.pyrun(local({ .pyL(subset_from_bitlist(c("a", "b", "c", "d"), "0011")$subset) })), "['c', 'd']")
})

test_that("subsets: idx", {
  expect_identical(.pyrun(local({ .pyL(subset_indices(c(1L, 3L), c(5L, 3L, 2L, 1L))) })), "[3, 1]")
  expect_identical(.pyrun(local({ .pyL(subset_indices(7L, c(5L, 3L, 2L, 1L))) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_indices(c("c", "a"), c("a", "b", "c"))) })), "[2, 0]")
})

test_that("subsets: ksubsets", {
  expect_identical(.pyrun(local({ lapply(ksubsets(1:5, 2), .pyT) })), "[(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(c("a", "b", "c"), 1), .pyT) })), "[('a',), ('b',), ('c',)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(1:2, 3), .pyT) })), "[]")
})

