skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_subsets.py);
# the rows come from the battery that compares this port with SymPy.

test_that("subsets: a matches SymPy", {
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_binary(a)$subset) })), "['b']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_prev_binary(a)$subset) })), "['c']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_lexicographic(a)$subset) })), "['d']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_prev_lexicographic(a)$subset) })), "['c']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_gray(a)$subset) })), "['c']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_prev_gray(a)$subset) })), "['d']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); subset_next_binary(a) == Subset("b", c("a", "b", "c", "d")) })), "True")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); subset_prev_gray(a) == Subset("d", c("a", "b", "c", "d")) })), "True")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_binary })), "3")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_lexicographic })), "14")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$rank_gray })), "2")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$cardinality })), "16")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$size })), "2")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); a$superset_size })), "4")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(a$superset) })), "['a', 'b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); bitlist_from_subset(a, c("a", "b", "c", "d")) })), "0011")
  expect_identical(.pyrun(local({ bitlist_from_subset(c("a", "c"), c("a", "b", "c", "d")) })), "1010")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_iterate_binary(a, -2)$subset) })), "['d']")
  expect_identical(.pyrun(local({ a <- Subset(c("a", "b", "c"), c("a", "b", "c", "d")); .pyL(subset_iterate_binary(a, 2)$subset) })), "[]")
  expect_identical(.pyrun(local({ a <- Subset(c("a", "b", "c", "d"), c("a", "b", "c", "d")); .pyL(subset_next_binary(a)$subset) })), "[]")
  expect_identical(.pyrun(local({ a <- Subset(character(), c("a", "b", "c", "d")); .pyL(subset_prev_binary(a)$subset) })), "['a', 'b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_iterate_graycode(a, 3)$subset) })), "['b', 'c', 'd']")
  expect_identical(.pyrun(local({ a <- Subset(c("c", "d"), c("a", "b", "c", "d")); .pyL(subset_iterate_graycode(a, -2)$subset) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_next_lexicographic(Subset(character(), c("a", "b", "c", "d")))$subset) })), "['a']")
  expect_identical(.pyrun(local({ .pyL(subset_next_lexicographic(Subset("d", c("a", "b", "c", "d")))$subset) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_prev_lexicographic(Subset(character(), c("a", "b", "c", "d")))$subset) })), "['d']")
  expect_identical(.pyrun(local({ .pyL(subset_prev_lexicographic(Subset("a", c("a", "b", "c", "d")))$subset) })), "[]")
  expect_identical(.pyrun(local({ Subset(c("a", "b", "c", "d"), c("a", "b", "c", "d"))$rank_lexicographic })), "4")
  expect_identical(.pyrun(local({ Subset(character(), c("a", "b", "c", "d"))$rank_lexicographic })), "0")
})

test_that("subsets: b matches SymPy", {
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_next_binary(a)$subset) })), "[2, 5, 6]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_prev_binary(a)$subset) })), "[2, 5]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_next_lexicographic(a)$subset) })), "[2, 6]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_prev_lexicographic(a)$subset) })), "[2, 5, 6, 7]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_next_gray(a)$subset) })), "[2, 5, 6, 7]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); .pyL(subset_prev_gray(a)$subset) })), "[2, 5]")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$rank_binary })), "37")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$rank_lexicographic })), "93")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$rank_gray })), "57")
  expect_identical(.pyrun(local({ a <- Subset(c(2L, 5L, 7L), 1:7); a$cardinality })), "128")
  expect_identical(.pyrun(local({ a <- Subset(c(7L, 2L), 1:7); .pyL(subset_next_lexicographic(a)$subset) })), "[3]")
  expect_identical(.pyrun(local({ a <- Subset(c(7L, 2L), 1:7); .pyL(subset_prev_lexicographic(a)$subset) })), "[2, 6, 7]")
  expect_identical(.pyrun(local({ a <- Subset(c(7L, 2L), 1:7); a$rank_lexicographic })), "96")
})

test_that("subsets: unrank matches SymPy", {
  expect_identical(.pyrun(local({ subset_unrank_binary(4, c("a", "b", "c", "d"))$rank_binary })), "4")
  expect_identical(.pyrun(local({ subset_unrank_gray(10, c("a", "b", "c", "d"))$rank_gray })), "10")
  expect_identical(.pyrun(local({ subset_unrank_binary(33, 1:9)$rank_binary })), "33")
  expect_identical(.pyrun(local({ subset_unrank_gray(25, 1:9)$rank_gray })), "25")
  expect_identical(.pyrun(local({ .pyL(subset_unrank_binary(4, c("a", "b", "c", "d"))$subset) })), "['b']")
  expect_identical(.pyrun(local({ .pyL(subset_unrank_gray(10, c("a", "b", "c", "d"))$subset) })), "['a', 'b', 'c', 'd']")
  expect_identical(.pyrun(local({ .pyL(subset_unrank_gray(25, 1:9)$subset) })), "[5, 7, 9]")
})

test_that("subsets: lexloop matches SymPy", {
  expect_identical(.pyrun(local({ a <- Subset(character(), c("a", "b", "c", "d")); out <- list(); for (j in 1:16) { out[[j]] <- .pyL(a$subset); a <- subset_next_lexicographic(a) }; out })), "[[], ['a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd'], ['a', 'b', 'd'], ['a', 'c'], ['a', 'c', 'd'], ['a', 'd'], ['b'], ['b', 'c'], ['b', 'c', 'd'], ['b', 'd'], ['c'], ['c', 'd'], ['d']]")
  expect_identical(.pyrun(local({ a <- Subset("d", c("a", "b", "c", "d")); out <- list(); for (j in 1:16) { out[[j]] <- .pyL(a$subset); a <- subset_prev_lexicographic(a) }; out })), "[['d'], ['c', 'd'], ['c'], ['b', 'd'], ['b', 'c', 'd'], ['b', 'c'], ['b'], ['a', 'd'], ['a', 'c', 'd'], ['a', 'c'], ['a', 'b', 'd'], ['a', 'b', 'c', 'd'], ['a', 'b', 'c'], ['a', 'b'], ['a'], []]")
})

test_that("subsets: grayloop matches SymPy", {
  expect_identical(.pyrun(local({ a <- Subset(integer(), 1:4); out <- list(); for (j in 1:16) { out[[j]] <- .pyL(a$subset); a <- subset_next_gray(a) }; out })), "[[], [4], [3, 4], [3], [2, 3], [2, 3, 4], [2, 4], [2], [1, 2], [1, 2, 4], [1, 2, 3, 4], [1, 2, 3], [1, 3], [1, 3, 4], [1, 4], [1]]")
})

test_that("subsets: binloop matches SymPy", {
  expect_identical(.pyrun(local({ a <- Subset(integer(), 1:3); out <- list(); for (j in 1:8) { out[[j]] <- .pyL(a$subset); a <- subset_next_binary(a) }; out })), "[[], [3], [2], [2, 3], [1], [1, 3], [1, 2], [1, 2, 3]]")
})

test_that("subsets: err matches SymPy", {
  expect_identical(.pyrun(local({ Subset(c("a", "b"), "a") })), "ERR ValueError: Invalid arguments have been provided. The superset must be larger than the subset.")
  expect_identical(.pyrun(local({ Subset("a", c("b", "c")) })), "ERR ValueError: The superset provided is invalid as it does not contain the element a")
  expect_identical(.pyrun(local({ subset_from_bitlist(c("a", "b"), "010") })), "ERR ValueError: The sizes of the lists are not equal")
  expect_identical(.pyrun(local({ Subset(3L, 1:2) })), "ERR ValueError: The superset provided is invalid as it does not contain the element 3")
})

test_that("subsets: eq matches SymPy", {
  expect_identical(.pyrun(local({ Subset("a", c("a", "b")) != Subset("b", c("a", "b")) })), "True")
  expect_identical(.pyrun(local({ Subset("a", c("a", "b")) != Subset("a", c("a", "c")) })), "True")
  expect_identical(.pyrun(local({ Subset("a", c("a", "b")) == Subset("a", c("a", "b")) })), "True")
})

test_that("subsets: fbl matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(subset_from_bitlist(c("a", "b", "c", "d"), "0011")$subset) })), "['c', 'd']")
  expect_identical(.pyrun(local({ .pyL(subset_from_bitlist(c("a", "b", "c", "d"), "1010")$subset) })), "['a', 'c']")
})

test_that("subsets: idx matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(subset_indices(c(1L, 3L), c(5L, 3L, 2L, 1L))) })), "[3, 1]")
  expect_identical(.pyrun(local({ .pyL(subset_indices(7L, c(5L, 3L, 2L, 1L))) })), "[]")
  expect_identical(.pyrun(local({ .pyL(subset_indices(c("c", "a"), c("a", "b", "c"))) })), "[2, 0]")
})

test_that("subsets: ksubsets matches SymPy", {
  expect_identical(.pyrun(local({ lapply(ksubsets(1:3, 2), .pyT) })), "[(1, 2), (1, 3), (2, 3)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(1:5, 2), .pyT) })), "[(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(1:4, 3), .pyT) })), "[(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(c("a", "b", "c"), 1), .pyT) })), "[('a',), ('b',), ('c',)]")
  expect_identical(.pyrun(local({ lapply(ksubsets(1:2, 3), .pyT) })), "[]")
})

