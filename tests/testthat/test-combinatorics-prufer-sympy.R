skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_prufer.py);
# the rows come from the battery that compares this port with SymPy.

test_that("prufer: nodes matches SymPy", {
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)), 5)$nodes })), "5")
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)))$nodes })), "5")
})

test_that("prufer: a matches SymPy", {
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)))$rank })), "0")
  expect_identical(.pyrun(local({ .pyL(Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)))$prufer_repr) })), "[0, 0, 0]")
  expect_identical(.pyrun(local({ Prufer(list(c(2, 4), c(1, 4), c(1, 3), c(0, 5), c(0, 4)))$rank })), "924")
  expect_identical(.pyrun(local({ Prufer(list(c(2, 4), c(1, 4), c(1, 3), c(0, 5), c(0, 4)))$nodes })), "6")
  expect_identical(.pyrun(local({ lapply(Prufer(list(c(2, 4), c(1, 4), c(1, 3), c(0, 5), c(0, 4)))$tree_repr, .pyL) })), "[[2, 4], [1, 4], [1, 3], [0, 5], [0, 4]]")
  expect_identical(.pyrun(local({ .pyL(Prufer(list(c(2, 4), c(1, 4), c(1, 3), c(0, 5), c(0, 4)))$prufer_repr) })), "[4, 1, 4, 0]")
})

test_that("prufer: edges matches SymPy", {
  expect_identical(.pyrun(local({ e <- prufer_edges(c(0, 1, 2, 3), c(1, 4, 5), c(1, 4, 6)); .pyT(list(lapply(e[[1]], .pyL), e[[2]])) })), "([[0, 1], [1, 2], [1, 4], [2, 3], [4, 5], [4, 6]], 7)")
})

test_that("prufer: size matches SymPy", {
  expect_identical(.pyrun(local({ Prufer(rep(0, 4))$size })), "1296")
  expect_identical(.pyrun(local({ Prufer(rep(6, 4))$size })), "1296")
  expect_identical(.pyrun(local({ Prufer(c(1, 2))$size })), "16")
})

test_that("prufer: tree matches SymPy", {
  expect_identical(.pyrun(local({ lapply(Prufer(list(c(0, 1), c(1, 5), c(0, 3), c(0, 2), c(2, 6), c(4, 7), c(2, 4)))$tree_repr, .pyL) })), "[[0, 1], [1, 5], [0, 3], [0, 2], [2, 6], [4, 7], [2, 4]]")
  expect_identical(.pyrun(local({ .pyL(Prufer(list(c(0, 1), c(1, 5), c(0, 3), c(0, 2), c(2, 6), c(4, 7), c(2, 4)))$prufer_repr) })), "[0, 1, 0, 2, 2, 4]")
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(1, 5), c(0, 3), c(0, 2), c(2, 6), c(4, 7), c(2, 4)))$rank })), "4244")
})

test_that("prufer: err matches SymPy", {
  expect_identical(.pyrun(local({ Prufer(list(c(1, 2), c(3, 4))) })), "ERR ValueError: Node 0 is missing.")
  expect_identical(.pyrun(local({ Prufer(list(c(2, 3), c(3, 4))) })), "ERR ValueError: Nodes [0, 1] are missing.")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 2), c(3, 4)); .pyL(Prufer(e[[1]], e[[2]])$prufer_repr) })), "[1, 3]")
  expect_identical(.pyrun(local({ prufer_edges(c(1, 3), c(3, 4)) })), "ERR ValueError: Node 3 is missing.")
  expect_identical(.pyrun(local({ prufer_edges(c(1, 2), c(5, 6)) })), "ERR ValueError: Nodes [4, 5] are missing.")
  expect_identical(.pyrun(local({ prufer_edges(c(1, 2), c(5, 7)) })), "ERR ValueError: Nodes [4, 5, 7] are missing.")
  expect_identical(.pyrun(local({ Prufer(list(integer())) })), "ERR ValueError: Prufer expects at least one edge in the tree.")
})

test_that("prufer: next matches SymPy", {
  expect_identical(.pyrun(local({ lapply(prufer_next(Prufer(list(c(0, 1), c(0, 2), c(0, 3))))$tree_repr, .pyL) })), "[[0, 2], [0, 1], [1, 3]]")
  expect_identical(.pyrun(local({ prufer_next(Prufer(list(c(0, 1), c(0, 2), c(0, 3))))$rank })), "1")
  expect_identical(.pyrun(local({ .pyL(prufer_next(Prufer(list(c(0, 1), c(0, 2), c(0, 3))), 5)$prufer_repr) })), "[1, 1]")
  expect_identical(.pyrun(local({ .pyL(prufer_prev(Prufer(list(c(0, 1), c(0, 2), c(0, 3))))$prufer_repr) })), "[3, 3]")
})

test_that("prufer: alias matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(Prufer(list(c(0, 1), c(1, 2), c(1, 3)))$prufer_repr) })), "[1, 1]")
  expect_identical(.pyrun(local({ lapply(Prufer(c(1, 1))$tree_repr, .pyL) })), "[[0, 1], [1, 2], [1, 3]]")
})

test_that("prufer: str matches SymPy", {
  expect_identical(.pyrun(local({ Prufer(c(1, 1)) })), "Prufer([1, 1])")
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(1, 2), c(1, 3))) })), "Prufer([[0, 1], [1, 2], [1, 3]])")
  expect_identical(.pyrun(local({ Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)), 5) })), "Prufer([[0, 1], [0, 2], [0, 3], [0, 4]], 5)")
  expect_identical(.pyrun(local({ prufer_unrank(6, 5) })), "Prufer([0, 1, 1])")
  expect_identical(.pyrun(local({ .pyL(prufer_unrank(6, 5)$prufer_repr) })), "[0, 1, 1]")
})

test_that("prufer: to matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(prufer_to_prufer(list(c(0, 1), c(0, 2), c(0, 3)), 4)) })), "[0, 0]")
  expect_identical(.pyrun(local({ lapply(prufer_to_tree(c(0, 0)), .pyL) })), "[[0, 1], [0, 2], [0, 3]]")
  expect_identical(.pyrun(local({ lapply(prufer_to_tree(c(3, 1, 2)), .pyL) })), "[[0, 3], [1, 3], [1, 2], [2, 4]]")
  expect_identical(.pyrun(local({ lapply(prufer_to_tree(integer()), .pyL) })), "[[0, 1]]")
})

test_that("prufer: exist matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(prufer_from_tree(list(c(0, 3), c(1, 3), c(2, 3), c(3, 4)), 5)) })), "[3, 3, 3]")
  expect_identical(.pyrun(local({ lapply(tree_from_prufer(c(3, 3, 3)), .pyL) })), "[[0, 3], [1, 3], [2, 3], [3, 4]]")
})

test_that("prufer: roundtrip matches SymPy", {
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 2)); t <- Prufer(e[[1]], e[[2]]); b <- c() - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([], True, [], 0)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 1, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(1) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0], True, [0], 0)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 3, 2)); t <- Prufer(e[[1]], e[[2]]); b <- c(3) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([2], True, [2], 2)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 2, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(2) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([1], True, [1], 1)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 1, 4), c(1, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(1, 1) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0, 0], True, [0, 0], 0)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(3, 2, 1, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(2, 1) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([1, 0], True, [1, 0], 4)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(3, 2, 1), c(2, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(2, 2) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([1, 1], True, [1, 1], 5)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 3, 2, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(3, 2) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([2, 1], True, [2, 1], 9)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 4, 2, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(4, 2) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([3, 1], True, [3, 1], 13)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(3, 1, 4, 2)); t <- Prufer(e[[1]], e[[2]]); b <- c(4, 1) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([3, 0], True, [3, 0], 12)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(4, 2, 1, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(1, 2) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0, 1], True, [0, 1], 1)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 2, 4, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(2, 4) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([1, 3], True, [1, 3], 7)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 3, 4, 2)); t <- Prufer(e[[1]], e[[2]]); b <- c(3, 4) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([2, 3], True, [2, 3], 11)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 4, 1), c(4, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(4, 4) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([3, 3], True, [3, 3], 15)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 2, 3, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(2, 3) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([1, 2], True, [1, 2], 6)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 3, 1), c(3, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(3, 3) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([2, 2], True, [2, 2], 10)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(1, 4, 3, 2)); t <- Prufer(e[[1]], e[[2]]); b <- c(4, 3) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([3, 2], True, [3, 2], 14)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 1, 4, 3)); t <- Prufer(e[[1]], e[[2]]); b <- c(1, 4) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0, 3], True, [0, 3], 3)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(2, 1, 3, 4)); t <- Prufer(e[[1]], e[[2]]); b <- c(1, 3) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0, 2], True, [0, 2], 2)")
  expect_identical(.pyrun(local({ e <- prufer_edges(c(6, 2, 1, 4), c(1, 3, 5, 8), c(3, 7)); t <- Prufer(e[[1]], e[[2]]); b <- c(1, 2, 1, 3, 3, 5) - 1; srt <- function(z) z[order(vapply(z, `[`, 0, 1), vapply(z, `[`, 0, 2))]; .pyT(list(.pyL(t$prufer_repr), identical(lapply(srt(Prufer(b)$tree_repr), as.numeric), lapply(srt(t$tree_repr), as.numeric)), .pyL(prufer_unrank(t$rank, e[[2]])$prufer_repr), t$rank)) })), "([0, 1, 0, 2, 2, 4], True, [0, 1, 0, 2, 2, 4], 4244)")
})

