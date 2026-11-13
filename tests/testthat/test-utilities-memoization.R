# Expected values are SymPy's outputs for the same calls (sympy/utilities/memoization.py);
# the rows come from the battery that compares this port with SymPy.

test_that("memoization: rec", {
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); f(0) })), "1")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); f(5) })), "120")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); f(20) })), "2432902008176640000")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); lapply(0:7, f) })), "[1, 1, 2, 6, 24, 120, 720, 5040]")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); .pyT(list(f(10), f(3))) })), "(3628800, 6)")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); f(6); attr(f, 'cache_length')() })), "7")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); f(6); attr(f, 'fetch_item')(4) })), "24")
  expect_identical(.pyrun(local({ f <- recurrence_memo(list(Integer(1)))(function(n, prev) n * prev[[length(prev)]]); attr(f, 'cache_length')() })), "1")
  expect_identical(.pyrun(local({ fib <- recurrence_memo(list(Integer(0), Integer(1)))(function(n, prev) prev[[length(prev)]] + prev[[length(prev) - 1L]]); fib(1) })), "1")
  expect_identical(.pyrun(local({ fib <- recurrence_memo(list(Integer(0), Integer(1)))(function(n, prev) prev[[length(prev)]] + prev[[length(prev) - 1L]]); fib(30) })), "832040")
  expect_identical(.pyrun(local({ fib <- recurrence_memo(list(Integer(0), Integer(1)))(function(n, prev) prev[[length(prev)]] + prev[[length(prev) - 1L]]); fib(100) })), "354224848179261915075")
  expect_identical(.pyrun(local({ fib <- recurrence_memo(list(Integer(0), Integer(1)))(function(n, prev) prev[[length(prev)]] + prev[[length(prev) - 1L]]); lapply(0:11, fib) })), "[0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]")
  expect_identical(.pyrun(local({ fib <- recurrence_memo(list(Integer(0), Integer(1)))(function(n, prev) prev[[length(prev)]] + prev[[length(prev) - 1L]]); fib(50); attr(fib, 'cache_length')() })), "51")
  expect_identical(.pyrun(local({ x <- Symbol("x"); h <- recurrence_memo(list(Integer(1), 2*x))(function(n, prev) 2*x*prev[[length(prev)]] - 2*(n - 1)*prev[[length(prev) - 1L]]); expand(h(6)) })), "64*x^6 - 480*x^4 + 720*x^2 - 120")
  expect_identical(.pyrun(local({ x <- Symbol("x"); h <- recurrence_memo(list(Integer(1), 2*x))(function(n, prev) 2*x*prev[[length(prev)]] - 2*(n - 1)*prev[[length(prev) - 1L]]); h(2) })), "4*x^2 - 2")
})

test_that("memoization: assoc", {
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); s(0, 0) })), "1")
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); s(4, 2) })), "11")
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); s(6, 3) })), "225")
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); lapply(0:5, function(k) s(5, k)) })), "[0, 24, 50, 35, 10, 1]")
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); s(7, 1) })), "720")
  expect_identical(.pyrun(local({ s <- assoc_recurrence_memo(function(i) if (i == 0) Integer(1) else Integer(0))(function(i, j, prev) if (j < i) prev[[i]][[j]] + (i - 1) * prev[[i]][[j + 1L]] else Integer(1)); s(8, 4); s(3, 3) })), "1")
  expect_identical(.pyrun(local({ b <- assoc_recurrence_memo(function(i) Integer(1))(function(i, j, prev) prev[[i]][[j]] + (if (j < i) prev[[i]][[j + 1L]] else Integer(0))); b(5, 2) })), "10")
  expect_identical(.pyrun(local({ b <- assoc_recurrence_memo(function(i) Integer(1))(function(i, j, prev) prev[[i]][[j]] + (if (j < i) prev[[i]][[j + 1L]] else Integer(0))); lapply(0:6, function(k) b(6, k)) })), "[1, 6, 15, 20, 15, 6, 1]")
  expect_identical(.pyrun(local({ b <- assoc_recurrence_memo(function(i) Integer(1))(function(i, j, prev) prev[[i]][[j]] + (if (j < i) prev[[i]][[j + 1L]] else Integer(0))); b(10, 5) })), "252")
  expect_identical(.pyrun(local({ b <- assoc_recurrence_memo(function(i) Integer(1))(function(i, j, prev) prev[[i]][[j]] + (if (j < i) prev[[i]][[j + 1L]] else Integer(0))); b(4, 4) })), "1")
})

