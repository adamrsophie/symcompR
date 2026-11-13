skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/strategies/tests/test_core.py, test_rl.py (debug, new));
# the rows come from the battery that compares this port with SymPy.

test_that("strategies: debug matches SymPy", {
  expect_identical(.pyrun(local({ posdec <- function(x) if (x > 0) x - 1L else x; tf <- tempfile(); rl <- strat_debug(posdec, tf); rl(5L); readChar(tf, 10000L) })), "Rule: posdec\\nIn:   5\\nOut:  4\\n\\n")
  expect_identical(.pyrun(local({ posdec <- function(x) if (x > 0) x - 1L else x; tf <- tempfile(); rl <- strat_debug(posdec, tf); rl(5L) })), "4")
  expect_identical(.pyrun(local({ posdec <- function(x) if (x > 0) x - 1L else x; tf <- tempfile(); rl <- strat_debug(posdec, tf); rl(0L); if (file.exists(tf)) readChar(tf, 10000L) else "" })), "")
  expect_identical(.pyrun(local({ posdec <- function(x) if (x > 0) x - 1L else x; tf <- tempfile(); rl <- strat_debug(posdec, tf); rl(3L); rl(1L); readChar(tf, 10000L) })), "Rule: posdec\\nIn:   3\\nOut:  2\\n\\nRule: posdec\\nIn:   1\\nOut:  0\\n\\n")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); sub1 <- function(e) subs(e, x, y); tf <- tempfile(); rl <- strat_debug(sub1, tf); rl(x + 1); readChar(tf, 10000L) })), "Rule: sub1\\nIn:   x + 1\\nOut:  y + 1\\n\\n")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); sub1 <- function(e) subs(e, x, y); tf <- tempfile(); rl <- strat_debug(sub1, tf); rl(x^2*y + x) })), "y^3 + y")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); sub1 <- function(e) subs(e, x, y); tf <- tempfile(); rl <- strat_debug(sub1, tf); rl(x^2*y + x); readChar(tf, 10000L) })), "Rule: sub1\\nIn:   x^2*y + x\\nOut:  y^3 + y\\n\\n")
})

test_that("strategies: new matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); strat_new(list(type = "Add"), Integer(1), Integer(2)) })), "1 + 2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); rebuild(strat_new(list(type = "Add"), Integer(1), Integer(2))) })), "3")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); strat_new(list(type = "Add"), x, x) })), "x + x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); rebuild(strat_new(list(type = "Add"), x, x)) })), "2*x")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); rebuild(strat_new(list(type = "Mul"), x, Integer(2), x)) })), "2*x^2")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); identical(strat_new(list(type = "Add"), Integer(1), Integer(2))$key, Integer(3)$key) })), "False")
})

