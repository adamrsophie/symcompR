# Expected values are SymPy's outputs for the same calls (sympy/solvers/tests/test_simplex.py);
# the rows come from the battery that compares this port with SymPy.

test_that("simplex: test_lp_err", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); x3 <- Symbol("x3"); x4 <- Symbol("x4"); y1 <- Symbol("y1"); y2 <- Symbol("y2"); .pyT(.lp("max", x, list(Gt(x, 0)))) })), "ERR TypeError: \\nonly equalities like Eq(x, y) or non-strict inequalities like x >= y\\nare allowed in lp, not x > 0")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); x3 <- Symbol("x3"); x4 <- Symbol("x4"); y1 <- Symbol("y1"); y2 <- Symbol("y2"); .pyT(.lp("max", x, list(Ne(x, 0)))) })), "ERR TypeError: \\nonly equalities like Eq(x, y) or non-strict inequalities like x >= y\\nare allowed in lp, not Ne(x, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); x3 <- Symbol("x3"); x4 <- Symbol("x4"); y1 <- Symbol("y1"); y2 <- Symbol("y2"); .pyT(lpmax(x, list(Le(x, oo)))) })), "ERR ValueError: only finite bounds are permitted")
})

test_that("simplex: test_simplex", {
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 1), list(-1, 1), list(0, 1), list(-1, 0)), list(5, 1, 2, -1), list(list(-1, -1)), list(1), dual = TRUE)) })), "(-6, [1, 0, 0, 0], [5, 0])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list()), list(), list(list(1)), list(0))) })), "(0, [0], [])")
})

test_that("simplex: simplex_extra", {
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 2), list(3, 1)), list(4, 5), list(list(-1, -1)), list(0), dual = FALSE)) })), "(-13/5, [6/5, 7/5], [2/5, 1/5])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 2), list(3, 1)), list(4, 5), list(list(-1, -1)), list(0), dual = TRUE)) })), "(-13/5, [2/5, 1/5], [6/5, 7/5])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 0), list(0, 1)), list(1, 1), list(list(-1, -2)), list(3), dual = FALSE)) })), "(-6, [1, 1], [1, 2])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, -1)), list(-2), list(list(1, 1)), list(0), dual = FALSE)) })), "(2, [0, 2], [1])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, -1)), list(-2), list(list(1, 1)), list(0), dual = TRUE)) })), "(2, [1], [0, 2])")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(-1, 0)), list(-1), list(list(-1, 0)), list(0), dual = FALSE)) })), "ERR UnboundedLPError: \\nObjective function can assume arbitrarily large values!")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(-1, 0)), list(-1), list(list(-1, 0)), list(0), dual = TRUE)) })), "ERR InfeasibleLPError: \\nThe constraint set is empty!")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 1)), list(-1), list(list(1, 1)), list(0), dual = FALSE)) })), "ERR InfeasibleLPError: \\nThe constraint set is empty!")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 1)), list(-1), list(list(1, 1)), list(0), dual = TRUE)) })), "ERR UnboundedLPError: \\nObjective function can assume arbitrarily large values!")
})

test_that("simplex: simplex_err", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); x3 <- Symbol("x3"); x4 <- Symbol("x4"); y1 <- Symbol("y1"); y2 <- Symbol("y2"); .pyT(.simplex(list(list(1, x)), list(1), list(list(1, 1)), list(0))) })), "ERR TypeError: \\nOnly rationals and floats are allowed.")
  expect_identical(.pyrun(local({ .pyT(.simplex(list(list(1, 1)), list(), list(list(1, 1)), list(0))) })), "ERR ValueError: must give A and B")
})

test_that("simplex: test_linprog", {
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1), list(list(1, 1)), list(5), bounds = list(`0` = list(NULL, NULL), `1` = list(NULL, NULL)))) })), "(0, [0, 0])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(-1, -1), list(list(1, 2), list(3, 1)), list(4, 5))) })), "(-13/5, [6/5, 7/5])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(-1, -1), list(list(1, 2), list(3, 1)), list(4, 5), bounds = list(0, NULL))) })), "(-13/5, [6/5, 7/5])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(Rational(1, 2), Rational(1, 3)), list(list(-1, -1)), list(-1))) })), "(1/3, [0, 1])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(list(1, 1), list(1, 1)))) })), "ERR ValueError: C must be a single row.")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1), list(list(1, 1, 1)), list(2))) })), "ERR ValueError: number of columns in A and C must match")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1), list(), list(2))) })), "ERR ValueError: A and b must both be given")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1), A_eq = NULL, b_eq = list(2))) })), "ERR ValueError: A_eq and b_eq must both be given")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1), list(list(1, 1)), list(-1))) })), "ERR InfeasibleLPError: \\nThe constraint set is empty!")
  expect_identical(.pyrun(local({ .pyT(linprog(list(Float(1.5), 1), list(list(-1, -1)), list(-2))) })), "(2.00000000000000, [0, 2])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1))) })), "(0, [0])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1, 1))) })), "(0, [0, 0])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(-1))) })), "ERR UnboundedLPError: \\nObjective function can assume arbitrarily large values!")
  expect_identical(.pyrun(local({ .pyT(linprog(list(1), bounds = list(0, NULL))) })), "(0, [0])")
  expect_identical(.pyrun(local({ .pyT(linprog(list(-1), bounds = list(0, NULL))) })), "ERR UnboundedLPError: \\nObjective function can assume arbitrarily large values!")
  expect_identical(.pyrun(local({ .pyT(linprog(list(-1), bounds = list(NULL, 4))) })), "(-4, [4])")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, -1), list(list(1, 1)), list(2), bounds = list(`1` = list(NULL, NULL)))) })), "(x1 - x2, [x1 + x2 <= 2, x1 >= 0])")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, 2))) })), "(x1 + 2*x2, [x1 >= 0, x2 >= 0])")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, 2), list(list(1, 1)))) })), "ERR ValueError: A and b must both be given")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, 2), list(), list(3))) })), "ERR ValueError: A and b must both be given")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, 2), list(list(1, 1, 1)), list(3))) })), "ERR ValueError: number of columns in A and C must match")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(1, 2), A_eq = list(list(1, 1)))) })), "ERR TypeError: \\nData type not understood; expecting list of lists or lists of values.")
  expect_identical(.pyrun(local({ .pyT(show_linprog(list(list(1, 2), list(3, 4)))) })), "ERR ValueError: C must be a single row.")
})

test_that("simplex: extra", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); x3 <- Symbol("x3"); x4 <- Symbol("x4"); y1 <- Symbol("y1"); y2 <- Symbol("y2"); .pyT(lpmax(x, list(Le(x*y, 3)))) })), "ERR NonlinearError: \\nnonlinear cross-term: (_z1 - x)*(_z2 - y)")
})

