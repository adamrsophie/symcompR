# Expected values are SymPy's outputs for the same calls (sympy/logic/tests/test_inference.py and test_dimacs.py);
# the rows come from the battery that compares this port with SymPy.

test_that("inference: test_literal", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); literal_symbol(Not(And(A, B))) })), "ERR ValueError: Argument must be a boolean literal.")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); literal_symbol(Not(Not(A))) })), "A")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); literal_symbol(S_true) })), "ERR ValueError: Argument must be a boolean literal.")
})

test_that("inference: test_find_pure_symbol", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol(list(A), list(A))) })), "(A, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol(list(A, B, C), list(Or(Not(A), Not(B)), Or(Not(B), Not(C)), Or(C, A)))) })), "(B, False)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol(list(A, B), list())) })), "(None, None)")
})

test_that("inference: test_find_pure_symbol_int_repr", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2, 3), list(c(1, -2), c(-2, -3), c(3, 1)))) })), "(1, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2, 3), list(c(-1, 2), c(2, -3), c(3, 1)))) })), "(2, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2, 3), list(c(-1, -2), c(-2, -3), c(3, 1)))) })), "(2, False)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2, 3), list(c(-1, 2), c(-2, -3), c(3, 1)))) })), "(None, None)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(3, 2, 1), list(c(-1, -2), c(-3, -2)))) })), "(3, False)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2, 3, 4), list(c(-1, 2), c(-2, -4), c(3, -1)))) })), "(3, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_pure_symbol_int_repr(list(1, 2), list())) })), "(None, None)")
})

test_that("inference: test_unit_clause", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(A), list())) })), "(A, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B)), list(A = TRUE))) })), "(B, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B)), list(B = TRUE))) })), "(A, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B, C), Or(B, Not(C)), Or(A, Not(B))), list(A = TRUE))) })), "(B, False)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B, C), Or(B, Not(C)), Or(A, B)), list(A = TRUE))) })), "(B, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B, C), Or(B, Not(C)), A), list())) })), "(A, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B), Or(B, C)), list())) })), "(None, None)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause(list(Or(A, B), Or(Not(B), C)), list(A = TRUE, B = FALSE))) })), "(C, True)")
})

test_that("inference: test_unit_clause_int_repr", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1)), list())) })), "(1, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1, 2)), list(`1` = TRUE))) })), "(2, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1, 2)), list(`2` = TRUE))) })), "(1, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1, 2, 3), c(2, -3), c(1, -2)), list(`1` = TRUE))) })), "(2, False)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1, 2, 3), c(3, -3), c(1, 2)), list(`1` = TRUE))) })), "(2, True)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(1, 2), c(2, 3)), list())) })), "(None, None)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyT(find_unit_clause_int_repr(list(c(-1, 2), c(-2, 3)), list(`1` = TRUE, `2` = TRUE))) })), "(3, True)")
})

test_that("inference: test_unit_propagate", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); unit_propagate(list(And(A, B), Or(Not(A), B)), A) })), "[A & B, B]")
})

test_that("inference: test_unit_propagate_int_repr", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); lapply(unit_propagate_int_repr(list(c(-1, 2, 3), c(1, -2)), 2), function(s) .pyL(sort(s))) })), "[[1]]")
})

test_that("inference: pl_true_int_repr", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); pl_true_int_repr(c(1, -2), list(`1` = FALSE, `2` = TRUE)) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); pl_true_int_repr(c(-1), list(`1` = FALSE)) })), "True")
})

test_that("inference: test_satisfiable", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); satisfiable(And(A, B), algorithm = 'nonexistent') })), "ERR NotImplementedError: ")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); satisfiable(And(A, B), algorithm = 'minisat22', use_lra_theory = TRUE) })), "ERR ValueError: Currently only dpll2 can handle using lra theory. minisat22 is not handled.")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); satisfiable(FALSE, all_models = TRUE) })), "[False]")
})

test_that("inference: test_valid", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); valid(S_false) })), "False")
})

test_that("inference: test_pl_true", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); pl_true(FALSE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); pl_true(A) })), "None")
})

test_that("inference: test_KB", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); KB()$ask(A) })), "ERR NotImplementedError: ")
})

test_that("dimacs: doc_load", {
  expect_identical(.pyrun(local({ dimacs_load('1') })), "cnf_1")
})

