# Expected values are SymPy's outputs for the same calls (sympy/logic/tests/test_boolalg.py (ANF, Xnor, Exclusive, gateinputcount));
# the rows come from the battery that compares this port with SymPy.

test_that("boolalg: test_Xnor", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor() })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A) })), "~A")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A, A) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(TRUE, A, A) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A, A, A, A, A) })), "~A")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(TRUE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(FALSE) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(TRUE, TRUE) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(TRUE, FALSE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(FALSE, FALSE) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(FALSE, A) })), "~A")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(TRUE, FALSE, FALSE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(FALSE, FALSE, A) })), "~A")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A, B) })), "~(A ^ B)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A, B, C) })), "~(A ^ B ^ C)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Xnor(A, Not(B)) })), "~(A ^ ~B)")
})

test_that("boolalg: test_Exclusive", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(FALSE, FALSE, FALSE) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(TRUE, TRUE, FALSE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(TRUE, TRUE, TRUE) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(A, B) })), "~(A & B)")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive() })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(A) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); Exclusive(A, TRUE) })), "~A")
})

test_that("boolalg: test_distribute", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); distribute_and_over_or(A) })), "A")
})

test_that("boolalg: test_to_anf", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); to_anf(And(x, y)) })), "x & y")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); to_anf(Not(x)) })), "x ^ True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); to_anf(S_true) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); to_anf(x) })), "x")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); to_anf(ITE(x, y, z)) })), "ITE(x, y, z)")
})

test_that("boolalg: test_is_anf", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); is_anf(S_true) })), "True")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); is_anf(Not(x)) })), "False")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); is_anf(Lt(x, 3)) })), "True")
})

test_that("boolalg: test_anf_coeffs", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyL(anf_coeffs(c(1, 1, 0, 1))) })), "[1, 0, 1, 1]")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyL(anf_coeffs(c(1))) })), "[1]")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); .pyL(anf_coeffs(c(1, 0, 1))) })), "ERR ValueError: The number of truth values must be a power of two, got 3")
})

test_that("boolalg: test_bool_minterm", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); bool_monomial(0, list(x, y)) })), "True")
})

test_that("boolalg: test_gateinputcount", {
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); gateinputcount(And(a, TRUE)) })), "0")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); gateinputcount(a) })), "0")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); gateinputcount(Implies(a, b)) })), "0")
  expect_identical(.pyrun(local({ A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C"); D <- Symbol("D"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); gateinputcount(Xnor(a, b)) })), "3")
})

