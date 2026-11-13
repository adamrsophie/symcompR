# Expected values are SymPy's outputs for the same calls (sympy/parsing/tests/test_smtlib.py);
# the rows come from the battery that compares this port with SymPy.

test_that("smtlib: test assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun y () Int)\n    (reset)\n    ")$assertions })), "[]")
  expect_identical(.pyrun(local({ parse_smtlib("(some-unknown-command 1 2)")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: some-unknown-command")
})

test_that("smtlib: test symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const a-1 Int)")$symbols) })), "{'a-1': a-1}")
})

test_that("smtlib: test syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-fun |unclosed () S)")) })), "ERR SMTLibSyntaxError")
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-sort B 1 2)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 0 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("")$assertions })), "[]")
})

test_that("smtlib: extra 0 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("")$symbols) })), "{}")
})

test_that("smtlib: extra 1 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("; only a comment")$assertions })), "[]")
})

test_that("smtlib: extra 1 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("; only a comment")$symbols) })), "{}")
})

test_that("smtlib: extra 2 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(check-sat)")$assertions })), "[]")
})

test_that("smtlib: extra 2 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(check-sat)")$symbols) })), "{}")
})

test_that("smtlib: extra 3 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(exit)")$assertions })), "[]")
})

test_that("smtlib: extra 3 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(exit)")$symbols) })), "{}")
})

test_that("smtlib: extra 4 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-model)")$assertions })), "[]")
})

test_that("smtlib: extra 4 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-model)")$symbols) })), "{}")
})

test_that("smtlib: extra 5 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Int)")$assertions })), "[Q.integer(x)]")
})

test_that("smtlib: extra 6 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)")$assertions })), "[Q.real(x)]")
})

test_that("smtlib: extra 6 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 7 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Bool)")$assertions })), "[]")
})

test_that("smtlib: extra 7 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Bool)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 14 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const -5 Int)")$assertions })), "[Q.integer(-5)]")
})

test_that("smtlib: extra 14 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const -5 Int)")$symbols) })), "{'-5': -5}")
})

test_that("smtlib: extra 31 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (- 5))")$assertions })), "[-5]")
})

test_that("smtlib: extra 33 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (/ 1 2))")$assertions })), "[1/2]")
})

test_that("smtlib: extra 51 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (abs x))")$symbols) })), "{}")
})

test_that("smtlib: extra 55 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (is_int 2))")$assertions })), "[True]")
})

test_that("smtlib: extra 55 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (is_int 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 57 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (to_int x))")$assertions })), "[floor(x)]")
})

test_that("smtlib: extra 66 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (not false))")$assertions })), "[True]")
})

test_that("smtlib: extra 66 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (not false))")$symbols) })), "{}")
})

test_that("smtlib: extra 68 syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (and))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 71 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (or a))")$symbols) })), "{}")
})

test_that("smtlib: extra 77 assertions", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (xor a true))")$assertions })), "[~a]")
})

test_that("smtlib: extra 78 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 81 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> a))")$symbols) })), "{}")
})

test_that("smtlib: extra 151 symbols", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 208 syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 210 syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 215 syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(push)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 244 syntax", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(check-sat-assuming ((and a)))")) })), "ERR SMTLibSyntaxError")
})

