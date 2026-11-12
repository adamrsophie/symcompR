# Expected values are SymPy's outputs for the same calls (unify).

skip_unless_full_tests()

for (v in c("a", "b", "c", "d", "n", "p", "q", "x", "y", "z")) assign(v, Symbol(v))
bi <- Symbol("b", integer = TRUE)
fmtd <- function(d) { ks <- sort(names(d), method = "radix"); paste0("{", paste(paste0(ks, ": ", vapply(ks, function(k) sstr(d[[k]]), "")), collapse = ", "), "}") }
fmt <- function(r) {
  if (is_sym(r)) return(sstr(r))
  if (is.list(r) && !is.null(names(r)) && length(r)) return(fmtd(r))
  if (is.list(r)) return(paste0("[", paste(vapply(r, fmt, ""), collapse = ", "), "]"))
  format(r)
}

test_that("unify match SymPy", {
  expect_equal(fmt(unify(new_sym("Add", list(Integer(1), Integer(2), Integer(3))), new_sym("Add", list(a, b, c)), variables = list(a, b, c))), "[{a: 1, b: 2, c: 3}, {a: 1, b: 3, c: 2}, {a: 2, b: 1, c: 3}, {a: 3, b: 1, c: 2}, {a: 2, b: 3, c: 1}, {a: 3, b: 2, c: 1}]")
  expect_equal(fmt(unify(new_sym("Add", list(Integer(1), Integer(2), Integer(3))), new_sym("Add", list(a, c)), variables = list(a, c))), "[{a: 1, c: 2 + 3}, {a: 1 + 2, c: 3}, {a: 1, c: 2 + 3}, {a: 1 + 3, c: 2}, {a: 2, c: 1 + 3}, {a: 1 + 2, c: 3}, {a: 2, c: 1 + 3}, {a: 2 + 3, c: 1}, {a: 3, c: 1 + 2}, {a: 1 + 3, c: 2}, {a: 3, c: 1 + 2}, {a: 2 + 3, c: 1}]")
  expect_equal(fmt(unify(sin(x) + cos(x)^2, sin(p) + cos(p)^2, variables = list(p, q))), "[{p: x}]")
  expect_equal(fmt(unify(FiniteSet(a, b, c), FiniteSet(x, y), variables = list(x, y))), "[{x: a, y: {b, c}}, {x: {a, b}, y: c}, {x: a, y: {b, c}}, {x: {a, c}, y: b}, {x: b, y: {a, c}}, {x: {a, b}, y: c}, {x: b, y: {a, c}}, {x: {b, c}, y: a}, {x: c, y: {a, b}}, {x: {a, c}, y: b}, {x: c, y: {a, b}}, {x: {b, c}, y: a}]")
  expect_equal(fmt(unify(And(Gt(z, 0), Lt(n, 3)), new_sym("And", list(x, y)), variables = list(x, y))), "[{x: z > 0, y: n < 3}]")
  expect_equal(fmt(unify(Tuple(1, 2), Tuple(1, x), variables = list(x))), "[{x: 2}]")
  expect_equal(fmt(unify(Tuple(1, 2), Tuple(a, b), list(list(a, 5)), variables = list(a, b))), "[]")
  expect_equal(fmt(unify(Tuple(1, 2), Tuple(a, b), variables = list(a, b))), "[{a: 1, b: 2}]")
  expect_equal(fmt(unify(x*y*z, a*b, variables = list(a, b))), "[{a: x, b: y*z}, {a: x*y, b: z}, {a: x, b: y*z}, {a: x*z, b: y}, {a: y, b: x*z}, {a: x*y, b: z}, {a: y, b: x*z}, {a: y*z, b: x}, {a: z, b: x*y}, {a: x*z, b: y}, {a: z, b: x*y}, {a: y*z, b: x}]")
  expect_equal(fmt(unify(2*x + 3, a*x + b, variables = list(a, b))), "[{a: 2, b: 3}]")
  expect_equal(fmt(unify(exp(x)*sin(y), exp(a)*b, variables = list(a, b))), "[{a: x, b: sin(y)}]")
  expect_equal(fmt(unify(sin(Integer(3))*sin(Integer(4))*sin(Integer(5)) + 4*cos(Integer(3))*cos(Integer(4)), a*cos(b)*cos(c) + d*sin(b)*sin(c), variables = list(a, b, c, d))[[1]]), "{a: 4, b: 3, c: 4, d: sin(5)}")
  expect_equal(fmt(rewriterule(p^2, p^3, variables = list(p))(x^2)), "[x^3]")
  expect_equal(fmt(rewriterule(p^2 + q^3, (p*q)^4, list(p, q))(x^2 + y^3)), "[x^4*y^4]")
  expect_equal(fmt(rewriterule(sin(p)^2 + sin(p)^2, 1, list(p, q))(sin(x)^2 + sin(x)^2)), "[1]")
  expect_equal(fmt(rewriterule(x, x + 1, list(x), function(x) .is_true(Lt(x, 10)))(Integer(15))), "[]")
  expect_equal(fmt(rewriterule(x, x + 1, list(x), function(x) .is_true(Lt(x, 10)))(Integer(5))), "[1 + 5]")
  expect_equal(fmt(rewriterule(x + y, x^y, list(x, y), function(x, y) isTRUE(ask_fact(x, "integer")))(a + bi)), "[b^a]")
  expect_equal(fmt(rewriterule(p + q, q + p, list(p, q))(x + y)), "[x + y, x + y]")
})
