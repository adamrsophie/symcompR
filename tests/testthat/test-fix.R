test_that("special functions follow SymPy", {
  x <- Symbol("x")
  expect_identical(sstr(subfactorial(5)), "44")
  expect_identical(sstr(subfactorial(10)), "1334961")
  expect_identical(sstr(Heaviside(oo)), "1")
  expect_identical(sstr(elliptic_f(x, 0)), "x")
  expect_identical(sstr(elliptic_pi(0, x)), "elliptic_k(x)")
  expect_identical(sstr(diff(asech(x), x)), "-sqrt(1/(x + 1))/(x*sqrt(1 - x))")
  expect_identical(sstr(Li(x)), "Li(x)")
})

test_that("number theory and logic follow SymPy", {
  p <- Symbol("p"); q <- Symbol("q")
  r <- integer_nthroot(gmp::as.bigz("1000000000000000000000000000000"), 3)
  expect_true(inherits(r[[1L]], "bigz"))
  expect_identical(as.character(r[[1L]]), "10000000000")
  expect_true(isTRUE(r[[2L]]))
  expect_identical(sstr(Or(p, Not(p))), "p | ~p")
  expect_identical(sstr(eliminate_implications(Equivalent(p, q))), "(p | ~q) & (q | ~p)")
  m <- satisfiable(And(Or(p, q), Not(And(p, q))), all_models = TRUE)
  expect_length(m, 2L)
})

test_that("traversal and sort keys follow SymPy", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(vapply(sym_args(x + y*z), sstr, ""), c("x", "y*z"))
  ks <- vapply(list(y, S(2), x^2, x), default_sort_key, "")
  expect_identical(order(ks), c(2L, 4L, 3L, 1L))
  expect_identical(sstr(bottom_up(x + 1, function(e) e)), "x + 1")
})

test_that("calculus helpers follow SymPy", {
  x <- Symbol("x"); t <- Symbol("t")
  expect_identical(sstr(pade(exp(x), x, 2, 2)), "(x^2/4 + 3*x/2 + 3)/(x^2/4 - 3*x/2 + 3)")
  expect_false(is_monotonic(x^3))
  expect_true(is_monotonic(exp(x)))
  expect_identical(sstr(line_integrate(x*Symbol("y"), Curve(list(t, 2*t), list(t, 0, 1)), list(x, Symbol("y")))), "2*sqrt(5)/3")
  out <- utils::capture.output(textplot(x^2, 0, 2, W = 30, H = 8))
  expect_identical(out[4L], "2.28571 |---------------------...------")
  expect_identical(out[9L], "         0              1.03448       2")
})

test_that("printing and drawing follow SymPy", {
  x <- Symbol("x"); n <- Symbol("n")
  expect_identical(latex(-oo), "-\\infty")
  expect_identical(latex(besselj(n, x)), "J_{n}\\left(x\\right)")
  expect_identical(format(ArrowStringDescription("mm", NULL, NULL, NULL, NULL, "d", "r", "_", "f")), "\\ar[dr]_{f}")
})
