# Regression cases for core arithmetic, evaluation, noncommutative products, traversal, subs and caches.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("expressions with keys over 10000 bytes are cached", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y")
  e <- powsimp(expand((x + y + 1)^25))
  expect_length(e$args, 351L)
  expect_equal(substr(sstr(e), 1, 51), "x^25 + 25*x^24*y + 25*x^24 + 300*x^23*y^2 + 600*x^2")
})

test_that("traversal and sort keys follow SymPy", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(vapply(sym_args(x + y*z), sstr, ""), c("x", "y*z"))
  ks <- vapply(list(y, S(2), x^2, x), default_sort_key, "")
  expect_identical(order(ks), c(2L, 4L, 3L, 1L))
  expect_identical(sstr(bottom_up(x + 1, function(e) e)), "x + 1")
})

test_that("subs replaces a function that is an integration variable", {
  x <- Symbol("x"); yy <- Symbol("y"); f <- Function("f")
  expect_equal(sstr(subs(Integral(x * f(x), f(x)), f(x), yy)), "Integral(x*y, y)")
})

test_that("Add, Mul and Pow honor evaluate = FALSE", {
  skip_on_cran()
  x <- Symbol("x")
  expect_equal(latex(Mul(2, x, evaluate = FALSE)), "2 x")
  expect_equal(sstr(Add(x, 1, x, evaluate = FALSE)), "x + x + 1")
  expect_equal(sstr(Pow(2, 2, evaluate = FALSE)), "2^2")
  expect_equal(sstr(Mul(x, 2, evaluate = FALSE)), "x*2")
})

test_that("noncommutative products keep their order", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y")
  A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE)
  X <- Symbol("X", commutative = FALSE); Y <- Symbol("Y", commutative = FALSE)
  expect_equal(sstr(expand((1 + X)^(-1) * (1 + Y)^(-1))), "(1 + X)^(-1)*(1 + Y)^(-1)")
  expect_equal(sstr(expand(X^(-1) * (Y + 1)^(-1))), "X^(-1)*(1 + Y)^(-1)")
  expect_equal(sstr(TR0(A*B - B*A)), "A*B - B*A")
  expect_equal(sstr(TR9(cos(x)*A*B + cos(y)*B*A)), "cos(x)*A*B + cos(y)*B*A")
  expect_equal(sstr(x*A + B), "x*A + B")
})

test_that("as_numer_denom, count_ops and TableForm follow SymPy", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_equal(sstr(TR0(x^oo)), "x^oo")
  expect_equal(count_ops(Max(x, y, z)), 2)
  expect_equal(sstr(count_ops(Min(x, y, z), visual = TRUE)), "2*MIN")
  expect_equal(count_ops(Piecewise(list(1, Lt(x, 0)), list(2, TRUE))), 1)
  v <- gmp::as.bigz("12345678901234567890")
  expect_equal(.py_format_int("%.3d", v), "12345678901234567890")
  expect_equal(.py_format_int("%.5d", gmp::as.bigz(12)), "00012")
  expect_equal(.py_format_int("id=%d", v), "id=12345678901234567890")
  expect_equal(.py_format_int("%05d", gmp::as.bigz(-42)), "-0042")
})

test_that("exact rational powers follow Integer._eval_power", {
  skip_on_cran()
  cases <- list(
    list(Integer(12), Rational(1, 4), "sqrt(2)*3^(1/4)"),
    list(Integer(504), Rational(1, 6), "sqrt(2)*3^(1/3)*7^(1/6)"),
    list(Integer(18), Rational(5, 6), "3*2^(5/6)*3^(2/3)"),
    list(Integer(6480), Rational(1, 8), "5^(1/8)*sqrt(6)"),
    list(Integer(-4), Rational(1073741825, 2147483647), "2*(-1)^(1073741825/2147483647)*2^(3/2147483647)"),
    list(Integer(-12), Rational(1, 4), "(-3)^(1/4)*sqrt(2)"),
    list(Rational(3, 4), Rational(2, 3), "6^(2/3)/4"),
    list(Integer(4), Rational(-1, 2), "1/2"),
    list(Rational(-1, 2), Rational(-1, 3), "-(-1)^(2/3)*2^(1/3)"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]^cs[[2L]]), cs[[3L]])
  expect_equal(sstr(Integer(-1)^Integer(10000001)), "-1")
})

test_that("Dict order, WildFunction and Poly domains", {
  skip_on_cran()
  x <- Symbol("x")
  a <- Symbol("a")
  k <- Integer("9007199254740992"); j <- k + 1
  expect_true(Dict(list(k, j), list(0, 0)) == Dict(list(j, k), list(0, 0)))
  F <- WildFunction("F")
  m <- sym_match(F, F)
  expect_equal(sstr(m$F), "F_")
  expect_equal(capture.output(print(Poly(x / 2 + a, x))), "Poly(1/2*x + a, x, domain='QQ[a]')")
})

test_that("expand takes the fraction of a product as Mul._eval_expand_mul", {
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  z <- Symbol("z")
  a <- Symbol("a")
  expect_equal(sstr(expand((sqrt(z) / 2)^(2 - 2 * a))), "2^(2*a)*z^(1 - a)/4")
  expect_equal(sstr(expand(exp(-x) * (y + 1))), "y*exp(-x) + exp(-x)")
  fr <- fraction(exp(-x) * y)
  expect_equal(sstr(fr[[2L]]), "exp(x)")
  fr <- fraction(y / x^2, exact = TRUE)
  expect_equal(sstr(fr[[2L]]), "x^2")
})

test_that("products rebuild factor exponents as Mul.flatten does", {
  skip_on_cran()
  y <- Symbol("y")
  z <- Symbol("z")
  a <- Symbol("a")
  expect_equal(sstr(powsimp(Rational(1, 2)^(2 - 2 * a) * y, combine = "base")), "(1/2)^(2 - 2*a)*y")
  expect_equal(sstr(powdenest((sqrt(z) / 2)^(2 - 2 * a), polar = TRUE)), "(1/2)^(2 - 2*a)*z^(1 - a)")
  expect_equal(sstr(simplify((sqrt(z) / 2)^(2 - 2 * a))), "2^(2*a - 2)*z^(1 - a)")
})

test_that("ordered() counts a Float as half a node", {
  skip_on_cran()
  z <- Symbol("z")
  expect_equal(sstr(hyper(list(), list(Float(2) / 3, Rational(1, 3)), z)), "hyper((), (0.666666666666667, 1/3), z)")
})

test_that("powers of powers, signs of sums and commutativity follow SymPy", {
  skip_on_cran()
  y <- Symbol("y")
  a <- Symbol("a")
  k <- Symbol("k", integer = TRUE); p <- Symbol("p", positive = TRUE)
  expect_equal(sstr((Integer(-2)^y)^(k + 1)), "(-2)^(y*(k + 1))")
  expect_equal(sstr((p^2)^y), "p^(2*y)")
  d <- Dummy("d", integer = TRUE, positive = TRUE)
  expect_true(isTRUE(ask_fact(d - 1, "nonnegative")))
  expect_equal(sstr(Ge(d, 1)), "True")
  B <- Symbol("B", commutative = FALSE)
  expect_equal(sstr(conjugate(B) + conjugate(B)^2 + conjugate(B)^3), "conjugate(B) + conjugate(B)^2 + conjugate(B)^3")
  f <- Function("f")
  expect_equal(sstr(f(Integer(10)) + f(Integer(2)) + f(Integer(3))), "f(2) + f(3) + f(10)")
  m <- Symbol("m")
  expect_equal(sstr(apart(-2 * a^2 * m + 2 * a * m + 4 * m, m)), "-m*(2*a^2 - 2*a - 4)")
})

test_that("equations between noncommutative objects stay undecided", {
  B <- Symbol("B", commutative = FALSE)
  expect_equal(sstr(Eq(B, 1)), "Eq(B, 1)")
  expect_equal(sstr(Ne(B, 1)), "Ne(B, 1)")
  expect_equal(sstr(Eq(B^2, 1)), "Eq(B^2, 1)")
})

test_that("invalid input is refused", {
  expect_error(Integer(NA), "ValueError")
  expect_error(Integer(Inf), "ValueError")
  expect_error(normalized(Matrix(list(3, 4))), "TRUE or FALSE")
  x <- Symbol("x")
  expect_warning(m <- as.matrix(Matrix(list(list(1, 2), list(3, x)))), "not numbers")
  expect_equal(m[1, 2], 2)
})

test_that("exact numbers compare exactly when arguments are sorted", {
  a <- Integer("9007199254740992"); b <- Integer("9007199254740993")
  expect_equal(symcompR:::.basic_compare(a, b), -1L)
  expect_equal(symcompR:::.basic_compare(b, a), 1L)
  expect_equal(symcompR:::.basic_compare(Rational("9007199254740993/2"), Rational("9007199254740991/2")), 1L)
})

test_that("the stable sort agrees with insertion order on ties and long inputs", {
  v <- c(5, 3, 5, 1, 3, 9, 2, 2, 8, 1, 7, 3, 6, 4, 0, 5, 5)
  o <- symcompR:::.stable_order(length(v), function(i, j) sign(v[i] - v[j]))
  expect_equal(o, order(v, method = "radix"))
  xs <- lapply(sprintf("v%03d", 200:1), Symbol)
  expect_equal(sstr(do.call(Add, xs[1:3])), "v198 + v199 + v200")
})

test_that("cse accepts expressions with keys over 10000 bytes", {
  skip_on_cran()
  xs <- lapply(sprintf("x%04d", seq_len(1500)), Symbol)
  r <- cse(do.call(Add, xs))
  expect_length(r$replacements, 0L)
  expect_equal(r$reduced[[1L]]$key, do.call(Add, xs)$key)
})

test_that("a sign query on a large integer does not need primality", {
  n <- Integer(gmp::as.bigz(2)^2203 - 1)
  expect_true(ask_fact(n, "positive"))
  expect_false(ask_fact(n, "even"))
  expect_true(ask_fact(Integer(7), "prime"))
  expect_false(ask_fact(Integer(1), "prime"))
})
