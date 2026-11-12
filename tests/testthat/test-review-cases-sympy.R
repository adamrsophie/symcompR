# Expected values are SymPy's outputs for the same calls: rows whose earlier
# expressions dropped the assumptions of SymPy's own tests, restored here.

skip_unless_full_tests()

skip_on_cran()

.bt <- new.env(parent = asNamespace("symcompR"))
eval(parse(text = r"---(x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); a <- Symbol("a"); b <- Symbol("b"); w <- Symbol("w")
n <- Symbol("n", integer = TRUE)
fmt <- function(v) {
  if (is.null(v)) return("None")
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (inherits(v, c("HolonomicFunction", "DifferentialOperator", "HolonomicSequence"))) return(format(v))
  if (inherits(v, "sym_generator")) return("<generator object Expr.series.<locals>.<genexpr>>")
  if (is.list(v) && !is_sym(v)) return(paste0(if (isTRUE(attr(v, "tuple"))) "(" else "[", paste(vapply(v, fmt, ""), collapse = ", "), if (isTRUE(attr(v, "tuple"))) ")" else "]"))
  if (is_sym(v)) return(gsub("\n", " ", sstr(v)))
  paste(format(v), collapse = " ")
}
tup <- function(v) structure(v, tuple = TRUE))---"), envir = .bt)
.run <- function(code) get("fmt", envir = .bt)(eval(parse(text = code), envir = .bt))

test_that("corrected cases match SymPy", {
  expect_equal(.run("local({ x <- Symbol(\"x\", real = TRUE); k <- Symbol(\"k\", positive = TRUE); n <- Symbol(\"n\", integer = TRUE, positive = TRUE); limit(x^n - x^(n - k - 1), x, oo) })"), "oo")  # p37:273
  expect_equal(.run("local({ k <- Symbol(\"k\", positive = TRUE); n <- Symbol(\"n\", positive = TRUE); limit((n + 1)^k / ((n + 1)^(k + 1) - n^(k + 1)), n, oo) })"), "1/(k + 1)")  # p37:388
  expect_equal(.run("local({ k <- Symbol(\"k\", positive = TRUE); n <- Symbol(\"n\", positive = TRUE); limit(expand((n + 1)^k / ((n + 1)^(k + 1) - n^(k + 1))), n, oo) })"), "1/(k + 1)")  # p37:389
  expect_equal(.run("local({ k <- Symbol(\"k\", positive = TRUE); n <- Symbol(\"n\", positive = TRUE); limit((n + 1)^k / (n * (-n^k + (n + 1)^k) + (n + 1)^k), n, oo) })"), "1/(k + 1)")  # p37:390
  expect_equal(.run("local({ b <- Symbol(\"b\", positive = TRUE); n <- Symbol(\"n\", integer = TRUE); limit(Mod(x, b), x, n * b, \"+\") })"), "0")  # p37:393
  expect_equal(.run("local({ b <- Symbol(\"b\", positive = TRUE); n <- Symbol(\"n\", integer = TRUE); limit(Mod(x, b), x, n * b, \"-\") })"), "b")  # p37:394
  expect_equal(.run("local({ n <- Symbol(\"n\", integer = TRUE); m <- Symbol(\"m\", integer = TRUE); fps_truncate(fps(x^(n + m) * sin(x), x)) })"), "x^(m + n + 5)/120 - x^(m + n + 3)/6 + x^(m + n + 1) + O(x^(m + n + 6))")  # p40:58
  expect_equal(.run("local({ i <- Symbol(\"i\", integer = TRUE); r <- IndexedBase(\"r\"); doit(Sum(r(i), list(i, 0, 3))) })"), "r[0] + r[1] + r[2] + r[3]")  # p21:102
  expect_equal(.run("local({ i <- Symbol(\"i\", integer = TRUE); j <- Symbol(\"j\", integer = TRUE); r <- IndexedBase(\"r\"); doit(Sum(r(i), list(i, j, j + 2))) })"), "r[j + 1] + r[j + 2] + r[j]")  # p21:103
  expect_equal(.run("local({ k <- Symbol(\"k\", integer = TRUE); A <- IndexedBase(\"A\"); doit(Sum(A(k), list(k, 1, 3))) })"), "A[1] + A[2] + A[3]")  # p21:104
  expect_equal(.run("local({ z <- Symbol(\"z\", positive = TRUE); n <- Symbol(\"n\", integer = TRUE, positive = TRUE); limit(binomial(n + z, n) * n^(-z), n, oo) })"), "1/gamma(z + 1)")  # p37:335
  expect_equal(.run("take(series(cos(x), x, 1, n = NULL), 6)"), "[cos(1), -(x - 1)*sin(1), -(x - 1)^2*cos(1)/2, (x - 1)^3*sin(1)/6, (x - 1)^4*cos(1)/24, -(x - 1)^5*sin(1)/120]")  # p35:5t
  expect_equal(.run("take(series(cos(x), x, 1, n = NULL, dir = \"-\"), 6)"), "[cos(1), (1 - x)*sin(1), -(1 - x)^2*cos(1)/2, -(1 - x)^3*sin(1)/6, (1 - x)^4*cos(1)/24, (1 - x)^5*sin(1)/120]")  # p35:6t
  expect_equal(.run("take(series(exp(x) / x, x, n = NULL), 6)"), "[1/x, 1, x/2, x^2/6, x^3/24, x^4/120]")  # p35:7t
})
