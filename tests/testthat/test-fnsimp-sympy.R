skip_unless_full_tests()

# Expected values are SymPy outputs (str() with ** shown as ^) for the
# same calls, many taken from sympy/simplify/tests/test_cse.py and the
# function test files.

x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); w <- Symbol("w")
a <- Symbol("a"); b <- Symbol("b")
n <- Symbol("n", integer = TRUE); p <- Symbol("p", positive = TRUE); xr <- Symbol("xr", real = TRUE)
f <- Function("f")
.fs_T <- function(v) {
  v <- if (is.list(v) && !is_sym(v)) unname(v) else as.list(v)
  do.call(Tuple, lapply(v, function(e) if (is.character(e)) sympify(e) else as_sym(e)))
}
.cse_str <- function(r) {
  rp <- paste0("[", paste(vapply(r$replacements, function(q) paste0("(", sstr(q[[1]]), ", ", sstr(q[[2]]), ")"), ""), collapse = ", "), "]")
  rd <- if (is_sym(r$reduced)) sstr(r$reduced) else paste0("[", paste(vapply(r$reduced, sstr, ""), collapse = ", "), "]")
  paste0("(", rp, ", ", rd, ")")
}
.lst <- function(v) paste0("['", paste(vapply(v, sstr, ""), collapse = "', '"), "']")
.fs_fmt <- function(v) {
  if (inherits(v, "Quaternion")) return(format(v))
  if (is_sym(v) || inherits(v, "symmatrix")) return(sstr(v))
  if (is.logical(v) && length(v) == 1L) return(if (v) "True" else "False")
  paste(as.character(v), collapse = " ")
}
.fs_err <- function(ex) tryCatch({ force(ex); "no error" }, error = function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition"))
  paste0("ERR ", if (length(cl)) cl[1] else "Error", ": ", conditionMessage(e))
})

test_that("expand(diracdelta = TRUE) splits DiracDelta over simple roots", {
  expect_identical(.fs_fmt(expand(DiracDelta(2*x), diracdelta = TRUE, wrt = x)), "DiracDelta(x)/2")
  expect_identical(.fs_fmt(expand(DiracDelta(x^2 - 1), diracdelta = TRUE, wrt = x)), "DiracDelta(x - 1)/2 + DiracDelta(x + 1)/2")
  expect_identical(.fs_fmt(expand(DiracDelta(x*y), diracdelta = TRUE, wrt = x)), "DiracDelta(x)/Abs(y)")
  expect_identical(.fs_fmt(expand(DiracDelta(x^2 + 1), diracdelta = TRUE)), "DiracDelta(x^2 + 1)")
  expect_identical(.fs_fmt(expand(DiracDelta(x^3 - x), diracdelta = TRUE)), "DiracDelta(x) + DiracDelta(x - 1)/2 + DiracDelta(x + 1)/2")
  expect_identical(.fs_fmt(expand(DiracDelta(x^2 - 2*x + 1), diracdelta = TRUE)), "DiracDelta(x^2 - 2*x + 1)")
  expect_identical(.fs_err(expand(DiracDelta(x*y), diracdelta = TRUE)), "ERR TypeError: \nWhen there is more than 1 free symbol or variable in the expression,\nthe 'wrt' keyword is required as a hint to expand when using the\nDiracDelta hint.")
})

test_that("rewrite knows the erfc, sign, hyper, SingularityFunction, asin and tractable targets", {
  expect_identical(.fs_fmt(rewrite(erf(x), "erfc")), "1 - erfc(x)")
  expect_identical(.fs_fmt(rewrite(Heaviside(x), "sign")), "Heaviside(x)")
  expect_identical(.fs_fmt(rewrite(Heaviside(xr), "sign")), "sign(xr)/2 + 1/2")
  expect_identical(.fs_fmt(rewrite(Abs(x), "sign")), "x/sign(x)")
  expect_identical(.fs_fmt(rewrite(erf(x), "hyper")), "2*x*hyper((1/2,), (3/2,), -x^2)/sqrt(pi)")
  expect_identical(.fs_fmt(rewrite(DiracDelta(x), "SingularityFunction")), "SingularityFunction(x, 0, -1)")
  expect_identical(.fs_fmt(rewrite(acos(x), "asin")), "-asin(x) + pi/2")
  expect_identical(.fs_fmt(rewrite(cosh(x), "tractable")), "exp(x)/2 + exp(-x)/2")
  expect_identical(.fs_fmt(rewrite(sin(x), "nonsense")), "sin(x)")
  expect_identical(.fs_fmt(rewrite(Heaviside(x-5), "SingularityFunction")), "SingularityFunction(x, 5, 0)")
  expect_identical(.fs_err(rewrite(DiracDelta(x-y), "SingularityFunction")), "ERR TypeError: \nrewrite(SingularityFunction) does not support arguments with more that\none variable.")
  expect_identical(.fs_fmt(rewrite(atan(x), "asin")), "(-asin(1/sqrt(x^2 + 1)) + pi/2)*sqrt(x^2)/x")
  expect_identical(.fs_fmt(rewrite(erfi(x), "erfc")), "I*erfc(I*x) - I")
  expect_identical(.fs_fmt(rewrite(elliptic_k(x), "hyper")), "pi*hyper((1/2, 1/2), (1,), x)/2")
  expect_identical(.fs_fmt(rewrite(gamma(x), "tractable")), "exp(loggamma(x))")
  expect_identical(.fs_fmt(rewrite(Heaviside(x), sign)), "Heaviside(x)")
})

test_that("simplify of And uses the relational patterns", {
  expect_identical(.fs_fmt(simplify(And(Eq(x, 2), x > 1))), "Eq(x, 2)")
  expect_identical(.fs_fmt(simplify(And(x >= 1, x <= 1))), "Eq(x, 1)")
  expect_identical(.fs_fmt(simplify(And(x > 1, x >= 2))), "x >= 2")
  expect_identical(.fs_fmt(simplify(And(x < 1, x > 3))), "False")
  expect_identical(.fs_fmt(simplify(And(Eq(x, 1), y > x))), "Eq(x, 1) & (y > 1)")
  expect_identical(.fs_fmt(simplify(And(x <= 1, Ne(x, 1)))), "x < 1")
})

test_that("nsolve follows mpmath findroot at the requested precision", {
  expect_identical(.fs_fmt(nsolve(x^2 - 2, x, 1, prec = 30)), "1.41421356237309504880168872421")
  expect_identical(.fs_fmt(nsolve(cos(x) - x, x, 1)), "0.739085133215161")
  expect_identical(.fs_fmt(nsolve(exp(x) - 3, x, 1, prec = 50)), "1.0986122886681096913952452369225257046474905578227")
  expect_identical(.fs_fmt(nsolve(x^2 + 1, x, sym_I)), "1.0*I")
  expect_identical(.fs_fmt(.lst(nsolve(list(x^2 + y^2 - 4, x - y), list(x, y), list(1, 1)))), "['1.41421356237310', '1.41421356237310']")
  expect_identical(.fs_err(nsolve(x^2 + 1, x, 1)), "ERR ValueError: Could not find root within given tolerance. (5.22912334436213473229 > 2.16840434497100886801e-19)\nTry another starting point or tweak arguments.")
})

test_that("dsolve keeps both branches of a separable Riccati equation", {
  expect_identical(.fs_fmt(.fs_T(dsolve(diff(f(x), x) - x*f(x) - x*f(x)^2, f(x)))), "(Eq(f(x), (sqrt(exp(2*C1 + x^2)) - exp(2*C1 + x^2))/(exp(2*C1 + x^2) - 1)), Eq(f(x), -(sqrt(exp(2*C1 + x^2)) + exp(2*C1 + x^2))/(exp(2*C1 + x^2) - 1)))")
})

test_that("heurisch gives field normal coefficients", {
  expect_identical(.fs_fmt(integrate(1/(sin(x) + cos(x)), x)),
                   "sqrt(2)*log(tan(x/2) - 1 + sqrt(2))/2 - sqrt(2)*log(tan(x/2) - sqrt(2) - 1)/2")
})

test_that("parametric integrals start with the degenerate case", {
  # the later pieces of SymPy's result depend on its hash seed
  expect_identical(sstr(integrate(exp(a*x)*sin(b*x), x)$args[[1L]]), "ExprCondPair(0, Eq(a, 0) & Eq(b, 0))")
  expect_identical(sstr(integrate(sin(a*x)*cos(b*x), x)$args[[1L]]), "ExprCondPair(0, Eq(a, 0) & Eq(b, 0))")
})
