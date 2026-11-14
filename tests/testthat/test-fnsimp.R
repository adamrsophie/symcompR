# Slow on one core: skipped on CRAN, run by devtools::test() and CI.
testthat::skip_on_cran()

# Expected values are SymPy outputs (str() with ** shown as ^) for the
# same calls.

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

test_that("products combine powers of minus one with the imaginary unit", {
  expect_identical(.fs_fmt(sym_I * (-1)^Rational(1, 3)), "(-1)^(5/6)")
  expect_identical(.fs_fmt(sym_I * (-1)^Rational(1, 3) * 2^Rational(5, 6)), "(-2)^(5/6)")
  expect_identical(.fs_fmt((-1)^Rational(1, 3) * (-1)^Rational(1, 4)), "(-1)^(7/12)")
  expect_identical(.fs_fmt(sym_I*(-1)^Rational(2,3)), "-(-1)^(1/6)")
  expect_identical(.fs_fmt(sym_I*sym_I*(-1)^Rational(1,3)), "-(-1)^(1/3)")
  expect_identical(.fs_fmt(3*x*sym_I*(-1)^Rational(3,4)), "-3*(-1)^(1/4)*x")
  expect_identical(.fs_fmt(sym_I*(-2)^Rational(1,3)), "(-1)^(5/6)*2^(1/3)")
})

test_that("nested floor and ceiling collapse", {
  expect_identical(.fs_fmt(ceiling(floor(x))), "floor(x)")
  expect_identical(.fs_fmt(floor(ceiling(x))), "ceiling(x)")
  expect_identical(.fs_fmt(ceiling(-floor(x))), "-floor(x)")
  expect_identical(.fs_fmt(floor(-ceiling(x)+1)), "floor(-ceiling(x)) + 1")
  expect_identical(.fs_fmt(floor(floor(x)+1)), "floor(x) + 1")
})

test_that("inverse trigonometric functions of trigonometric values evaluate", {
  expect_identical(.fs_fmt(acos(cos(Integer(1)))), "1")
  expect_identical(.fs_fmt(atan(tan(Integer(2)))), "2 - pi")
  expect_identical(.fs_fmt(atan(tan(Rational(1,2)))), "1/2")
  expect_identical(.fs_fmt(acot(cot(Integer(1)))), "1")
  expect_identical(.fs_fmt(acos(-cos(Integer(1)))), "-1 + pi")
  expect_identical(.fs_fmt(acos(sin(Integer(1)))), "-1 + pi/2")
  expect_identical(.fs_fmt(acos(cos(Integer(5)))), "-5 + 2*pi")
  expect_identical(.fs_fmt(atan(cot(Integer(1)))), "-1 + pi/2")
  expect_identical(.fs_fmt(acot(tan(Integer(3)))), "-3 + pi/2")
  expect_identical(.fs_fmt(acos(-sin(Integer(2)))), "-2 + 3*pi/2")
})

test_that("fibonacci and lucas follow SymPy for negative indices and polynomials", {
  expect_identical(.fs_fmt(fibonacci(Integer(5), x)), "x^4 + 3*x^2 + 1")
  expect_identical(.fs_fmt(fibonacci(Integer(-3))), "2")
  expect_identical(.fs_fmt(lucas(Integer(-3))), "-4")
  expect_identical(.fs_fmt(fibonacci(Integer(3), 2*y)), "4*y^2 + 1")
  expect_identical(.fs_fmt(lucas(Integer(10))), "123")
  expect_identical(.fs_err(fibonacci(Integer(0), x)), "ERR ValueError: Fibonacci polynomials are defined only for positive integer indices.")
})

test_that("Union and Intersection take evaluate = FALSE", {
  expect_identical(.fs_fmt(Union(Interval(1, 2), Interval(0, 1), evaluate = FALSE)), "Union(Interval(0, 1), Interval(1, 2))")
  expect_identical(.fs_fmt(Intersection(Interval(1, 2), Interval(0, 1), evaluate = FALSE)), "Intersection(Interval(0, 1), Interval(1, 2))")
})

test_that("Quaternion powers and scalar arithmetic follow SymPy", {
  expect_identical(.fs_fmt(Quaternion(1, 2, 3, 4)^2), "(-28) + 4*i + 6*j + 8*k")
  expect_identical(.fs_fmt(Quaternion(1, 2, 3, 4) * sym_I), "(-2) + 1*i + 4*j + (-3)*k")
  expect_identical(.fs_fmt(Quaternion(1, 2, 3, 4)^-1), "1/30 + (-1/15)*i + (-1/10)*j + (-2/15)*k")
  expect_identical(.fs_fmt(2/Quaternion(1, 2, 3, 4)), "1/15 + (-2/15)*i + (-1/5)*j + (-4/15)*k")
  expect_identical(.fs_fmt(Quaternion(1, 2, 3, 4) + 5), "6 + 2*i + 3*j + 4*k")
})

test_that("Mul and Add order symbols with assumptions as Basic.compare", {
  expect_identical(.fs_fmt(.fs_T(sym_args(x*f(x)*n))), "(x, n, f(x))")
  expect_identical(.fs_fmt(.fs_T(sym_args(p*n*x))), "(x, n, p)")
  expect_identical(.fs_fmt(.fs_T(sym_args(x + n + p))), "(x, n, p)")
})

test_that("sympify parses method calls and lambda", {
  expect_identical(.fs_fmt(sympify("Matrix([[1, 2], [3, 4]]).det()")), "-2")
  expect_identical(.fs_fmt(sympify("lambda x: x + 1")), "Lambda(x, x + 1)")
  expect_identical(.fs_fmt(sympify("(x**2).diff(x)")), "2*x")
  expect_identical(.fs_fmt(sympify("Matrix([[1, 2], [3, 4]]).inv()")), "Matrix([[-2, 1], [3/2, -1/2]])")
})

test_that("has and sym_replace take function types", {
  expect_identical(.fs_fmt(has(x^2 + sin(y), sin)), "True")
  expect_identical(.fs_fmt(has(x^2 + sin(y), cos)), "False")
  expect_identical(.fs_fmt(sym_replace(sin(x) + sin(y), sin, cos)), "cos(x) + cos(y)")
  expect_identical(.fs_fmt(has(erf(x) + 1, erf)), "True")
})

test_that("sqf_list multiplicities are integers", {
  expect_identical(.fs_fmt(.fs_T(sapply(sqf_list(x^5 - x^4 - 2*x^3 + 2*x^2 + x - 1)[[2]], function(p) class(p[[2]])))), "(integer, integer)")
})

test_that("factorint and factorrat keys follow SymPy order", {
  expect_identical(.fs_fmt(.fs_T(names(factorint(-12)))), "(2, 3, -1)")
  expect_identical(.fs_fmt(.fs_T(factorint(-12, multiple = TRUE))), "(-1, 2, 2, 3)")
  expect_identical(.fs_fmt(.fs_T(names(factorrat(Rational(-8, 9))))), "(2, -1, 3)")
})

test_that("cse finds common arguments, negations and negative powers", {
  expect_identical(.fs_fmt(.cse_str(cse(list(x*y*z, x*y*t)))), "([(x0, x*y)], [x0*z, t*x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y + z, x + y + t)))), "([(x0, x + y)], [x0 + z, t + x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x - y)*(z - y) + sqrt((x - y)*(z - y)))))), "([(x0, (x - y)*(-y + z))], [sqrt(x0) + x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x*y*z, x*y*t), symbols = numbered_symbols("t")))), "([(t0, x*y)], [t0*z, t*t0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x^2 + y^2, sqrt(x^2 + y^2), (x^2 + y^2)^2)))), "([(x0, x^2 + y^2)], [x0, sqrt(x0), x0^2])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x+y)^2 + sin(x+y), (x+y)*3)))), "([(x0, x + y)], [x0^2 + sin(x0), 3*x + 3*y])")
  expect_identical(.fs_fmt(.cse_str(cse(list(1/x^2, 1/x^3)))), "([], [x^(-2), x^(-3)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y + z, w + x + y, t + x + y + z)))), "([(x0, x + y), (x1, x0 + z)], [x1, w + x0, t + x1])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x*y*z*w, x*y*z, x*y, w*x*y*z*t)))), "([(x0, x*y), (x1, x0*z), (x2, w*x1)], [x2, x1, x0, t*x2])")
  expect_identical(.fs_fmt(.cse_str(cse(list(-x - y, x + y + 1)))), "([(x0, x + y)], [-x0, x0 + 1])")
  expect_identical(.fs_fmt(.cse_str(cse(list(sin(x+y)+cos(x+y), exp(x+y)), ignore = list(y)))), "([], [sin(x + y) + cos(x + y), exp(x + y)])")
  expect_identical(.fs_fmt(.cse_str(cse(x*y + x*y*z, list = FALSE))), "([(x0, x*y)], x0*z + x0)")
  expect_identical(.fs_fmt(.cse_str(cse(list(exp(-x)*exp(-y), exp(-x))))), "([(x0, exp(-x))], [x0*exp(-y), x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x^-2 + y^-2, 1/(x^2))))), "([(x0, x^(-2))], [x0 + y^(-2), x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x+y, x+y), order = "none"))), "([(x0, x + y)], [x0, x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(2*x*y + 2*x*z, x*y)))), "([(x0, x*y)], [2*x*z + 2*x0, x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x + y)/(x + y + 1), (x+y+1)^2)))), "([(x0, x + y), (x1, x0 + 1)], [x0/x1, x1^2])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x+y)^2 + sqrt(x+y))))), "([(x0, x + y)], [sqrt(x0) + x0^2])")
  expect_identical(.fs_fmt(.cse_str(cse(list((w*x + y)^2 + sqrt(w*x + y))))), "([(x0, w*x + y)], [sqrt(x0) + x0^2])")
  expect_identical(.fs_fmt(.cse_str(cse(list(((w + x + y + z)*(w - y - z))/(w + x)^3)))), "([(x0, w + x), (x1, y + z)], [(w - x1)*(x0 + x1)/x0^3])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x + y)*z, (x + y)*w)))), "([(x0, x + y)], [x0*z, w*x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(w*x*y + z, w*y)))), "([(x0, w*y)], [x*x0 + z, x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(w*x*y, w*x*y + z, w*y)))), "([(x0, w*y), (x1, x*x0)], [x1, x1 + z, x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list((x - z)*(y - z), x - z, y - z)))), "([(x0, -z), (x1, x + x0), (x2, x0 + y)], [x1*x2, x1, x2])")
  expect_identical(.fs_fmt(.cse_str(cse(list(y - z, x - z, (x - z)*(y - z))))), "([(x0, -z), (x1, x0 + y), (x2, x + x0)], [x1, x2, x1*x2])")
  expect_identical(.fs_fmt(.cse_str(cse(list(w*y + w + x + y + z, w*x*y)))), "([(x0, w*y)], [w + x + x0 + y + z, x*x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y, x + y + z)))), "([(x0, x + y)], [x0, x0 + z])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y, x + z)))), "([], [x + y, x + z])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x*y, z + x*y, x*y*z + 3)))), "([(x0, x*y)], [x0, x0 + z, x0*z + 3])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y)))), "([], [x + y])")
  expect_identical(.fs_fmt(.cse_str(cse(list(Integer(42))))), "([], [42])")
  expect_identical(.fs_fmt(.cse_str(cse(list(exp(x)*exp(-x), exp(2*x)*exp(-x))))), "([], [1, exp(x)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(1/x + 1/y, 1/(x*y))))), "([(x0, 1/x), (x1, 1/y)], [x0 + x1, x0*x1])")
  expect_identical(.fs_fmt(.cse_str(cse(list(1/x^2 + x^2)))), "([(x0, x^2)], [x0 + 1/x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x^2 + (1 + 1/x^2)/x^2)))), "([(x0, x^2), (x1, 1/x0)], [x0 + x1*(x1 + 1)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(1/x^2 + (1 + 1/x^2)*x^2)))), "([(x0, x^2), (x1, 1/x0)], [x0*(x1 + 1) + x1])")
  expect_identical(.fs_fmt(.cse_str(cse(list(cos(1/x^2) + sin(1/x^2))))), "([(x0, x^(-2))], [sin(x0) + cos(x0)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(y/(2 + x^2) + z/x^2/y)))), "([(x0, x^2)], [y/(x0 + 2) + z/(x0*y)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(exp(x^2) + x^2*cos(1/x^2))))), "([(x0, x^2)], [x0*cos(1/x0) + exp(x0)])")
  expect_identical(.fs_fmt(.cse_str(cse(list((1 + 1/x^2)/x^2)))), "([(x0, x^(-2))], [x0*(x0 + 1)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x^(2*y) + x^(-2*y))))), "([(x0, x^(2*y))], [x0 + 1/x0])")
  expect_identical(.fs_fmt(.cse_str(cse(list(exp(y)*(3*y + 3*sqrt(x+1)), exp(y)*(5*y + 5*sqrt(x+1)))))), "([(x0, exp(y)), (x1, sqrt(x + 1))], [x0*(3*x1 + 3*y), x0*(5*x1 + 5*y)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(exp(y)*(3*y + 3*sqrt(x+1)), exp(y)*(5*y + 5*sqrt(x+1))), ignore = list(y)))), "([(x0, sqrt(x + 1))], [(3*x0 + 3*y)*exp(y), (5*x0 + 5*y)*exp(y)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(w*exp(x)*exp(-z), exp(y)*exp(x)*exp(-z)), ignore = list(x)))), "([(x0, exp(-z))], [w*x0*exp(x), x0*exp(x)*exp(y)])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x + y, 2 + x + y, x + y + z, 3 + x + y + z)))), "([(x0, x + y), (x1, x0 + z)], [x0, x0 + 2, x1, x1 + 3])")
  expect_identical(.fs_fmt(.cse_str(cse(list(x/(-4*x^2 + y^2))))), "([], [x/(-4*x^2 + y^2)])")
})

test_that("satisfiable lists the model in SymPy's order", {
  # SymPy iterates a set of literals here; the order shown is the one of
  # CPython's set iteration, which may vary with the hash seed
  m <- satisfiable(And(Or(x, y), Not(z), Or(Not(x), y)))
  expect_identical(names(m), c("y", "z", "x"))
  expect_identical(unname(unlist(m)), c(TRUE, FALSE, FALSE))
})

test_that("cse numbers more than ten replacements", {
  r <- cse(lapply(1:12, function(k) sin(x + k*y)^2 + sin(x + k*y)))
  expect_length(r$replacements, 12L)
  expect_identical(sstr(r$replacements[[12]][[1]]), "x11")
  expect_identical(sstr(r$reduced[[12]]), "x11^2 + x11")
})

test_that("nsolve returns a start point that is already a root", {
  expect_identical(sstr(nsolve(list(x^2 - 4), list(x), list(2))$x), "2.00000000000000")
  expect_identical(.fs_err(nsolve(x*y - 1, x, 1)), "ERR ValueError: expected a one-dimensional and numerical function")
})

test_that("parametric integrals start with the degenerate case", {
  r <- integrate(exp(a*x)*cos(b*x), x)
  expect_identical(sstr(r$args[[1L]]), "ExprCondPair(x, Eq(a, 0) & Eq(b, 0))")
})

test_that("cse keeps a matrix as a matrix", {
  r <- cse(list(x + y, Matrix(list(list(x + y, 0), list(0, 0)))))
  expect_identical(sstr(r$replacements[[1]][[2]]), "x + y")
  expect_true(inherits(r$reduced[[2]], "symmatrix"))
  expect_identical(sstr(r$reduced[[2]]), "Matrix([[x0, 0], [0, 0]])")
})
