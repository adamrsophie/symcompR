# Expected values are SymPy's outputs for the same calls (Parsing (latex2)).

skip_unless_full_tests()

skip_on_cran()

.bt <- new.env(parent = asNamespace("symcompR"))
eval(parse(text = r"---(x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
.fmt <- function(v) {
  f <- get0(".pe_str", envir = asNamespace("symcompR"))
  if (!is.null(f)) return(f(v))
  if (is.null(v)) "None" else if (inherits(v, "lark_tree")) format(v) else sstr(v)
}
.norm_err <- function(cls, msg) {
  lines <- strsplit(msg, "\n", fixed = TRUE)[[1L]]
  if (endsWith(msg, "\n")) lines <- c(lines, "")
  if (cls == "UnexpectedCharacters") {
    b <- startsWith(lines, "\t* ")
    rest <- lines[!b]
    lines <- c(rest[-length(rest)], sort(lines[b], method = "radix"), rest[length(rest)])
  }
  paste0("ERR ", cls, ": ", paste(lines, collapse = "\\n"))
}
.fmt1 <- function(v) paste(.fmt(v), collapse = "\n")
.err_msg <- function(e) {
  cls <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  m <- .norm_err(if (length(cls)) cls[1L] else "Error", conditionMessage(e))
  gsub("\\n", "\n", sub("^ERR [^:]*: ", "", m), fixed = TRUE)
}
)---"), envir = .bt)
.run <- function(code) {
  e <- new.env(parent = .bt)
  NULL
  get(".fmt1", envir = .bt)(eval(parse(text = code), envir = e))
}
.setup <- function(code) invisible(eval(parse(text = code), envir = .bt))
.err <- function(code) tryCatch({ .run(code); "no error" }, error = function(e) get(".err_msg", envir = .bt)(e))

test_that("Parsing (latex2) match SymPy", {
  expect_equal(.run("parse_latex_lark_all(\"x + y\")"), "[x + y]")
  expect_equal(.run("parse_latex_lark_all(\"2x \\\\cdot y\")"), "[2*x*y]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xy\")"), "[sin(x*y), y*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"f^{-1}(x)\")"), "[f^{-1}(x), x/f]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\cos ab + \\\\sin xy\")"), "[sin(x*y) + cos(a*b), b*cos(a) + y*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\cos ab + \\\\sin xyz\")"), "[sin(x*y*z) + cos(a*b), b*cos(a) + y*z*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xy \\\\cos ab\")"), "[sin(x*y)*cos(a*b), b*y*sin(x)*cos(a)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\ln xy + \\\\sin^2 ab\")"), "[log(x*y) + sin(a*b)^2, b*sin(a)^2 + y*log(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xy = \\\\cos ab\")"), "[Eq(sin(x*y), cos(a*b)), Eq(y*sin(x), b*cos(a))]")
  expect_equal(.run("parse_latex_lark_all(\"f^{-1}(x) + g^{-1}(y)\")"), "[f^{-1}(x) + g^{-1}(y), y/g + x/f]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\tan hk\")"), "[tan(h*k), k*tan(h)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xyz\")"), "[sin(x*y*z), y*z*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin x y z\")"), "[sin(x*y*z), y*z*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xyza\")"), "[sin(a*x*y*z), a*y*z*sin(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin^2 xyz\")"), "[sin(x*y*z)^2, y*z*sin(x)^2]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xyz \\\\cos t\")"), "[sin(x*y*z)*cos(t), y*z*sin(x)*cos(t)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\ln xy\")"), "[log(x*y), y*log(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\sin xy + f^{-1}(x)\")"), "[f^{-1}(x) + sin(x*y), sin(x*y) + x/f, y*sin(x) + f^{-1}(x), y*sin(x) + x/f]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\operatorname{f} x y\")"), "[f(x*y), y*f(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\int_{f(a)}^{f(b)} f(z) dz\")"), "[Integral(f(z), (z, f(a), f(b))), Integral(f*z, (z, a*f, b*f))]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\frac{d f(x)}{dx}\")"), "[Derivative(f(x), x), Derivative(f*x, x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\frac{d\\\\theta(x)}{dx}\")"), "[Derivative(theta(x), x), Derivative(theta*x, x)]")
  expect_equal(.run("parse_latex_lark_all(\"f(x)\")"), "[f(x), f*x]")
  expect_equal(.run("parse_latex_lark_all(\"f'_1(x)\")"), "[f'_{1}(x), f'_{1}*x]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\ln x y\")"), "[log(x*y), y*log(x)]")
  expect_equal(.run("parse_latex_lark_all(\"\\\\log xy\")"), "[log(x*y), y*log(x)]")
  expect_equal(.run("parse_latex_lark(\"2\\\\pi r\", locals = list(\"\\\\pi\" = Symbol(\"pi\")))"), "2*pi*r")
  expect_equal(.run("parse_latex_lark(\"\\\\mathrm{e}^{\\\\pi}\", locals = list(\"\\\\mathrm{e}\" = Symbol(\"e\")))"), "e^pi")
  expect_equal(.run("parse_latex_lark(\"\\\\imaginaryunit x\", locals = list(\"\\\\imaginaryunit\" = Symbol(\"j\")))"), "j*x")
  expect_equal(.run("parse_latex_lark(\"\\\\sin x + \\\\cos x\", locals = list(\"\\\\sin\" = Function(\"g\")))"), "g(x) + cos(x)")
  expect_equal(.run("parse_latex_lark(\"\\\\sin^2 x\", locals = list(\"\\\\sin\" = Function(\"g\")))"), "g(x)^2")
  expect_equal(.run("parse_latex_lark(\"\\\\log_2 x\", locals = list(\"\\\\log\" = Function(\"g\")))"), "g(x, 2)")
  expect_equal(.run("parse_latex_lark(\"\\\\operatorname{erf}(x)\", locals = list(\"\\\\operatorname{erf}\" = Function(\"g\")))"), "g(x)")
  expect_equal(.run("parse_latex_lark(\"\\\\operatorname{sin}(x, y)\", locals = list(\"\\\\operatorname{sin}\" = Function(\"g\")))"), "g(x, y)")
  expect_equal(.run("parse_latex_lark(\"\\\\operatorname{sin} x + \\\\sin x\", locals = list(\"\\\\operatorname{sin}\" = Function(\"g\")))"), "g(x) + sin(x)")
  expect_equal(.err("parse_latex_lark(\"\\\\sqrt{x}\", locals = list(\"\\\\sqrt\" = Function(\"g\")))"), "The reading of \\sqrt cannot be overridden.")
  expect_equal(.run("parse_latex_lark(\"2\\\\pi r\", locals = list())"), "2*pi*r")
  expect_equal(.run("parse_latex_lark_all(\"\\\\ln x y\", locals = list(\"\\\\ln\" = Function(\"g\")))"), "[g(x*y), y*g(x)]")
  expect_equal(.err("parse_latex(\"\\\\begin{cases} x=1 \\\\\\\\ y=2 \\\\end{cases}\", backend = \"lark\")"), "The 'cases' environment is not currently supported by parse_latex. \nIf you're trying to parse a piecewise function or a block of expressions, \nconsider breaking them into separate parse_latex calls.")
  expect_equal(.err("parse_latex(\"\\\\begin{pmatrix}a & b\", backend = \"lark\")"), "Extra '\\begin{pmatrix}' at index 0 or missing corresponding '\\end{pmatrix}' in LaTeX string: \\begin{pmatrix}a & b")
  expect_equal(.err("parse_latex(\"a \\\\end{pmatrix}\", backend = \"lark\")"), "Extra '\\end{pmatrix}' at index 2 or missing corresponding '\\begin{pmatrix}' in LaTeX string: a \\end{pmatrix}")
  expect_equal(.err("parse_latex(\"\\\\left(\\\\begin{matrix}a\\\\end{matrix}\\\\right]\", backend = \"lark\")"), "Expected '\\end{matrix}\\right)' to close the '\\left(\\begin{matrix}' at index 0 but found '\\end{matrix}\\right]' at index 21 of LaTeX string instead: \\left(\\begin{matrix}a\\end{matrix}\\right]")
  expect_equal(.err("parse_latex(\"x\", backend = \"other\")"), "Using the 'other' backend in the LaTeX parser is not supported, backend must be one of ('antlr', 'lark')")
})
