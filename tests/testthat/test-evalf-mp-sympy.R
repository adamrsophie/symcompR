skip_unless_full_tests()
# N() against SymPy's str(N(expr, n)) for the functions evaluated through the
# mpmath port: hypergeometric and Meijer G functions, the review examples and
# the special functions (Bessel, Airy, error functions, exponential
# integrals, incomplete gamma, elliptic integrals, orthogonal polynomials,
# zeta, polylog, Dirichlet eta, Lambert W) at 15 and 30 digits. Each row of
# evalf-mp-sympy-battery.tsv holds the R expression, the digits and SymPy's
# output; errors compare as "ERR <class>: <message>".

P <- function(s) parse_expr(s)
rows <- utils::read.delim(test_path("evalf-mp-sympy-battery.tsv"), header = FALSE, quote = "",
                          colClasses = "character", col.names = c("code", "digits", "want"))
err_string <- function(err) {
  cl <- setdiff(class(err), c("simpleError", "error", "condition", "mpmath_error", "sympy_error"))
  paste0("ERR ", if (length(cl)) cl[1] else "Error", ": ", conditionMessage(err))
}

test_that("N() reproduces SymPy's output for the mpmath battery", {
  for (i in seq_len(nrow(rows))) {
    got <- tryCatch(suppressWarnings(sstr(N(eval(parse(text = rows$code[i])), as.integer(rows$digits[i])))),
                    error = err_string)
    expect_identical(got, rows$want[i], label = paste0(rows$code[i], " at ", rows$digits[i], " digits"))
  }
})
