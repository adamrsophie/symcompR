skip_unless_full_tests()
# N() against SymPy's evalf for the rows of sympy/core/tests/test_evalf.py,
# the evalf rows of sympy/functions/*/tests and the N and evalf rows of
# sympy/polys/tests/test_rootoftools.py. Each row of
# evalf-core-sympy-battery.tsv holds the R call and SymPy's str() of the
# result with ** shown as ^; errors compare as "ERR <class>: <message>".

rows <- utils::read.delim(test_path("evalf-core-sympy-battery.tsv"), header = FALSE, quote = "",
                          colClasses = "character", col.names = c("code", "want"))
err_string <- function(err) {
  cl <- setdiff(class(err), c("simpleError", "error", "condition", "sympy_error", "mpmath_error"))
  paste0("ERR ", if (length(cl)) cl[1L] else "Error", ": ", conditionMessage(err))
}

test_that("N() reproduces SymPy's evalf battery", {
  for (i in seq_len(nrow(rows))) {
    got <- tryCatch({
      v <- suppressWarnings(eval(parse(text = rows$code[i]), envir = new.env(parent = globalenv())))
      gsub("**", "^", sstr(v), fixed = TRUE)
    }, error = err_string)
    expect_identical(gsub("\n", "\\n", got, fixed = TRUE), rows$want[i], label = rows$code[i])
  }
})
