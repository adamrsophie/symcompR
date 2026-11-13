# Acceptance target: 84 polynomial checks in 42 unknowns taken from a C
# program. Skipped when the source file is not present.

test_that("the 84 equation key system solves to a unique solution", {
  path <- path.expand("~/Downloads/main.c")
  skip_if_not(file.exists(path))
  skip_on_cran()
  pattern <- "^\\s*if \\((.*) != (-?[0-9]+)\\)$"
  checks <- grep(pattern, readLines(path), value = TRUE)
  expect_length(checks, 84L)
  polys <- sub(pattern, "(\\1) - (\\2)", checks)
  polys <- gsub("key\\[([0-9]+)\\]", "x\\1", polys)
  eqs <- lapply(polys, sympify)
  xs <- lapply(0:41, function(i) Symbol(paste0("x", i)))
  sol <- solve(eqs, xs)
  # SymPy returns a dict for a linear system and a list of solutions otherwise
  if (!is.null(names(sol))) sol <- list(sol)
  expect_length(sol, 1L)
  vals <- lapply(sol[[1L]], as_sym)
  expect_true(all(vapply(vals, function(v) v$type == "Integer", TRUE)))
  for (e in eqs) expect_true(is_zero_num(subs(e, sol[[1L]])))
  key <- intToUtf8(vapply(vals, as.integer, 0L))
  expect_equal(nchar(key), 42L)
})
