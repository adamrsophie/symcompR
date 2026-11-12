# Package files share one namespace, so a helper defined at top level in two
# files is silently replaced by whichever file is collated last. This check
# reads the sources, so it runs from a source checkout and is skipped on an
# installed package.

test_that("no top-level name is defined in two source files", {
  # An installed package (covr, R CMD check) has an R/ folder without .R files.
  files <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(length(files) == 0L, "package sources not available")
  owner <- list()
  for (f in files) {
    for (e in as.list(parse(f, keep.source = FALSE))) {
      if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") && is.symbol(e[[2L]])) {
        nm <- as.character(e[[2L]])
        owner[[nm]] <- union(owner[[nm]], basename(f))
      }
    }
  }
  shared <- Filter(function(files) length(files) > 1L, owner)
  expect_equal(names(shared), character())
})
