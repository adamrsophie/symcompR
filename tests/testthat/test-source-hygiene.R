# Package files share one namespace, so a helper defined at top level in two
# files is silently replaced by whichever file is collated last. This check
# reads the sources, so it runs from a source checkout and is skipped on an
# installed package.

test_that("no top-level name is defined in two source files", {
  rdir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "package sources not available")
  owner <- list()
  for (f in list.files(rdir, pattern = "[.]R$", full.names = TRUE)) {
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
