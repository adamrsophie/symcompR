# Test tiers. The SymPy comparison suites take more than an hour on one
# core, so a file that calls skip_unless_full_tests() at its top runs only
# when the environment variable SYMCOMPR_FULL_TESTS is "true": the weekly
# full check (.github/workflows/full-check.yaml) sets it, and so can a local
# run, for example Sys.setenv(SYMCOMPR_FULL_TESTS = "true") before
# devtools::test(). Pushes, pull requests and CRAN run the fast tier.
skip_unless_full_tests <- function() {
  if (!identical(Sys.getenv("SYMCOMPR_FULL_TESTS"), "true"))
    testthat::skip("full test tier; set SYMCOMPR_FULL_TESTS=true to run it")
}
