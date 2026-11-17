## Resubmission

The incoming pretests of 2026-11-16 reported an overall check time of 28
minutes, then 24 minutes, on r-devel-windows-x86_64.

* All 39 vignettes are now precomputed (the sources are kept as `.Rmd.orig`
  and excluded from the build), so rebuilding them runs no package code
  (189 s to 40 s on Windows).
* The slowest tests now skip on CRAN (174 s to 81 s on Windows); the full
  comparison with SymPy's test suite runs in continuous integration.
* One slow example no longer runs.
* About 13 minutes on Windows is the step "checking R code for possible
  problems". That step scales with the amount of R code (about 187,000
  lines), and it cannot be shortened without removing functionality. We would
  be grateful if the check time could be accepted for this reason.
* Following the review of 2026-11-17, `Authors@R` lists the authors and
  copyright holders of all ported code as contributors and copyright holders.
  The SymPy authors are too many to list individually, so they appear as the
  "SymPy Development Team", whose entry points to `inst/SYMPY-AUTHORS`
  (SymPy's own AUTHORS file at the ported commit). `inst/COPYRIGHTS` refers
  to it as well.

The words flagged as possibly misspelled, "combinatorics" and
"multiprecision", are spelled correctly; they are standard terms in
mathematics and numerical computing.

## Submission

This is the first submission of symcompR, a native R port of the Python
computer algebra system SymPy. It has no Python dependency and no compiled
code; it imports only gmp and base R packages.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

The check also reports an installed size of about 22.6 Mb (R code 16.5 Mb).
The package ports a large library (about 187,000 lines of R covering the
SymPy API), so the R code itself is large; there is no bundled data. Installing
from source takes about three minutes because of the amount of R code.

## Test environments

* macOS (Apple silicon), R release, local.
* GitHub Actions: ubuntu-latest (R devel and release), windows-latest
  (R release).

## Tests and examples

The tests that run on CRAN take about 40 seconds. Slower comparisons with
SymPy's own test suite are skipped on CRAN with `skip_on_cran()` and run in
continuous integration. Examples that take longer than a few seconds are
wrapped in `\donttest{}`.

## Copyright

SymPy, mpmath and the other projects whose code is ported are credited as
copyright holders in `Authors@R`, and their license notices are kept in
`inst/COPYRIGHTS`.
