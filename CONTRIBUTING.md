# Contributing to symcompR

Thank you for considering a contribution to `symcompR`. The package
ports SymPy to R, so a change is judged on two things: whether the R
code is sound, and whether it gives the answer SymPy gives.

## Before You Start

- Use GitHub issues for bug reports, feature requests, documentation
  gaps, and questions about porting a SymPy module.
- For larger changes, such as porting a new module, open an issue first
  so the design can be discussed before implementation.
- Keep pull requests focused. A new module, a bug fix, and a refactor
  should usually be separate pull requests.

## Reporting a Wrong Answer

A wrong result is the most valuable bug report. Please include:

- the R input that produces it,
- what symcompR returns,
- what SymPy returns for the same input, if you can run it,
- [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

## Development Setup

Install the development dependencies from the repository root:

``` r

install.packages(c("remotes", "testthat", "devtools", "roxygen2"))
remotes::install_deps(dependencies = TRUE)
```

symcompR has no compiled code, so no C or C++ toolchain is needed.

## Porting Rules

- Follow the SymPy source at the pinned commit recorded in
  `inst/COPYRIGHTS`. Port the algorithm, not just the output of one
  example.
- Take expected values in tests from SymPy itself, run on the same
  input. Do not write an expected value by hand when SymPy can produce
  it.
- Attaching the package must not mask any base R function. A SymPy name
  that clashes with base R gets a `sym_` prefix, or an S3 method when
  base R has a generic. The test in `tests/testthat/test-aliases.R`
  enforces this.
- Keep R code ASCII only; write other characters as `\u` escapes.
- Use exact arithmetic (`gmp` big integers and rationals) wherever SymPy
  is exact. Never round through doubles in an exact code path.
- Keep roxygen tags directly above the function they document; a helper
  placed between them takes over the export and the help page.

## Branches and Pull Requests

- Use one branch per logical change, named after the change
  (`port-fu-trigsimp`, `fix-range-intersection`).
- Write commit subjects in the imperative (“Port the Fu trigonometric
  simplification routines”).
- Keep the implementation, its tests, the roxygen source, the
  regenerated documentation, and the `NEWS.md` bullet together in the
  same pull request.

## Code Style

- Follow the style already used in the surrounding R files.
- Prefer clear validation and explicit errors over silent coercion.
- Do not add new package dependencies unless they are necessary and
  justified.
- Regenerate the `.Rd` files with `devtools::document()` after changing
  roxygen comments.

## Testing

Run the tests while developing:

``` r

devtools::test()
```

Before a pull request, run the full check:

``` sh
R CMD build .
R CMD check --as-cran symcompR_*.tar.gz
```

## Pull Request Checklist

Before opening a pull request, confirm:

- Tests relevant to the change pass, with expected values from SymPy.
- No base R function is masked.
- Documentation is updated.
- `NEWS.md` has an entry for user-facing changes.
- No generated build artifacts, local caches, `.Rcheck` directories, or
  tarballs are committed.

## Code of Conduct

All contributors are expected to follow the project [Code of
Conduct](https://adamrsophie.github.io/symcompR/CODE_OF_CONDUCT.md).
