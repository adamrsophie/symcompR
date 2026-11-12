# Changelog

## symcompR 0.1.0

- First version: a native R port of SymPy (pinned to commit `86150a81`),
  with no Python dependency. It covers 867 of the 901 names in SymPy’s
  top level API: core arithmetic and assumptions, elementary and special
  functions, calculus, integration and integral transforms, polynomials,
  algebraic number theory and Galois groups, simplification, solvers
  (algebraic, inequalities, ODE, PDE, recurrences, Diophantine),
  matrices and matrix expressions, number theory, combinatorics and
  group theory, logic, sets, geometry, statistics, sequences and Fourier
  series, vectors, tensors, physical units, parsers (SymPy, LaTeX and
  Mathematica syntax), and printers (LaTeX, MathML, two dimensional
  text, and code in thirteen languages).
- Attaching the package masks no base R function. SymPy names that clash
  with base R carry a `sym_` prefix
  ([`sym_factor()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md),
  [`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md),
  [`sym_prod()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md),
  [`sym_reshape()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md),
  …) or are S3 methods of the base generic.
- Machine readable citation and software metadata: `CITATION.cff`,
  `codemeta.json`, `.zenodo.json`, `ro-crate-metadata.json`, and
  `inst/CITATION`.
