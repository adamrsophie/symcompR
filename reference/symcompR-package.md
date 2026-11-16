# symcompR: Symbolic Computation in R

A native R computer algebra system ported from the Python library SymPy,
using 'gmp' for exact integers and rationals. The SymPy Development Team
is gratefully acknowledged.

## Details

Attaching the package masks no function of base R. SymPy names that
clash with base R carry a `sym_` prefix, such as
[`sym_factor()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md),
[`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md),
[`sym_symbols()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md),
[`sym_var()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md)
and
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md);
the constants are `sym_pi`, `sym_E` and `sym_I`.

## See also

Useful links:

- <https://github.com/adamrsophie/symcompR>

- <https://adamrsophie.github.io/symcompR/>

- Report bugs at <https://github.com/adamrsophie/symcompR/issues>

## Author

**Maintainer**: Adam R Sophie <adamrsophie@gmail.com> \[copyright
holder\]

Authors:

- Adam R Sophie <adamrsophie@gmail.com> \[copyright holder\]
