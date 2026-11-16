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

Other contributors:

- SymPy Development Team (Author of 'SymPy', from which the algorithms
  are ported) \[copyright holder\]

- Fredrik Johansson (Author of 'mpmath', from which the multiprecision
  code is ported) \[copyright holder\]

- mpmath contributors \[copyright holder\]

- Erez Shinan (Author of 'lark', from which the Earley parser is ported)
  \[copyright holder\]

- Matthew Rocklin (Author of 'multipledispatch', from which the dispatch
  code is ported) \[copyright holder\]

- latex2sympy authors (Authors of 'latex2sympy', from which SymPy's
  LaTeX parser derives) \[copyright holder\]

- Sergey B Kirpichev (Author of patches from 'Diofant' included in
  'SymPy') \[copyright holder\]
