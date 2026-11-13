# Create symbols

`Symbol()` creates one symbol. `sym_symbols()` splits a string on spaces
or commas and returns one symbol or a list of them; ranges such as
`"x0:3"` expand to `x0, x1, x2`. `sym_var()` does the same and assigns
the symbols in the calling environment. Assumptions such as
`positive = TRUE` or `integer = TRUE` are recorded on the symbol and
drive simplification, as in SymPy. Symbols with the same name and
different assumptions are distinct.

## Usage

``` r
Symbol(name, ...)

sym_var(names, ..., envir = parent.frame())

sym_symbols(names, ..., cls = NULL)
```

## Arguments

- name:

  Symbol name.

- ...:

  Assumptions given as `name = TRUE` or `name = FALSE`. Valid names are
  `commutative`, `complex`, `real`, `rational`, `irrational`, `integer`,
  `noninteger`, `even`, `odd`, `prime`, `composite`, `positive`,
  `negative`, `zero`, `nonzero`, `nonnegative`, `nonpositive`, `finite`,
  `infinite`, `imaginary`, `algebraic`, `transcendental`,
  `extended_real`, `hermitian`, `antihermitian` and `polar`, and the
  Fermi level assumptions `above_fermi` and `below_fermi` used by
  [`wicks()`](https://adamrsophie.github.io/symcompR/reference/secondquant-wicks.md).

- names:

  A string of names separated by spaces or commas.

- envir:

  Environment in which `sym_var()` assigns.

- cls:

  For `sym_symbols()`, the constructor to use instead of `Symbol()`,
  such as `Dummy` or `"Dummy"`.

## Value

A `sym` symbol, or a list of symbols.

## Examples

``` r
x <- Symbol("x")
xyz <- sym_symbols("x y z")
n <- Symbol("n", integer = TRUE, positive = TRUE)
```
