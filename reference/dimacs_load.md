# Read DIMACS CNF

`dimacs_load(s)` reads a boolean expression in the DIMACS CNF format
(variable `n` becomes the symbol `cnf_n`), and
`dimacs_load_file(location)` reads it from a file, as
[`load()`](https://rdrr.io/r/base/load.html) and `load_file()` in
`sympy.logic.utilities.dimacs`.

## Usage

``` r
dimacs_load(s)

dimacs_load_file(location)
```

## Arguments

- s:

  A string.

- location:

  A file path.

## Value

A boolean expression.

## Examples

``` r
dimacs_load("1 2 \n 3")
#> cnf_3 & (cnf_1 | cnf_2)
```
