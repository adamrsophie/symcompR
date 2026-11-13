# C math macros

`get_math_macros()` returns the constants that the C printers print as
the macros of `math.h`, as SymPy's `sympy.printing.c.get_math_macros`: a
dictionary from expressions to macro names, which `py_getitem()`-style
lookups read with `cgdict_get()`.

## Usage

``` r
get_math_macros()

cgdict_get(d, key)
```

## Arguments

- d:

  A dictionary from `get_math_macros()`.

- key:

  An expression.

## Value

`get_math_macros()` returns a dictionary; `cgdict_get()` the value for
`key`, or `NULL`.

## Examples

``` r
m <- get_math_macros()
cgdict_get(m, exp(Integer(1)))
#> [1] "M_E"
```
