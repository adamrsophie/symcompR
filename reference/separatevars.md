# Separate variables

`separatevars(e)` factors `e` into a product of factors that each depend
on fewer symbols, when possible.

## Usage

``` r
separatevars(e, syms = NULL, dict = FALSE)
```

## Arguments

- e:

  An expression.

- syms:

  Optional symbols; with `dict = TRUE` returns a named list.

- dict:

  Return a list of factors per symbol.

## Value

An expression or a list.
