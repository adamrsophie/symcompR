# Dynamic symbols

Functions of time `t`, as SymPy's `dynamicsymbols()`: `q(t)`, or its
`level`-th time derivative. Names follow
[`sym_symbols()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md):
`"q1 q2"`, `"q1, q2"` or `"q1:3"` give a list.

## Usage

``` r
dynamicsymbols(names, level = 0, ...)
```

## Arguments

- names:

  Names separated by spaces or commas, or a range.

- level:

  Order of the time derivative.

- ...:

  Assumptions such as `positive = TRUE`, given to the functions.

## Value

A `sym` expression or a list of them.

## Examples

``` r
q <- dynamicsymbols("q")
dynamicsymbols("q", 1)
#> Derivative(q(t), t)
dynamicsymbols("q1:4")
#> [[1]]
#> q1(t)
#> 
#> [[2]]
#> q2(t)
#> 
#> [[3]]
#> q3(t)
#> 
```
