# Convert an expression to a string

`sstr()` gives the SymPy style string form with R operators, so the
result can be parsed again with
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md).
`srepr()` shows the expression tree.

## Usage

``` r
sstr(e)

srepr(e)
```

## Arguments

- e:

  A `sym` expression.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x")
sstr(x^2 + 2*x + 1)
#> [1] "x^2 + 2*x + 1"
srepr(x^2 + 1)
#> [1] "Add(Integer(1), Pow(Symbol('x'), Integer(2)))"
```
