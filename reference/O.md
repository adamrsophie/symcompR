# Order terms

`O(expr, ...)` (alias `Order()`) is SymPy's big O term. The extra
arguments are the variables (each a symbol, or `list(x, x0)` for an
expansion point other than 0); by default the free symbols of `expr` at
0. Added to an expression it absorbs the terms it contains; `removeO()`
drops it and `getO()` returns it (or `NULL`).

## Usage

``` r
O(expr, ...)

Order(expr, ...)

removeO(e)

getO(e)
```

## Arguments

- expr:

  An expression.

- ...:

  Variables, or `list(variable, point)` pairs.

- e:

  An expression with an order term.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
O(x^3)
#> O(x^3)
x + x^2 + O(x^2)
#> x + O(x^2)
O(1/x, list(x, Inf))
#> O(1/x, (x, oo))
```
