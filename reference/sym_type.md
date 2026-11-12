# Type of a symbolic expression

Returns the head of the expression tree, such as `"Add"`, `"Mul"`,
`"Pow"`, `"Symbol"`, `"Integer"` or `"Function"`. This plays the role of
SymPy's `func` and class checks such as `expr.is_Add`.

## Usage

``` r
sym_type(e)
```

## Arguments

- e:

  A `sym` expression.

## Value

A character scalar.
