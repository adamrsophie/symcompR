# Whether a leading minus sign can be pulled out of an expression

Mirrors SymPy's `could_extract_minus_sign()`, which decides canonical
forms such as `sin(-x) = -sin(x)`.

## Usage

``` r
could_extract_minus_sign(e)
```

## Arguments

- e:

  A `sym` expression.

## Value

`TRUE` or `FALSE`.
