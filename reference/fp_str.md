# Python style string of a group theory object

`fp_str()` returns the string SymPy's
[`str()`](https://rdrr.io/r/utils/str.html) gives for the objects of the
free group, finitely presented group, coset table, homomorphism and
polyhedron functions, with `**` shown as `^`. R lists print as Python
lists, lists made by these functions as tuples where SymPy returns
tuples, `NA` as `None` and `TRUE` as `True`.

## Usage

``` r
fp_str(x)
```

## Arguments

- x:

  An object.

## Value

A character scalar.

## Examples

``` r
G <- free_group("x, y")
x <- G[[2]]; y <- G[[3]]
fp_str(x^2 * y^-1)
#> [1] "x^2*y^-1"
fp_str(word_array_form(x^2 * y^-1))
#> [1] "((x, 2), (y, -1))"
```
