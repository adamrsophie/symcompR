# Convert R values to symbolic expressions

`as_sym()` (alias `S()`) converts numbers, strings, logical values,
`gmp` big numbers, R complex numbers, formulas and calls into `sym`
expressions. Integer valued doubles become exact integers and other
doubles become
[`Float()`](https://adamrsophie.github.io/symcompR/reference/Integer.md)
values. Strings are parsed with
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md).

## Usage

``` r
as_sym(x)

S(x)
```

## Arguments

- x:

  Value to convert.

## Value

A `sym` expression.

## Examples

``` r
S(2) / 3
#> 2/3
S("x^2 + 1")
#> x^2 + 1
```
