# Simplify a noncommutative expression

`nc_simplify(expr)` simplifies a noncommutative product by collecting
repeated subterms into powers and inverses, `a*b*a*b*c` becoming
`(a*b)^2*c`, as SymPy's `nc_simplify()`. It works on products of
noncommutative symbols and on matrix expressions.

## Usage

``` r
nc_simplify(expr, deep = TRUE)
```

## Arguments

- expr:

  An expression.

- deep:

  Also simplify the arguments of the products found.

## Value

An expression.

## Examples

``` r
a <- Symbol("a", commutative = FALSE); b <- Symbol("b", commutative = FALSE)
nc_simplify(a * b * a * b * a * b)
#> (a*b)^3
```
