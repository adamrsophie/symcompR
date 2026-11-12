# Pretty printing

`pretty(e)` lays an expression out in two dimensions, as SymPy's
`pretty()`: raised exponents, fractions over a bar, root signs and tall
parentheses, in Unicode or ASCII. `pprint()` and `pretty_print()` print
it. `pprint_use_unicode(flag)` sets the default and returns the previous
value; `pprint_try_use_unicode()` turns Unicode on when the session
encoding is UTF-8. `pager_print()` shows the result in the pager.

## Usage

``` r
# S3 method for class 'sym'
pretty(x, use_unicode = NULL, ...)

pprint(e, use_unicode = NULL)

pretty_print(e, use_unicode = NULL)

pprint_use_unicode(flag = NULL)

pprint_try_use_unicode()

pager_print(e, use_unicode = NULL)
```

## Arguments

- use_unicode:

  Unicode or ASCII; the session default when `NULL`.

- ...:

  Unused.

- e, x:

  An expression.

- flag:

  `TRUE`, `FALSE`, or `NULL` to only query.

## Value

A character scalar; the printing functions return `NULL` invisibly.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
cat(pretty((x + 1)/(y - 2), use_unicode = FALSE), "\n")
#> x + 1
#> -----
#> y - 2 
pprint(sqrt(x^2 + 1), use_unicode = FALSE)
#>    ________
#>   /  2     
#> \/  x  + 1 
```
