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

# S3 method for class 'symmatrix'
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

  Printer settings of SymPy's `PrettyPrinter` for `pretty()` of an
  expression: `order` (`"lex"`, `"rev-lex"` or `"none"`), `full_prec`,
  `use_unicode_sqrt_char`, `root_notation`, `mat_symbol_style`
  (`"plain"` or `"bold"`), `imaginary_unit` and `perm_cyclic`; unused
  otherwise.

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
cat(pretty(Symbol("theta") + Symbol("x1"), use_unicode = TRUE), "\n")
#> θ + x₁ 
cat(pretty(Integral(x^2, list(x, 0, 1)), use_unicode = FALSE), "\n")
#>   1      
#>   /      
#>  |       
#>  |   2   
#>  |  x  dx
#>  |       
#> /        
#> 0         
```
