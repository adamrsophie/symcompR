# Common subexpression elimination

`cse(exprs)` finds common subexpressions as SymPy's `cse()`: repeated
subtrees, negated terms, negative powers, and the common arguments of
sums and products. It returns a list with `replacements` (a list of
`list(symbol, expression)` pairs) and `reduced` (the rewritten
expressions).

## Usage

``` r
cse(
  exprs,
  prefix = "x",
  symbols = NULL,
  optimizations = NULL,
  postprocess = NULL,
  order = "canonical",
  ignore = list(),
  list = TRUE
)
```

## Arguments

- exprs:

  An expression, a matrix, or a list of them.

- prefix:

  Prefix of the new symbols when `symbols` is `NULL`.

- symbols:

  A list of symbols to use for the replacements, such as
  `numbered_symbols("t")`, or `NULL` for `prefix` followed by a number.

- optimizations:

  A list of `list(pre, post)` pairs of functions applied to each
  expression before and after the elimination, or `NULL`.

- postprocess:

  A function of the replacements and the reduced expressions whose value
  is returned, or `NULL`.

- order:

  `"canonical"` to visit the arguments of sums and products in sorted
  order, or `"none"`.

- ignore:

  Symbols; subexpressions containing any of them are not replaced.

- list:

  If `FALSE`, a single expression gives a single reduced expression
  rather than a list.

## Value

A list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
cse(list((x + y)^2 + sin(x + y), (x + y) * 3))
#> $replacements
#> $replacements[[1]]
#> $replacements[[1]][[1]]
#> x0
#> 
#> $replacements[[1]][[2]]
#> x + y
#> 
#> 
#> 
#> $reduced
#> $reduced[[1]]
#> x0^2 + sin(x0)
#> 
#> $reduced[[2]]
#> 3*x + 3*y
#> 
#> 
cse(list(x*y*z, x*y*t))
#> $replacements
#> $replacements[[1]]
#> $replacements[[1]][[1]]
#> x0
#> 
#> $replacements[[1]][[2]]
#> x*y
#> 
#> 
#> 
#> $reduced
#> $reduced[[1]]
#> x0*z
#> 
#> $reduced[[2]]
#> t*x0
#> 
#> 
```
