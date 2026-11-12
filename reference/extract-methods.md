# Extracting factors and terms

`extract_multiplicatively(e, c)` writes `e` as `c` times something "in a
nice way" (keeping the properties of the arguments) and returns that
something, or `NULL`, as SymPy's `Expr.extract_multiplicatively()`.
`extract_additively(e, c)` does the same for a term.
`as_coefficient(e, expr)` is the coefficient of `expr` in `e` when `e`
is a multiple of `expr`, else `NULL`. `extract_branch_factor(e)` writes
`e` as `exp_polar(2*pi*I*n)*z` and returns `list(z, n)`.
`add_primitive(e)` returns the rational content of a sum and the
primitive remainder, as SymPy's `Add.primitive()`.

## Usage

``` r
add_primitive(e)

extract_multiplicatively(e, c)

extract_additively(e, c)

as_coefficient(e, expr)

extract_branch_factor(e, allow_half = FALSE)
```

## Arguments

- e:

  An expression.

- c, expr:

  The factor or term to extract.

- allow_half:

  Also extract `exp_polar(I*pi)`, as half a branch.

## Value

An expression, `NULL`, or a list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
extract_multiplicatively(2 * x, 2)
#> x
extract_multiplicatively((x * y)^3, x^2 * y)
#> x*y^2
extract_additively(x + 3, 2)
#> x + 1
as_coefficient(2 * sym_I * x, sym_I)
#> 2*x
extract_branch_factor(exp_polar(3 * sym_pi * sym_I + x))
#> [[1]]
#> exp_polar(x + I*pi)
#> 
#> [[2]]
#> 1
#> 
```
