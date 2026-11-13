# Dynamic symbols of an expression

`find_dynamicsymbols()` returns the functions of time and their
derivatives found in `expression`, as SymPy's `find_dynamicsymbols()`. A
vector needs the frame in which to express it.

## Usage

``` r
find_dynamicsymbols(expression, exclude = NULL, reference_frame = NULL)
```

## Arguments

- expression:

  An expression, matrix or vector.

- exclude:

  Dynamic symbols to leave out.

- reference_frame:

  The
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
  for a vector expression.

## Value

A list of expressions (a set in SymPy).

## Examples

``` r
q <- dynamicsymbols("q")
find_dynamicsymbols(q^2 + dynamicsymbols("q", 1))
#> [[1]]
#> Derivative(q(t), t)
#> 
#> [[2]]
#> q(t)
#> 
#> attr(,"pyset")
#> [1] TRUE
```
