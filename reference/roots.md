# Roots of a polynomial

`roots(f, x)` returns the roots of a univariate polynomial, solved by
factorization over the rationals and by radical formulas up to degree
four (and for binomials), like SymPy's `roots()`. The result is a list
of distinct roots with a `multiplicity` attribute; with
`multiple = TRUE` each root is repeated by its multiplicity. Factors
that are not solvable by these formulas are left out unless
`rootof = TRUE`, in which case they appear as
[`CRootOf()`](https://adamrsophie.github.io/symcompR/reference/CRootOf.md)
objects.

## Usage

``` r
roots(
  f,
  x = NULL,
  multiple = FALSE,
  rootof = FALSE,
  filter = NULL,
  cubics = TRUE,
  trig = FALSE,
  quartics = TRUE,
  quintics = FALSE,
  predicate = NULL,
  strict = FALSE,
  auto = TRUE
)
```

## Arguments

- f:

  A polynomial expression or
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

- x:

  The variable; when `k` is missing, `x` is the index and the variable
  is the single free symbol of `f`.

- multiple:

  Repeat roots by multiplicity.

- rootof:

  Include `CRootOf` roots for unsolvable factors.

- filter:

  Keep only `"real"`, `"complex"`, `"rational"`, or all roots (`NULL`).

- cubics, quartics, quintics:

  Use the formulas for degrees three to five.

- trig:

  Use the trigonometric form for cubics with three real roots.

- predicate:

  A function selecting roots.

- strict:

  Signal an error when some roots are not found.

- auto:

  Convert ring coefficients to their field.

## Value

A list of expressions.

## Examples

``` r
x <- Symbol("x")
roots(x^2 - 3*x + 2, x)
#> [[1]]
#> 2
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 1
roots(x^3 - 1, x)
#> [[1]]
#> 1
#> 
#> [[2]]
#> -1/2 - sqrt(3)*I/2
#> 
#> [[3]]
#> -1/2 + sqrt(3)*I/2
#> 
#> attr(,"multiplicity")
#> [1] 1 1 1
roots((x - 1)^2 * (x + 3), x)
#> [[1]]
#> -3
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 2
```
