# Methods of Poly objects

`p$cancel(g, include = FALSE)` cancels the common factors of the
polynomials `p` and `g`, as SymPy's `Poly.cancel()`: it returns the
constant factor and the two reduced polynomials, or with
`include = TRUE` the two polynomials with the constant multiplied in.
`p$factor_list()` is SymPy's `Poly.factor_list()`: the constant and a
list of (factor, multiplicity) pairs with the factors as polynomials.
`p$domain` is the coefficient domain. Other names give the components of
the object.

## Usage

``` r
# S3 method for class 'Poly'
x$name
```

## Arguments

- x:

  A polynomial made by
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

- name:

  A method or component name.

## Value

The method as a function, or the component.

## Examples

``` r
x <- Symbol("x")
Poly(4*x^2 - 4, x)$cancel(Poly(2*x - 2, x))
#> (1, 2*x + 2, 1) 
```
