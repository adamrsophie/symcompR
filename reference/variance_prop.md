# Propagation of variance

Symbolically propagates the variance of the symbols in `expr` (SymPy's
`variance_prop()`); every symbol is treated as a random variable unless
listed in `consts`.

## Usage

``` r
variance_prop(expr, consts = list(), include_covar = FALSE)
```

## Arguments

- expr:

  An expression.

- consts:

  Symbols that are constants.

- include_covar:

  Add covariance terms.

## Value

An expression in
[`Variance()`](https://adamrsophie.github.io/symcompR/reference/symbolic_probability.md)
and
[`Covariance()`](https://adamrsophie.github.io/symcompR/reference/symbolic_probability.md)
terms.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
variance_prop(x * y)
#> x^2*Variance(y) + y^2*Variance(x)
```
