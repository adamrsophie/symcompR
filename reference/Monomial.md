# Monomials

`Monomial(exponents, gens)` is a monomial given by its exponent vector,
as SymPy's `Monomial`. `*`, `/` and `^` combine monomials, and
`as_expr(m, ...)` writes it in the generators.

## Usage

``` r
Monomial(exponents, gens = NULL)

# S3 method for class 'sym_monomial'
format(x, ...)

# S3 method for class 'sym_monomial'
print(x, ...)

# S3 method for class 'sym_monomial'
Ops(e1, e2)
```

## Arguments

- exponents:

  Nonnegative integer exponents.

- gens:

  Optional generators.

- x, e1, e2:

  Monomials.

- ...:

  Generators for
  [`as_expr()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

## Value

A monomial, or an expression from
[`as_expr()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

## Examples

``` r
m <- Monomial(c(1, 2))
m
#> Monomial((1, 2))
m * Monomial(c(2, 1))
#> Monomial((3, 3))
```
