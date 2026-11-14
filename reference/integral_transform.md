# Change of variables in an integral

`integral_transform(I, x, u)` performs the substitution `x = u` in the
integral `I`, as SymPy's `Integral.transform()`: when `x` is a symbol,
`u` is an expression in a new variable (give `list(u, uvar)` when it has
several free symbols); when `u` is a symbol, `x` is an expression in the
integration variable that becomes `u`. Limits are mapped, and swapped
with a sign change when they come out in decreasing order.

## Usage

``` r
integral_transform(I, x, u)
```

## Arguments

- I:

  An unevaluated
  [`Integral()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md).

- x:

  An integration variable, or an expression in one.

- u:

  The expression replacing it, or `list(expr, symbol)`.

## Value

An `Integral`.

## Examples

``` r
x <- Symbol("x"); u <- Symbol("u")
integral_transform(Integral(x, list(x, 0, 1)), x, 2*u)
#> Integral(4*u, (u, 0, 1/2))
integral_transform(Integral(x * exp(x^2), list(x, 0, 1)), x^2, u)
#> Integral(exp(u)/2, (u, 0, 1))
```
