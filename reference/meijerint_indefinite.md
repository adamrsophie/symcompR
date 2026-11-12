# Integration by Meijer G functions

Ports of SymPy's `meijerint_indefinite`, `meijerint_definite` and
`meijerint_inversion` (sympy/integrals/meijerint.py). The integrand is
rewritten as a product of Meijer G functions, which are integrated with
the formulas of Prudnikov, Brychkov and Marichev and then expanded back
to elementary and special functions with
[`hyperexpand()`](https://adamrsophie.github.io/symcompR/reference/hyperexpand.md).

## Usage

``` r
meijerint_indefinite(f, x)

meijerint_definite(f, x, a, b)

meijerint_inversion(f, x, t)
```

## Arguments

- f:

  The integrand.

- x:

  The integration variable.

- a, b:

  The limits.

- t:

  The variable of the inverse Laplace transform.

## Value

An expression, a list `(value, condition)`, or `NULL`.

## Details

`meijerint_indefinite(f, x)` returns an antiderivative (possibly a
Piecewise) or `NULL`. `meijerint_definite(f, x, a, b)` returns a list of
the value and its convergence condition, or `NULL`.
`meijerint_inversion(f, x, t)` computes the inverse Laplace transform of
`f`, or `NULL`.

## Examples

``` r
x <- Symbol("x")
meijerint_indefinite(sin(x), x)
#> -cos(x)
meijerint_definite(exp(-x^2), x, 0, oo)
#> [[1]]
#> sqrt(pi)/2
#> 
#> [[2]]
#> True
#> 
```
