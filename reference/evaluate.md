# Turn evaluation off

SymPy's `with evaluate(False):` context. `evaluate(FALSE, expr)`
evaluates the R expression `expr` with automatic simplification off: the
arithmetic operators build unevaluated `Add`, `Mul` and `Pow` nodes,
functions stay unevaluated (`sin(0)`), and the constructors
[`Add()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md),
[`Mul()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md),
[`Pow()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md),
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) and the
other relations default to `evaluate = FALSE`.

## Usage

``` r
evaluate(x, expr)
```

## Arguments

- x:

  `FALSE` to turn evaluation off, `TRUE` to turn it on.

- expr:

  An R expression.

## Value

The value of `expr`.

## Examples

``` r
x <- Symbol("x")
evaluate(FALSE, x + x)
#> x + x
evaluate(FALSE, sin(Integer(0)))
#> sin(0)
x + x
#> 2*x
```
