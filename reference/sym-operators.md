# Arithmetic, comparison and logic operators for symbolic expressions

`+`, `-`, `*`, `/`, `^`, `%%` and `%/%` build canonical symbolic
expressions. `==` and `!=` test structural equality, as in SymPy; use
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) for a
symbolic equation. `<`, `<=`, `>` and `>=` return `TRUE` or `FALSE` when
the answer is decidable and an unevaluated relational otherwise. `&`,
`|` and `!` build logical expressions.

## Usage

``` r
# S3 method for class 'sym'
Ops(e1, e2)

# S3 method for class 'sym'
Math(x, ...)

# S3 method for class 'sym'
Summary(..., na.rm = FALSE)

# S3 method for class 'sym'
Complex(z)
```

## Arguments

- e1, e2:

  Operands; at least one is a `sym`.

- x:

  A `sym` expression.

- ...:

  Further arguments, such as `base` for
  [`log()`](https://rdrr.io/r/base/Log.html).

- na.rm:

  Ignored.

- z:

  A `sym` expression.

## Value

A `sym` expression or a logical value.

## Examples

``` r
x <- Symbol("x")
(x + 1)^2 - 1
#> (x + 1)^2 - 1
x + x == 2 * x
#> [1] TRUE
Symbol("p", positive = TRUE) > 0
#> [1] TRUE
```
