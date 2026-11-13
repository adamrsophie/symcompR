# C99 math functions

Unevaluated functions for the C99 math library, after SymPy's
`sympy.codegen.cfunctions`: `sym_expm1(x)` is `exp(x) - 1`,
`sym_log1p(x)` is `log(x + 1)`, `exp2(x)` is `2^x`, `sym_log2(x)` and
`sym_log10(x)` are the base 2 and base 10 logarithms, `fma(x, y, z)` is
`x*y + z`, `Sqrt(x)` and `Cbrt(x)` are the square and cube roots,
`hypot(x, y)` is `sqrt(x^2 + y^2)`, and `isnan(x)` and `isinf(x)` test
for a NaN and an infinity. They differentiate, and
`expand(e, func = TRUE)` or
[`rewrite()`](https://adamrsophie.github.io/symcompR/reference/rewrite.md)
turns them into elementary functions; the code printers print them as
the matching C calls.

## Usage

``` r
sym_expm1(x)

sym_log1p(x)

exp2(x)

sym_log2(x)

fma(x, y, z)

sym_log10(x)

Sqrt(x)

Cbrt(x)

hypot(x, y)

isnan(x)

isinf(x)
```

## Arguments

- x, y, z:

  Expressions.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
sym_expm1(x)
#> expm1(x)
expand(sym_log1p(x), func = TRUE)
#> log(x + 1)
diff(hypot(x, 2), x)
#> x/hypot(x, 2)
```
