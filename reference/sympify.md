# Parse a string into a symbolic expression

`sympify()` reads a mathematical expression written with R or Python
operators (`^` or `**` for powers) and returns a `sym` expression.
Number literals keep their exact digits, so
`"123456789012345678901234567890"` and `"1/3"` stay exact, while
decimals become
[`Float()`](https://adamrsophie.github.io/symcompR/reference/Integer.md)
values. Names become symbols except for known constants (`pi`, `E`, `I`,
`oo`) and functions; unknown function names become undefined functions.
The input is parsed, never evaluated as R code. `==` in a string means
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md).

## Usage

``` r
sympify(s, locals = list(), rational = FALSE)

sym_from_call(x, locals = list(), rational = FALSE, from_string = FALSE)
```

## Arguments

- s:

  A character string, or a value accepted by
  [`as_sym()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md).

- locals:

  Optional named list mapping names to `sym` objects or R functions,
  used before the built-in names.

- rational:

  If `TRUE`, decimal literals become exact rationals.

- x:

  An R call, formula, or name.

- from_string:

  Internal flag; `TRUE` when called from `sympify()`.

## Value

A `sym` expression.

## Examples

``` r
sympify("x**2 + 2*x + 1")
#> x^2 + 2*x + 1
sympify("sin(x)^2 + cos(x)^2")
#> sin(x)^2 + cos(x)^2
sympify("1/3 + 0.5")
#> 0.833333333333333
```
