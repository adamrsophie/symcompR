# Add a numerical implementation to an undefined function

`implemented_function(symfunc, implementation)` returns an undefined
function named after `symfunc` whose applications carry
`implementation`, as SymPy's `implemented_function()`.
[`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
evaluates the applications with the implementation, and the code
printers inline an implementation that is a
[`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md).

## Usage

``` r
implemented_function(symfunc, implementation)
```

## Arguments

- symfunc:

  A function name, or an undefined function made by
  [`Function()`](https://adamrsophie.github.io/symcompR/reference/Function.md).

- implementation:

  A
  [`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md)
  or an R function.

## Value

An undefined function (class `sym_function`).

## Examples

``` r
x <- Symbol("x")
f <- implemented_function("f", Lambda(x, 2 * x))
ccode(f(x))
#> [1] "2*x"
```
