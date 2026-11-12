# Undefined functions

`Function("f")` returns an R function that builds unevaluated
applications such as `f(x)`. Derivatives of these stay as
[`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md)
objects, which is how SymPy represents unknown functions in differential
equations.

## Usage

``` r
Function(name)
```

## Arguments

- name:

  Function name.

## Value

An R function of any number of arguments.

## Examples

``` r
f <- Function("f")
x <- Symbol("x")
f(x) + f(x)
#> 2*f(x)
```
