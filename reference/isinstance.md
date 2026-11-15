# Class membership in SymPy's hierarchy

`isinstance(e, cls)` is Python's `isinstance()` for SymPy classes:
whether `cls` (a class name such as `"Expr"`, `"Atom"`, `"Number"`,
`"Set"`, `"MatrixBase"`, or a vector of names, as a tuple of classes) is
among the classes of `e`. `AppliedUndef` is the class of an applied
undefined function, `FunctionClass` the class of a function such as
`Function("f")`. `RealNumber()` is SymPy's alias of
[`Float()`](https://adamrsophie.github.io/symcompR/reference/Integer.md).

## Usage

``` r
isinstance(e, cls)

RealNumber(...)
```

## Arguments

- e:

  An object.

- cls:

  A character vector of class names.

- ...:

  Arguments of
  [`Float()`](https://adamrsophie.github.io/symcompR/reference/Integer.md).

## Value

`TRUE` or `FALSE`.

## Examples

``` r
x <- Symbol("x")
isinstance(x, "Atom")
#> [1] TRUE
isinstance(x + 1, c("Number", "Add"))
#> [1] TRUE
isinstance(Function("f")(x), "AppliedUndef")
#> [1] TRUE
```
