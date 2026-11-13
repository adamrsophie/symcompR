# More code printers

`rust_code()` prints Rust (method call syntax such as `x.powi(2)`),
`glsl_code()` GLSL, `smtlib_code()` SMT-LIB s-expressions (relations
become `assert`s, with `declare-const` lines for the symbols when
`auto_declare` is `TRUE`) and `python()` a Python snippet that rebuilds
the expression with SymPy. `multiline_latex(lhs, rhs, terms_per_line)`
breaks a long right side over lines of an `align*` environment.
`print_glsl()`, `print_python()` and `print_maple_code()` print the
strings.

## Usage

``` r
rust_code(e, assign_to = NULL, ...)

glsl_code(e, assign_to = NULL)

print_glsl(e, assign_to = NULL)

print_maple_code(e, assign_to = NULL)

python(e)

print_python(e)

smtlib_code(e, auto_declare = TRUE)

multiline_latex(lhs, rhs, terms_per_line = 1, environment = "align*")
```

## Arguments

- e:

  An expression, or for `smtlib_code()` a list of them.

- assign_to:

  Optional variable name.

- ...:

  Printer settings of `rust_code()`, as in SymPy (for example
  `user_functions` or `contract`).

- auto_declare:

  Declare the symbols as reals.

- lhs, rhs:

  Sides of the equation.

- terms_per_line:

  Terms per line.

- environment:

  The LaTeX environment.

## Value

A character scalar; the `print_` functions print it.

## Examples

``` r
x <- Symbol("x")
rust_code(x^2 + sin(x)/3)
#> [1] "x.powi(2) + (1_f64/3.0)*x.sin()"
glsl_code(sqrt(x + 1))
#> [1] "sqrt(x + 1)"
smtlib_code(Eq(x + 1, 3))
#> [1] "(declare-const x Real)\n(assert (= (+ 1 x) 3))"
```
