# Parsing

Expressions can be parsed from strings in SymPy syntax, from LaTeX, and
from Mathematica syntax.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); k <- Symbol("k")
f <- Function("f")
```

## SymPy syntax

[`parse_expr()`](https://adamrsophie.github.io/symcompR/reference/parse_expr.md)
(and
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md))
parse a string in Python and SymPy syntax. Undefined names become
symbols.

``` r

parse_expr("x**2 + 3*x - 1/2")
#> x^2 + 3*x - 1/2
parse_expr("0.[3] + x!")
#> factorial(x) + 1/3
parse_expr("1 + 2", evaluate = FALSE)
#> 1 + 2
simplify(parse_expr("sin(x)**2 + cos(x)**2"))
#> 1
```

Transformations extend the grammar. `"implicit"` allows implicit
multiplication and function application, as typed by hand:

``` r

parse_expr("2x sin x", transformations = "implicit")
#> 2*x*sin(x)
parse_expr("2x", transformations = c(standard_transformations, implicit_multiplication))
#> 2*x
parse_expr("x^2", transformations = c(standard_transformations, convert_xor))
#> x^2
parse_expr("3x^2 + 2xy", transformations = c(standard_transformations,
           implicit_multiplication_application, convert_xor))
#> 3*x^2 + 2*x*y
parse_expr("sin**2 x", transformations = c(standard_transformations,
           implicit_multiplication_application))
#> sin(x)^2
parse_expr("x = 2", transformations = c(standard_transformations, convert_equals_signs))
#> Eq(x, 2)
parse_expr("0.1 + 2/3", transformations = c(standard_transformations, t_rationalize))
#> 23/30
```

`convert_xor` makes `^` mean exponentiation (in Python it is exclusive
or), and `t_rationalize` reads decimal numbers as exact rationals.

### Token transformations

The transformations are functions on the token stream. Besides the
predefined collections, the individual steps can be listed:
[`auto_symbol()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
turns undefined names into symbols,
[`auto_number()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
wraps numbers,
[`factorial_notation()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
handles `!`,
[`repeated_decimals()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
reads `0.[3]`,
[`lambda_notation()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
handles `lambda`,
[`function_exponentiation()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
reads `sin**2 x`,
[`implicit_application()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
reads `sin x`, and
[`split_symbols()`](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md)
splits `xy` into `x*y`.
[`parser_T()`](https://adamrsophie.github.io/symcompR/reference/parser_T.md)
selects transformations by number, as `T[...]` does in SymPy.

``` r

parse_expr("x!", transformations = list(auto_symbol, auto_number, factorial_notation))
#> factorial(x)
parse_expr("0.[3]", transformations = list(auto_symbol, repeated_decimals, auto_number))
#> 1/3
parse_expr("lambda x: x + 1", transformations = list(lambda_notation, auto_symbol, auto_number))
#> Lambda(x, x + 1)
parse_expr("sin**2(x)", transformations = c(standard_transformations, list(function_exponentiation)))
#> sin(x)^2
parse_expr("sin x", transformations = c(standard_transformations, list(implicit_application)))
#> sin(x)
parse_expr("xy", transformations = c(standard_transformations, list(split_symbols, implicit_multiplication)))
#> x*y
parse_expr("2x", transformations = parser_T(0:4, 8))
#> 2*x
keep_mu <- split_symbols_custom(function(s) s != "mu")
parse_expr("mu xy", transformations = c(standard_transformations,
           list(keep_mu, implicit_multiplication)))
#> mu*x*y
```

## Converting R calls

[`sym_from_call()`](https://adamrsophie.github.io/symcompR/reference/sympify.md)
converts an unevaluated R call into an expression, and
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md)
converts strings and R values.

``` r

sym_from_call(quote(x^2 + sin(y)))
#> x^2 + sin(y)
sympify("x + 1/2")
#> x + 1/2
sympify(0.5, rational = TRUE)
#> 0.500000000000000
```

## LaTeX

[`parse_latex()`](https://adamrsophie.github.io/symcompR/reference/parse_latex.md)
reads LaTeX math. Backslashes must be doubled inside R strings (or use
raw strings, `r"(...)"`).

``` r

parse_latex("\\frac{1 + x}{\\sqrt{y}}")
#> (x + 1)/sqrt(y)
parse_latex("x^2 + 2x + 1 = 0")
#> Eq(x^2 + 2*x + 1, 0)
parse_latex("\\sin^2(x) + \\cos^2(x)")
#> sin(x)^2 + cos(x)^2
parse_latex("\\int_0^1 x^2 dx")
#> Integral(x^2, (x, 0, 1))
parse_latex("\\sum_{k=1}^{n} k^2")
#> Sum(k^2, (k, 1, n))
parse_latex(r"(\frac{d}{dx} x^3)")
#> Derivative(x^3, x)
parse_latex(r"(\binom{n}{k})")
#> binomial(n, k)
```

Integrals, sums, and derivatives are returned unevaluated; call
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md) to
evaluate them.

``` r

doit(parse_latex("\\int_0^1 x^2 dx"))
#> 1/3
```

[`parse_latex_lark()`](https://adamrsophie.github.io/symcompR/reference/parse_latex_lark.md)
is the Lark grammar based LaTeX parser of SymPy, and
[`parse_latex_lark_all()`](https://adamrsophie.github.io/symcompR/reference/parse_latex_lark_all.md)
returns every parse of an ambiguous input.

``` r

parse_latex_lark("2\\pi r")
#> 2*pi*r
parse_latex_lark("\\frac{a}{b} + c^2")
#> a/b + c^2
parse_latex_lark_all("x + y")
#> [[1]]
#> x + y
```

## Mathematica

[`parse_mathematica()`](https://adamrsophie.github.io/symcompR/reference/parse_mathematica.md)
reads Mathematica syntax, including implicit multiplication and
bracketed function calls.

``` r

parse_mathematica("Sin[x]^2 Tan[y]")
#> sin(x)^2*tan(y)
parse_mathematica("Log[2, 8] + 1/2")
#> 7/2
parse_mathematica("x*(a + b)")
#> x*(a + b)
```

[`parse_mathematica_to_fullformlist()`](https://adamrsophie.github.io/symcompR/reference/parse_mathematica.md)
stops at the FullForm stage and returns nested lists, with the head of
each expression first.

``` r

parse_mathematica_to_fullformlist("x + Sin[y]")
#> [[1]]
#> [1] "Plus"
#> 
#> [[2]]
#> [1] "x"
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] "Sin"
#> 
#> [[3]][[2]]
#> [1] "y"
```
