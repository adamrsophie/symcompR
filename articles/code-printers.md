# Code printers

The code printers of `sympy.printing` translate an expression into
source code for another programming language.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); k <- Symbol("k")
f <- Function("f")
```

## Expressions and assignments

The code printers translate an expression into an expression in the
target language. With `assign_to`, they produce a complete assignment
statement.

``` r

e <- x^2 + sin(x) / 3
ccode(e)
#> [1] "pow(x, 2) + (1.0/3.0)*sin(x)"
ccode(sin(x)^2 + cos(x)^2, assign_to = "z")
#> [1] "z = pow(sin(x), 2) + pow(cos(x), 2);"
cxxcode(e)
#> [1] "std::pow(x, 2) + (1.0/3.0)*std::sin(x)"
jscode(e)
#> [1] "Math.pow(x, 2) + (1/3)*Math.sin(x)"
pycode(e)
#> [1] "x**2 + (1/3)*math.sin(x)"
rust_code(e)
#> [1] "x.powi(2) + (1_f64/3.0)*x.sin()"
maple_code(e)
#> [1] "x^2 + (1/3)*sin(x)"
rcode(x^2 + exp(x))
#> [1] "x^2 + exp(x)"
glsl_code(sqrt(x + 1))
#> [1] "sqrt(x + 1)"
mathematica_code(sin(x)^2 + exp(x))
#> [1] "Exp[x] + Sin[x]^2"
```

[`smtlib_code()`](https://adamrsophie.github.io/symcompR/reference/more-code-printers.md)
writes SMT-LIB, the input language of SMT solvers such as Z3:

``` r

cat(smtlib_code(Eq(x + 1, 3)), "\n")
#> (declare-const x Real)
#> (assert (= (+ 1 x) 3))
```

For evaluating expressions in R itself, use
[`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
(see the getting started vignette), which builds an R function directly.

[`fcode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md)
writes free form Fortran (SymPy’s `source_format = "free"`), and
[`octave_code()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md)
writes Octave and MATLAB code. The `print_*()` forms print the code
instead of returning it.

``` r

fcode(x^3 + Rational(1, 2))
#> [1] "      x**3 + 1.0d0/2.0d0"
fcode(sin(x) + 2 * y, assign_to = "r")
#> [1] "      r = 2*y + sin(x)"
octave_code(sin(x) + 2 * y)
#> [1] "2*y + sin(x)"
print_ccode(sin(x) + 2 * y)
#> 2*y + sin(x)
print_fcode(x^3 + Rational(1, 2))
#>       x**3 + 1.0d0/2.0d0
print_jscode(sin(x) + 2 * y)
#> 2*y + Math.sin(x)
print_rcode(sin(x) + 2 * y)
#> 2*y + sin(x)
print_glsl(sin(x) + 2 * y)
#> 2*y + sin(x)
print_maple_code(sin(x) + 2 * y)
#> 2*y + sin(x)
julia_code(sin(x) + 2 * y)
#> [1] "2 * y + sin(x)"
```
