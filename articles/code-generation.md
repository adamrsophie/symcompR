# Code generation

The printing vignette shows
[`ccode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md),
[`fcode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md),
[`pycode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md)
and the other functions that print an expression as code. This vignette
covers the rest of `sympy.codegen` and `sympy.utilities.codegen`: an
abstract syntax tree for whole programs, nodes specific to C, C++,
Fortran, NumPy and SciPy, special functions that map to fast library
calls, rewriting optimizations, the printer objects themselves, and
generators that write complete source files.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
i <- Symbol("i", integer = TRUE); n <- Symbol("n", integer = TRUE)
xr <- Symbol("x", real = TRUE)
```

## Statements and blocks

Assignments, augmented assignments and code blocks:

``` r

Assignment(x, y + 1)
#> Assignment(x, y + 1)
AddAugmentedAssignment(x, 1)
#> AddAugmentedAssignment(x, 1)
SubAugmentedAssignment(x, 1)
#> SubAugmentedAssignment(x, 1)
MulAugmentedAssignment(x, 2)
#> MulAugmentedAssignment(x, 2)
DivAugmentedAssignment(x, 2)
#> DivAugmentedAssignment(x, 2)
ModAugmentedAssignment(x, 2)
#> ModAugmentedAssignment(x, 2)
aug_assign(x, "+", y)
#> AddAugmentedAssignment(x, y)
cat(ccode(CodeBlock(Assignment(x, 1), AddAugmentedAssignment(x, y))), "\n")
#> x = 1;
#> x += y;
ast_attr(CodeBlock(Assignment(x, 1), Assignment(y, x + 1)), "left_hand_sides")
#> (x, y)
```

[`codeblock_topological_sort()`](https://adamrsophie.github.io/symcompR/reference/codeblock-methods.md)
orders assignments so that every variable is assigned before it is used,
and
[`codeblock_cse()`](https://adamrsophie.github.io/symcompR/reference/codeblock-methods.md)
eliminates common subexpressions:

``` r

cat(ccode(codeblock_topological_sort(list(Assignment(x, y + z), Assignment(z, 1),
                                          Assignment(y, 2)))), "\n")
#> z = 1;
#> y = 2;
#> x = y + z;
cat(ccode(codeblock_cse(CodeBlock(Assignment(x, sin(y)^2 + sin(y)),
                                  Assignment(z, sin(y))))), "\n")
#> x0 = sin(y);
#> x = pow(x0, 2) + x0;
#> z = x0;
```

Loops, with `break_` and `continue_`:

``` r

cat(ccode(For(i, Range(0, 10), list(AddAugmentedAssignment(x, i)))), "\n")
#> for (i = 0; i < 10; i += 1) {
#>    x += i;
#> }
cat(ccode(While(x < 10, list(AddAugmentedAssignment(x, 1), break_, continue_))), "\n")
#> while (x < 10) {
#>    x += 1;
#>    break;
#>    continue;
#> }
cat(ccode(Scope(list(Assignment(x, 1)))), "\n")
#> {
#>    x = 1;
#> }
```

## Variables, types and declarations

[`Variable()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
and
[`Pointer()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
carry a type and attributes such as `value_const` and `pointer_const`;
[`Declaration()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
and
[`as_Declaration()`](https://adamrsophie.github.io/symcompR/reference/codegen-types.md)
declare them.

``` r

ccode(Declaration(Variable(x, type = real, value = 1.5)))
#> [1] "double x = 1.5"
ccode(as_Declaration(Variable(x, type = float32), value = 2))
#> [1] "float x = 2"
ccode(as_Declaration(Variable(x, type = int32, attrs = list(value_const)), value = 3))
#> [1] "const int32_t x = 3"
ccode(as_Declaration(Pointer(x, type = real)))
#> [1] "double * x"
ccode(Declaration(Pointer(x, type = real, attrs = list(pointer_const))))
#> [1] "double * const x"
ccode(ast_getitem(Pointer(x, type = real), i))
#> [1] "x[i]"
ccode(ast_getitem(Pointer(x, type = real), list(i, 2)))
#> [1] "x[i][2]"
ccode(Element("x", "ij"))
#> [1] "x[i][j]"
ccode(Element("x", list(i, 2)))
#> [1] "x[i][2]"
```

[`variable_deduced()`](https://adamrsophie.github.io/symcompR/reference/codegen-types.md)
and
[`pointer_deduced()`](https://adamrsophie.github.io/symcompR/reference/codegen-types.md)
choose the type from the assumptions on the symbol, as does
[`type_from_expr()`](https://adamrsophie.github.io/symcompR/reference/codegen-types.md):

``` r

variable_deduced(xr)
#> Variable(x, type=real)
pointer_deduced(xr)
#> Pointer(x, type=real)
type_from_expr(n)
#> integer
type_from_expr(xr)
#> real
```

The predefined types carry their limits, read with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md):

``` r

ast_attr(float16, "dig")
#> 3
ast_attr(float32, "max")
#> 340282346638528859811704183484516925440
ast_attr(float64, "eps")
#> 1/4503599627370496
ast_attr(float80, "nbits")
#> 80
ast_attr(float128, "nbits")
#> 128
ast_attr(float256, "nbits")
#> 256
ast_attr(int8, "min")
#> -128
ast_attr(int16, "max")
#> 32767
ast_attr(int32, "max")
#> 2147483647
ast_attr(int64, "min")
#> -9223372036854775808
ast_attr(uint8, "max")
#> 255
ast_attr(uint16, "max")
#> 65535
ast_attr(uint32, "max")
#> 4294967295
ast_attr(uint64, "max")
#> 18446744073709551615
ast_attr(complex128, "nbits")
#> 128
complex64
#> complex64
bool_
#> bool
complex_
#> complex
intc
#> intc
intp
#> intp
sym_integer
#> integer
untyped
#> untyped
type_cast_check(int8, 127)
#> 127
type_cast_nocheck(float32, S(1) / 3)
#> 0.333333333
```

New types are built from the type classes:

``` r

Type("float")
#> float
FloatType("f32", nbits = 32, nmant = 23, nexp = 8)
#> f32
ComplexType("c64", nbits = 64, nmant = 23, nexp = 8)
#> c64
SignedIntType("i8", 8)
#> i8
UnsignedIntType("u8", 8)
#> u8
IntBaseType("int")
#> int
FloatBaseType("float")
#> float
ComplexBaseType("complex")
#> complex
```

## Functions, strings and other nodes

``` r

ccode(FunctionCall("foo", list(x, y)))
#> [1] "foo(x, y)"
ccode(FunctionPrototype(real, "f", list(Variable(x, type = real))))
#> [1] "double f(double x)"
fdef <- FunctionDefinition(real, "f", list(Variable(x, type = real)), list(Return(x^2)))
cat(ccode(fdef), "\n")
#> double f(double x){
#>    return pow(x, 2);
#> }
ccode(from_FunctionDefinition(fdef))
#> [1] "double f(double x)"
cat(ccode(from_FunctionPrototype(FunctionPrototype(real, "f", list(Variable(x, type = real))),
                                 list(Return(x^2)))), "\n")
#> double f(double x){
#>    return pow(x, 2);
#> }
ccode(Return(x + 1))
#> [1] "return x + 1"
ccode(Comment("hello"))
#> [1] "/* hello */"
String("foo")
#> foo
QuotedString("foo")
#> foo
ccode(QuotedString("foo"))
#> [1] "\"foo\""
pycode(Print(list(x, y), "x=%s y=%s"))
#> [1] "print(\"x=%s y=%s\" % (x, y), end=\"\")"
pycode(Print(list(x), file = stderr_))
#> [1] "print(x, file=sys.stderr)"
stdout_
#> Stream(stdout)
Stream("stdout")
#> Stream(stdout)
pycode(Raise(RuntimeError_("boom")))
#> [1] "raise RuntimeError(boom)"
pycode(KeywordFunctionCall("f", list(x), list(k = 1)))
#> [1] "f(x)"
Attribute("volatile")
#> volatile
Token()
#> Token()
Node()
#> Node()
NoneToken()
#> NoneToken()
none
#> NoneToken()
BreakToken()
#> BreakToken()
ContinueToken()
#> ContinueToken()
```

Nodes are inspected with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md),
[`ast_kwargs()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md),
[`attr_params()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md)
and
[`ast_func()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md)
(which returns the constructor):

``` r

v <- Variable(x, type = real)
ast_attr(v, "type")
#> real
ast_kwargs(v, exclude = "type")
#> $symbol
#> x
#> 
#> $value
#> NoneToken()
#> 
#> $attrs
#> ()
attr_params(Variable(x, type = real, attrs = list(value_const)), "value_const")
#> ()
ast_func(Assignment(x, 1))(y, 2)
#> Assignment(y, 2)
```

## C and C++ nodes

``` r

ccode(sizeof("unsigned int"))
#> [1] "sizeof(unsigned int)"
ccode(alignof("unsigned int"))
#> [1] "alignof(unsigned int)"
ccode(CommaOperator(PreIncrement(x), 2 * x))
#> [1] "(++(x), 2*x)"
ccode(PostIncrement(x))
#> [1] "(x)++"
ccode(PreDecrement(x))
#> [1] "--(x)"
ccode(PostDecrement(x))
#> [1] "(x)--"
cat(ccode(Label("start", list(Assignment(x, 1)))), "\n")
#> start:
#> x = 1;
ccode(goto("start"))
#> [1] "goto start"
ccode(Declaration(Variable(x, type = real, attrs = list(volatile))))
#> [1] "double x"
ccode(Declaration(Variable(x, type = real, attrs = list(static))))
#> [1] "double x"
ccode(Declaration(Pointer(x, type = real, attrs = list(restrict))))
#> [1] "double * restrict x"
cat(ccode(struct("pt", list(Variable(x, type = real), Variable(y, type = real)))), "\n")
#> struct pt {
#>    double x;
#>    double y;
#> }
cat(ccode(c_union("u", list(Variable(x, type = real), Variable(y, type = int32)))), "\n")
#> union u {
#>    double x;
#>    int32_t y;
#> }
cxxcode(using("std::vector"))
#> [1] "using std::vector"
cxxcode(using("std::vector", "vec"))
#> [1] "using vec = std::vector"
constexpr
#> constexpr
```

## Fortran nodes

``` r

fc <- function(e, ...) cat(fcode(e, source_format = "free", ...), "\n")
fc(Do(list(AddAugmentedAssignment(x, i)), i, 1, n))
#> do i = 1, n
#>    x = x + i
#> end do
fc(isign(1, x))
#> isign(1, x)
fc(dsign(1, x))
#> dsign(1, x)
fc(cmplx(x, y))
#> cmplx(x, y)
fc(kind(x))
#> kind(x)
fc(literal_dp(0.5))
#> 5d-1
fc(literal_sp(0.5))
#> 5e-1
fc(fortran_merge(x, y, x > y), standard = 95)
#> merge(x, y, x > y)
a <- Symbol("a")
fc(size(a))
#> size(a)
fc(lbound(a, 1))
#> lbound(a, 1)
fc(ubound(a, 1))
#> ubound(a, 1)
fc(fortran_shape(a))
#> shape(a)
fc(fortran_reshape(a, list(2, 3)))
#> reshape(a, (/2, 3/))
fc(allocated(a))
#> allocated(a)
fc(sum_(a))
#> sum(a)
fc(product_(a))
#> product(a)
fc(fortran_use("signallib", only = list("snr")))
#> use signallib, only: snr
fc(fortran_use("signallib", rename = list(list("thingy", "convolution2d"))))
#> use signallib, thingy => convolution2d
fc(use_rename("thingy", "convolution2d"))
#> thingy => convolution2d
fc(Program("prog", list(Print(list(x)))))
#> program prog
#>    print *, x
#> end program
fc(Module("mod", list(), list(Declaration(Variable(x, type = real)))), standard = 2003)
#> module mod
#> 
#> 
#> contains
#> 
#> real*8 :: x
#> end module
fc(Subroutine("sub", list(Variable(x, type = real, attrs = list(intent_inout))),
              list(Assignment(x, 1))), standard = 2003)
#> subroutine sub(x)
#> real*8, intent(inout) :: x
#> x = 1
#> end subroutine
fc(SubroutineCall("sub", list(x)))
#> call sub(x)
fc(FortranReturn(x))
#> return x
fc(GoTo(list(10, 20), i))
#> go to (10, 20), i
fc(ArrayConstructor(list(1, 2, 3)))
#> (/1, 2, 3/)
fc(ArrayConstructor(list(ImpliedDoLoop(i^2, i, 1, 4))))
#> (/(i**2, i = 1, 4)/)
```

Arrays are declared with
[`fortran_array()`](https://adamrsophie.github.io/symcompR/reference/codegen-fnodes.md),
or with the `intent_*`, `allocatable` and
[`fortran_dimension()`](https://adamrsophie.github.io/symcompR/reference/codegen-fnodes.md)
attributes; `assumed_extent` marks a dimension of unknown size and
[`Extent()`](https://adamrsophie.github.io/symcompR/reference/codegen-fnodes.md)
a range.

``` r

fc(Declaration(fortran_array(a, list(4), "in", type = real)), standard = 2003)
#> real*8, dimension(4), intent(in) :: a
fc(Declaration(Variable(a, type = real, attrs = list(intent_in, fortran_dimension(3)))),
   standard = 2003)
#> real*8, dimension(3), intent(in) :: a
fc(Declaration(Variable(a, type = real, attrs = list(intent_out, fortran_dimension(assumed_extent)))),
   standard = 2003)
#> real*8, dimension(:), intent(out) :: a
fc(Declaration(Variable(a, type = real, attrs = list(allocatable, fortran_dimension(assumed_extent)))),
   standard = 2003)
#> real*8, dimension(:), allocatable :: a
Extent(1, 3)
#> 1:3
assumed_size
#> dimension(*)
intent_out
#> intent_out
fc(FunctionDefinition(real, "f", list(Variable(x, type = real, attrs = list(intent_in))),
                      list(Return(x^2)), attrs = list(bind_C("f"), elemental, pure)),
   standard = 2003)
#> elemental real*8 function f(x) bind(C, name="f")
#> real*8, intent(in) :: x
#> f = x**2
#> end function
void
#> void
```

## Library functions

`sympy.codegen.cfunctions` defines functions that C libraries compute
accurately, with derivatives and rewrites:

``` r

sym_expm1(x)
#> expm1(x)
diff(sym_expm1(x), x)
#> exp(x)
rewrite(sym_expm1(x), "exp")
#> exp(x) - 1
diff(sym_log1p(x), x)
#> 1/(x + 1)
expand(sym_log1p(x), func = TRUE)
#> log(x + 1)
diff(exp2(x), x)
#> log(2)*exp2(x)
expand(exp2(x), func = TRUE)
#> 2^x
expand(sym_log2(x), func = TRUE)
#> log(x)/log(2)
expand(sym_log10(x), func = TRUE)
#> log(x)/log(10)
expand(fma(x, y, z), func = TRUE)
#> x*y + z
diff(fma(x, y, z), x)
#> y
expand(Sqrt(x), func = TRUE)
#> sqrt(x)
expand(Cbrt(x), func = TRUE)
#> x^(1/3)
expand(hypot(x, y), func = TRUE)
#> sqrt(x^2 + y^2)
diff(hypot(x, 2), x)
#> x/hypot(x, 2)
ccode(isnan(x))
#> [1] "isnan(x)"
ccode(isinf(x))
#> [1] "isinf(x)"
ccode(sym_expm1(x), standard = "C99")
#> [1] "expm1(x)"
```

NumPy and SciPy nodes, and the matrix solve node:

``` r

diff(logaddexp(x, y), x)
#> 1/(exp(-x + y) + 1)
rewrite(logaddexp(x, y), "log")
#> log(exp(x) + exp(y))
rewrite(logaddexp2(x, y), "log")
#> log(2^x + 2^y)/log(2)
rewrite(cosm1(x), "cos")
#> cos(x) - 1
rewrite(powm1(x, y), "Pow")
#> x^y - 1
printer_doprint(NumPyPrinter(), amax(a))
#> [1] "numpy.amax(a, axis=None)"
printer_doprint(NumPyPrinter(), amin(a))
#> [1] "numpy.amin(a, axis=None)"
printer_doprint(NumPyPrinter(), logaddexp(x, y))
#> [1] "numpy.logaddexp(x, y)"
printer_doprint(SciPyPrinter(), cosm1(x))
#> [1] "scipy.special.cosm1(x)"
printer_doprint(SciPyPrinter(), powm1(x, y))
#> [1] "scipy.special.powm1(x, y)"
A <- MatrixSymbol("A", 3, 3); b <- MatrixSymbol("b", 3, 1)
printer_doprint(NumPyPrinter(), MatrixSolve(A, b))
#> [1] "numpy.linalg.solve(A, b)"
List(1, 2, 3)
#> [1, 2, 3]
NumExprEvaluate(x + y)
#> NumExprEvaluate(x + y)
```

## Printer objects

Each language has a printer object;
[`printer_doprint()`](https://adamrsophie.github.io/symcompR/reference/code-printer-objects.md)
prints an expression with it, optionally assigning the result.

``` r

printer_doprint(C89CodePrinter(), x^2 + sqrt(x))
#> [1] "sqrt(x) + pow(x, 2)"
printer_doprint(C99CodePrinter(), sym_expm1(x) + Max(x, y))
#> [1] "expm1(x) + fmax(x, y)"
printer_doprint(C99CodePrinter(), x^2, assign_to = y)
#> [1] "y = pow(x, 2);"
printer_doprint(C11CodePrinter(), alignof("int"))
#> [1] "alignof(int)"
printer_doprint(CXX98CodePrinter(), Max(x, y))
#> [1] "std::max(x, y)"
printer_doprint(CXX11CodePrinter(), Max(x, y) + sym_log1p(x))
#> [1] "std::log1p(x) + std::max(x, y)"
printer_doprint(CXX17CodePrinter(), gamma(x))
#> [1] "std::tgamma(x)"
printer_doprint(FCodePrinter(list(source_format = "free")), x^2)
#> [1] "x**2"
printer_doprint(JavascriptCodePrinter(), Max(x, y) + sin(x))
#> [1] "Math.sin(x) + Math.max(x, y)"
printer_doprint(JuliaCodePrinter(), x^2 + sin(x))
#> [1] "x .^ 2 + sin(x)"
printer_doprint(OctaveCodePrinter(), x^2 + sin(x))
#> [1] "x.^2 + sin(x)"
printer_doprint(PythonCodePrinter(), sin(x) + x^2)
#> [1] "x**2 + math.sin(x)"
printer_doprint(MpmathPrinter(), sin(x) + Rational(1, 3))
#> [1] "mpmath.sin(x) + mpmath.mpf(1)/mpmath.mpf(3)"
printer_doprint(NumPyPrinter(), sin(x) + x^2)
#> [1] "x**2 + numpy.sin(x)"
printer_doprint(SciPyPrinter(), gamma(x))
#> [1] "scipy.special.gamma(x)"
printer_doprint(LambdaPrinter(), And(x > 0, y > 0))
#> [1] "(((x > 0)) and ((y > 0)))"
printer_doprint(RCodePrinter(), x^2 + sin(x))
#> [1] "x^2 + sin(x)"
printer_doprint(RustCodePrinter(), x^2 + sin(x))
#> [1] "x.powi(2) + x.sin()"
cgdict_get(get_math_macros(), exp(Integer(1)))
#> [1] "M_E"
cgdict_get(get_math_macros(), sqrt(Integer(2)))
#> [1] "M_SQRT2"
```

## Rewriting optimizations

[`sym_optimize()`](https://adamrsophie.github.io/symcompR/reference/codegen-rewriting.md)
applies rewriting rules that make generated code faster or more
accurate:

``` r

sym_optimize(3 * exp(x) - 3, list(expm1_opt))
#> 3*expm1(x)
sym_optimize(log(x + 1), list(log1p_opt))
#> log1p(x)
sym_optimize(2^x, list(exp2_opt))
#> exp2(x)
sym_optimize(log(x) / log(Integer(2)), list(log2_opt))
#> log2(x)
sym_optimize(log(x) / log(Integer(2)), list(log2const_opt))
#> log(x)/log(2)
sym_optimize(log(exp(x) + exp(y)), list(logaddexp_opt))
#> logaddexp(x, y)
sym_optimize(log(exp(x) + exp(y)), list(logsumexp_2terms_opt))
#> log1p(exp(Min(x, y))) + Max(x, y)
sym_optimize(log(2^x + 2^y) / log(Integer(2)), list(logaddexp2_opt))
#> logaddexp2(x, y)
sym_optimize(cos(x) - 1, list(cosm1_opt))
#> cosm1(x)
sym_optimize(x^y - 1, list(powm1_opt))
#> powm1(x, y)
sym_optimize(sin(x) / x, list(sinc_opt1))
#> sinc(x)
sym_optimize(-sin(x) / x, list(sinc_opt2))
#> -sin(x)/x
sym_optimize(sin(x) / x, sinc_opts)
#> sinc(x)
sym_optimize(log(x) / log(Integer(2)) + log(x + 1), optims_c99)
#> log1p(x) + log2(x)
sym_optimize(log(exp(x) + exp(y)), optims_numpy)
#> logaddexp(x, y)
sym_optimize(cos(x) - 1, optims_scipy)
#> cosm1(x)
sym_optimize(MatMul(Inverse(A), b), list(matinv_opt))
#> A^(-1)*b
ccode(sym_optimize(x^4 + x^2, list(create_expand_pow_optimization(4))))
#> [1] "x*x + x*x*x*x"
sym_optimize(exp(x) - 1 + log(1 + x), list(FuncMinusOneOptim(exp, sym_expm1)))
#> log(x + 1) + expm1(x)
w <- Wild("w")
sym_optimize(log(x + 1) + log(y + 1), list(ReplaceOptim(log(w + 1), sym_log1p(w))))
#> log1p(x) + log1p(y)
```

[`SumApprox()`](https://adamrsophie.github.io/symcompR/reference/SumApprox.md)
drops terms that are small within given bounds, and
[`SeriesApprox()`](https://adamrsophie.github.io/symcompR/reference/SumApprox.md)
replaces functions by truncated series:

``` r

sa <- SumApprox(list(list(x, list(0, 0.001)), list(y, list(100, 1000))), reltol = 0.01)
sym_optimize(x + y * z, list(sa))
#> x + y*z
sym_optimize(sin(x), list(SeriesApprox(list(list(x, list(-0.1, 0.1))), reltol = 0.01)))
#> x
```

## Algorithms and rendering

[`newtons_method()`](https://adamrsophie.github.io/symcompR/reference/newtons_method.md)
builds a loop for Newton’s method, and
[`newtons_method_function()`](https://adamrsophie.github.io/symcompR/reference/newtons_method.md)
a complete function:

``` r

cat(ccode(newtons_method(cos(xr) - xr^3, xr)), "\n")
#> {
#>    double _Dummy_120 = INFINITY;
#>    while (fabs(_Dummy_120) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {
#>       _Dummy_120 = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));
#>       x += _Dummy_120;
#>    };
#> }
cat(ccode(newtons_method_function(cos(xr) - xr^3, xr)), "\n")
#> double newton(double x){
#>    double d_x = INFINITY;
#>    while (fabs(d_x) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {
#>       d_x = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));
#>       x += d_x;
#>    };
#>    return x;
#> }
```

The render functions write a complete source file or module:

``` r

cat(render_as_source_file(Assignment(xr, sin(xr))), "\n")
#> #include <math.h>
#> 
#> x = sin(x);
cat(py_render_as_module(Print(list(x), "%f")), "\n")
#> 
#> 
#> print("%f" % (x), end="")
cat(f_render_as_module(list(FunctionDefinition(real, "f",
      list(Variable(x, type = real, attrs = list(intent_in))), list(Return(x^2)))), "mymod"), "\n")
#> module mymod
#>    
#>    implicit none
#>    private
#>    public f
#> 
#> 
#> contains
#> 
#> real*8 function f(x)
#> real*8, intent(in) :: x
#> f = x**2
#> end function
#> 
#> end module
```

[`lambdastr()`](https://adamrsophie.github.io/symcompR/reference/lambdastr.md)
writes the source of a Python lambda, and
[`implemented_function()`](https://adamrsophie.github.io/symcompR/reference/implemented_function.md)
attaches a numerical implementation to a new function:

``` r

lambdastr(x, x^2)
#> [1] "lambda x: (x**2)"
lambdastr(list(x, list(y)), x + y)
#> [1] "lambda _0,_1: (lambda x,y: (x + y))(_0,_1[0])"
f <- implemented_function("f", Lambda(x, 2 * x))
ccode(f(x))
#> [1] "2*x"
```

## Routines and code generators

[`make_routine()`](https://adamrsophie.github.io/symcompR/reference/Routine.md)
collects the arguments and results of an expression; the argument and
result classes are read with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md).

``` r

r <- make_routine("f", x^2 + y)
lapply(ast_attr(r, "arguments"), ast_attr, "name")
#> [[1]]
#> x
#> 
#> [[2]]
#> y
ast_attr(ast_attr(r, "results")[[1]], "expr")
#> x^2 + y
ast_attr(get_default_datatype(x), "cname")
#> [1] "double"
ast_attr(get_default_datatype(n), "cname")
#> [1] "int"
ast_attr(DataType("double", "REAL*8", "float", "", "", "f64"), "cname")
#> [1] "double"
ast_attr(InputArgument(x), "name")
#> x
ast_attr(OutputArgument(y, y, x^2), "expr")
#> x^2
ast_attr(InOutArgument(x, x, x + 1), "expr")
#> x + 1
ast_attr(Result(x^2), "expr")
#> x^2
ast_attr(Routine("f", list(InputArgument(x)), list(Result(x^2, result_var = Symbol("r"))),
                 list(), list()), "name")
#> [1] "f"
```

[`codegen()`](https://adamrsophie.github.io/symcompR/reference/code-generators.md)
writes source and header files for named expressions; the generator
classes do the same step by step.

``` r

res <- codegen(list("f", x * y), "C99", header = FALSE, empty = FALSE)
cat(res[[1]][[2]])
#> #include "f.h"
#> #include <math.h>
#> double f(double x, double y) {
#>    double f_result;
#>    f_result = x*y;
#>    return f_result;
#> }
cat(res[[2]][[2]])
#> #ifndef PROJECT__F__H
#> #define PROJECT__F__H
#> double f(double x, double y);
#> #endif
cat(codegen(list("f", x * y), "F95", header = FALSE, empty = FALSE)[[1]][[2]])
#> REAL*8 function f(x, y)
#> implicit none
#> REAL*8, intent(in) :: x
#> REAL*8, intent(in) :: y
#> f = x*y
#> end function
write_with <- function(g, prefix = "test") {
  rt <- ast_attr(g, "routine")("f", x * y, NULL, NULL)
  cat(ast_attr(g, "write")(list(rt), prefix, header = FALSE, empty = FALSE)[[1]][[2]])
}
write_with(CCodeGen())
#> #include "test.h"
#> #include <math.h>
#> double f(double x, double y) {
#>    double f_result;
#>    f_result = x*y;
#>    return f_result;
#> }
write_with(C89CodeGen())
#> #include "test.h"
#> #include <math.h>
#> double f(double x, double y) {
#>    double f_result;
#>    f_result = x*y;
#>    return f_result;
#> }
write_with(C99CodeGen())
#> #include "test.h"
#> #include <math.h>
#> double f(double x, double y) {
#>    double f_result;
#>    f_result = x*y;
#>    return f_result;
#> }
write_with(FCodeGen())
#> REAL*8 function f(x, y)
#> implicit none
#> REAL*8, intent(in) :: x
#> REAL*8, intent(in) :: y
#> f = x*y
#> end function
write_with(JuliaCodeGen())
#> function f(x, y)
#>     out1 = x .* y
#>     return out1
#> end
write_with(OctaveCodeGen(), "f")
#> function out1 = f(x, y)
#>   out1 = x.*y;
#> end
write_with(RustCodeGen())
#> fn f(x: f64, y: f64) -> f64 {
#>     let out1 = x*y;
#>     out1
#> }
unname(class(get_code_generator("C99"))[1])
#> [1] "C99CodeGen"
unname(class(get_code_generator("F95"))[1])
#> [1] "FCodeGen"
```
