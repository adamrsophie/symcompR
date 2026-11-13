# Reference values from SymPy (sympy.codegen and sympy.utilities.codegen).

test_that("AST nodes print through the code printers", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(ccode(Assignment(x, y + 1)), "x = y + 1;")
  expect_equal(ccode(aug_assign(x, "+", y)), "x += y;")
  expect_equal(sstr(Assignment(x, y)), "Assignment(x, y)")
  n <- Symbol("n", integer = TRUE); i <- Idx("i", n)
  A <- IndexedBase("A"); B <- IndexedBase("B")
  expect_equal(fcode(Assignment(A(i), B(i)), source_format = "free"), "do i = 1, n\n   A(i) = B(i)\nend do")
  expect_equal(fcode(Do(list(Assignment(x, x + 1)), Symbol("i", integer = TRUE), 1, 10), source_format = "free"),
               "do i = 1, 10\n   x = x + 1\nend do")
  expect_equal(pycode(While(Lt(abs(x), 1), list(AddAugmentedAssignment(x, 1)))),
               "while (abs(x) < 1):\n    x += 1")
})

test_that("codegen() returns file names and contents", {
  x <- Symbol("x"); y <- Symbol("y")
  r <- codegen(list("f", x + y), "C99", "test", header = FALSE, empty = FALSE)
  expect_equal(r[[1L]][[1L]], "test.c")
  expect_equal(r[[1L]][[2L]], "#include \"test.h\"\n#include <math.h>\ndouble f(double x, double y) {\n   double f_result;\n   f_result = x + y;\n   return f_result;\n}\n")
  expect_equal(r[[2L]][[2L]], "#ifndef PROJECT__TEST__H\n#define PROJECT__TEST__H\ndouble f(double x, double y);\n#endif\n")
})

test_that("rewriting optimizations", {
  x <- Symbol("x")
  expect_equal(sstr(sym_optimize(3*exp(x) - 3, list(expm1_opt))), "3*expm1(x)")
  expect_equal(sstr(sym_optimize(log(x)/log(Integer(2)) + log(x + 1), optims_c99)), "log1p(x) + log2(x)")
  expect_equal(ccode(sym_optimize(x^4 + x^2, list(create_expand_pow_optimization(4)))), "x*x + x*x*x*x")
})

test_that("Newton's method and source rendering", {
  x <- Symbol("x", real = TRUE)
  f <- newtons_method_function(cos(x) - x^3, x)
  expect_equal(pycode(f), paste0("def newton(x):\n    d_x = math.inf\n    while (abs(d_x) > 4.0e-16*abs(x) + 1.0e-12):\n",
                                 "        d_x = (x**3 - math.cos(x))/(-3*x**2 - math.sin(x))\n        x += d_x\n    return x"))
  expect_equal(render_as_source_file(Assignment(x, sin(x))), "#include <math.h>\n\nx = sin(x);")
  expect_error(newtons_method_function(Symbol("a") * x, x), "Missing symbols in params: a")
})

test_that("implemented_function and lambdastr", {
  x <- Symbol("x"); y <- Symbol("y")
  f <- implemented_function("f", function(t) 2 * t)
  expect_equal(lambdify(x, f(x))(3), 6)
  expect_equal(lambdastr(list(x, list(y)), x + y), "lambda _0,_1: (lambda x,y: (x + y))(_0,_1[0])")
  expect_equal(lambdastr(x, x^2), "lambda x: (x**2)")
})

test_that("language printers follow SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(rcode(sqrt(x + 1), assign_to = "y"), "y = sqrt(x + 1);")
  expect_equal(rcode(Piecewise(list(x, Lt(x, 1)), list(x^2, TRUE))), "ifelse(x < 1,x,x^2)")
  expect_equal(octave_code(x^2), "x.^2")
  expect_equal(jscode(Max(x, y)), "Math.max(x, y)")
  expect_equal(ccode(x^Rational(1, 3)), "cbrt(x)")
  expect_equal(julia_code(x^2 + sin(x)/3), "x .^ 2 + sin(x) / 3")
  expect_equal(mathematica_code(Integral(x^2, list(x, 0, 1))), "Hold[Integrate[x^2, {x, 0, 1}]]")
})
