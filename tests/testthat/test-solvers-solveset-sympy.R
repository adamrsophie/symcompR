skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/solvers/tests/test_solveset.py, sympy/sets/tests/test_conditionset.py, sympy/sets/tests/test_fancysets.py, sympy/sets/tests/test_sets.py);
# the rows come from the battery that compares this port with SymPy.

test_that("solveset: test_invert_real matches SymPy", {
  expect_identical(.pyrun(invert_real(exp(Symbol("x", real = TRUE)), Symbol("z"), Symbol("x", real = TRUE))), "(x, Intersection({log(z)}, Reals))")
  expect_identical(.pyrun(invert_real((Integer(3) + Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {y - 3})")
  expect_identical(.pyrun(invert_real((Integer(3)*Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {y/3})")
  expect_identical(.pyrun(invert_real(exp(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {log(y)})")
  expect_identical(.pyrun(invert_real(exp((Integer(3)*Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {log(y)/3})")
  expect_identical(.pyrun(invert_real(exp((Integer(3) + Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {log(y) - 3})")
  expect_identical(.pyrun(invert_real((Integer(3) + exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, Intersection({log(y - 3)}, Reals))")
  expect_identical(.pyrun(invert_real((Integer(3)*exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {log(y/3)})")
  expect_identical(.pyrun(invert_real(log(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {exp(y)})")
  expect_identical(.pyrun(invert_real(log((Integer(3)*Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {exp(y)/3})")
  expect_identical(.pyrun(invert_real(log((Integer(3) + Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {exp(y) - 3})")
  expect_identical(.pyrun(invert_real(Abs(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {-y, y})")
  expect_identical(.pyrun(invert_real((Integer(2))^(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {log(y)/log(2)})")
  expect_identical(.pyrun(invert_real((Integer(2))^(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, Intersection({log(log(y)/log(2))}, Reals))")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(Integer(2)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {-sqrt(y), sqrt(y)})")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(Rational(1, 2)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {y^2})")
  expect_identical(.pyrun(invert_real(Symbol("x", real = TRUE), Symbol("x", real = TRUE), Symbol("x", real = TRUE))), "ERR ValueError: y should be independent of x ")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(S_Pi), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, {y^(1/pi)})")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(S_Pi), (Integer(-1)*sym_E), Symbol("x", real = TRUE))), "(x, EmptySet)")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(Rational(3, 2)), Integer(1000), Symbol("x", real = TRUE))), "(x, {100})")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE))^(Float("1.00000000000000", precision = 53)), Integer(1), Symbol("x", real = TRUE))), "(x^1.0, {1})")
  expect_identical(.pyrun(invert_real(Integer(1), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "ERR ValueError: Inverse of constant function doesn't exist")
  expect_identical(.pyrun(invert_real((Symbol("x", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(31))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x^31 + x, {y})")
  expect_identical(.pyrun(invert_real(Abs((Integer(1) + Symbol("x", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(31)))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x^31 + x, {-y - 1, y - 1})")
  expect_identical(.pyrun(invert_real(sin(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, y <= 1, Union(ImageSet(Lambda(_n, 2*_n*pi + asin(y)), Integers), ImageSet(Lambda(_n, 2*_n*pi - asin(y) + pi), Integers))))")
  expect_identical(.pyrun(invert_real(sin(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, y <= 1, Union(ImageSet(Lambda(_n, log(2*_n*pi + asin(y))), Integers), ImageSet(Lambda(_n, log(2*_n*pi - asin(y) + pi)), Integers))))")
  expect_identical(.pyrun(invert_real(csc(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, 1 <= y, Union(ImageSet(Lambda(_n, 2*_n*pi + acsc(y)), Integers), ImageSet(Lambda(_n, 2*_n*pi - acsc(y) + pi), Integers))))")
  expect_identical(.pyrun(invert_real(csc(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, 1 <= y, Union(ImageSet(Lambda(_n, log(2*_n*pi + acsc(y))), Integers), ImageSet(Lambda(_n, log(2*_n*pi - acsc(y) + pi)), Integers))))")
  expect_identical(.pyrun(invert_real(cos(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, y <= 1, Union(ImageSet(Lambda(_n, 2*_n*pi - acos(y)), Integers), ImageSet(Lambda(_n, 2*_n*pi + acos(y)), Integers))))")
  expect_identical(.pyrun(invert_real(cos(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, y <= 1, Union(ImageSet(Lambda(_n, log(2*_n*pi - acos(y))), Integers), ImageSet(Lambda(_n, log(2*_n*pi + acos(y))), Integers))))")
  expect_identical(.pyrun(invert_real(sec(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, 1 <= y, Union(ImageSet(Lambda(_n, 2*_n*pi - asec(y)), Integers), ImageSet(Lambda(_n, 2*_n*pi + asec(y)), Integers))))")
  expect_identical(.pyrun(invert_real(sec(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, 1 <= y, Union(ImageSet(Lambda(_n, log(2*_n*pi - asec(y))), Integers), ImageSet(Lambda(_n, log(2*_n*pi + asec(y))), Integers))))")
  expect_identical(.pyrun(invert_real(tan(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ImageSet(Lambda(_n, _n*pi + atan(y)), Integers))")
  expect_identical(.pyrun(invert_real(tan(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ImageSet(Lambda(_n, log(_n*pi + atan(y))), Integers))")
  expect_identical(.pyrun(invert_real(cot(Symbol("x", real = TRUE)), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ImageSet(Lambda(_n, _n*pi + acot(y)), Integers))")
  expect_identical(.pyrun(invert_real(cot(exp(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ImageSet(Lambda(_n, log(_n*pi + acot(y))), Integers))")
  expect_identical(.pyrun(invert_real(tan(tan(Symbol("x", real = TRUE))), Symbol("y", positive = TRUE), Symbol("x", real = TRUE))), "(x, ConditionSet(x, Eq(tan(tan(x)), y), Reals))")
  expect_identical(.pyrun(invert_real((Symbol("x", positive = TRUE))^(S_Pi), Symbol("y", positive = TRUE), Symbol("x", positive = TRUE))), "(x, {y^(1/pi)})")
  expect_identical(.pyrun(invert_real(sinh(Symbol("x", positive = TRUE)), Symbol("r", real = TRUE), Symbol("x", positive = TRUE))), "(x, {asinh(r)})")
  expect_identical(.pyrun(invert_real(sinh(log(Symbol("x", positive = TRUE))), Symbol("p", positive = TRUE), Symbol("x", positive = TRUE))), "(x, {exp(asinh(p))})")
  expect_identical(.pyrun(invert_real(cosh(Symbol("x", positive = TRUE)), Symbol("r", real = TRUE), Symbol("x", positive = TRUE))), "(x, Intersection({-acosh(r), acosh(r)}, Reals))")
  expect_identical(.pyrun(invert_real(cosh(Symbol("x", positive = TRUE)), (Integer(1) + Symbol("p", positive = TRUE)), Symbol("x", positive = TRUE))), "(x, {-acosh(p + 1), acosh(p + 1)})")
  expect_identical(.pyrun(invert_real(tanh(Symbol("x", positive = TRUE)), Symbol("r", real = TRUE), Symbol("x", positive = TRUE))), "(x, Intersection({atanh(r)}, Reals))")
  expect_identical(.pyrun(invert_real(coth(Symbol("x", positive = TRUE)), (Integer(1) + Symbol("p", positive = TRUE)), Symbol("x", positive = TRUE))), "(x, {acoth(p + 1)})")
  expect_identical(.pyrun(invert_real(sech(Symbol("x", positive = TRUE)), Symbol("r", real = TRUE), Symbol("x", positive = TRUE))), "(x, Intersection({-asech(r), asech(r)}, Reals))")
  expect_identical(.pyrun(invert_real(csch(Symbol("x", positive = TRUE)), Symbol("p", positive = TRUE), Symbol("x", positive = TRUE))), "(x, {acsch(p)})")
  expect_identical(.pyrun(invert_real(tanh(sin(Symbol("x", positive = TRUE))), Symbol("r", real = TRUE), Symbol("x", positive = TRUE))), "(x, ConditionSet(x, (-1 <= atanh(r)) & (atanh(r) <= 1), Union(ImageSet(Lambda(_n, 2*_n*pi + asin(atanh(r))), Integers), ImageSet(Lambda(_n, 2*_n*pi - asin(atanh(r)) + pi), Integers))))")
})

test_that("solveset: test_invert_trig_hyp_real matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_invert_trig_hyp_real(cosh(Symbol("x")), Range(Integer(-5), Integer(10), Integer(1)), Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, -acosh(_n)), Range(1, 10, 1)), ImageSet(Lambda(_n, acosh(_n)), Range(1, 10, 1))))")
  expect_identical(.pyrun(symcompR:::.sset_invert_trig_hyp_real(coth(Symbol("x")), Interval(Integer(-3), Integer(2), FALSE, FALSE), Symbol("x"))), "(x, Union(Interval(-oo, -acoth(3)), Interval(acoth(2), oo)))")
  expect_identical(.pyrun(symcompR:::.sset_invert_trig_hyp_real(tanh(Symbol("x")), Interval(Rational(-1, 2), Integer(1), FALSE, FALSE), Symbol("x"))), "(x, Interval(-atanh(1/2), oo))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n"); symcompR:::.sset_invert_trig_hyp_real(sech(Symbol("x")), ImageSet(Lambda(.D1, (Rational(1, 2) + (Rational(1, 3)*.D1))), Naturals0), Symbol("x"))})), "(x, {-asech(1/2), asech(1/2), -asech(5/6), asech(5/6)})")
  expect_identical(.pyrun(symcompR:::.sset_invert_trig_hyp_real(csch(Symbol("x")), Reals, Symbol("x"))), "(x, Union(Interval.open(-oo, 0), Interval.open(0, oo)))")
})

test_that("solveset: test_invert_complex matches SymPy", {
  expect_identical(.pyrun(invert_complex((Integer(3) + Symbol("x")), Symbol("y"), Symbol("x"))), "(x, {y - 3})")
  expect_identical(.pyrun(invert_complex((Integer(3)*Symbol("x")), Symbol("y"), Symbol("x"))), "(x, {y/3})")
  expect_identical(.pyrun(invert_complex(((Integer(-1) + Symbol("x")))^(Integer(3)), Integer(0), Symbol("x"))), "(x, {1})")
  expect_identical(.pyrun(invert_complex(exp(Symbol("x")), Symbol("y"), Symbol("x"))), "(x, ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers))")
  expect_identical(.pyrun(invert_complex(log(Symbol("x")), Symbol("y"), Symbol("x"))), "(x, {exp(y)})")
  expect_identical(.pyrun(invert_real(Integer(1), Symbol("y"), Symbol("x"))), "ERR ValueError: Inverse of constant function doesn't exist")
  expect_identical(.pyrun(invert_complex(Symbol("x"), Symbol("x"), Symbol("x"))), "ERR ValueError: y should be independent of x ")
  expect_identical(.pyrun(invert_complex(Symbol("x"), Symbol("x"), Integer(1))), "ERR ValueError: x must be a symbol")
  expect_identical(.pyrun(invert_complex(sin(Symbol("x")), sym_I, Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, 2*_n*pi + I*log(1 + sqrt(2))), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi - I*log(1 + sqrt(2))), Integers)))")
  expect_identical(.pyrun(invert_complex(cos(Symbol("x")), (Integer(1) + sym_I), Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, 2*_n*pi - acos(1 + I)), Integers), ImageSet(Lambda(_n, 2*_n*pi + acos(1 + I)), Integers)))")
  expect_identical(.pyrun(invert_complex(tan((Integer(2)*Symbol("x"))), Integer(1), Symbol("x"))), "(x, ImageSet(Lambda(_n, _n*pi/2 + pi/8), Integers))")
  expect_identical(.pyrun(invert_complex(cot(Symbol("x")), (Integer(2)*sym_I), Symbol("x"))), "(x, ImageSet(Lambda(_n, _n*pi - I*acoth(2)), Integers))")
  expect_identical(.pyrun(invert_complex(sinh(Symbol("x")), Integer(0), Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + I*pi), Integers)))")
  expect_identical(.pyrun(invert_complex(cosh(Symbol("x")), Integer(0), Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, 2*_n*I*pi + I*pi/2), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + 3*I*pi/2), Integers)))")
  expect_identical(.pyrun(invert_complex(tanh(Symbol("x")), Integer(1), Symbol("x"))), "(x, EmptySet)")
  expect_identical(.pyrun(invert_complex(tanh(Symbol("x")), Symbol("a"), Symbol("x"))), "(x, ConditionSet(x, Ne(a, -1) & Ne(a, 1), ImageSet(Lambda(_n, _n*I*pi + atanh(a)), Integers)))")
  expect_identical(.pyrun(invert_complex(coth(Symbol("x")), Integer(1), Symbol("x"))), "(x, EmptySet)")
  expect_identical(.pyrun(invert_complex(coth(Symbol("x")), Symbol("a"), Symbol("x"))), "(x, ConditionSet(x, Ne(a, -1) & Ne(a, 1), ImageSet(Lambda(_n, _n*I*pi + acoth(a)), Integers)))")
  expect_identical(.pyrun(invert_complex(sech(Symbol("x")), Integer(2), Symbol("x"))), "(x, Union(ImageSet(Lambda(_n, 2*_n*I*pi + I*pi/3), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + 5*I*pi/3), Integers)))")
})

test_that("solveset: test_domain_check matches SymPy", {
  expect_identical(.pyrun(domain_check(((Integer(1) + ((Integer(1) + Symbol("x")))^(Integer(-2))))^(Integer(-1)), Symbol("x"), Integer(-1))), "False")
  expect_identical(.pyrun(domain_check((Symbol("x"))^(Integer(2)), Symbol("x"), Integer(0))), "True")
  expect_identical(.pyrun(domain_check(Symbol("x"), Symbol("x"), oo)), "False")
  expect_identical(.pyrun(domain_check(Integer(0), Symbol("x"), oo)), "False")
})

test_that("solveset: test_issue_11536 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-100) + (Integer(0))^(Symbol("x"))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(-1) + (Integer(0))^(Symbol("x"))), Symbol("x"), Reals)), "{0}")
})

test_that("solveset: test_issue_17479 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-8)*Symbol("x")) + (Integer(4)*Symbol("x")*((Symbol("x"))^(Integer(2)) + (Symbol("y"))^(Integer(2)))) + (Integer(4)*Symbol("x")*((Symbol("x"))^(Integer(2)) + (Symbol("z"))^(Integer(2))))), ((Integer(-4)*Symbol("y")) + (Integer(4)*Symbol("y")*((Symbol("x"))^(Integer(2)) + (Symbol("y"))^(Integer(2))))), ((Integer(-4)*Symbol("z")) + (Integer(4)*Symbol("z")*((Symbol("x"))^(Integer(2)) + (Symbol("z"))^(Integer(2)))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(0, y, 0), (0, y, z), (0, -z, z), (0, z, z), (0, -sqrt(z^2), z), (0, sqrt(z^2), z), (-sqrt(2)*sqrt(2 - y^2)/2, y, 0), (sqrt(2)*sqrt(2 - y^2)/2, y, 0), (-sqrt((-z^3 + z)/z), -z, z), (sqrt((-z^3 + z)/z), -z, z), (-sqrt(-(z^3 - z)/z), z, z), (sqrt(-(z^3 - z)/z), z, z), (-sqrt(-(z - 1)*(z + 1)), -z, z), (-sqrt(-(z - 1)*(z + 1)), z, z), (-sqrt(-(z - 1)*(z + 1)), -sqrt(z^2), z), (-sqrt(-(z - 1)*(z + 1)), sqrt(z^2), z), (sqrt(-(z - 1)*(z + 1)), -z, z), (sqrt(-(z - 1)*(z + 1)), z, z), (sqrt(-(z - 1)*(z + 1)), -sqrt(z^2), z), (sqrt(-(z - 1)*(z + 1)), sqrt(z^2), z)}")
})

test_that("solveset: test_issue_21047 matches SymPy", {
  expect_identical(.pyrun(solveset((((Integer(-1) + ((Integer(-1) + Symbol("x")))^(Rational(1, 2))))^(Integer(6)) + ((Integer(2) + (Integer(-1)*Symbol("x"))))^(Integer(2))), Symbol("x"), Reals)), "{2}")
  expect_identical(.pyrun(solveset(((Integer(2))^(Rational(1, 2)) + ((Integer(1) + (Symbol("x"))^(Rational(1, 2))))^(Integer(2)) + ((Integer(-1) + (Symbol("x"))^(Rational(1, 2))))^(Integer(2)) + (Integer(-2)*(Symbol("x"))^(Integer(2)))), Symbol("x"), Reals)), "{1/2 - sqrt(2*sqrt(2) + 5)/2, 1/2 + sqrt(2*sqrt(2) + 5)/2}")
})

test_that("solveset: test_is_function_class_equation matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", tan(Symbol("x")), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (Integer(-1) + tan(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (sin(Symbol("x")) + tan(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", ((Integer(-1)*Symbol("a")) + sin(Symbol("x")) + tan(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", ((sin(Symbol("x"))*tan(Symbol("x"))) + sin(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", ((sin(Symbol("x"))*tan((Symbol("a") + Symbol("x")))) + sin(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", ((sin(Symbol("x"))*tan((Symbol("a")*Symbol("x")))) + sin(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (Integer(-1) + (Symbol("a")*tan(Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (Integer(-1) + (tan(Symbol("x")))^(Integer(2)) + sin(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (Symbol("x") + tan(Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", tan((Symbol("x"))^(Integer(2))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (sin(Symbol("x")) + tan((Symbol("x"))^(Integer(2)))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (tan(Symbol("x")))^(sin(Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("TrigonometricFunction", (sin(Symbol("x")) + tan(sin(Symbol("x")))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", tanh(Symbol("x")), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (Integer(-1) + tanh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (sinh(Symbol("x")) + tanh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", ((Integer(-1)*Symbol("a")) + sinh(Symbol("x")) + tanh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", ((sinh(Symbol("x"))*tanh(Symbol("x"))) + sinh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", ((sinh(Symbol("x"))*tanh((Symbol("a") + Symbol("x")))) + sinh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", ((sinh(Symbol("x"))*tanh((Symbol("a")*Symbol("x")))) + sinh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (Integer(-1) + (Symbol("a")*tanh(Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (Integer(-1) + (tanh(Symbol("x")))^(Integer(2)) + sinh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (Symbol("x") + tanh(Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", tanh((Symbol("x"))^(Integer(2))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (sinh(Symbol("x")) + tanh((Symbol("x"))^(Integer(2)))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (tanh(Symbol("x")))^(sinh(Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_function_class_equation("HyperbolicFunction", (sinh(Symbol("x")) + tanh(sinh(Symbol("x")))), Symbol("x"))), "False")
})

test_that("solveset: test_garbage_input matches SymPy", {
  expect_identical(.pyrun(solveset_real(list(Symbol("y")), Symbol("y"))), "ERR ValueError: [y] is not a valid SymPy expression")
  expect_identical(.pyrun(solveset_real(Symbol("x", real = TRUE), Integer(1))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(-1) + Symbol("x", real = TRUE)), Integer(1))), "{x}")
  expect_identical(.pyrun(solveset_real(Symbol("x", real = TRUE), S_Pi)), "EmptySet")
  expect_identical(.pyrun(solveset_real(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(list(Symbol("x", real = TRUE)), Symbol("x", real = TRUE))), "ERR ValueError: [x] is not a valid SymPy expression")
  expect_identical(.pyrun(solveset_complex(Symbol("x", real = TRUE), S_Pi)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(1) + Symbol("x", real = TRUE)), Reals)), "ERR ValueError: Reals is not a valid SymPy symbol")
  expect_identical(.pyrun(solveset((Integer(1) + Symbol("x", real = TRUE)), Symbol("x", real = TRUE), Integer(2))), "ERR ValueError: 2 is not a valid domain")
})

test_that("solveset: test_solve_mul matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Integer(-3) + exp(Symbol("x")))*(Symbol("b") + (Symbol("a")*Symbol("x")))), Symbol("x"))), "Union({log(3)}, Intersection({-b/a}, Reals))")
  expect_identical(.pyrun(solveset_real(((Integer(-3) + exp(Symbol("x")))*(Symbol("bb", real = TRUE) + (Symbol("x")*Symbol("anz", nonzero = TRUE)))), Symbol("x"))), "{-bb/anz, log(3)}")
  expect_identical(.pyrun(solveset_real(((Integer(8) + (Integer(2)*Symbol("x")))*(Integer(8) + exp(Symbol("x")))), Symbol("x"))), "{-4}")
  expect_identical(.pyrun(solveset_real((Symbol("x")*(log(Symbol("x")))^(Integer(-1))), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_solve_invert matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-3) + exp(Symbol("x"))), Symbol("x"))), "{log(3)}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + log(Symbol("x"))), Symbol("x"))), "{exp(3)}")
  expect_identical(.pyrun(solveset_real((Integer(3))^((Integer(2) + Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(3))^((Integer(2) + (Integer(-1)*Symbol("x")))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Symbol("y") + (Integer(-1)*Symbol("b")*exp((Symbol("a")*(Symbol("x"))^(Integer(-1)))))), Symbol("x"))), "Intersection({a/log(y/b)}, Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-10) + (Integer(2))^(Symbol("x"))), Symbol("x"))), "{1 + log(5)/log(2)}")
})

test_that("solveset: test_issue_25768 matches SymPy", {
  expect_identical(.pyrun(solveset_real((Rational(-1, 2) + sin(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers))")
  expect_identical(.pyrun(solveset_real((Float("-0.500000000000000", precision = 53) + sin(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi + 0.523598775598299), Integers), ImageSet(Lambda(_n, 2*_n*pi - 0.523598775598299 + pi), Integers))")
})

test_that("solveset: test_errorinverses matches SymPy", {
  expect_identical(.pyrun(solveset_real((Rational(-1, 2) + erf(Symbol("x"))), Symbol("x"))), "{erfinv(1/2)}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + erfinv(Symbol("x"))), Symbol("x"))), "{erf(2)}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + erfc(Symbol("x"))), Symbol("x"))), "{0}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + erfcinv(Symbol("x"))), Symbol("x"))), "{erfc(2)}")
})

test_that("solveset: test_solve_polynomial matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Integer(3)*Symbol("x", real = TRUE))), Symbol("x", real = TRUE))), "{2/3}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Symbol("x", real = TRUE))^(Integer(2))), Symbol("x", real = TRUE))), "{-1, 1}")
  expect_identical(.pyrun(solveset_real((Symbol("x", real = TRUE) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(3)))), Symbol("x", real = TRUE))), "{y^3}")
  expect_identical(.pyrun(solveset_real((Integer(-4) + (Symbol("x", real = TRUE))^(Integer(3)) + (Integer(-15)*Symbol("x", real = TRUE))), Symbol("x", real = TRUE))), "{4, -2 - sqrt(3), -2 + sqrt(3)}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Symbol("x", real = TRUE))^(Rational(1, 2))), Symbol("x", real = TRUE))), "{1}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Symbol("x", real = TRUE))^(Rational(1, 2))), Symbol("x", real = TRUE))), "{4}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Symbol("x", real = TRUE))^(Rational(1, 4))), Symbol("x", real = TRUE))), "{16}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + (Symbol("x", real = TRUE))^(Rational(1, 3))), Symbol("x", real = TRUE))), "{27}")
  expect_identical(.pyrun(solveset_real((Integer(1) + (Symbol("x", real = TRUE))^(Integer(3)) + (Symbol("x", real = TRUE))^(Integer(5))), Symbol("x", real = TRUE))), "{CRootOf(x^5 + x^3 + 1, 0)}")
  expect_identical(.pyrun(solveset_real((Integer(6) + (Integer(-2)*Symbol("x", real = TRUE)) + (Integer(-2)*(Symbol("x", real = TRUE))^(Integer(3))) + (Integer(4)*(Symbol("x", real = TRUE))^(Integer(2)))), Symbol("x", real = TRUE))), "{1/(9*(sqrt(77)/6 + 79/54)^(1/3)) + 2/3 + (sqrt(77)/6 + 79/54)^(1/3)}")
})

test_that("solveset: test_return_root_of matches SymPy", {
  expect_identical(.pyrun(solveset_complex((Integer(20) + (Symbol("x"))^(Integer(5)) + (Integer(-15)*(Symbol("x"))^(Integer(3))) + (Integer(-5)*(Symbol("x"))^(Integer(2))) + (Integer(10)*Symbol("x"))), Symbol("x"))), "{CRootOf(x^5 - 15*x^3 - 5*x^2 + 10*x + 20, 0), CRootOf(x^5 - 15*x^3 - 5*x^2 + 10*x + 20, 1), CRootOf(x^5 - 15*x^3 - 5*x^2 + 10*x + 20, 2), CRootOf(x^5 - 15*x^3 - 5*x^2 + 10*x + 20, 3), CRootOf(x^5 - 15*x^3 - 5*x^2 + 10*x + 20, 4)}")
  expect_identical(.pyrun(solveset_complex((Integer(7) + (Symbol("x"))^(Integer(5)) + (Integer(3)*(Symbol("x"))^(Integer(3)))), Symbol("x"))), "{CRootOf(x^5 + 3*x^3 + 7, 0), CRootOf(x^5 + 3*x^3 + 7, 1), CRootOf(x^5 + 3*x^3 + 7, 2), CRootOf(x^5 + 3*x^3 + 7, 3), CRootOf(x^5 + 3*x^3 + 7, 4)}")
  expect_identical(.pyrun(solveset_complex((Integer(2) + (Symbol("x"))^(Integer(6)) + (Integer(-2)*Symbol("x"))), Symbol("x"))), "{CRootOf(x^6 - 2*x + 2, 0), CRootOf(x^6 - 2*x + 2, 1), CRootOf(x^6 - 2*x + 2, 2), CRootOf(x^6 - 2*x + 2, 3), CRootOf(x^6 - 2*x + 2, 4), CRootOf(x^6 - 2*x + 2, 5)}")
  expect_identical(.pyrun(solveset_complex((Rational(7, 4) + (Symbol("x"))^(Integer(5)) + (Integer(3)*(Symbol("x"))^(Integer(2))) + (Integer(4)*(Symbol("x"))^(Integer(3)))), Symbol("x"))), "{CRootOf(4*x^5 + 16*x^3 + 12*x^2 + 7, 0), CRootOf(4*x^5 + 16*x^3 + 12*x^2 + 7, 1), CRootOf(4*x^5 + 16*x^3 + 12*x^2 + 7, 2), CRootOf(4*x^5 + 16*x^3 + 12*x^2 + 7, 3), CRootOf(4*x^5 + 16*x^3 + 12*x^2 + 7, 4)}")
  expect_identical(.pyrun(solveset_complex((Symbol("x")*((Integer(-1) + Symbol("x")))^(Integer(2))*(Integer(1) + Symbol("x"))*(Integer(1) + (Symbol("x"))^(Integer(6)) + (Integer(-1)*Symbol("x")))), Symbol("x"))), "{-1, 0, 1, CRootOf(x^6 - x + 1, 0), CRootOf(x^6 - x + 1, 1), CRootOf(x^6 - x + 1, 2), CRootOf(x^6 - x + 1, 3), CRootOf(x^6 - x + 1, 4), CRootOf(x^6 - x + 1, 5)}")
})

test_that("solveset: test_solveset_sqrt_1 matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-2) + ((Integer(6) + (Integer(5)*Symbol("x"))))^(Rational(1, 2)) + (Integer(-1)*Symbol("x"))), Symbol("x"))), "{-1, 2}")
  expect_identical(.pyrun(solveset_real((Integer(7) + ((Integer(-1) + Symbol("x")))^(Rational(1, 2)) + (Integer(-1)*Symbol("x"))), Symbol("x"))), "{10}")
  expect_identical(.pyrun(solveset_real((Integer(-5) + ((Integer(-2) + Symbol("x")))^(Rational(1, 2))), Symbol("x"))), "{27}")
  expect_identical(.pyrun(solveset_real((Integer(-7) + (Symbol("x"))^(Rational(1, 2))), Symbol("x"))), "{49}")
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Integer(3)))^(Rational(1, 2)), Symbol("x"))), "{0}")
  expect_identical(.pyrun(solveset_real(((Integer(-1) + Symbol("x")))^(Rational(1, 2)), Symbol("x"))), "{1}")
  expect_identical(.pyrun(solveset_real((((Symbol("x"))^(Integer(-1))*(Integer(-3) + Symbol("x"))))^(Rational(1, 2)), Symbol("x"))), "{3}")
  expect_identical(.pyrun(solveset_real((Rational(-1, 2) + (((Symbol("x"))^(Integer(-1))*(Integer(-3) + Symbol("x"))))^(Rational(1, 2))), Symbol("x"))), "{4}")
})

test_that("solveset: test_solveset_sqrt_2 matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-2) + ((Integer(-1) + (Integer(2)*Symbol("x", real = TRUE))))^(Rational(1, 2)) + (Integer(-1)*((Integer(-4) + Symbol("x", real = TRUE)))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{5, 13}")
  expect_identical(.pyrun(solveset_real((Integer(2) + ((Integer(7) + Symbol("x", real = TRUE)))^(Rational(1, 2)) + (Integer(-1)*((Integer(3) + (Integer(-1)*Symbol("x", real = TRUE))))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{-6}")
  expect_identical(.pyrun(solveset_real((Integer(-7) + (((Integer(-1)*((Integer(-5) + (Symbol("x", real = TRUE))^(Integer(2))))^(Rational(1, 2))) + (Integer(17)*Symbol("x", real = TRUE))))^(Rational(1, 2))), Symbol("x", real = TRUE))), "{3}")
  expect_identical(.pyrun(solveset_real((Integer(1) + Symbol("x", real = TRUE) + (Integer(-1)*(((Symbol("x", real = TRUE))^(Integer(4)) + (Integer(-1)*Symbol("x", real = TRUE)) + (Integer(4)*(Symbol("x", real = TRUE))^(Integer(3)))))^(Rational(1, 4)))), Symbol("x", real = TRUE))), "{-1/2, -1/3}")
  expect_identical(.pyrun(solveset_real((((Integer(9) + (Integer(2)*Symbol("x", real = TRUE))))^(Rational(1, 2)) + (Integer(-1)*((Integer(1) + Symbol("x", real = TRUE)))^(Rational(1, 2))) + (Integer(-1)*((Integer(4) + Symbol("x", real = TRUE)))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{0}")
  expect_identical(.pyrun(solveset_real((((Integer(-1) + (Integer(2)*Symbol("x", real = TRUE))))^(Rational(1, 2)) + ((Integer(4) + Symbol("x", real = TRUE)))^(Rational(1, 2)) + (Integer(-3)*((Integer(-1) + Symbol("x", real = TRUE)))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{5}")
  expect_identical(.pyrun(solveset_real((Integer(-12) + ((Symbol("x", real = TRUE))^(Rational(1, 2))*((Integer(-7) + Symbol("x", real = TRUE)))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{16}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + (Symbol("x", real = TRUE))^(Rational(1, 2)) + ((Integer(-3) + Symbol("x", real = TRUE)))^(Rational(1, 2))), Symbol("x", real = TRUE))), "{4}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + Symbol("x", real = TRUE) + ((Integer(-7) + (Integer(2)*(Symbol("x", real = TRUE))^(Integer(2)))))^(Rational(1, 2))), Symbol("x", real = TRUE))), "{-8, 2}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + ((Integer(4) + (Integer(9)*(Symbol("x", real = TRUE))^(Integer(2)))))^(Rational(1, 2)) + (Integer(-3)*Symbol("x", real = TRUE))), Symbol("x", real = TRUE))), "{0}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + ((Integer(-3) + Symbol("x", real = TRUE)))^(Rational(1, 2)) + (Integer(-1)*(Symbol("x", real = TRUE))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(-3) + ((Integer(-5) + (Integer(2)*Symbol("x", real = TRUE))))^(Rational(1, 3))), Symbol("x", real = TRUE))), "{16}")
  expect_identical(.pyrun(solveset_real((Integer(-4) + (Symbol("x", real = TRUE))^(Rational(1, 2)) + (Symbol("x", real = TRUE))^(Rational(1, 4))), Symbol("x", real = TRUE))), "{(-1/2 + sqrt(17)/2)^4}")
  expect_identical(.pyrun(solveset_real(((Symbol("x", real = TRUE))^(Rational(1, 2)) + (Symbol("x", real = TRUE))^(Rational(1, 4)) + (Integer(-1)*((Integer(-1) + Symbol("x", real = TRUE)))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "EmptySet")
  expect_identical(.pyrun(solveset_real((((Integer(-4) + Symbol("x", real = TRUE)))^(Integer(2)) + ((Integer(-2) + (Symbol("x", real = TRUE))^(Rational(1, 2))))^(Integer(4))), Symbol("x", real = TRUE))), "{-4, 4}")
  expect_identical(.pyrun(solveset_real(((Symbol("x", real = TRUE))^(Rational(1, 2)) + ((Integer(1) + Symbol("x", real = TRUE)))^(Rational(1, 2)) + ((Integer(1) + (Integer(-1)*Symbol("x", real = TRUE))))^(Rational(1, 2)) + (Rational(-6, 5)*(Integer(5))^(Rational(1, 2)))), Symbol("x", real = TRUE))), "{4/5, -1484/375 + 172564/(140625*(114*sqrt(12657)/78125 + 12459439/52734375)^(1/3)) + 4*(114*sqrt(12657)/78125 + 12459439/52734375)^(1/3)}")
  expect_identical(.pyrun(solveset_real(((Symbol("x", real = TRUE))^(Rational(1, 2)) + (Symbol("x", real = TRUE))^(Rational(1, 3)) + (Symbol("x", real = TRUE))^(Rational(1, 4))), Symbol("x", real = TRUE))), "{0}")
  expect_identical(.pyrun(solveset_real((Symbol("x", real = TRUE)*((Integer(1) + (Symbol("x", real = TRUE))^(Integer(2))))^(Rational(-1, 2))), Symbol("x", real = TRUE))), "{0}")
  expect_identical(.pyrun(solveset_real(((Symbol("y", real = TRUE))^(Integer(-2))*((Integer(1) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2)))))^(Rational(-1, 2))*(Symbol("x", real = TRUE) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(3))))), Symbol("x", real = TRUE))), "{y^3}")
  expect_identical(.pyrun(solveset_real((Integer(-9) + ((Integer(5) + Symbol("x", real = TRUE)))^(Rational(-1, 5))), Symbol("x", real = TRUE))), "{-295244/59049}")
})

test_that("solveset: test_solve_sqrt_fail matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(1) + (((Symbol("x"))^(Integer(3)) + (Integer(-3)*(Symbol("x"))^(Integer(2)))))^(Rational(1, 3)) + (Integer(-1)*Symbol("x"))), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_solve_sqrt_3 matches SymPy", {
  expect_identical(.pyrun(solveset_real((((((Rational(-1, 4) + (Rational(-1, 2)*(Symbol("m"))^(Integer(2))) + (Rational(-1, 4)*((Integer(1) + (Integer(-4)*(Symbol("m"))^(Integer(2))) + (Integer(4)*(Symbol("m"))^(Integer(4))) + (Integer(8)*Symbol("m"))))^(Rational(1, 2)))))^(Integer(2)) + ((Rational(-1, 4) + (Rational(1, 2)*(Symbol("m"))^(Integer(2))) + (Integer(-1)*Symbol("m")) + (Rational(-1, 4)*((Integer(1) + (Integer(-4)*(Symbol("m"))^(Integer(2))) + (Integer(4)*(Symbol("m"))^(Integer(4))) + (Integer(8)*Symbol("m"))))^(Rational(1, 2)))))^(Integer(2))))^(Rational(1, 2)) + (Integer(-1)*((((Rational(1, 2) + (Rational(-1, 2)*Symbol("m")*(Symbol("q"))^(Integer(-1)))))^(Integer(2)) + ((Symbol("m") + (Integer(-1)*Symbol("q"))))^(Integer(2))))^(Rational(1, 2)))), Symbol("q"))), "ConditionSet(q, Eq(sqrt((m - q)^2 + (-m/(2*q) + 1/2)^2) - sqrt((-m^2/2 - sqrt(4*m^4 - 4*m^2 + 8*m + 1)/4 - 1/4)^2 + (m^2/2 - m - sqrt(4*m^4 - 4*m^2 + 8*m + 1)/4 - 1/4)^2), 0), Reals)")
})

test_that("solveset: test_solve_polynomial_symbolic_param matches SymPy", {
  expect_identical(.pyrun(solveset_complex((((Integer(-1) + (Symbol("x"))^(Integer(2))))^(Integer(2)) + (Integer(-1)*Symbol("a"))), Symbol("x"))), "{-sqrt(1 - sqrt(a)), sqrt(1 - sqrt(a)), -sqrt(sqrt(a) + 1), sqrt(sqrt(a) + 1)}")
  expect_identical(.pyrun(solveset_complex((Symbol("y") + (Integer(-1)*Symbol("b")*((Integer(1) + (Symbol("a")*Symbol("x"))))^(Integer(-1)))), Symbol("x"))), "Complement({(b/y - 1)/a}, {-1/a})")
  expect_identical(.pyrun(solveset_complex((Symbol("y") + (Integer(-1)*Symbol("b")*Symbol("x")*((Symbol("a") + Symbol("x")))^(Integer(-1)))), Symbol("x"))), "Complement({-a*y/(-b + y)}, {-a})")
})

test_that("solveset: test_solve_rational matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(1) + (Symbol("x"))^(Integer(-1))), Symbol("x"))), "{-1}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + exp((Integer(-1)*Symbol("x")))), Symbol("x"))), "{0}")
  expect_identical(.pyrun(solveset_real((Symbol("x")*(Integer(1) + (Integer(-5)*(Symbol("x"))^(Integer(-1))))), Symbol("x"))), "{5}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Integer(2)*Symbol("x")*((Integer(2) + Symbol("x")))^(Integer(-1)))), Symbol("x"))), "{2}")
  expect_identical(.pyrun(solveset_real((((Symbol("x"))^(Integer(2))*((Integer(7) + (Integer(-1)*Symbol("x"))))^(Integer(-2))) + (Integer(2)*Symbol("x")*((Integer(7) + (Integer(-1)*Symbol("x"))))^(Integer(-1)))), Symbol("x"))), "{0, 14}")
})

test_that("solveset: test_solveset_real_gen_is_pow matches SymPy", {
  expect_identical(.pyrun(solveset_real(Integer(2), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_no_sol matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + ((-oo)*Symbol("x"))))), "EmptySet")
  expect_identical(.pyrun(solveset((oo*Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset(((-oo) + (oo*Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(Integer(4), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(exp(Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(1) + (Symbol("x"))^(Integer(2))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(-3)*Symbol("a")*(Symbol("x"))^(Rational(-1, 2))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Symbol("x"))^(Integer(-1)), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((((Integer(2) + Symbol("x")))^(Integer(-1)) + (((Integer(2) + Symbol("x")))^(Integer(-2))*(Integer(-1) + (Integer(-1)*Symbol("x"))))), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_sol_zero_real matches SymPy", {
  expect_identical(.pyrun(solveset_real(Integer(0), Symbol("x"))), "Reals")
  expect_identical(.pyrun(solveset(Integer(0), Symbol("x"), Interval(Integer(1), Integer(2), FALSE, FALSE))), "Interval(1, 2)")
  expect_identical(.pyrun(solveset_real((Integer(-1) + ((Integer(1) + Symbol("x")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x"))), Symbol("x"))), "Reals")
})

test_that("solveset: test_no_sol_rational_extragenous matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Integer(3) + (Symbol("x")*((Integer(1) + Symbol("x")))^(Integer(-1)))))^(Integer(-2)), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((((Integer(1) + ((Integer(-1) + Symbol("x")))^(Integer(-1))))^(Integer(-1))*(Integer(-1) + Symbol("x"))), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_solve_polynomial_cv_1a matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Symbol("x"))^(Rational(1, 2))), Symbol("x"))), "{1}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Symbol("x"))^(Rational(1, 2))), Symbol("x"))), "{4}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Symbol("x"))^(Rational(1, 4))), Symbol("x"))), "{16}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + (Symbol("x"))^(Rational(1, 3))), Symbol("x"))), "{27}")
  expect_identical(.pyrun(solveset_real((Symbol("x")*(Integer(-3) + (Symbol("x"))^(Rational(1, 3)))), Symbol("x"))), "{0, 27}")
})

test_that("solveset: test_solveset_real_rational matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Integer(2)*Symbol("x", real = TRUE)*((Integer(2) + Symbol("x", real = TRUE)))^(Integer(-1)))), Symbol("x", real = TRUE))), "{2}")
})

test_that("solveset: test_solveset_real_log matches SymPy", {
  expect_identical(.pyrun(solveset_real(log(((Integer(1) + Symbol("x"))*(Integer(-1) + Symbol("x")))), Symbol("x"))), "{-sqrt(2), sqrt(2)}")
})

test_that("solveset: test_poly_gens matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-8) + (Integer(4))^(((Integer(2)*Symbol("x")) + (Integer(2)*(Symbol("x"))^(Integer(2)))))), Symbol("x"))), "{-3/2, 1/2}")
})

test_that("solveset: test_solve_abs matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-1) + Abs(Symbol("x"))), Symbol("x"))), "ERR ValueError: \\nAbsolute values cannot be inverted in the complex domain.")
  expect_identical(.pyrun(local({.D1 <- Dummy("n"); solveset(((Integer(-1)*.D1) + Abs(Symbol("x"))), Symbol("x"), Reals)})), "ConditionSet(x, Contains(_n, Interval(0, oo)), {-_n, _n})")
  expect_identical(.pyrun(solveset_real((Integer(-2) + Abs(Symbol("x"))), Symbol("x"))), "{-2, 2}")
  expect_identical(.pyrun(solveset_real((Integer(2) + Abs(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(((Integer(-2)*Abs((Integer(-3) + Symbol("x")))) + Abs((Integer(3) + Symbol("x")))), Symbol("x"))), "{1, 9}")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Abs((Integer(-1) + Symbol("x")))) + (Integer(2)*Abs(Symbol("x")))), Symbol("x"))), "{-1, 1/3}")
  expect_identical(.pyrun(invert_real(((Integer(-1)*Symbol("b")) + Abs((Symbol("a") + (Integer(-1)*Abs((Integer(3) + Symbol("x"))))))), Integer(0), Symbol("x"))), "(x, ConditionSet(x, Contains(b, Interval(0, oo)) & Contains(a + b, Interval(0, oo)) & Contains(a - b, Interval(0, oo)), {-a - b - 3, -a + b - 3, a - b - 3, a + b - 3}))")
  expect_identical(.pyrun(solveset(Eq(sin(Abs(Symbol("x"))), Integer(1)), Symbol("x"), domain = Reals)), "Union(ImageSet(Lambda(_n, -2*_n*pi - pi/2), Range(0, oo, 1)), ImageSet(Lambda(_n, 2*_n*pi + pi/2), Range(0, oo, 1)))")
})

test_that("solveset: test_issue_9824 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + (sin(Symbol("x")))^(Integer(2)) + (Integer(-2)*sin(Symbol("x")))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)")
  expect_identical(.pyrun(solveset((Integer(1) + (cos(Symbol("x")))^(Integer(2)) + (Integer(-2)*cos(Symbol("x")))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*pi), Integers)")
})

test_that("solveset: test_issue_9565 matches SymPy", {
  expect_identical(.pyrun(solveset_real(Le(Abs((((Integer(-5) + Symbol("x")))^(Integer(-1))*(Integer(-1) + Symbol("x")))), Rational(1, 3)), Symbol("x"))), "Interval(-1, 2)")
})

test_that("solveset: test_issue_10069 matches SymPy", {
  expect_identical(.pyrun(solveset_real(Gt((Integer(-1) + (Abs((Integer(-1) + Symbol("x"))))^(Integer(-1))), Integer(0)), Symbol("x"))), "Union(Interval.open(0, 1), Interval.open(1, 2))")
})

test_that("solveset: test_real_imag_splitting matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-3) + (((Symbol("a", real = TRUE))^(Integer(2)) + (Integer(-1)*(Symbol("b", real = TRUE))^(Integer(2)))))^(Rational(1, 2))), Symbol("a", real = TRUE))), "{-sqrt(b^2 + 9), sqrt(b^2 + 9)}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + (((Symbol("a", real = TRUE))^(Integer(2)) + (Symbol("b", real = TRUE))^(Integer(2))))^(Rational(1, 2))), Symbol("a", real = TRUE))), "Intersection({-sqrt(9 - b^2), sqrt(9 - b^2)}, Reals)")
})

test_that("solveset: test_solve_only_exp_1 matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Symbol("y", positive = TRUE)) + exp(Symbol("x"))), Symbol("x"))), "{log(y)}")
  expect_identical(.pyrun(solveset_real((Integer(-4) + exp(Symbol("x")) + exp((Integer(-1)*Symbol("x")))), Symbol("x"))), "{log(2 - sqrt(3)), log(sqrt(3) + 2)}")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Symbol("y", positive = TRUE)) + exp(Symbol("x")) + exp((Integer(-1)*Symbol("x")))), Symbol("x"))), "Union(Intersection({log(y/2 - sqrt(y - 2)*sqrt(y + 2)/2)}, Reals), Intersection({log(y/2 + sqrt(y - 2)*sqrt(y + 2)/2)}, Reals))")
})

test_that("solveset: test_atan2 matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Rational(-1, 3)*S_Pi) + sym_atan2(Symbol("x"), Integer(2))), Symbol("x"))), "{2*sqrt(3)}")
})

test_that("solveset: test_piecewise_solveset matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-3) + Piecewise(list((Integer(-2) + Symbol("x")), Gt(Symbol("x"), Integer(2))), list((Integer(2) + (Integer(-1)*Symbol("x"))), S_true))), Symbol("x"))), "{-1, 5}")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Symbol("y", positive = TRUE)) + Piecewise(list((Integer(-3) + Symbol("x")), Ge(Symbol("x"), Integer(3))), list((Integer(3) + (Integer(-1)*Symbol("x"))), S_true))), Symbol("x"))), "{3 - y, y + 3}")
  expect_identical(.pyrun(solveset(Piecewise(list(((Integer(-2) + Symbol("x")))^(Integer(2)), Ge(Symbol("x"), Integer(0))), list(Integer(0), S_true)), Symbol("x"), domain = Reals)), "Union({2}, Interval.open(-oo, 0))")
  expect_identical(.pyrun(solveset(((Integer(-1)*sym_I) + Piecewise(list((Integer(1) + Symbol("x")), Gt(Symbol("x"), Integer(0))), list(sym_I, S_true))), Symbol("x"), Reals)), "Interval(-oo, 0)")
  expect_identical(.pyrun(solveset(Piecewise(list((Integer(-1) + Symbol("x")), Ne(Symbol("x"), sym_I)), list(Symbol("x"), S_true)), Symbol("x"))), "{1}")
  expect_identical(.pyrun(solveset(Gt(Piecewise(list(Integer(1), Gt(Symbol("x"), Integer(10))), list(Integer(0), S_true)), Integer(0)), Symbol("x"), Reals)), "Interval.open(10, oo)")
  expect_identical(.pyrun(solveset(S_true, Symbol("x"), domain = Interval(Integer(-3), Integer(10), FALSE, FALSE))), "Interval(-3, 10)")
  expect_identical(.pyrun(solveset(Piecewise(list(Integer(0), Eq(Symbol("x"), Integer(0))), list(((Symbol("x"))^(Integer(2))*(Abs(Symbol("x")))^(Integer(-1))), S_true)), Symbol("x"), domain = Reals)), "{0}")
  expect_identical(.pyrun(solveset(Piecewise(list(Integer(0), Eq(Symbol("x"), S_Pi)), list(((sin(Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*S_Pi))), S_true)))), "{pi}")
})

test_that("solveset: test_solveset_complex_polynomial matches SymPy", {
  expect_identical(.pyrun(solveset_complex((Symbol("c") + (Symbol("a")*(Symbol("x"))^(Integer(2))) + (Symbol("b")*Symbol("x"))), Symbol("x"))), "{-b/(2*a) - sqrt(-4*a*c + b^2)/(2*a), -b/(2*a) + sqrt(-4*a*c + b^2)/(2*a)}")
  expect_identical(.pyrun(solveset_complex((Symbol("x") + (Integer(-1)*(Symbol("y"))^(Integer(3)))), Symbol("y"))), "{x^(1/3), -x^(1/3)/2 - sqrt(3)*I*x^(1/3)/2, -x^(1/3)/2 + sqrt(3)*I*x^(1/3)/2}")
  expect_identical(.pyrun(solveset_complex((Integer(-1) + Symbol("x") + (Symbol("x"))^(Integer(-1))), Symbol("x"))), "{1/2 - sqrt(3)*I/2, 1/2 + sqrt(3)*I/2}")
})

test_that("solveset: test_sol_zero_complex matches SymPy", {
  expect_identical(.pyrun(solveset_complex(Integer(0), Symbol("x"))), "Complexes")
})

test_that("solveset: test_solveset_complex_rational matches SymPy", {
  expect_identical(.pyrun(solveset_complex((((Integer(-3) + Symbol("x")))^(Integer(-1))*(Integer(-1) + Symbol("x"))*(Symbol("x") + (Integer(-1)*sym_I))), Symbol("x"))), "{1, I}")
  expect_identical(.pyrun(solveset_complex(((Symbol("y"))^(Integer(-2))*((Integer(1) + (Integer(-1)*(Symbol("y"))^(Integer(2)))))^(Rational(-1, 2))*(Symbol("x") + (Integer(-1)*(Symbol("y"))^(Integer(3))))), Symbol("x"))), "{y^3}")
  expect_identical(.pyrun(solveset_complex(((Integer(-1)*sym_I) + (Integer(-1)*(Symbol("x"))^(Integer(2)))), Symbol("x"))), "{-sqrt(2)/2 + sqrt(2)*I/2, sqrt(2)/2 - sqrt(2)*I/2}")
})

test_that("solveset: test_solveset_complex_exp matches SymPy", {
  expect_identical(.pyrun(solveset_complex((Integer(-1) + exp(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*I*pi), Integers)")
  expect_identical(.pyrun(solveset_complex(((Integer(-1)*sym_I) + exp(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, I*(2*_n*pi + pi/2)), Integers)")
  expect_identical(.pyrun(solveset_complex(exp((Integer(-1)*Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(((Rational(1, 2)*exp(Symbol("x"))) + (Rational(-1, 2)*exp((Integer(-1)*Symbol("x"))))), Symbol("x"))), "ImageSet(Lambda(_n, _n*I*pi), Integers)")
})

test_that("solveset: test_solveset_real_exp matches SymPy", {
  expect_identical(.pyrun(solveset(Eq((Integer(-2))^(Symbol("x")), Integer(4)), Symbol("x"), Reals)), "{2}")
  expect_identical(.pyrun(solveset(Eq((Integer(-1)*(Integer(2))^(Symbol("x"))), Integer(4)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Eq((Integer(-3))^(Symbol("x")), Integer(27)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Eq((Integer(-5))^((Integer(1) + Symbol("x"))), Integer(625)), Symbol("x"), Reals)), "{3}")
  expect_identical(.pyrun(solveset(Eq((Integer(2))^((Integer(-3) + Symbol("x"))), Integer(-16)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Eq((Integer(-3))^((Integer(-3) + Symbol("x"))), Integer("-4052555153018976267")), Symbol("x"), Reals)), "{42}")
  expect_identical(.pyrun(solveset(Eq((Integer(2))^(Symbol("x")), Symbol("y")), Symbol("x"), Reals)), "Intersection({log(y)/log(2)}, Reals)")
  expect_identical(.pyrun(invert_real((Integer(-16) + (Integer(-2))^((Integer(2)*Symbol("x")))), Integer(0), Symbol("x"))), "(x, {2})")
})

test_that("solveset: test_solve_complex_log matches SymPy", {
  expect_identical(.pyrun(solveset_complex(log(Symbol("x")), Symbol("x"))), "{1}")
  expect_identical(.pyrun(solveset_complex((Integer(1) + (Integer(-1)*log((Symbol("a") + (Integer(4)*(Symbol("x"))^(Integer(2))))))), Symbol("x"))), "{-sqrt(E - a)/2, sqrt(E - a)/2}")
})

test_that("solveset: test_solve_complex_sqrt matches SymPy", {
  expect_identical(.pyrun(solveset_complex((Integer(-2) + ((Integer(6) + (Integer(5)*Symbol("x"))))^(Rational(1, 2)) + (Integer(-1)*Symbol("x"))), Symbol("x"))), "{-1, 2}")
  expect_identical(.pyrun(solveset_complex((Integer(-2) + ((Integer(6) + (Integer(5)*Symbol("x"))))^(Rational(1, 2)) + (Integer(-1)*Symbol("x")) + (Integer(-2)*sym_I)), Symbol("x"))), "{-2, 3 - 4*I}")
  expect_identical(.pyrun(solveset_complex((Integer(4)*Symbol("x")*(Integer(1) + (Integer(-1)*Symbol("a")*(Symbol("x"))^(Rational(1, 2))))), Symbol("x"))), "{0, a^(-2)}")
})

test_that("solveset: test_solveset_complex_tan matches SymPy", {
  expect_identical(.pyrun(solveset_complex((sym_I*((exp((sym_I*Symbol("x"))) + exp((Integer(-1)*sym_I*Symbol("x")))))^(Integer(-1))*((Integer(-1)*exp((sym_I*Symbol("x")))) + exp((Integer(-1)*sym_I*Symbol("x"))))), Symbol("x"))), "Complement(ImageSet(Lambda(_n, _n*pi), Integers), ImageSet(Lambda(_n, _n*pi + pi/2), Integers))")
})

test_that("solveset: test_solve_trig matches SymPy", {
  expect_identical(.pyrun(solveset_real(sin(Symbol("x")), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))")
  expect_identical(.pyrun(solveset_real((Integer(-1) + sin(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)")
  expect_identical(.pyrun(solveset_real(cos(Symbol("x")), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers))")
  expect_identical(.pyrun(solveset_real((cos(Symbol("x")) + sin(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi + 3*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + 7*pi/4), Integers))")
  expect_identical(.pyrun(solveset_real(((cos(Symbol("x")))^(Integer(2)) + (sin(Symbol("x")))^(Integer(2))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_complex((Rational(-1, 2) + cos(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/3), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/3), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*sin(Symbol("y"))) + sin((Symbol("a") + Symbol("y")))), Symbol("a"), domain = Reals)), "ConditionSet(a, (-1 <= sin(y)) & (sin(y) <= 1), Union(ImageSet(Lambda(_n, 2*_n*pi - y + asin(sin(y))), Integers), ImageSet(Lambda(_n, 2*_n*pi - y - asin(sin(y)) + pi), Integers)))")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (cos(Symbol("x"))*sin((Integer(2)*Symbol("x")))) + (cos((Integer(2)*Symbol("x")))*sin(Symbol("x")))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*pi/3 + pi/6), Integers)")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (cos((Integer(2)*Symbol("x")))*cos((Integer(4)*Symbol("x"))))), Symbol("x"))), "ImageSet(Lambda(_n, _n*pi), Integers)")
  expect_identical(.pyrun(solveset((Rational(3, 4) + sin((Rational(1, 10)*Symbol("x")))))), "Union(ImageSet(Lambda(_n, 20*_n*pi - 10*asin(3/4) + 20*pi), Integers), ImageSet(Lambda(_n, 20*_n*pi + 10*asin(3/4) + 10*pi), Integers))")
  expect_identical(.pyrun(solveset((cos((Rational(1, 5)*Symbol("x"))) + cos((Rational(1, 15)*Symbol("x")))))), "Union(ImageSet(Lambda(_n, 30*_n*pi + 15*pi/2), Integers), ImageSet(Lambda(_n, 30*_n*pi + 45*pi/2), Integers), ImageSet(Lambda(_n, 30*_n*pi + 75*pi/4), Integers), ImageSet(Lambda(_n, 30*_n*pi + 45*pi/4), Integers), ImageSet(Lambda(_n, 30*_n*pi + 105*pi/4), Integers), ImageSet(Lambda(_n, 30*_n*pi + 15*pi/4), Integers))")
  expect_identical(.pyrun(solveset((Integer(5) + sec((Rational(1, 3)*Symbol("x")*(Integer(2))^(Rational(1, 2))))))), "Union(ImageSet(Lambda(_n, 3*sqrt(2)*(2*_n*pi - asec(-5))/2), Integers), ImageSet(Lambda(_n, 3*sqrt(2)*(2*_n*pi + asec(-5))/2), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*cot((Rational(1, 2)*S_Pi*Symbol("x")))) + tan((S_Pi*Symbol("x")))))), "Union(ImageSet(Lambda(_n, 2*(2*_n*pi + pi/2)/pi), Integers), ImageSet(Lambda(_n, 2*(2*_n*pi - pi/2)/pi), Integers), ImageSet(Lambda(_n, 2*(2*_n*pi - 5*pi/6)/pi), Integers), ImageSet(Lambda(_n, 2*(2*_n*pi + 5*pi/6)/pi), Integers), ImageSet(Lambda(_n, 2*(2*_n*pi - pi/6)/pi), Integers), ImageSet(Lambda(_n, 2*(2*_n*pi + pi/6)/pi), Integers))")
  expect_identical(.pyrun(solveset(cos((Integer(9)*Symbol("x"))))), "Union(ImageSet(Lambda(_n, 2*_n*pi/9 + pi/18), Integers), ImageSet(Lambda(_n, 2*_n*pi/9 + pi/6), Integers))")
  expect_identical(.pyrun(solveset((cot((Integer(12)*Symbol("x"))) + sin((Integer(8)*Symbol("x")))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, _n*pi/2 + pi/8), Integers), ImageSet(Lambda(_n, _n*pi/2 + 3*pi/8), Integers), ImageSet(Lambda(_n, _n*pi/2 + 5*pi/16), Integers), ImageSet(Lambda(_n, _n*pi/2 + 3*pi/16), Integers), ImageSet(Lambda(_n, _n*pi/2 + 7*pi/16), Integers), ImageSet(Lambda(_n, _n*pi/2 + pi/16), Integers))")
  expect_identical(.pyrun(solveset((Rational(-1, 2) + sin((Rational(1, 180)*S_Pi*Symbol("x")))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 180*(2*_n*pi + pi/6)/pi), Integers), ImageSet(Lambda(_n, 180*(2*_n*pi + 5*pi/6)/pi), Integers))")
})

test_that("solveset: test_solve_trig_hyp_by_inversion matches SymPy", {
  expect_identical(.pyrun(solveset_real((Rational(-1, 2) + sin((Integer(3) + (Integer(2)*Symbol("x"))))), Symbol("x"))), "Union(ImageSet(Lambda(_n, _n*pi - 3/2 + 13*pi/12), Integers), ImageSet(Lambda(_n, _n*pi - 3/2 + 17*pi/12), Integers))")
  expect_identical(.pyrun(solveset_complex((Rational(-1, 2) + sin((Integer(3) + (Integer(2)*Symbol("x"))))), Symbol("x"))), "Union(ImageSet(Lambda(_n, _n*pi - 3/2 + 13*pi/12), Integers), ImageSet(Lambda(_n, _n*pi - 3/2 + 17*pi/12), Integers))")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*((Integer(1) + (Rational(-2, 5)*(Integer(5))^(Rational(1, 2)))))^(Rational(1, 2))) + tan(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, _n*pi + pi/10), Integers)")
  expect_identical(.pyrun(solveset_complex(((Integer(-1)*((Integer(1) + (Rational(-2, 5)*(Integer(5))^(Rational(1, 2)))))^(Rational(1, 2))) + tan(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, _n*pi + pi/10), Integers)")
  expect_identical(.pyrun(solveset_real((Integer(-5) + (Integer(3)*cosh((Integer(2)*Symbol("x"))))), Symbol("x"))), "{-acosh(5/3)/2, acosh(5/3)/2}")
  expect_identical(.pyrun(solveset_complex((Integer(-5) + (Integer(3)*cosh((Integer(2)*Symbol("x"))))), Symbol("x"))), "Union(ImageSet(Lambda(_n, _n*I*pi - acosh(5/3)/2), Integers), ImageSet(Lambda(_n, _n*I*pi + acosh(5/3)/2), Integers))")
  expect_identical(.pyrun(solveset_real((Integer(-2) + sinh((Integer(-3) + Symbol("x")))), Symbol("x"))), "{asinh(2) + 3}")
  expect_identical(.pyrun(solveset_complex((Integer(-2) + sinh((Integer(-3) + Symbol("x")))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*I*pi + asinh(2) + 3), Integers), ImageSet(Lambda(_n, 2*_n*I*pi - asinh(2) + 3 + I*pi), Integers))")
  expect_identical(.pyrun(solveset_real(((Rational(-1, 4)*(Integer(2))^(Rational(1, 2))) + (Rational(-1, 4)*(Integer(6))^(Rational(1, 2))) + cos(sinh(Symbol("x")))), Symbol("x"))), "Union(ImageSet(Lambda(_n, asinh(2*_n*pi + pi/12)), Integers), ImageSet(Lambda(_n, asinh(2*_n*pi + 23*pi/12)), Integers))")
  expect_identical(.pyrun(solveset(((Rational(-1, 4)*(Integer(2))^(Rational(1, 2))) + (Rational(-1, 4)*(Integer(6))^(Rational(1, 2))) + cos(sinh(Symbol("x")))), Symbol("x"), Interval(Integer(2), Integer(3), FALSE, FALSE))), "Union(Intersection(ImageSet(Lambda(_n, asinh(2*_n*pi + pi/12)), Integers), Interval(2, 3)), Intersection(ImageSet(Lambda(_n, asinh(2*_n*pi + 23*pi/12)), Integers), Interval(2, 3)))")
  expect_identical(.pyrun(solveset_real((Integer(-2) + cosh((Integer(-1) + (Symbol("x"))^(Integer(2))))), Symbol("x"))), "{-sqrt(1 + acosh(2)), sqrt(1 + acosh(2))}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + sin(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(2) + cos(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(sec(Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(csc(Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(1) + cosh(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(coth(Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(-2) + sech(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real(sech(Symbol("x")), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(1) + tanh(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(tanh(Symbol("x")), Integer(1))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(coth(Symbol("x")), Integer(-1))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(sech(Symbol("x")), Integer(0))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(csch(Symbol("x")), Integer(0))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(-3) + Abs(csch(Symbol("x")))), Symbol("x"))), "{-acsch(3), acsch(3)}")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*exp(Integer(-9))) + tanh((Integer(-1) + (Symbol("x"))^(Integer(2))))), Symbol("x"))), "{-sqrt(atanh(exp(-9)) + 1), sqrt(atanh(exp(-9)) + 1)}")
  expect_identical(.pyrun(solveset_real((Integer(2) + coth(log(Symbol("x")))), Symbol("x"))), "{exp(-acoth(2))}")
  expect_identical(.pyrun(solveset_real((Integer(2) + coth(exp(Symbol("x")))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_complex(((Rational(-1, 2)*sym_I) + sinh(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*I*pi + I*pi/6), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + 5*I*pi/6), Integers))")
  expect_identical(.pyrun(solveset_complex((Rational(3, 4) + sinh((Rational(1, 10)*Symbol("x")))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 20*_n*I*pi - 10*asinh(3/4)), Integers), ImageSet(Lambda(_n, 20*_n*I*pi + 10*asinh(3/4) + 10*I*pi), Integers))")
  expect_identical(.pyrun(solveset_complex((Integer(5) + sech((Rational(1, 3)*Symbol("x")*(Integer(2))^(Rational(1, 2))))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 3*sqrt(2)*(2*_n*I*pi - asech(-5))/2), Integers), ImageSet(Lambda(_n, 3*sqrt(2)*(2*_n*I*pi + asech(-5))/2), Integers))")
  expect_identical(.pyrun(solveset_complex(cosh((Integer(9)*Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*I*pi/9 + I*pi/18), Integers), ImageSet(Lambda(_n, 2*_n*I*pi/9 + I*pi/6), Integers))")
  expect_identical(.pyrun(solveset((Integer(1) + (coth(Symbol("z")))^(Integer(5)) + (Integer(-4)*coth(Symbol("z")))), Symbol("z"), Complexes)), "Union(ImageSet(Lambda(_n, _n*I*pi + acoth(CRootOf(x^5 - 4*x + 1, 0))), Integers), ImageSet(Lambda(_n, _n*I*pi + acoth(CRootOf(x^5 - 4*x + 1, 1))), Integers), ImageSet(Lambda(_n, _n*I*pi + acoth(CRootOf(x^5 - 4*x + 1, 2))), Integers), ImageSet(Lambda(_n, _n*I*pi + acoth(CRootOf(x^5 - 4*x + 1, 3))), Integers), ImageSet(Lambda(_n, _n*I*pi + acoth(CRootOf(x^5 - 4*x + 1, 4))), Integers))")
  expect_identical(.pyrun(solveset((Integer(1) + (coth(Symbol("z")))^(Integer(5)) + (Integer(-4)*coth(Symbol("z")))), Symbol("z"), Reals)), "{acoth(CRootOf(x^5 - 4*x + 1, 0)), acoth(CRootOf(x^5 - 4*x + 1, 2))}")
  expect_identical(.pyrun(solveset(((cos(Symbol("x")))^(Integer(2)) + (Integer(-1)*(Integer(3))^(Rational(1, 2))) + (Integer(2)*cos(Symbol("x"))) + (Rational(-1, 2)*(Integer(3))^(Rational(1, 2))*cos(Symbol("x")))), Symbol("x"), Complexes)), "Union(ImageSet(Lambda(_n, 2*_n*pi - acos(-2)), Integers), ImageSet(Lambda(_n, 2*_n*pi + acos(-2)), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 11*pi/6), Integers))")
  expect_identical(.pyrun(solveset(((cos(Symbol("x")))^(Integer(2)) + (Integer(-1)*(Integer(3))^(Rational(1, 2))) + (Integer(2)*cos(Symbol("x"))) + (Rational(-1, 2)*(Integer(3))^(Rational(1, 2))*cos(Symbol("x")))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 11*pi/6), Integers))")
  expect_identical(.pyrun(solveset((((Integer(1) + (Integer(-1)*sec((Integer(4) + (Symbol("x")*(Integer(3))^(Rational(1, 2))))))))^(Integer(-1))*(Integer(1) + (sec((Integer(4) + (Symbol("x")*(Integer(3))^(Rational(1, 2))))))^(Integer(2)))))), "Union(ImageSet(Lambda(_n, sqrt(3)*(2*_n*pi - 4 - asec(I))/3), Integers), ImageSet(Lambda(_n, sqrt(3)*(2*_n*pi - 4 + asec(I))/3), Integers), ImageSet(Lambda(_n, sqrt(3)*(2*_n*pi - 4 - asec(-I))/3), Integers), ImageSet(Lambda(_n, sqrt(3)*(2*_n*pi - 4 + asec(-I))/3), Integers))")
  expect_identical(.pyrun(solveset((Float("-5.67800000000000", precision = 53) + (tan((Float("3.14000000000000", precision = 53)*Symbol("x"))))^(Rational(3, 2))), Symbol("x"), Interval(Integer(0), Integer(3), FALSE, FALSE))), "{0.403301114561067, 0.403301114561067 + 0.318471337579618*pi, 0.403301114561067 + 0.636942675159236*pi}")
})

test_that("solveset: test_old_trig_issues matches SymPy", {
  expect_identical(.pyrun(solveset(sinh(Symbol("x")), Symbol("x"), Reals)), "{0}")
  expect_identical(.pyrun(solveset(sinh(Symbol("x")), Symbol("x"), Complexes)), "Union(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + I*pi), Integers))")
  expect_identical(.pyrun(solveset(sin((S_Pi*Symbol("x"))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n), Integers), ImageSet(Lambda(_n, (2*_n*pi + pi)/pi), Integers))")
  expect_identical(.pyrun(solveset(sin((S_Pi*Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n), Integers), ImageSet(Lambda(_n, (2*_n*pi + pi)/pi), Integers))")
  expect_identical(.pyrun(solveset((sym_I*cot(((Integer(-8)*sym_E) + (Integer(8)*Symbol("x"))))), Symbol("x"))), "ImageSet(Lambda(_n, _n*pi/8 - 13*pi/16 + E), Integers)")
  expect_identical(.pyrun(solveset((Float("-0.500000000000000", precision = 53) + cos((Integer(2)*Symbol("x")))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, _n*pi - 0.523598775598299), Integers), ImageSet(Lambda(_n, _n*pi + 0.523598775598299), Integers))")
  expect_identical(.pyrun(solveset((Float("-0.500000000000000", precision = 53) + cos((Integer(2)*Symbol("x")))), Symbol("x"), Complexes)), "Union(ImageSet(Lambda(_n, _n*pi - 0.523598775598299), Integers), ImageSet(Lambda(_n, _n*pi + 0.523598775598299), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*(Integer(2))^(Rational(1, 2))) + tan(Symbol("x"))), Symbol("x"), Interval(Integer(0), (Rational(1, 2)*S_Pi), FALSE, FALSE))), "{atan(sqrt(2))}")
  expect_identical(.pyrun(solveset(((Integer(-1)*S_Pi) + tan(Symbol("x"))), Symbol("x"), Interval(Integer(0), (Rational(1, 2)*S_Pi), FALSE, FALSE))), "{atan(pi)}")
  expect_identical(.pyrun(solveset(((Integer(-1)*Symbol("y")) + cos(Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, (-1 <= y) & (y <= 1), Union(ImageSet(Lambda(_n, 2*_n*pi - acos(y)), Integers), ImageSet(Lambda(_n, 2*_n*pi + acos(y)), Integers)))")
  expect_identical(.pyrun(solveset(sin(log(Symbol("x"))), Symbol("x"), Interval(Integer(0), Integer(1), TRUE, FALSE))), "Union(Intersection(ImageSet(Lambda(_n, exp(2*_n*pi)), Integers), Interval.Lopen(0, 1)), Intersection(ImageSet(Lambda(_n, exp(2*_n*pi + pi)), Integers), Interval.Lopen(0, 1)))")
  expect_identical(.pyrun(solveset(((Integer(-1)*sin(Integer(1))) + sin(Symbol("x"))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi + 1), Integers), ImageSet(Lambda(_n, 2*_n*pi - 1 + pi), Integers))")
  expect_identical(.pyrun(solveset(((Rational(-1, 3)*(Integer(5))^(Rational(1, 2))) + sin(Symbol("x"))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi + asin(sqrt(5)/3)), Integers), ImageSet(Lambda(_n, 2*_n*pi - asin(sqrt(5)/3) + pi), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*cosh(Integer(2))) + sinh(Symbol("x"))), Symbol("x"), Reals)), "{asinh(cosh(2))}")
  expect_identical(.pyrun(solveset(Eq(tan(Symbol("x")), Symbol("y")), Symbol("x"), domain = Reals)), "ConditionSet(x, (-oo < y) & (y < oo), ImageSet(Lambda(_n, _n*pi + atan(y)), Integers))")
  expect_identical(.pyrun(solveset(Eq(tan(Symbol("x")), Symbol("r", real = TRUE)), Symbol("x"), domain = Reals)), "ImageSet(Lambda(_n, _n*pi + atan(r)), Integers)")
})

test_that("solveset: test_solve_hyperbolic matches SymPy", {
  expect_identical(.pyrun(solveset((cosh(Symbol("x")) + sinh(Symbol("x"))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset((cos(Symbol("x")) + sinh(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(cos(x) + sinh(x), 0), Complexes)")
  expect_identical(.pyrun(solveset_real((sech(Symbol("x")) + sinh(Symbol("x"))), Symbol("x"))), "{log(sqrt(-2 + sqrt(5)))}")
  expect_identical(.pyrun(solveset_real((Integer(-5) + (Integer(2)*sinh(Symbol("x"))) + cosh((Integer(2)*Symbol("x")))), Symbol("x"))), "{log(-2 + sqrt(5)), log(1 + sqrt(2))}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + ((cosh(Symbol("x")))^(Integer(-1))*(coth(Symbol("x")) + sinh((Integer(2)*Symbol("x")))))), Symbol("x"))), "{log(1/2 + sqrt(5)/2), log(1 + sqrt(2))}")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (cosh(Symbol("x"))*sinh(Symbol("x")))), Symbol("x"))), "{log(4 + sqrt(17))/2}")
  expect_identical(.pyrun(solveset_real((Integer(-1) + sinh(Symbol("x")) + tanh(Symbol("x"))), Symbol("x"))), "{log(sqrt(2)/2 + sqrt(-1/2 + sqrt(2)))}")
  expect_identical(.pyrun(solveset_complex((sech(Symbol("x")) + sinh(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*I*pi + log(sqrt(-2 + sqrt(5)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi/2) + log(sqrt(2 + sqrt(5)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(sqrt(-2 + sqrt(5)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi - pi/2) + log(sqrt(2 + sqrt(5)))), Integers))")
  expect_identical(.pyrun(solveset((cosh((Rational(1, 5)*Symbol("x"))) + cosh((Rational(1, 15)*Symbol("x")))))), "Union(ImageSet(Lambda(_n, 15*I*(2*_n*pi + pi/2)), Integers), ImageSet(Lambda(_n, 15*I*(2*_n*pi - pi/2)), Integers), ImageSet(Lambda(_n, 15*I*(2*_n*pi - 3*pi/4)), Integers), ImageSet(Lambda(_n, 15*I*(2*_n*pi + 3*pi/4)), Integers), ImageSet(Lambda(_n, 15*I*(2*_n*pi - pi/4)), Integers), ImageSet(Lambda(_n, 15*I*(2*_n*pi + pi/4)), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*coth((Rational(1, 2)*S_Pi*Symbol("x")))) + tanh((S_Pi*Symbol("x")))))), "Union(ImageSet(Lambda(_n, 2*I*(2*_n*pi + pi/2)/pi), Integers), ImageSet(Lambda(_n, 2*I*(2*_n*pi - pi/2)/pi), Integers))")
  expect_identical(.pyrun(solveset(((Integer(-1)*cosh((Integer(5)*Symbol("x")))) + cosh(Symbol("x")) + cosh((Integer(3)*Symbol("x")))), Symbol("x"), Reals)), "ConditionSet(x, Eq(cosh(x) + cosh(3*x) - cosh(5*x), 0), Reals)")
  expect_identical(.pyrun(solveset((coth((Integer(12)*Symbol("x"))) + sinh((Integer(8)*Symbol("x")))))), "ConditionSet(x, Eq(sinh(8*x) + coth(12*x), 0), Complexes)")
})

test_that("solveset: test_solve_trig_hyp_symbolic matches SymPy", {
  expect_identical(.pyrun(solveset(sin((Symbol("a")*Symbol("x"))), Symbol("x"))), "ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi + pi)/a), Integers), ImageSet(Lambda(_n, 2*_n*pi/a), Integers)))")
  expect_identical(.pyrun(solveset(cosh((Symbol("x")*(Symbol("a"))^(Integer(-1)))), Symbol("x"))), "ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, a*(2*_n*I*pi + I*pi/2)), Integers), ImageSet(Lambda(_n, a*(2*_n*I*pi + 3*I*pi/2)), Integers)))")
  expect_identical(.pyrun(solveset((cos((Rational(4, 3)*Symbol("x")*(Integer(3))^(Rational(1, 2))*(S_Pi)^(Integer(-1))*(Symbol("a"))^(Integer(2))*(Symbol("b"))^(Integer(-1)))) + sin((Rational(2, 3)*Symbol("x")*(Integer(3))^(Rational(1, 2))*(S_Pi)^(Integer(-1))*(Symbol("a"))^(Integer(2))*(Symbol("b"))^(Integer(-1))))), Symbol("x"))), "ConditionSet(x, Ne(b, 0) & Ne(a^2, 0), Union(ImageSet(Lambda(_n, sqrt(3)*pi*b*(2*_n*pi + pi/2)/(2*a^2)), Integers), ImageSet(Lambda(_n, sqrt(3)*pi*b*(2*_n*pi - 5*pi/6)/(2*a^2)), Integers), ImageSet(Lambda(_n, sqrt(3)*pi*b*(2*_n*pi - pi/6)/(2*a^2)), Integers)))")
  expect_identical(.pyrun(solveset((Integer(-3) + cosh((Symbol("x")*(Integer(1) + (Symbol("a"))^(Integer(2)))))), Symbol("x"))), "ConditionSet(x, Ne(a^2 + 1, 0), Union(ImageSet(Lambda(_n, (2*_n*I*pi - acosh(3))/(a^2 + 1)), Integers), ImageSet(Lambda(_n, (2*_n*I*pi + acosh(3))/(a^2 + 1)), Integers)))")
  expect_identical(.pyrun(solveset((Integer(-2) + cosh((Symbol("x")*(Integer(1) + (Symbol("ar", real = TRUE))^(Integer(2)))))), Symbol("x"), Reals)), "{-acosh(2)/(ar^2 + 1), acosh(2)/(ar^2 + 1)}")
  expect_identical(.pyrun(solveset(((Integer(-1)*cot((Symbol("x")*(Integer(3) + (Integer(3)*sym_I))))) + cot((Symbol("x")*(Integer(1) + sym_I)))), Symbol("x"))), "Union(ImageSet(Lambda(_n, (1 - I)*(2*_n*pi + pi/2)/2), Integers), ImageSet(Lambda(_n, (1 - I)*(2*_n*pi - pi/2)/2), Integers))")
})

test_that("solveset: test_issue_9616 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-1) + sinh(Symbol("x")) + tanh(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*I*pi + log(sqrt(2)/2 + sqrt(-1/2 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi - atan(sqrt(2)*sqrt(1/2 + sqrt(2))) + pi) + log(sqrt(1 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(-sqrt(2)/2 + sqrt(-1/2 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi - pi + atan(sqrt(2)*sqrt(1/2 + sqrt(2)))) + log(sqrt(1 + sqrt(2)))), Integers))")
})

test_that("solveset: test_solve_invalid_sol matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Integer(-1))*sin(Symbol("x"))), Symbol("x"))), "Union(Complement(ImageSet(Lambda(_n, 2*_n*pi), Integers), {0}), Complement(ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), {0}))")
  expect_identical(.pyrun(solveset_complex(((Symbol("x"))^(Integer(-1))*(Integer(-1) + exp(Symbol("x")))), Symbol("x"))), "Complement(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), {0})")
})

test_that("solveset: test_solve_lambert matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Symbol("x")*exp(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(x*exp(x) - 1, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Symbol("x") + exp(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(x + exp(x), 0), Reals)")
  expect_identical(.pyrun(solveset_real((Symbol("x") + (Integer(2))^(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(2^x + x, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(5) + (Integer(2))^((Integer(3) + (Integer(-5)*Symbol("x")))) + (Integer(3)*Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(2^(3 - 5*x) + 3*x + 5, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(-6)*(Integer(7))^((Integer(9) + (Integer(3)*Symbol("x"))))) + (Integer(2)*((Integer(4) + (Integer(3)*Symbol("x"))))^(Integer(5)))), Symbol("x"))), "ConditionSet(x, Eq(-6*7^(3*x + 9) + 2*(3*x + 4)^5, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(2048) + (Integer(-242121642)*(Integer(7))^((Integer(3)*Symbol("x")))) + (Integer(486)*(Symbol("x"))^(Integer(5))) + (Integer(3240)*(Symbol("x"))^(Integer(4))) + (Integer(7680)*Symbol("x")) + (Integer(8640)*(Symbol("x"))^(Integer(3))) + (Integer(11520)*(Symbol("x"))^(Integer(2)))), Symbol("x"))), "ConditionSet(x, Eq(-242121642*7^(3*x) + 486*x^5 + 3240*x^4 + 8640*x^3 + 11520*x^2 + 7680*x + 2048, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (Integer(3)*exp((Integer(2) + (Integer(-7)*Symbol("x"))))) + (Integer(5)*Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(5*x + 3*exp(2 - 7*x) - 1, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(5) + (Integer(2)*Symbol("x")) + log((Integer(-2) + (Integer(3)*Symbol("x"))))), Symbol("x"))), "ConditionSet(x, Eq(2*x + log(3*x - 2) + 5, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(3)*Symbol("x")) + log((Integer(4)*Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(3*x + log(4*x), 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Symbol("x"))^(Symbol("x"))))), "ERR TypeError: solveset_real() missing 1 required positional argument: 'symbol'")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Symbol("a")*Symbol("x")) + (Integer(2)*Symbol("x")*log(Symbol("x")))), Symbol("x"))), "Union({0}, Intersection({exp(a/2)}, Reals))")
  expect_identical(.pyrun(solveset_real(((Symbol("a", real = TRUE)*(Symbol("x"))^(Integer(-1))) + exp((Rational(1, 2)*Symbol("x")))), Symbol("x"))), "Complement(ConditionSet(x, Eq(a + x*exp(x/2), 0), Reals), {0})")
  expect_identical(.pyrun(solveset_real(((Rational(1, 2)*exp((Rational(1, 2)*Symbol("x")))) + (Integer(-1)*Symbol("a", real = TRUE)*(Symbol("x"))^(Integer(-2)))), Symbol("x"))), "Complement(ConditionSet(x, Eq(-2*a + x^2*exp(x/2), 0), Reals), {0})")
  expect_identical(.pyrun(solveset_real((Integer(-1) + (tanh((Integer(-3) + Symbol("x")))*tanh((Integer(3) + Symbol("x"))))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset_real((Integer(1) + (((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2)) + (Integer(-6)*Symbol("x")) + (Integer(-2)*log(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(-6*x + (3*x + log(x))^2 - 2*log(x) + 1, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(3) + ((Integer(-1) + (((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2))))^(Integer(2)) + (Integer(-2)*(((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2)))), Symbol("x"))), "ConditionSet(x, Eq(-2*(3*x + log(x))^2 + ((3*x + log(x))^2 - 1)^2 + 3, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2)) + (Integer(-6)*Symbol("x")) + (Integer(-2)*log(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(-6*x + (3*x + log(x))^2 - 2*log(x) - 2, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(1) + (Integer(3)*Symbol("x")) + (Symbol("x")*log(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(x*log(x) + 3*x + 1, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-3) + (Symbol("x")*exp(Symbol("x"))*exp((Symbol("x")*exp(Symbol("x")))))), Symbol("x"))), "ConditionSet(x, Eq(x*exp(x)*exp(x*exp(x)) - 3, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Symbol("a", real = TRUE))^((Integer(5) + (Integer(3)*Symbol("x")))) + (Integer(3)*log((Symbol("a", real = TRUE))^((Integer(5) + (Integer(3)*Symbol("x"))))))), Symbol("x"))), "ConditionSet(x, Eq(a^(3*x + 5) + 3*log(a^(3*x + 5)), 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Symbol("p", positive = TRUE))^((Integer(5) + (Integer(3)*Symbol("x")))) + (Integer(3)*log((Symbol("p", positive = TRUE))^((Integer(5) + (Integer(3)*Symbol("x"))))))), Symbol("x"))), "ConditionSet(x, Eq(p^(3*x + 5) + 3*log(p^(3*x + 5)), 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Symbol("a", real = TRUE))^((Integer(5) + (Integer(3)*Symbol("x")))) + (Integer(3)*log((Symbol("a", real = TRUE))^((Integer(5) + (Integer(3)*Symbol("x")))))) + (Symbol("b")*log((Symbol("a", real = TRUE))^((Integer(5) + (Integer(3)*Symbol("x"))))))), Symbol("x"))), "ConditionSet(x, Eq(a^(3*x + 5) + b*log(a^(3*x + 5)) + 3*log(a^(3*x + 5)), 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Rational(1, 4)*exp((Rational(1, 2)*Symbol("x")))) + (Integer(2)*Symbol("a", real = TRUE)*(Symbol("x"))^(Integer(-3)))), Symbol("x"))), "Complement(ConditionSet(x, Eq(8*a + x^3*exp(x/2), 0), Reals), {0})")
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Integer(3)) + (Integer(-1)*(Integer(3))^(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(-3^x + x^3, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(3))^(cos(Symbol("x"))) + (Integer(-1)*(cos(Symbol("x")))^(Integer(3)))))), "ERR TypeError: solveset_real() missing 1 required positional argument: 'symbol'")
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Integer(2))^(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(-2^x + x^2, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(2))^(Symbol("x")) + (Integer(-1)*(Symbol("x"))^(Integer(2)))), Symbol("x"))), "ConditionSet(x, Eq(2^x - x^2, 0), Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(3)*log(Symbol("x"))) + (Integer(-1)*Symbol("x")*log(Integer(3)))))), "ERR TypeError: solveset_real() missing 1 required positional argument: 'symbol'")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*Symbol("y")) + LambertW((Integer(2)*Symbol("x")))))), "ERR TypeError: solveset_real() missing 1 required positional argument: 'symbol'")
})

test_that("solveset: test_other_lambert matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Rational(6, 5)) + (Integer(-1)*(Rational(6, 5))^(Symbol("x")))), Symbol("x"))), "ConditionSet(x, Eq(5^(5*x)*x^6 - 6^(5*x), 0), Reals)")
})

test_that("solveset: test_solveset matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("x") + Symbol("y")))), "ERR ValueError: \\nThe independent variable must be specified for a multivariate\\nequation.")
  expect_identical(.pyrun(solveset(Symbol("x"), Integer(1))), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(1) + Symbol("y") + (Function("f")(Integer(1)))^(Integer(2))), Function("f")(Integer(1)))), "{-sqrt(-y - 1), sqrt(-y - 1)}")
  expect_identical(.pyrun(solveset((Integer(-1) + (Function("f")(Integer(1)))^(Integer(2))), Function("f")(Integer(1)), Reals)), "{-1, 1}")
  expect_identical(.pyrun(solveset((Integer(1) + (Function("f")(Integer(1)))^(Integer(2))), Function("f")(Integer(1)))), "{-I, I}")
  expect_identical(.pyrun(solveset((Integer(-1) + Symbol("x")), Integer(1))), "{x}")
  expect_identical(.pyrun(solveset(((Integer(-1)*cos(Symbol("x"))) + sin(Symbol("x"))), sin(Symbol("x")))), "{cos(x)}")
  expect_identical(.pyrun(solveset(Integer(0), domain = Reals)), "Reals")
  expect_identical(.pyrun(solveset(Integer(1))), "EmptySet")
  expect_identical(.pyrun(solveset(TRUE, domain = Reals)), "Reals")
  expect_identical(.pyrun(solveset(FALSE, domain = Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(-1) + exp(Symbol("x"))), domain = Reals)), "{0}")
  expect_identical(.pyrun(solveset((Integer(-1) + exp(Symbol("x"))), Symbol("x"), Reals)), "{0}")
  expect_identical(.pyrun(solveset(Eq(exp(Symbol("x")), Integer(1)), Symbol("x"), Reals)), "{0}")
  expect_identical(.pyrun(solveset((Integer(-1) + exp(Symbol("x"))), exp(Symbol("x")), Reals)), "{1}")
  expect_identical(.pyrun(solveset(Ge((Integer(-1) + Symbol("x")), Integer(0)), Symbol("x"), Reals)), "Interval(1, oo)")
  expect_identical(.pyrun(solveset(Ge((Integer(-1) + exp(Symbol("x"))), Integer(0)), Symbol("x"), Reals)), "Interval(0, oo)")
  expect_identical(.pyrun(solveset((Integer(-1) + exp(Symbol("x"))), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*I*pi), Integers)")
  expect_identical(.pyrun(solveset(Eq(exp(Symbol("x")), Integer(1)), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*I*pi), Integers)")
  expect_identical(.pyrun(solveset((Integer(1) + (Symbol("x"))^(Integer(2)) + Function("f")(Integer(0))), Symbol("x"))), "{-sqrt(-f(0) - 1), sqrt(-f(0) - 1)}")
  expect_identical(.pyrun(solveset(Gt(atan(log(Symbol("x"))), Integer(0)), Symbol("x"), domain = Interval(Integer(0), oo, TRUE, TRUE))), "Interval.open(1, oo)")
})

test_that("solveset: test_multi_exp matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-5) + exp(exp(Symbol("x")))), Symbol("x"))), "ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(2*_n*I*pi + log(5))) + log(Abs(2*_n*I*pi + log(5)))), ProductSet(Integers, Integers))")
  expect_identical(.pyrun(solveset((Symbol("c") + (Symbol("d")*exp(exp((Symbol("b") + (Symbol("a")*Symbol("x"))))))), Symbol("x"))), "ImageSet(Lambda(x, (-b + x)/a), ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d))))), ProductSet(Integers, Integers)))")
  expect_identical(.pyrun(solveset((Symbol("c") + (Symbol("d")*exp(exp(exp((Symbol("b") + (Symbol("a")*Symbol("x")))))))), Symbol("x"))), "ImageSet(Lambda(x, (-b + x)/a), ImageSet(Lambda(((_k2, _k1, _n),), I*(2*_k2*pi + arg(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))))) + log(Abs(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d))))))), ProductSet(Integers, Integers, Integers)))")
  expect_identical(.pyrun(solveset((Symbol("c") + (Symbol("d")*exp(exp(exp(exp((Symbol("b") + (Symbol("a")*Symbol("x"))))))))), Symbol("x"))), "ImageSet(Lambda(x, (-b + x)/a), ImageSet(Lambda(((_k3, _k2, _k1, _n),), I*(2*_k3*pi + arg(I*(2*_k2*pi + arg(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))))) + log(Abs(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))))))) + log(Abs(I*(2*_k2*pi + arg(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))))) + log(Abs(I*(2*_k1*pi + arg(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d)))) + log(Abs(I*(2*_n*pi + arg(-c/d)) + log(Abs(c/d))))))))), ProductSet(Integers, Integers, Integers, Integers)))")
})

test_that("solveset: test__solveset_multi matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list(Symbol("x")), list(Reals))), "{(-1,), (1,)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Integer(1) + Symbol("x"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "{(-1, 1)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Integer(1) + Symbol("x"))), list(Symbol("y"), Symbol("x")), list(Reals, Reals))), "{(1, -1)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Integer(-1) + Symbol("x") + (Integer(-1)*Symbol("y")))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "{(1/2, -1/2)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + Symbol("x")), (Integer(-2) + Symbol("y"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "{(1, 2)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "Union(ImageSet(Lambda(((x,),), (x, -x)), ProductSet(Reals)), ImageSet(Lambda(((y,),), (-y, y)), ProductSet(Reals)))")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Integer(1) + Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "EmptySet")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Symbol("x") + (Integer(-1)*Symbol("y"))), (Integer(-1) + Symbol("x"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "EmptySet")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Symbol("x") + Symbol("y")), (Symbol("x") + (Integer(-1)*Symbol("y"))), (Integer(-1) + Symbol("x"))), list(Symbol("y"), Symbol("x")), list(Reals, Reals))), "EmptySet")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + Symbol("x") + Symbol("y") + Symbol("z")), (Integer(-2) + Symbol("x") + Symbol("y") + (Integer(-1)*Symbol("z"))), (Integer(-3) + Symbol("x") + (Integer(-1)*Symbol("y")) + (Integer(-1)*Symbol("z")))), list(Symbol("x"), Symbol("y"), Symbol("z")), list(Reals, Reals, Reals))), "{(2, -1/2, -1/2)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-2) + (Symbol("x"))^(Integer(2)) + (Symbol("y"))^(Integer(2))), (Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "{(-1, 1), (1, -1)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + (Symbol("x"))^(Integer(2))), Symbol("y")), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "{(-1, 0), (1, 0)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list(((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2))))), list(Symbol("x"), Symbol("y")), list(Reals, Reals))), "Union(ImageSet(Lambda(((x,),), (x, -Abs(x))), ProductSet(Reals)), ImageSet(Lambda(((x,),), (x, Abs(x))), ProductSet(Reals)), ImageSet(Lambda(((y,),), (-Abs(y), y)), ProductSet(Reals)), ImageSet(Lambda(((y,),), (Abs(y), y)), ProductSet(Reals)))")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + (Symbol("r")*cos(Symbol("theta")))), (Symbol("r")*sin(Symbol("theta")))), list(Symbol("theta"), Symbol("r")), list(Interval(Integer(0), S_Pi, FALSE, FALSE), Interval(Integer(-1), Integer(1), FALSE, FALSE)))), "{(0, 1), (pi, -1)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list((Integer(-1) + (Symbol("r")*cos(Symbol("theta")))), (Symbol("r")*sin(Symbol("theta")))), list(Symbol("r"), Symbol("theta")), list(Interval(Integer(0), Integer(1), FALSE, FALSE), Interval(Integer(0), S_Pi, FALSE, FALSE)))), "{(1, 0)}")
  expect_identical(.pyrun(symcompR:::.sset_solveset_multi(list(((Integer(-1)*Symbol("r")) + (Symbol("r")*cos(Symbol("theta")))), (Symbol("r")*sin(Symbol("theta")))), list(Symbol("r"), Symbol("theta")), list(Interval(Integer(0), Integer(1), FALSE, FALSE), Interval(Integer(0), S_Pi, FALSE, FALSE)))), "Union(ImageSet(Lambda(((r,),), (r, 0)), ImageSet(Lambda(r, (r,)), Interval(0, 1))), ImageSet(Lambda(((theta,),), (0, theta)), ImageSet(Lambda(theta, (theta,)), Interval(0, pi))))")
})

test_that("solveset: test_conditionset matches SymPy", {
  expect_identical(.pyrun(solveset(Eq(((cos(Symbol("x")))^(Integer(2)) + (sin(Symbol("x")))^(Integer(2))), Integer(1)), Symbol("x"), domain = Reals)), "Reals")
  expect_identical(.pyrun(solveset(Eq(((Symbol("x"))^(Integer(2)) + (Symbol("x")*sin(Symbol("x")))), Integer(1)), Symbol("x"), domain = Reals)), "ConditionSet(x, Eq(x^2 + x*sin(x) - 1, 0), Reals)")
  expect_identical(.pyrun(solveset(Eq((Rational(-1, 2)*sym_I*((Integer(-1)*exp((Integer(-1)*sym_I*Symbol("x")))) + exp((sym_I*Symbol("x"))))), Integer(1)), Symbol("x"))), "ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)")
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + sin(Symbol("x"))), Integer(1)), Symbol("x"), domain = Reals)), "ConditionSet(x, x + sin(x) > 1, Reals)")
  expect_identical(.pyrun(solveset(Eq(sin(Abs(Symbol("x"))), Symbol("x")), Symbol("x"), domain = Reals)), "ConditionSet(x, Eq(-x + sin(Abs(x)), 0), Reals)")
  expect_identical(.pyrun(solveset(((Symbol("y"))^(Symbol("x")) + (Integer(-1)*Symbol("z"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(y^x - z, 0), Reals)")
})

test_that("solveset: test_conditionset_equality matches SymPy", {
  expect_identical(.pyrun(solveset(Eq(tan(Symbol("x")), Symbol("y")), Symbol("x"))), "ImageSet(Lambda(_n, _n*pi + atan(y)), Integers)")
})

test_that("solveset: test_solveset_domain matches SymPy", {
  expect_identical(.pyrun(solveset(Eq((Symbol("x") + Symbol("y")), (sym_I + (Integer(2)*Symbol("y")))), Symbol("x"), Reals)), "Intersection({y + I}, Reals)")
  expect_identical(.pyrun(solveset((Integer(-6) + (Symbol("x"))^(Integer(2)) + (Integer(-1)*Symbol("x"))), Symbol("x"), Interval(Integer(0), oo, FALSE, TRUE))), "{3}")
  expect_identical(.pyrun(solveset((Integer(-1) + (Symbol("x"))^(Integer(2))), Symbol("x"), Interval(Integer(0), oo, FALSE, TRUE))), "{1}")
  expect_identical(.pyrun(solveset((Integer(-16) + (Symbol("x"))^(Integer(4))), Symbol("x"), Interval(Integer(0), Integer(10), FALSE, FALSE))), "{2}")
})

test_that("solveset: test_improve_coverage matches SymPy", {
  expect_identical(.pyrun(solveset((exp(Symbol("x")) + sin(Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(exp(x) + sin(x), 0), Reals)")
})

test_that("solveset: test_issue_9522 matches SymPy", {
  expect_identical(.pyrun(solveset(Eq((Symbol("x") + ((Integer(-4) + (Symbol("x"))^(Integer(2))))^(Integer(-1))), (Integer(2) + ((Integer(-4) + (Symbol("x"))^(Integer(2))))^(Integer(-1)))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Eq((Symbol("x") + (Symbol("x"))^(Integer(-1))), (Symbol("x"))^(Integer(-1))), Symbol("x"), Reals)), "EmptySet")
})

test_that("solveset: test_solvify matches SymPy", {
  expect_identical(.pyrun(solvify((Integer(10) + (Symbol("x"))^(Integer(2))), Symbol("x"), Reals)), "[]")
  expect_identical(.pyrun(solvify((Integer(1) + (Symbol("x"))^(Integer(3))), Symbol("x"), Complexes)), "[-1, 1/2 - sqrt(3)*I/2, 1/2 + sqrt(3)*I/2]")
  expect_identical(.pyrun(solvify(log(Symbol("x")), Symbol("x"), Reals)), "[1]")
  expect_identical(.pyrun(solvify(cos(Symbol("x")), Symbol("x"), Reals)), "[pi/2, 3*pi/2]")
  expect_identical(.pyrun(solvify((Integer(1) + sin(Symbol("x"))), Symbol("x"), Reals)), "[3*pi/2]")
  expect_identical(.pyrun(solvify(sin(exp(Symbol("x"))), Symbol("x"), Complexes)), "ERR NotImplementedError: solveset is unable to solve this equation.")
})

test_that("solveset: test_solvify_piecewise matches SymPy", {
  expect_identical(.pyrun(solvify(Piecewise(list(Integer(0), Lt(Symbol("x"), Integer(-1))), list((Symbol("x"))^(Integer(2)), Le(Symbol("x"), Integer(1))), list(log(Symbol("x")), S_true)), Symbol("x"), Reals)), "[0]")
  expect_identical(.pyrun(solvify(Piecewise(list(Integer(0), Lt(Symbol("x"), Integer(-10))), list((Integer(-6) + (Symbol("x"))^(Integer(2)) + (Integer(5)*Symbol("x"))), Ge(Symbol("x"), Integer(-9)))), Symbol("x"), Reals)), "[-6, 1]")
  expect_identical(.pyrun(solvify(Piecewise(list(Integer(0), Eq(Symbol("x"), Integer(0))), list(((Symbol("x"))^(Integer(2))*(Abs(Symbol("x")))^(Integer(-1))), S_true)), Symbol("x"), Reals)), "[0]")
  expect_identical(.pyrun(solvify(Piecewise(list(Integer(0), Eq(Symbol("x"), S_Pi)), list(((sin(Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*S_Pi))), S_true)), Symbol("x"), Reals)), "[pi]")
})

test_that("solveset: test_abs_invert_solvify matches SymPy", {
  expect_identical(.pyrun(solvify(sin(Symbol("x", positive = TRUE)), Symbol("x", positive = TRUE), Reals)), "[0, pi]")
  expect_identical(.pyrun(solvify(sin(Abs(Symbol("x"))), Symbol("x"), Reals)), "None")
})

test_that("solveset: test_linear_eq_to_matrix matches SymPy", {
  expect_identical(.pyrun(linear_eq_to_matrix(Integer(0), Symbol("x"))), "(Matrix([[0]]), Matrix([[0]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(Integer(1), Symbol("x"))), "(Matrix([[0]]), Matrix([[-1]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Integer(-3) + Symbol("y") + (Integer(-2)*Symbol("z")) + (Integer(2)*Symbol("x"))), (Symbol("x") + (Integer(-1)*Symbol("y")) + (Integer(-1)*Symbol("z"))), (Integer(-12) + Symbol("x") + Symbol("y") + (Integer(3)*Symbol("z")))), Symbol("x"), Symbol("y"), Symbol("z"))), "(Matrix([\\n[2,  1, -2],\\n[1, -1, -1],\\n[1,  1,  3]]), Matrix([\\n[ 3],\\n[ 0],\\n[12]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list(Eq(((Integer(-1)*Symbol("z")) + (Integer(2)*Symbol("y")) + (Integer(3)*Symbol("x"))), Integer(1)), Eq(((Integer(-2)*Symbol("y")) + (Integer(2)*Symbol("x")) + (Integer(4)*Symbol("z"))), Integer(-2)), (Symbol("y") + (Integer(-2)*Symbol("x")) + (Integer(-2)*Symbol("z")))), Symbol("x"), Symbol("y"), Symbol("z"))), "(Matrix([\\n[ 3,  2, -1],\\n[ 2, -2,  4],\\n[-2,  1, -2]]), Matrix([\\n[ 1],\\n[-2],\\n[ 0]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list(((Integer(-1)*Symbol("d")) + (Symbol("b")*Symbol("y")) + (Symbol("c")*Symbol("z")) + (Symbol("a")*Symbol("b")*Symbol("x"))), ((Integer(-1)*Symbol("h")) + (Symbol("d")*Symbol("x")) + (Symbol("e")*Symbol("x")) + (Symbol("f")*Symbol("y")) + (Symbol("g")*Symbol("z"))), ((Integer(-1)*Symbol("l")) + (Symbol("i")*Symbol("x")) + (Symbol("j")*Symbol("y")) + (Symbol("k")*Symbol("z")))), Symbol("x"), Symbol("y"), Symbol("z"))), "(Matrix([\\n[  a*b, b, c],\\n[d + e, f, g],\\n[    i, j, k]]), Matrix([\\n[d],\\n[h],\\n[l]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list(((Integer(-1)*Symbol("d")) + (Symbol("b")*Symbol("y")) + (Symbol("c")*Symbol("z")) + (Symbol("a")*Symbol("b")*Symbol("x"))), ((Integer(-1)*Symbol("h")) + (Symbol("d")*Symbol("x")) + (Symbol("e")*Symbol("x")) + (Symbol("f")*Symbol("y")) + (Symbol("g")*Symbol("z"))), ((Integer(-1)*Symbol("l")) + (Symbol("i")*Symbol("x")) + (Symbol("j")*Symbol("y")) + (Symbol("k")*Symbol("z")))))), "ERR ValueError: \\nSymbols must be given, for which coefficients are to be found.")
  expect_identical(.pyrun(linear_eq_to_matrix(list(((Integer(-1)*Symbol("d")) + (Symbol("b")*Symbol("y")) + (Symbol("c")*Symbol("z")) + (Symbol("a")*Symbol("b")*Symbol("x"))), ((Integer(-1)*Symbol("h")) + (Symbol("d")*Symbol("x")) + (Symbol("e")*Symbol("x")) + (Symbol("f")*Symbol("y")) + (Symbol("g")*Symbol("z"))), ((Integer(-1)*Symbol("l")) + (Symbol("i")*Symbol("x")) + (Symbol("j")*Symbol("y")) + (Symbol("k")*Symbol("z")))), list(Symbol("x"), Symbol("x"), Symbol("y")))), "ERR ValueError: Symbols must be unique")
  expect_identical(.pyrun(linear_eq_to_matrix(Eq((Symbol("x") + (Symbol("x"))^(Integer(-1))), (Symbol("x"))^(Integer(-1))), list(Symbol("x")))), "ERR NonlinearError: nonlinear term: 1/x")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("x"))^(Integer(2))), list(Symbol("x")))), "ERR NonlinearError: nonlinear term: x^2")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("x")*Symbol("y"))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: \\nnonlinear cross-term: x*y")
  expect_identical(.pyrun(linear_eq_to_matrix(S_true, Symbol("x"))), "ERR ValueError: \\nEquation(s) must be given as a sequence, Expr, Eq or Matrix.")
  expect_identical(.pyrun(linear_eq_to_matrix(S_false, Symbol("x"))), "ERR ValueError: \\nEquation(s) must be given as a sequence, Expr, Eq or Matrix.")
  expect_identical(.pyrun(linear_eq_to_matrix(Eq(Integer(1), Integer(2)), Symbol("x"))), "ERR ValueError: \\nEquation(s) must be given as a sequence, Expr, Eq or Matrix.")
  expect_identical(.pyrun(linear_eq_to_matrix(list(Symbol("x")), list((Symbol("x"))^(Integer(-1))))), "(Matrix([[0]]), Matrix([[-x]]))")
  expect_identical(.pyrun(linear_eq_to_matrix((Symbol("x") + (Symbol("y")*(Integer(3) + (Symbol("z")*(Integer(2) + (Integer(3)*Symbol("x"))))))), Symbol("x"))), "(Matrix([[3*y*z + 1]]), Matrix([[-y*(2*z + 3)]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(Matrix(list(list((Integer(-7) + (Symbol("a")*Symbol("x")) + (Symbol("b")*Symbol("y")))), list(((Integer(-1)*Symbol("c")) + (Integer(5)*Symbol("x")) + (Integer(6)*Symbol("y")))))), Symbol("x"), Symbol("y"))), "(Matrix([\\n[a, b],\\n[5, 6]]), Matrix([\\n[7],\\n[c]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(Eq((Integer(2) + Symbol("x")), Integer(1)), Symbol("x"))), "(Matrix([[1]]), Matrix([[-1]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list(), list(Symbol("x"), Symbol("y")))), "(Matrix(0, 2, []), Matrix(0, 1, []))")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "(Matrix([[1, 1]]), Matrix([[0]]))")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "(Matrix([[1, 1]]), Matrix([[0]]))")
})

test_that("solveset: test_issue_16577 matches SymPy", {
  expect_identical(.pyrun(linear_eq_to_matrix(Eq(((Integer(4)*Symbol("y")) + (Symbol("a")*((Integer(2)*Symbol("x")) + (Integer(3)*Symbol("y"))))), Integer(5)), Symbol("x"), Symbol("y"))), "(Matrix([[2*a, 3*a + 4]]), Matrix([[5]]))")
})

test_that("solveset: test_issue_10085 matches SymPy", {
  expect_identical(.pyrun(invert_real(exp(Symbol("x")), Integer(0), Symbol("x"))), "(x, EmptySet)")
})

test_that("solveset: test_linsolve matches SymPy", {
  expect_identical(.pyrun(linsolve(list((Integer(-7) + Symbol("x1") + Symbol("x3") + Symbol("x4") + (Integer(2)*Symbol("x2"))), (Integer(-12) + Symbol("x1") + (Integer(-1)*Symbol("x4")) + (Integer(2)*Symbol("x2")) + (Integer(2)*Symbol("x3"))), (Integer(-4) + (Integer(2)*Symbol("x1")) + (Integer(4)*Symbol("x2")) + (Integer(6)*Symbol("x4")))), list(Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4")))), "{(-2*x2 - 3*x4 + 2, x2, 2*x4 + 5, x4)}")
  expect_identical(.pyrun(linsolve(list((Integer(-7) + Symbol("x1") + Symbol("x3") + Symbol("x4") + (Integer(2)*Symbol("x2"))), (Integer(-12) + Symbol("x1") + (Integer(-1)*Symbol("x4")) + (Integer(2)*Symbol("x2")) + (Integer(2)*Symbol("x3"))), (Integer(-4) + (Integer(2)*Symbol("x1")) + (Integer(4)*Symbol("x2")) + (Integer(6)*Symbol("x4")))), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"))), "{(-2*x2 - 3*x4 + 2, x2, 2*x4 + 5, x4)}")
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Integer(1), Integer(2), Integer(1), Integer(1)), list(Integer(1), Integer(2), Integer(2), Integer(-1)), list(Integer(2), Integer(4), Integer(0), Integer(6)))), Matrix(list(list(Integer(7)), list(Integer(12)), list(Integer(4))))), list(Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4")))), "{(-2*x2 - 3*x4 + 2, x2, 2*x4 + 5, x4)}")
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Integer(1), Integer(2), Integer(1), Integer(1)), list(Integer(1), Integer(2), Integer(2), Integer(-1)), list(Integer(2), Integer(4), Integer(0), Integer(6)))), Matrix(list(list(Integer(7)), list(Integer(12)), list(Integer(4))))), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"))), "{(-2*x2 - 3*x4 + 2, x2, 2*x4 + 5, x4)}")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); .D2 <- Dummy("x1"); .D3 <- Dummy("x2"); .D4 <- Dummy("x3"); linsolve(list(Matrix(list(list(Integer(1), Integer(2), Integer(1), Integer(1)), list(Integer(1), Integer(2), Integer(2), Integer(-1)), list(Integer(2), Integer(4), Integer(0), Integer(6)))), Matrix(list(list(Integer(7)), list(Integer(12)), list(Integer(4))))), .D1, .D2, .D3, .D4)})), "{(-2*_x1 - 3*_x3 + 2, _x1, 2*_x3 + 5, _x3)}")
  expect_identical(.pyrun(linsolve(list((Integer(-7) + Symbol("x1") + Symbol("x3") + Symbol("x4") + (Integer(2)*Symbol("x2"))), (Integer(-12) + Symbol("x1") + (Integer(-1)*Symbol("x4")) + (Integer(2)*Symbol("x2")) + (Integer(2)*Symbol("x3"))), (Integer(-4) + (Integer(2)*Symbol("x1")) + (Integer(4)*Symbol("x2")) + (Integer(6)*Symbol("x4")))))), "ERR ValueError: \\nWhen passing a system of equations, the explicit symbols for which a\\nsolution is being sought must be given as a sequence, too.")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); linsolve(.D1)})), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); .D2 <- Dummy("x1"); linsolve(.D1, .D2)})), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); .D2 <- Dummy("x1"); linsolve(list(Matrix(list(list(Integer(1), Integer(2), Integer(1), Integer(1)), list(Integer(1), Integer(2), Integer(2), Integer(-1)), list(Integer(2), Integer(4), Integer(0), Integer(6))))), .D1, .D2)})), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); .D2 <- Dummy("x1"); linsolve(Matrix(list(list(Integer(1), Integer(2), Integer(1), Integer(1)), list(Integer(1), Integer(2), Integer(2), Integer(-1)), list(Integer(2), Integer(4), Integer(0), Integer(6)))), Matrix(list(list(Integer(7)), list(Integer(12)), list(Integer(4)))), .D1, .D2)})), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({.D1 <- Dummy("x0"); linsolve(list(.D1), .D1, .D1)})), "ERR ValueError: duplicate symbols given")
  expect_identical(.pyrun(linsolve(list((Integer(-1) + Symbol("x") + Symbol("y")), (Integer(-3) + Symbol("y") + (Symbol("x"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: x^2")
  expect_identical(.pyrun(linsolve(list((Symbol("y") + cos(Symbol("x"))), (Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: cos(x)")
  expect_identical(.pyrun(linsolve(list((Integer(-1) + Symbol("x") + Symbol("z")), (Integer(-3) + Symbol("y") + (Symbol("x"))^(Integer(2)))), list(Symbol("z"), Symbol("y")))), "{(1 - x, 3 - x^2)}")
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Symbol("a"), Symbol("b")), list(Symbol("c"), Symbol("d")))), Matrix(list(list(Symbol("e")), list(Symbol("g"))))), list(Symbol("x"), Symbol("y")))), "{((-b*g + d*e)/(a*d - b*c), (a*g - c*e)/(a*d - b*c))}")
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Integer(1), Integer(2), Integer(3)), list(Integer(2), Integer(4), Integer(6)), list(Integer(3), Integer(6), Integer(9)))), Matrix(list(list(Integer(0)), list(Integer(0)), list(Integer(1))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "EmptySet")
  expect_identical(.pyrun(linsolve(Matrix(list(list((Integer(2)*sym_I*Symbol("J1")), (Integer(2)*sym_I*Symbol("J2")), (Integer(-2)*(Symbol("J1"))^(Integer(-1)))), list((Integer(-2)*sym_I*Symbol("J2")), (Integer(-2)*sym_I*Symbol("J1")), (Integer(2)*(Symbol("J2"))^(Integer(-1)))), list(Integer(0), Integer(2), (Integer(2)*sym_I*(Symbol("J1"))^(Integer(-1))*(Symbol("J2"))^(Integer(-1)))), list(Integer(2), Integer(0), Integer(0)))), Symbol("A"), Symbol("B"))), "{(0, I/(J1*J2))}")
  expect_identical(.pyrun(linsolve(Matrix(list(list(Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0)))), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d"), Symbol("e"))), "{(a, 0, c, 0, e)}")
  expect_identical(.pyrun(linsolve(Matrix(list(list(Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Symbol("_tau0")))))), "{(tau0, 0, tau1, _tau0, tau2)}")
  expect_identical(.pyrun(linsolve(Matrix(list(list(Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Symbol("tau0")))))), "{(tau00, 0, tau01, tau0, tau02)}")
  expect_identical(.pyrun(linsolve(Matrix(list(list(Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Symbol("tau1")))))), "{(tau00, 0, tau01, tau1, tau02)}")
  expect_identical(.pyrun(linsolve(list((Integer(-5) + Function("f")(Symbol("x")))), Function("f")(Symbol("x")))), "{(5,)}")
  expect_identical(.pyrun(linsolve(list(Eq(Symbol("x"), (Symbol("x") + Symbol("y")))), list(Symbol("x"), Symbol("y")))), "{(x, 0)}")
  expect_identical(.pyrun(linsolve(list(Eq((Symbol("x") + (Symbol("x")*Symbol("y"))), (Integer(1) + Symbol("y")))), list(Symbol("x")))), "{(1,)}")
  expect_identical(.pyrun(linsolve(list(Eq((Integer(1) + Symbol("y")), (Symbol("x") + (Symbol("x")*Symbol("y"))))), list(Symbol("x")))), "{(1,)}")
  expect_identical(.pyrun(linsolve(list(Eq((Symbol("x"))^(Integer(2)), (Symbol("y") + (Symbol("x"))^(Integer(2))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: x^2")
  expect_identical(.pyrun(linsolve(list(), list(Symbol("x")))), "EmptySet")
  expect_identical(.pyrun(linsolve(list(Integer(0)), list(Symbol("x")))), "{(x,)}")
  expect_identical(.pyrun(linsolve(list(Symbol("x")), list(Symbol("x"), Symbol("y")))), "{(0, y)}")
  expect_identical(.pyrun(linsolve(list(Symbol("x"), Integer(0)), list(Symbol("x"), Symbol("y")))), "{(0, y)}")
})

test_that("solveset: test_linsolve_large_sparse matches SymPy", {
  expect_identical(.pyrun(linsolve(list((Symbol("x0") + Symbol("y0")), (Integer(1) + Symbol("x0") + (Integer(-1)*Symbol("y0"))), (Symbol("x1") + Symbol("y1")), (Integer(1) + Symbol("x1") + (Integer(-1)*Symbol("y1"))), (Symbol("x2") + Symbol("y2")), (Integer(1) + Symbol("x2") + (Integer(-1)*Symbol("y2"))), (Symbol("x3") + Symbol("y3")), (Integer(1) + Symbol("x3") + (Integer(-1)*Symbol("y3"))), (Symbol("x4") + Symbol("y4")), (Integer(1) + Symbol("x4") + (Integer(-1)*Symbol("y4"))), (Symbol("x5") + Symbol("y5")), (Integer(1) + Symbol("x5") + (Integer(-1)*Symbol("y5"))), (Symbol("x6") + Symbol("y6")), (Integer(1) + Symbol("x6") + (Integer(-1)*Symbol("y6"))), (Symbol("x7") + Symbol("y7")), (Integer(1) + Symbol("x7") + (Integer(-1)*Symbol("y7"))), (Symbol("x8") + Symbol("y8")), (Integer(1) + Symbol("x8") + (Integer(-1)*Symbol("y8"))), (Symbol("x9") + Symbol("y9")), (Integer(1) + Symbol("x9") + (Integer(-1)*Symbol("y9"))), (Symbol("x10") + Symbol("y10")), (Integer(1) + Symbol("x10") + (Integer(-1)*Symbol("y10"))), (Symbol("x11") + Symbol("y11")), (Integer(1) + Symbol("x11") + (Integer(-1)*Symbol("y11"))), (Symbol("x12") + Symbol("y12")), (Integer(1) + Symbol("x12") + (Integer(-1)*Symbol("y12"))), (Symbol("x13") + Symbol("y13")), (Integer(1) + Symbol("x13") + (Integer(-1)*Symbol("y13"))), (Symbol("x14") + Symbol("y14")), (Integer(1) + Symbol("x14") + (Integer(-1)*Symbol("y14"))), (Symbol("x15") + Symbol("y15")), (Integer(1) + Symbol("x15") + (Integer(-1)*Symbol("y15"))), (Symbol("x16") + Symbol("y16")), (Integer(1) + Symbol("x16") + (Integer(-1)*Symbol("y16"))), (Symbol("x17") + Symbol("y17")), (Integer(1) + Symbol("x17") + (Integer(-1)*Symbol("y17"))), (Symbol("x18") + Symbol("y18")), (Integer(1) + Symbol("x18") + (Integer(-1)*Symbol("y18"))), (Symbol("x19") + Symbol("y19")), (Integer(1) + Symbol("x19") + (Integer(-1)*Symbol("y19"))), (Symbol("x20") + Symbol("y20")), (Integer(1) + Symbol("x20") + (Integer(-1)*Symbol("y20"))), (Symbol("x21") + Symbol("y21")), (Integer(1) + Symbol("x21") + (Integer(-1)*Symbol("y21"))), (Symbol("x22") + Symbol("y22")), (Integer(1) + Symbol("x22") + (Integer(-1)*Symbol("y22"))), (Symbol("x23") + Symbol("y23")), (Integer(1) + Symbol("x23") + (Integer(-1)*Symbol("y23"))), (Symbol("x24") + Symbol("y24")), (Integer(1) + Symbol("x24") + (Integer(-1)*Symbol("y24"))), (Symbol("x25") + Symbol("y25")), (Integer(1) + Symbol("x25") + (Integer(-1)*Symbol("y25"))), (Symbol("x26") + Symbol("y26")), (Integer(1) + Symbol("x26") + (Integer(-1)*Symbol("y26"))), (Symbol("x27") + Symbol("y27")), (Integer(1) + Symbol("x27") + (Integer(-1)*Symbol("y27"))), (Symbol("x28") + Symbol("y28")), (Integer(1) + Symbol("x28") + (Integer(-1)*Symbol("y28"))), (Symbol("x29") + Symbol("y29")), (Integer(1) + Symbol("x29") + (Integer(-1)*Symbol("y29"))), (Symbol("x30") + Symbol("y30")), (Integer(1) + Symbol("x30") + (Integer(-1)*Symbol("y30"))), (Symbol("x31") + Symbol("y31")), (Integer(1) + Symbol("x31") + (Integer(-1)*Symbol("y31"))), (Symbol("x32") + Symbol("y32")), (Integer(1) + Symbol("x32") + (Integer(-1)*Symbol("y32"))), (Symbol("x33") + Symbol("y33")), (Integer(1) + Symbol("x33") + (Integer(-1)*Symbol("y33"))), (Symbol("x34") + Symbol("y34")), (Integer(1) + Symbol("x34") + (Integer(-1)*Symbol("y34"))), (Symbol("x35") + Symbol("y35")), (Integer(1) + Symbol("x35") + (Integer(-1)*Symbol("y35"))), (Symbol("x36") + Symbol("y36")), (Integer(1) + Symbol("x36") + (Integer(-1)*Symbol("y36"))), (Symbol("x37") + Symbol("y37")), (Integer(1) + Symbol("x37") + (Integer(-1)*Symbol("y37"))), (Symbol("x38") + Symbol("y38")), (Integer(1) + Symbol("x38") + (Integer(-1)*Symbol("y38"))), (Symbol("x39") + Symbol("y39")), (Integer(1) + Symbol("x39") + (Integer(-1)*Symbol("y39"))), (Symbol("x40") + Symbol("y40")), (Integer(1) + Symbol("x40") + (Integer(-1)*Symbol("y40"))), (Symbol("x41") + Symbol("y41")), (Integer(1) + Symbol("x41") + (Integer(-1)*Symbol("y41"))), (Symbol("x42") + Symbol("y42")), (Integer(1) + Symbol("x42") + (Integer(-1)*Symbol("y42"))), (Symbol("x43") + Symbol("y43")), (Integer(1) + Symbol("x43") + (Integer(-1)*Symbol("y43"))), (Symbol("x44") + Symbol("y44")), (Integer(1) + Symbol("x44") + (Integer(-1)*Symbol("y44"))), (Symbol("x45") + Symbol("y45")), (Integer(1) + Symbol("x45") + (Integer(-1)*Symbol("y45"))), (Symbol("x46") + Symbol("y46")), (Integer(1) + Symbol("x46") + (Integer(-1)*Symbol("y46"))), (Symbol("x47") + Symbol("y47")), (Integer(1) + Symbol("x47") + (Integer(-1)*Symbol("y47"))), (Symbol("x48") + Symbol("y48")), (Integer(1) + Symbol("x48") + (Integer(-1)*Symbol("y48"))), (Symbol("x49") + Symbol("y49")), (Integer(1) + Symbol("x49") + (Integer(-1)*Symbol("y49"))), (Symbol("x50") + Symbol("y50")), (Integer(1) + Symbol("x50") + (Integer(-1)*Symbol("y50"))), (Symbol("x51") + Symbol("y51")), (Integer(1) + Symbol("x51") + (Integer(-1)*Symbol("y51"))), (Symbol("x52") + Symbol("y52")), (Integer(1) + Symbol("x52") + (Integer(-1)*Symbol("y52"))), (Symbol("x53") + Symbol("y53")), (Integer(1) + Symbol("x53") + (Integer(-1)*Symbol("y53"))), (Symbol("x54") + Symbol("y54")), (Integer(1) + Symbol("x54") + (Integer(-1)*Symbol("y54"))), (Symbol("x55") + Symbol("y55")), (Integer(1) + Symbol("x55") + (Integer(-1)*Symbol("y55"))), (Symbol("x56") + Symbol("y56")), (Integer(1) + Symbol("x56") + (Integer(-1)*Symbol("y56"))), (Symbol("x57") + Symbol("y57")), (Integer(1) + Symbol("x57") + (Integer(-1)*Symbol("y57"))), (Symbol("x58") + Symbol("y58")), (Integer(1) + Symbol("x58") + (Integer(-1)*Symbol("y58"))), (Symbol("x59") + Symbol("y59")), (Integer(1) + Symbol("x59") + (Integer(-1)*Symbol("y59"))), (Symbol("x60") + Symbol("y60")), (Integer(1) + Symbol("x60") + (Integer(-1)*Symbol("y60"))), (Symbol("x61") + Symbol("y61")), (Integer(1) + Symbol("x61") + (Integer(-1)*Symbol("y61"))), (Symbol("x62") + Symbol("y62")), (Integer(1) + Symbol("x62") + (Integer(-1)*Symbol("y62"))), (Symbol("x63") + Symbol("y63")), (Integer(1) + Symbol("x63") + (Integer(-1)*Symbol("y63"))), (Symbol("x64") + Symbol("y64")), (Integer(1) + Symbol("x64") + (Integer(-1)*Symbol("y64"))), (Symbol("x65") + Symbol("y65")), (Integer(1) + Symbol("x65") + (Integer(-1)*Symbol("y65"))), (Symbol("x66") + Symbol("y66")), (Integer(1) + Symbol("x66") + (Integer(-1)*Symbol("y66"))), (Symbol("x67") + Symbol("y67")), (Integer(1) + Symbol("x67") + (Integer(-1)*Symbol("y67"))), (Symbol("x68") + Symbol("y68")), (Integer(1) + Symbol("x68") + (Integer(-1)*Symbol("y68"))), (Symbol("x69") + Symbol("y69")), (Integer(1) + Symbol("x69") + (Integer(-1)*Symbol("y69"))), (Symbol("x70") + Symbol("y70")), (Integer(1) + Symbol("x70") + (Integer(-1)*Symbol("y70"))), (Symbol("x71") + Symbol("y71")), (Integer(1) + Symbol("x71") + (Integer(-1)*Symbol("y71"))), (Symbol("x72") + Symbol("y72")), (Integer(1) + Symbol("x72") + (Integer(-1)*Symbol("y72"))), (Symbol("x73") + Symbol("y73")), (Integer(1) + Symbol("x73") + (Integer(-1)*Symbol("y73"))), (Symbol("x74") + Symbol("y74")), (Integer(1) + Symbol("x74") + (Integer(-1)*Symbol("y74"))), (Symbol("x75") + Symbol("y75")), (Integer(1) + Symbol("x75") + (Integer(-1)*Symbol("y75"))), (Symbol("x76") + Symbol("y76")), (Integer(1) + Symbol("x76") + (Integer(-1)*Symbol("y76"))), (Symbol("x77") + Symbol("y77")), (Integer(1) + Symbol("x77") + (Integer(-1)*Symbol("y77"))), (Symbol("x78") + Symbol("y78")), (Integer(1) + Symbol("x78") + (Integer(-1)*Symbol("y78"))), (Symbol("x79") + Symbol("y79")), (Integer(1) + Symbol("x79") + (Integer(-1)*Symbol("y79"))), (Symbol("x80") + Symbol("y80")), (Integer(1) + Symbol("x80") + (Integer(-1)*Symbol("y80"))), (Symbol("x81") + Symbol("y81")), (Integer(1) + Symbol("x81") + (Integer(-1)*Symbol("y81"))), (Symbol("x82") + Symbol("y82")), (Integer(1) + Symbol("x82") + (Integer(-1)*Symbol("y82"))), (Symbol("x83") + Symbol("y83")), (Integer(1) + Symbol("x83") + (Integer(-1)*Symbol("y83"))), (Symbol("x84") + Symbol("y84")), (Integer(1) + Symbol("x84") + (Integer(-1)*Symbol("y84"))), (Symbol("x85") + Symbol("y85")), (Integer(1) + Symbol("x85") + (Integer(-1)*Symbol("y85"))), (Symbol("x86") + Symbol("y86")), (Integer(1) + Symbol("x86") + (Integer(-1)*Symbol("y86"))), (Symbol("x87") + Symbol("y87")), (Integer(1) + Symbol("x87") + (Integer(-1)*Symbol("y87"))), (Symbol("x88") + Symbol("y88")), (Integer(1) + Symbol("x88") + (Integer(-1)*Symbol("y88"))), (Symbol("x89") + Symbol("y89")), (Integer(1) + Symbol("x89") + (Integer(-1)*Symbol("y89"))), (Symbol("x90") + Symbol("y90")), (Integer(1) + Symbol("x90") + (Integer(-1)*Symbol("y90"))), (Symbol("x91") + Symbol("y91")), (Integer(1) + Symbol("x91") + (Integer(-1)*Symbol("y91"))), (Symbol("x92") + Symbol("y92")), (Integer(1) + Symbol("x92") + (Integer(-1)*Symbol("y92"))), (Symbol("x93") + Symbol("y93")), (Integer(1) + Symbol("x93") + (Integer(-1)*Symbol("y93"))), (Symbol("x94") + Symbol("y94")), (Integer(1) + Symbol("x94") + (Integer(-1)*Symbol("y94"))), (Symbol("x95") + Symbol("y95")), (Integer(1) + Symbol("x95") + (Integer(-1)*Symbol("y95"))), (Symbol("x96") + Symbol("y96")), (Integer(1) + Symbol("x96") + (Integer(-1)*Symbol("y96"))), (Symbol("x97") + Symbol("y97")), (Integer(1) + Symbol("x97") + (Integer(-1)*Symbol("y97"))), (Symbol("x98") + Symbol("y98")), (Integer(1) + Symbol("x98") + (Integer(-1)*Symbol("y98"))), (Symbol("x99") + Symbol("y99")), (Integer(1) + Symbol("x99") + (Integer(-1)*Symbol("y99"))), (Symbol("x100") + Symbol("y100")), (Integer(1) + Symbol("x100") + (Integer(-1)*Symbol("y100"))), (Symbol("x101") + Symbol("y101")), (Integer(1) + Symbol("x101") + (Integer(-1)*Symbol("y101"))), (Symbol("x102") + Symbol("y102")), (Integer(1) + Symbol("x102") + (Integer(-1)*Symbol("y102"))), (Symbol("x103") + Symbol("y103")), (Integer(1) + Symbol("x103") + (Integer(-1)*Symbol("y103"))), (Symbol("x104") + Symbol("y104")), (Integer(1) + Symbol("x104") + (Integer(-1)*Symbol("y104"))), (Symbol("x105") + Symbol("y105")), (Integer(1) + Symbol("x105") + (Integer(-1)*Symbol("y105"))), (Symbol("x106") + Symbol("y106")), (Integer(1) + Symbol("x106") + (Integer(-1)*Symbol("y106"))), (Symbol("x107") + Symbol("y107")), (Integer(1) + Symbol("x107") + (Integer(-1)*Symbol("y107"))), (Symbol("x108") + Symbol("y108")), (Integer(1) + Symbol("x108") + (Integer(-1)*Symbol("y108"))), (Symbol("x109") + Symbol("y109")), (Integer(1) + Symbol("x109") + (Integer(-1)*Symbol("y109"))), (Symbol("x110") + Symbol("y110")), (Integer(1) + Symbol("x110") + (Integer(-1)*Symbol("y110"))), (Symbol("x111") + Symbol("y111")), (Integer(1) + Symbol("x111") + (Integer(-1)*Symbol("y111"))), (Symbol("x112") + Symbol("y112")), (Integer(1) + Symbol("x112") + (Integer(-1)*Symbol("y112"))), (Symbol("x113") + Symbol("y113")), (Integer(1) + Symbol("x113") + (Integer(-1)*Symbol("y113"))), (Symbol("x114") + Symbol("y114")), (Integer(1) + Symbol("x114") + (Integer(-1)*Symbol("y114"))), (Symbol("x115") + Symbol("y115")), (Integer(1) + Symbol("x115") + (Integer(-1)*Symbol("y115"))), (Symbol("x116") + Symbol("y116")), (Integer(1) + Symbol("x116") + (Integer(-1)*Symbol("y116"))), (Symbol("x117") + Symbol("y117")), (Integer(1) + Symbol("x117") + (Integer(-1)*Symbol("y117"))), (Symbol("x118") + Symbol("y118")), (Integer(1) + Symbol("x118") + (Integer(-1)*Symbol("y118"))), (Symbol("x119") + Symbol("y119")), (Integer(1) + Symbol("x119") + (Integer(-1)*Symbol("y119"))), (Symbol("x120") + Symbol("y120")), (Integer(1) + Symbol("x120") + (Integer(-1)*Symbol("y120"))), (Symbol("x121") + Symbol("y121")), (Integer(1) + Symbol("x121") + (Integer(-1)*Symbol("y121"))), (Symbol("x122") + Symbol("y122")), (Integer(1) + Symbol("x122") + (Integer(-1)*Symbol("y122"))), (Symbol("x123") + Symbol("y123")), (Integer(1) + Symbol("x123") + (Integer(-1)*Symbol("y123"))), (Symbol("x124") + Symbol("y124")), (Integer(1) + Symbol("x124") + (Integer(-1)*Symbol("y124"))), (Symbol("x125") + Symbol("y125")), (Integer(1) + Symbol("x125") + (Integer(-1)*Symbol("y125"))), (Symbol("x126") + Symbol("y126")), (Integer(1) + Symbol("x126") + (Integer(-1)*Symbol("y126"))), (Symbol("x127") + Symbol("y127")), (Integer(1) + Symbol("x127") + (Integer(-1)*Symbol("y127"))), (Symbol("x128") + Symbol("y128")), (Integer(1) + Symbol("x128") + (Integer(-1)*Symbol("y128"))), (Symbol("x129") + Symbol("y129")), (Integer(1) + Symbol("x129") + (Integer(-1)*Symbol("y129"))), (Symbol("x130") + Symbol("y130")), (Integer(1) + Symbol("x130") + (Integer(-1)*Symbol("y130"))), (Symbol("x131") + Symbol("y131")), (Integer(1) + Symbol("x131") + (Integer(-1)*Symbol("y131"))), (Symbol("x132") + Symbol("y132")), (Integer(1) + Symbol("x132") + (Integer(-1)*Symbol("y132"))), (Symbol("x133") + Symbol("y133")), (Integer(1) + Symbol("x133") + (Integer(-1)*Symbol("y133"))), (Symbol("x134") + Symbol("y134")), (Integer(1) + Symbol("x134") + (Integer(-1)*Symbol("y134"))), (Symbol("x135") + Symbol("y135")), (Integer(1) + Symbol("x135") + (Integer(-1)*Symbol("y135"))), (Symbol("x136") + Symbol("y136")), (Integer(1) + Symbol("x136") + (Integer(-1)*Symbol("y136"))), (Symbol("x137") + Symbol("y137")), (Integer(1) + Symbol("x137") + (Integer(-1)*Symbol("y137"))), (Symbol("x138") + Symbol("y138")), (Integer(1) + Symbol("x138") + (Integer(-1)*Symbol("y138"))), (Symbol("x139") + Symbol("y139")), (Integer(1) + Symbol("x139") + (Integer(-1)*Symbol("y139"))), (Symbol("x140") + Symbol("y140")), (Integer(1) + Symbol("x140") + (Integer(-1)*Symbol("y140"))), (Symbol("x141") + Symbol("y141")), (Integer(1) + Symbol("x141") + (Integer(-1)*Symbol("y141"))), (Symbol("x142") + Symbol("y142")), (Integer(1) + Symbol("x142") + (Integer(-1)*Symbol("y142"))), (Symbol("x143") + Symbol("y143")), (Integer(1) + Symbol("x143") + (Integer(-1)*Symbol("y143"))), (Symbol("x144") + Symbol("y144")), (Integer(1) + Symbol("x144") + (Integer(-1)*Symbol("y144"))), (Symbol("x145") + Symbol("y145")), (Integer(1) + Symbol("x145") + (Integer(-1)*Symbol("y145"))), (Symbol("x146") + Symbol("y146")), (Integer(1) + Symbol("x146") + (Integer(-1)*Symbol("y146"))), (Symbol("x147") + Symbol("y147")), (Integer(1) + Symbol("x147") + (Integer(-1)*Symbol("y147"))), (Symbol("x148") + Symbol("y148")), (Integer(1) + Symbol("x148") + (Integer(-1)*Symbol("y148"))), (Symbol("x149") + Symbol("y149")), (Integer(1) + Symbol("x149") + (Integer(-1)*Symbol("y149"))), (Symbol("x150") + Symbol("y150")), (Integer(1) + Symbol("x150") + (Integer(-1)*Symbol("y150"))), (Symbol("x151") + Symbol("y151")), (Integer(1) + Symbol("x151") + (Integer(-1)*Symbol("y151"))), (Symbol("x152") + Symbol("y152")), (Integer(1) + Symbol("x152") + (Integer(-1)*Symbol("y152"))), (Symbol("x153") + Symbol("y153")), (Integer(1) + Symbol("x153") + (Integer(-1)*Symbol("y153"))), (Symbol("x154") + Symbol("y154")), (Integer(1) + Symbol("x154") + (Integer(-1)*Symbol("y154"))), (Symbol("x155") + Symbol("y155")), (Integer(1) + Symbol("x155") + (Integer(-1)*Symbol("y155"))), (Symbol("x156") + Symbol("y156")), (Integer(1) + Symbol("x156") + (Integer(-1)*Symbol("y156"))), (Symbol("x157") + Symbol("y157")), (Integer(1) + Symbol("x157") + (Integer(-1)*Symbol("y157"))), (Symbol("x158") + Symbol("y158")), (Integer(1) + Symbol("x158") + (Integer(-1)*Symbol("y158"))), (Symbol("x159") + Symbol("y159")), (Integer(1) + Symbol("x159") + (Integer(-1)*Symbol("y159"))), (Symbol("x160") + Symbol("y160")), (Integer(1) + Symbol("x160") + (Integer(-1)*Symbol("y160"))), (Symbol("x161") + Symbol("y161")), (Integer(1) + Symbol("x161") + (Integer(-1)*Symbol("y161"))), (Symbol("x162") + Symbol("y162")), (Integer(1) + Symbol("x162") + (Integer(-1)*Symbol("y162"))), (Symbol("x163") + Symbol("y163")), (Integer(1) + Symbol("x163") + (Integer(-1)*Symbol("y163"))), (Symbol("x164") + Symbol("y164")), (Integer(1) + Symbol("x164") + (Integer(-1)*Symbol("y164"))), (Symbol("x165") + Symbol("y165")), (Integer(1) + Symbol("x165") + (Integer(-1)*Symbol("y165"))), (Symbol("x166") + Symbol("y166")), (Integer(1) + Symbol("x166") + (Integer(-1)*Symbol("y166"))), (Symbol("x167") + Symbol("y167")), (Integer(1) + Symbol("x167") + (Integer(-1)*Symbol("y167"))), (Symbol("x168") + Symbol("y168")), (Integer(1) + Symbol("x168") + (Integer(-1)*Symbol("y168"))), (Symbol("x169") + Symbol("y169")), (Integer(1) + Symbol("x169") + (Integer(-1)*Symbol("y169"))), (Symbol("x170") + Symbol("y170")), (Integer(1) + Symbol("x170") + (Integer(-1)*Symbol("y170"))), (Symbol("x171") + Symbol("y171")), (Integer(1) + Symbol("x171") + (Integer(-1)*Symbol("y171"))), (Symbol("x172") + Symbol("y172")), (Integer(1) + Symbol("x172") + (Integer(-1)*Symbol("y172"))), (Symbol("x173") + Symbol("y173")), (Integer(1) + Symbol("x173") + (Integer(-1)*Symbol("y173"))), (Symbol("x174") + Symbol("y174")), (Integer(1) + Symbol("x174") + (Integer(-1)*Symbol("y174"))), (Symbol("x175") + Symbol("y175")), (Integer(1) + Symbol("x175") + (Integer(-1)*Symbol("y175"))), (Symbol("x176") + Symbol("y176")), (Integer(1) + Symbol("x176") + (Integer(-1)*Symbol("y176"))), (Symbol("x177") + Symbol("y177")), (Integer(1) + Symbol("x177") + (Integer(-1)*Symbol("y177"))), (Symbol("x178") + Symbol("y178")), (Integer(1) + Symbol("x178") + (Integer(-1)*Symbol("y178"))), (Symbol("x179") + Symbol("y179")), (Integer(1) + Symbol("x179") + (Integer(-1)*Symbol("y179"))), (Symbol("x180") + Symbol("y180")), (Integer(1) + Symbol("x180") + (Integer(-1)*Symbol("y180"))), (Symbol("x181") + Symbol("y181")), (Integer(1) + Symbol("x181") + (Integer(-1)*Symbol("y181"))), (Symbol("x182") + Symbol("y182")), (Integer(1) + Symbol("x182") + (Integer(-1)*Symbol("y182"))), (Symbol("x183") + Symbol("y183")), (Integer(1) + Symbol("x183") + (Integer(-1)*Symbol("y183"))), (Symbol("x184") + Symbol("y184")), (Integer(1) + Symbol("x184") + (Integer(-1)*Symbol("y184"))), (Symbol("x185") + Symbol("y185")), (Integer(1) + Symbol("x185") + (Integer(-1)*Symbol("y185"))), (Symbol("x186") + Symbol("y186")), (Integer(1) + Symbol("x186") + (Integer(-1)*Symbol("y186"))), (Symbol("x187") + Symbol("y187")), (Integer(1) + Symbol("x187") + (Integer(-1)*Symbol("y187"))), (Symbol("x188") + Symbol("y188")), (Integer(1) + Symbol("x188") + (Integer(-1)*Symbol("y188"))), (Symbol("x189") + Symbol("y189")), (Integer(1) + Symbol("x189") + (Integer(-1)*Symbol("y189"))), (Symbol("x190") + Symbol("y190")), (Integer(1) + Symbol("x190") + (Integer(-1)*Symbol("y190"))), (Symbol("x191") + Symbol("y191")), (Integer(1) + Symbol("x191") + (Integer(-1)*Symbol("y191"))), (Symbol("x192") + Symbol("y192")), (Integer(1) + Symbol("x192") + (Integer(-1)*Symbol("y192"))), (Symbol("x193") + Symbol("y193")), (Integer(1) + Symbol("x193") + (Integer(-1)*Symbol("y193"))), (Symbol("x194") + Symbol("y194")), (Integer(1) + Symbol("x194") + (Integer(-1)*Symbol("y194"))), (Symbol("x195") + Symbol("y195")), (Integer(1) + Symbol("x195") + (Integer(-1)*Symbol("y195"))), (Symbol("x196") + Symbol("y196")), (Integer(1) + Symbol("x196") + (Integer(-1)*Symbol("y196"))), (Symbol("x197") + Symbol("y197")), (Integer(1) + Symbol("x197") + (Integer(-1)*Symbol("y197"))), (Symbol("x198") + Symbol("y198")), (Integer(1) + Symbol("x198") + (Integer(-1)*Symbol("y198"))), (Symbol("x199") + Symbol("y199")), (Integer(1) + Symbol("x199") + (Integer(-1)*Symbol("y199"))), (Symbol("x200") + Symbol("y200")), (Integer(1) + Symbol("x200") + (Integer(-1)*Symbol("y200"))), (Symbol("x201") + Symbol("y201")), (Integer(1) + Symbol("x201") + (Integer(-1)*Symbol("y201"))), (Symbol("x202") + Symbol("y202")), (Integer(1) + Symbol("x202") + (Integer(-1)*Symbol("y202"))), (Symbol("x203") + Symbol("y203")), (Integer(1) + Symbol("x203") + (Integer(-1)*Symbol("y203"))), (Symbol("x204") + Symbol("y204")), (Integer(1) + Symbol("x204") + (Integer(-1)*Symbol("y204"))), (Symbol("x205") + Symbol("y205")), (Integer(1) + Symbol("x205") + (Integer(-1)*Symbol("y205"))), (Symbol("x206") + Symbol("y206")), (Integer(1) + Symbol("x206") + (Integer(-1)*Symbol("y206"))), (Symbol("x207") + Symbol("y207")), (Integer(1) + Symbol("x207") + (Integer(-1)*Symbol("y207"))), (Symbol("x208") + Symbol("y208")), (Integer(1) + Symbol("x208") + (Integer(-1)*Symbol("y208"))), (Symbol("x209") + Symbol("y209")), (Integer(1) + Symbol("x209") + (Integer(-1)*Symbol("y209"))), (Symbol("x210") + Symbol("y210")), (Integer(1) + Symbol("x210") + (Integer(-1)*Symbol("y210"))), (Symbol("x211") + Symbol("y211")), (Integer(1) + Symbol("x211") + (Integer(-1)*Symbol("y211"))), (Symbol("x212") + Symbol("y212")), (Integer(1) + Symbol("x212") + (Integer(-1)*Symbol("y212"))), (Symbol("x213") + Symbol("y213")), (Integer(1) + Symbol("x213") + (Integer(-1)*Symbol("y213"))), (Symbol("x214") + Symbol("y214")), (Integer(1) + Symbol("x214") + (Integer(-1)*Symbol("y214"))), (Symbol("x215") + Symbol("y215")), (Integer(1) + Symbol("x215") + (Integer(-1)*Symbol("y215"))), (Symbol("x216") + Symbol("y216")), (Integer(1) + Symbol("x216") + (Integer(-1)*Symbol("y216"))), (Symbol("x217") + Symbol("y217")), (Integer(1) + Symbol("x217") + (Integer(-1)*Symbol("y217"))), (Symbol("x218") + Symbol("y218")), (Integer(1) + Symbol("x218") + (Integer(-1)*Symbol("y218"))), (Symbol("x219") + Symbol("y219")), (Integer(1) + Symbol("x219") + (Integer(-1)*Symbol("y219"))), (Symbol("x220") + Symbol("y220")), (Integer(1) + Symbol("x220") + (Integer(-1)*Symbol("y220"))), (Symbol("x221") + Symbol("y221")), (Integer(1) + Symbol("x221") + (Integer(-1)*Symbol("y221"))), (Symbol("x222") + Symbol("y222")), (Integer(1) + Symbol("x222") + (Integer(-1)*Symbol("y222"))), (Symbol("x223") + Symbol("y223")), (Integer(1) + Symbol("x223") + (Integer(-1)*Symbol("y223"))), (Symbol("x224") + Symbol("y224")), (Integer(1) + Symbol("x224") + (Integer(-1)*Symbol("y224"))), (Symbol("x225") + Symbol("y225")), (Integer(1) + Symbol("x225") + (Integer(-1)*Symbol("y225"))), (Symbol("x226") + Symbol("y226")), (Integer(1) + Symbol("x226") + (Integer(-1)*Symbol("y226"))), (Symbol("x227") + Symbol("y227")), (Integer(1) + Symbol("x227") + (Integer(-1)*Symbol("y227"))), (Symbol("x228") + Symbol("y228")), (Integer(1) + Symbol("x228") + (Integer(-1)*Symbol("y228"))), (Symbol("x229") + Symbol("y229")), (Integer(1) + Symbol("x229") + (Integer(-1)*Symbol("y229"))), (Symbol("x230") + Symbol("y230")), (Integer(1) + Symbol("x230") + (Integer(-1)*Symbol("y230"))), (Symbol("x231") + Symbol("y231")), (Integer(1) + Symbol("x231") + (Integer(-1)*Symbol("y231"))), (Symbol("x232") + Symbol("y232")), (Integer(1) + Symbol("x232") + (Integer(-1)*Symbol("y232"))), (Symbol("x233") + Symbol("y233")), (Integer(1) + Symbol("x233") + (Integer(-1)*Symbol("y233"))), (Symbol("x234") + Symbol("y234")), (Integer(1) + Symbol("x234") + (Integer(-1)*Symbol("y234"))), (Symbol("x235") + Symbol("y235")), (Integer(1) + Symbol("x235") + (Integer(-1)*Symbol("y235"))), (Symbol("x236") + Symbol("y236")), (Integer(1) + Symbol("x236") + (Integer(-1)*Symbol("y236"))), (Symbol("x237") + Symbol("y237")), (Integer(1) + Symbol("x237") + (Integer(-1)*Symbol("y237"))), (Symbol("x238") + Symbol("y238")), (Integer(1) + Symbol("x238") + (Integer(-1)*Symbol("y238"))), (Symbol("x239") + Symbol("y239")), (Integer(1) + Symbol("x239") + (Integer(-1)*Symbol("y239"))), (Symbol("x240") + Symbol("y240")), (Integer(1) + Symbol("x240") + (Integer(-1)*Symbol("y240"))), (Symbol("x241") + Symbol("y241")), (Integer(1) + Symbol("x241") + (Integer(-1)*Symbol("y241"))), (Symbol("x242") + Symbol("y242")), (Integer(1) + Symbol("x242") + (Integer(-1)*Symbol("y242"))), (Symbol("x243") + Symbol("y243")), (Integer(1) + Symbol("x243") + (Integer(-1)*Symbol("y243"))), (Symbol("x244") + Symbol("y244")), (Integer(1) + Symbol("x244") + (Integer(-1)*Symbol("y244"))), (Symbol("x245") + Symbol("y245")), (Integer(1) + Symbol("x245") + (Integer(-1)*Symbol("y245"))), (Symbol("x246") + Symbol("y246")), (Integer(1) + Symbol("x246") + (Integer(-1)*Symbol("y246"))), (Symbol("x247") + Symbol("y247")), (Integer(1) + Symbol("x247") + (Integer(-1)*Symbol("y247"))), (Symbol("x248") + Symbol("y248")), (Integer(1) + Symbol("x248") + (Integer(-1)*Symbol("y248"))), (Symbol("x249") + Symbol("y249")), (Integer(1) + Symbol("x249") + (Integer(-1)*Symbol("y249"))), (Symbol("x250") + Symbol("y250")), (Integer(1) + Symbol("x250") + (Integer(-1)*Symbol("y250"))), (Symbol("x251") + Symbol("y251")), (Integer(1) + Symbol("x251") + (Integer(-1)*Symbol("y251"))), (Symbol("x252") + Symbol("y252")), (Integer(1) + Symbol("x252") + (Integer(-1)*Symbol("y252"))), (Symbol("x253") + Symbol("y253")), (Integer(1) + Symbol("x253") + (Integer(-1)*Symbol("y253"))), (Symbol("x254") + Symbol("y254")), (Integer(1) + Symbol("x254") + (Integer(-1)*Symbol("y254"))), (Symbol("x255") + Symbol("y255")), (Integer(1) + Symbol("x255") + (Integer(-1)*Symbol("y255"))), (Symbol("x256") + Symbol("y256")), (Integer(1) + Symbol("x256") + (Integer(-1)*Symbol("y256"))), (Symbol("x257") + Symbol("y257")), (Integer(1) + Symbol("x257") + (Integer(-1)*Symbol("y257"))), (Symbol("x258") + Symbol("y258")), (Integer(1) + Symbol("x258") + (Integer(-1)*Symbol("y258"))), (Symbol("x259") + Symbol("y259")), (Integer(1) + Symbol("x259") + (Integer(-1)*Symbol("y259"))), (Symbol("x260") + Symbol("y260")), (Integer(1) + Symbol("x260") + (Integer(-1)*Symbol("y260"))), (Symbol("x261") + Symbol("y261")), (Integer(1) + Symbol("x261") + (Integer(-1)*Symbol("y261"))), (Symbol("x262") + Symbol("y262")), (Integer(1) + Symbol("x262") + (Integer(-1)*Symbol("y262"))), (Symbol("x263") + Symbol("y263")), (Integer(1) + Symbol("x263") + (Integer(-1)*Symbol("y263"))), (Symbol("x264") + Symbol("y264")), (Integer(1) + Symbol("x264") + (Integer(-1)*Symbol("y264"))), (Symbol("x265") + Symbol("y265")), (Integer(1) + Symbol("x265") + (Integer(-1)*Symbol("y265"))), (Symbol("x266") + Symbol("y266")), (Integer(1) + Symbol("x266") + (Integer(-1)*Symbol("y266"))), (Symbol("x267") + Symbol("y267")), (Integer(1) + Symbol("x267") + (Integer(-1)*Symbol("y267"))), (Symbol("x268") + Symbol("y268")), (Integer(1) + Symbol("x268") + (Integer(-1)*Symbol("y268"))), (Symbol("x269") + Symbol("y269")), (Integer(1) + Symbol("x269") + (Integer(-1)*Symbol("y269"))), (Symbol("x270") + Symbol("y270")), (Integer(1) + Symbol("x270") + (Integer(-1)*Symbol("y270"))), (Symbol("x271") + Symbol("y271")), (Integer(1) + Symbol("x271") + (Integer(-1)*Symbol("y271"))), (Symbol("x272") + Symbol("y272")), (Integer(1) + Symbol("x272") + (Integer(-1)*Symbol("y272"))), (Symbol("x273") + Symbol("y273")), (Integer(1) + Symbol("x273") + (Integer(-1)*Symbol("y273"))), (Symbol("x274") + Symbol("y274")), (Integer(1) + Symbol("x274") + (Integer(-1)*Symbol("y274"))), (Symbol("x275") + Symbol("y275")), (Integer(1) + Symbol("x275") + (Integer(-1)*Symbol("y275"))), (Symbol("x276") + Symbol("y276")), (Integer(1) + Symbol("x276") + (Integer(-1)*Symbol("y276"))), (Symbol("x277") + Symbol("y277")), (Integer(1) + Symbol("x277") + (Integer(-1)*Symbol("y277"))), (Symbol("x278") + Symbol("y278")), (Integer(1) + Symbol("x278") + (Integer(-1)*Symbol("y278"))), (Symbol("x279") + Symbol("y279")), (Integer(1) + Symbol("x279") + (Integer(-1)*Symbol("y279"))), (Symbol("x280") + Symbol("y280")), (Integer(1) + Symbol("x280") + (Integer(-1)*Symbol("y280"))), (Symbol("x281") + Symbol("y281")), (Integer(1) + Symbol("x281") + (Integer(-1)*Symbol("y281"))), (Symbol("x282") + Symbol("y282")), (Integer(1) + Symbol("x282") + (Integer(-1)*Symbol("y282"))), (Symbol("x283") + Symbol("y283")), (Integer(1) + Symbol("x283") + (Integer(-1)*Symbol("y283"))), (Symbol("x284") + Symbol("y284")), (Integer(1) + Symbol("x284") + (Integer(-1)*Symbol("y284"))), (Symbol("x285") + Symbol("y285")), (Integer(1) + Symbol("x285") + (Integer(-1)*Symbol("y285"))), (Symbol("x286") + Symbol("y286")), (Integer(1) + Symbol("x286") + (Integer(-1)*Symbol("y286"))), (Symbol("x287") + Symbol("y287")), (Integer(1) + Symbol("x287") + (Integer(-1)*Symbol("y287"))), (Symbol("x288") + Symbol("y288")), (Integer(1) + Symbol("x288") + (Integer(-1)*Symbol("y288"))), (Symbol("x289") + Symbol("y289")), (Integer(1) + Symbol("x289") + (Integer(-1)*Symbol("y289"))), (Symbol("x290") + Symbol("y290")), (Integer(1) + Symbol("x290") + (Integer(-1)*Symbol("y290"))), (Symbol("x291") + Symbol("y291")), (Integer(1) + Symbol("x291") + (Integer(-1)*Symbol("y291"))), (Symbol("x292") + Symbol("y292")), (Integer(1) + Symbol("x292") + (Integer(-1)*Symbol("y292"))), (Symbol("x293") + Symbol("y293")), (Integer(1) + Symbol("x293") + (Integer(-1)*Symbol("y293"))), (Symbol("x294") + Symbol("y294")), (Integer(1) + Symbol("x294") + (Integer(-1)*Symbol("y294"))), (Symbol("x295") + Symbol("y295")), (Integer(1) + Symbol("x295") + (Integer(-1)*Symbol("y295"))), (Symbol("x296") + Symbol("y296")), (Integer(1) + Symbol("x296") + (Integer(-1)*Symbol("y296"))), (Symbol("x297") + Symbol("y297")), (Integer(1) + Symbol("x297") + (Integer(-1)*Symbol("y297"))), (Symbol("x298") + Symbol("y298")), (Integer(1) + Symbol("x298") + (Integer(-1)*Symbol("y298"))), (Symbol("x299") + Symbol("y299")), (Integer(1) + Symbol("x299") + (Integer(-1)*Symbol("y299"))), (Symbol("x300") + Symbol("y300")), (Integer(1) + Symbol("x300") + (Integer(-1)*Symbol("y300"))), (Symbol("x301") + Symbol("y301")), (Integer(1) + Symbol("x301") + (Integer(-1)*Symbol("y301"))), (Symbol("x302") + Symbol("y302")), (Integer(1) + Symbol("x302") + (Integer(-1)*Symbol("y302"))), (Symbol("x303") + Symbol("y303")), (Integer(1) + Symbol("x303") + (Integer(-1)*Symbol("y303"))), (Symbol("x304") + Symbol("y304")), (Integer(1) + Symbol("x304") + (Integer(-1)*Symbol("y304"))), (Symbol("x305") + Symbol("y305")), (Integer(1) + Symbol("x305") + (Integer(-1)*Symbol("y305"))), (Symbol("x306") + Symbol("y306")), (Integer(1) + Symbol("x306") + (Integer(-1)*Symbol("y306"))), (Symbol("x307") + Symbol("y307")), (Integer(1) + Symbol("x307") + (Integer(-1)*Symbol("y307"))), (Symbol("x308") + Symbol("y308")), (Integer(1) + Symbol("x308") + (Integer(-1)*Symbol("y308"))), (Symbol("x309") + Symbol("y309")), (Integer(1) + Symbol("x309") + (Integer(-1)*Symbol("y309"))), (Symbol("x310") + Symbol("y310")), (Integer(1) + Symbol("x310") + (Integer(-1)*Symbol("y310"))), (Symbol("x311") + Symbol("y311")), (Integer(1) + Symbol("x311") + (Integer(-1)*Symbol("y311"))), (Symbol("x312") + Symbol("y312")), (Integer(1) + Symbol("x312") + (Integer(-1)*Symbol("y312"))), (Symbol("x313") + Symbol("y313")), (Integer(1) + Symbol("x313") + (Integer(-1)*Symbol("y313"))), (Symbol("x314") + Symbol("y314")), (Integer(1) + Symbol("x314") + (Integer(-1)*Symbol("y314"))), (Symbol("x315") + Symbol("y315")), (Integer(1) + Symbol("x315") + (Integer(-1)*Symbol("y315"))), (Symbol("x316") + Symbol("y316")), (Integer(1) + Symbol("x316") + (Integer(-1)*Symbol("y316"))), (Symbol("x317") + Symbol("y317")), (Integer(1) + Symbol("x317") + (Integer(-1)*Symbol("y317"))), (Symbol("x318") + Symbol("y318")), (Integer(1) + Symbol("x318") + (Integer(-1)*Symbol("y318"))), (Symbol("x319") + Symbol("y319")), (Integer(1) + Symbol("x319") + (Integer(-1)*Symbol("y319"))), (Symbol("x320") + Symbol("y320")), (Integer(1) + Symbol("x320") + (Integer(-1)*Symbol("y320"))), (Symbol("x321") + Symbol("y321")), (Integer(1) + Symbol("x321") + (Integer(-1)*Symbol("y321"))), (Symbol("x322") + Symbol("y322")), (Integer(1) + Symbol("x322") + (Integer(-1)*Symbol("y322"))), (Symbol("x323") + Symbol("y323")), (Integer(1) + Symbol("x323") + (Integer(-1)*Symbol("y323"))), (Symbol("x324") + Symbol("y324")), (Integer(1) + Symbol("x324") + (Integer(-1)*Symbol("y324"))), (Symbol("x325") + Symbol("y325")), (Integer(1) + Symbol("x325") + (Integer(-1)*Symbol("y325"))), (Symbol("x326") + Symbol("y326")), (Integer(1) + Symbol("x326") + (Integer(-1)*Symbol("y326"))), (Symbol("x327") + Symbol("y327")), (Integer(1) + Symbol("x327") + (Integer(-1)*Symbol("y327"))), (Symbol("x328") + Symbol("y328")), (Integer(1) + Symbol("x328") + (Integer(-1)*Symbol("y328"))), (Symbol("x329") + Symbol("y329")), (Integer(1) + Symbol("x329") + (Integer(-1)*Symbol("y329"))), (Symbol("x330") + Symbol("y330")), (Integer(1) + Symbol("x330") + (Integer(-1)*Symbol("y330"))), (Symbol("x331") + Symbol("y331")), (Integer(1) + Symbol("x331") + (Integer(-1)*Symbol("y331"))), (Symbol("x332") + Symbol("y332")), (Integer(1) + Symbol("x332") + (Integer(-1)*Symbol("y332"))), (Symbol("x333") + Symbol("y333")), (Integer(1) + Symbol("x333") + (Integer(-1)*Symbol("y333"))), (Symbol("x334") + Symbol("y334")), (Integer(1) + Symbol("x334") + (Integer(-1)*Symbol("y334"))), (Symbol("x335") + Symbol("y335")), (Integer(1) + Symbol("x335") + (Integer(-1)*Symbol("y335"))), (Symbol("x336") + Symbol("y336")), (Integer(1) + Symbol("x336") + (Integer(-1)*Symbol("y336"))), (Symbol("x337") + Symbol("y337")), (Integer(1) + Symbol("x337") + (Integer(-1)*Symbol("y337"))), (Symbol("x338") + Symbol("y338")), (Integer(1) + Symbol("x338") + (Integer(-1)*Symbol("y338"))), (Symbol("x339") + Symbol("y339")), (Integer(1) + Symbol("x339") + (Integer(-1)*Symbol("y339"))), (Symbol("x340") + Symbol("y340")), (Integer(1) + Symbol("x340") + (Integer(-1)*Symbol("y340"))), (Symbol("x341") + Symbol("y341")), (Integer(1) + Symbol("x341") + (Integer(-1)*Symbol("y341"))), (Symbol("x342") + Symbol("y342")), (Integer(1) + Symbol("x342") + (Integer(-1)*Symbol("y342"))), (Symbol("x343") + Symbol("y343")), (Integer(1) + Symbol("x343") + (Integer(-1)*Symbol("y343"))), (Symbol("x344") + Symbol("y344")), (Integer(1) + Symbol("x344") + (Integer(-1)*Symbol("y344"))), (Symbol("x345") + Symbol("y345")), (Integer(1) + Symbol("x345") + (Integer(-1)*Symbol("y345"))), (Symbol("x346") + Symbol("y346")), (Integer(1) + Symbol("x346") + (Integer(-1)*Symbol("y346"))), (Symbol("x347") + Symbol("y347")), (Integer(1) + Symbol("x347") + (Integer(-1)*Symbol("y347"))), (Symbol("x348") + Symbol("y348")), (Integer(1) + Symbol("x348") + (Integer(-1)*Symbol("y348"))), (Symbol("x349") + Symbol("y349")), (Integer(1) + Symbol("x349") + (Integer(-1)*Symbol("y349"))), (Symbol("x350") + Symbol("y350")), (Integer(1) + Symbol("x350") + (Integer(-1)*Symbol("y350"))), (Symbol("x351") + Symbol("y351")), (Integer(1) + Symbol("x351") + (Integer(-1)*Symbol("y351"))), (Symbol("x352") + Symbol("y352")), (Integer(1) + Symbol("x352") + (Integer(-1)*Symbol("y352"))), (Symbol("x353") + Symbol("y353")), (Integer(1) + Symbol("x353") + (Integer(-1)*Symbol("y353"))), (Symbol("x354") + Symbol("y354")), (Integer(1) + Symbol("x354") + (Integer(-1)*Symbol("y354"))), (Symbol("x355") + Symbol("y355")), (Integer(1) + Symbol("x355") + (Integer(-1)*Symbol("y355"))), (Symbol("x356") + Symbol("y356")), (Integer(1) + Symbol("x356") + (Integer(-1)*Symbol("y356"))), (Symbol("x357") + Symbol("y357")), (Integer(1) + Symbol("x357") + (Integer(-1)*Symbol("y357"))), (Symbol("x358") + Symbol("y358")), (Integer(1) + Symbol("x358") + (Integer(-1)*Symbol("y358"))), (Symbol("x359") + Symbol("y359")), (Integer(1) + Symbol("x359") + (Integer(-1)*Symbol("y359"))), (Symbol("x360") + Symbol("y360")), (Integer(1) + Symbol("x360") + (Integer(-1)*Symbol("y360"))), (Symbol("x361") + Symbol("y361")), (Integer(1) + Symbol("x361") + (Integer(-1)*Symbol("y361"))), (Symbol("x362") + Symbol("y362")), (Integer(1) + Symbol("x362") + (Integer(-1)*Symbol("y362"))), (Symbol("x363") + Symbol("y363")), (Integer(1) + Symbol("x363") + (Integer(-1)*Symbol("y363"))), (Symbol("x364") + Symbol("y364")), (Integer(1) + Symbol("x364") + (Integer(-1)*Symbol("y364"))), (Symbol("x365") + Symbol("y365")), (Integer(1) + Symbol("x365") + (Integer(-1)*Symbol("y365"))), (Symbol("x366") + Symbol("y366")), (Integer(1) + Symbol("x366") + (Integer(-1)*Symbol("y366"))), (Symbol("x367") + Symbol("y367")), (Integer(1) + Symbol("x367") + (Integer(-1)*Symbol("y367"))), (Symbol("x368") + Symbol("y368")), (Integer(1) + Symbol("x368") + (Integer(-1)*Symbol("y368"))), (Symbol("x369") + Symbol("y369")), (Integer(1) + Symbol("x369") + (Integer(-1)*Symbol("y369"))), (Symbol("x370") + Symbol("y370")), (Integer(1) + Symbol("x370") + (Integer(-1)*Symbol("y370"))), (Symbol("x371") + Symbol("y371")), (Integer(1) + Symbol("x371") + (Integer(-1)*Symbol("y371"))), (Symbol("x372") + Symbol("y372")), (Integer(1) + Symbol("x372") + (Integer(-1)*Symbol("y372"))), (Symbol("x373") + Symbol("y373")), (Integer(1) + Symbol("x373") + (Integer(-1)*Symbol("y373"))), (Symbol("x374") + Symbol("y374")), (Integer(1) + Symbol("x374") + (Integer(-1)*Symbol("y374"))), (Symbol("x375") + Symbol("y375")), (Integer(1) + Symbol("x375") + (Integer(-1)*Symbol("y375"))), (Symbol("x376") + Symbol("y376")), (Integer(1) + Symbol("x376") + (Integer(-1)*Symbol("y376"))), (Symbol("x377") + Symbol("y377")), (Integer(1) + Symbol("x377") + (Integer(-1)*Symbol("y377"))), (Symbol("x378") + Symbol("y378")), (Integer(1) + Symbol("x378") + (Integer(-1)*Symbol("y378"))), (Symbol("x379") + Symbol("y379")), (Integer(1) + Symbol("x379") + (Integer(-1)*Symbol("y379"))), (Symbol("x380") + Symbol("y380")), (Integer(1) + Symbol("x380") + (Integer(-1)*Symbol("y380"))), (Symbol("x381") + Symbol("y381")), (Integer(1) + Symbol("x381") + (Integer(-1)*Symbol("y381"))), (Symbol("x382") + Symbol("y382")), (Integer(1) + Symbol("x382") + (Integer(-1)*Symbol("y382"))), (Symbol("x383") + Symbol("y383")), (Integer(1) + Symbol("x383") + (Integer(-1)*Symbol("y383"))), (Symbol("x384") + Symbol("y384")), (Integer(1) + Symbol("x384") + (Integer(-1)*Symbol("y384"))), (Symbol("x385") + Symbol("y385")), (Integer(1) + Symbol("x385") + (Integer(-1)*Symbol("y385"))), (Symbol("x386") + Symbol("y386")), (Integer(1) + Symbol("x386") + (Integer(-1)*Symbol("y386"))), (Symbol("x387") + Symbol("y387")), (Integer(1) + Symbol("x387") + (Integer(-1)*Symbol("y387"))), (Symbol("x388") + Symbol("y388")), (Integer(1) + Symbol("x388") + (Integer(-1)*Symbol("y388"))), (Symbol("x389") + Symbol("y389")), (Integer(1) + Symbol("x389") + (Integer(-1)*Symbol("y389"))), (Symbol("x390") + Symbol("y390")), (Integer(1) + Symbol("x390") + (Integer(-1)*Symbol("y390"))), (Symbol("x391") + Symbol("y391")), (Integer(1) + Symbol("x391") + (Integer(-1)*Symbol("y391"))), (Symbol("x392") + Symbol("y392")), (Integer(1) + Symbol("x392") + (Integer(-1)*Symbol("y392"))), (Symbol("x393") + Symbol("y393")), (Integer(1) + Symbol("x393") + (Integer(-1)*Symbol("y393"))), (Symbol("x394") + Symbol("y394")), (Integer(1) + Symbol("x394") + (Integer(-1)*Symbol("y394"))), (Symbol("x395") + Symbol("y395")), (Integer(1) + Symbol("x395") + (Integer(-1)*Symbol("y395"))), (Symbol("x396") + Symbol("y396")), (Integer(1) + Symbol("x396") + (Integer(-1)*Symbol("y396"))), (Symbol("x397") + Symbol("y397")), (Integer(1) + Symbol("x397") + (Integer(-1)*Symbol("y397"))), (Symbol("x398") + Symbol("y398")), (Integer(1) + Symbol("x398") + (Integer(-1)*Symbol("y398"))), (Symbol("x399") + Symbol("y399")), (Integer(1) + Symbol("x399") + (Integer(-1)*Symbol("y399"))), (Symbol("x400") + Symbol("y400")), (Integer(1) + Symbol("x400") + (Integer(-1)*Symbol("y400"))), (Symbol("x401") + Symbol("y401")), (Integer(1) + Symbol("x401") + (Integer(-1)*Symbol("y401"))), (Symbol("x402") + Symbol("y402")), (Integer(1) + Symbol("x402") + (Integer(-1)*Symbol("y402"))), (Symbol("x403") + Symbol("y403")), (Integer(1) + Symbol("x403") + (Integer(-1)*Symbol("y403"))), (Symbol("x404") + Symbol("y404")), (Integer(1) + Symbol("x404") + (Integer(-1)*Symbol("y404"))), (Symbol("x405") + Symbol("y405")), (Integer(1) + Symbol("x405") + (Integer(-1)*Symbol("y405"))), (Symbol("x406") + Symbol("y406")), (Integer(1) + Symbol("x406") + (Integer(-1)*Symbol("y406"))), (Symbol("x407") + Symbol("y407")), (Integer(1) + Symbol("x407") + (Integer(-1)*Symbol("y407"))), (Symbol("x408") + Symbol("y408")), (Integer(1) + Symbol("x408") + (Integer(-1)*Symbol("y408"))), (Symbol("x409") + Symbol("y409")), (Integer(1) + Symbol("x409") + (Integer(-1)*Symbol("y409"))), (Symbol("x410") + Symbol("y410")), (Integer(1) + Symbol("x410") + (Integer(-1)*Symbol("y410"))), (Symbol("x411") + Symbol("y411")), (Integer(1) + Symbol("x411") + (Integer(-1)*Symbol("y411"))), (Symbol("x412") + Symbol("y412")), (Integer(1) + Symbol("x412") + (Integer(-1)*Symbol("y412"))), (Symbol("x413") + Symbol("y413")), (Integer(1) + Symbol("x413") + (Integer(-1)*Symbol("y413"))), (Symbol("x414") + Symbol("y414")), (Integer(1) + Symbol("x414") + (Integer(-1)*Symbol("y414"))), (Symbol("x415") + Symbol("y415")), (Integer(1) + Symbol("x415") + (Integer(-1)*Symbol("y415"))), (Symbol("x416") + Symbol("y416")), (Integer(1) + Symbol("x416") + (Integer(-1)*Symbol("y416"))), (Symbol("x417") + Symbol("y417")), (Integer(1) + Symbol("x417") + (Integer(-1)*Symbol("y417"))), (Symbol("x418") + Symbol("y418")), (Integer(1) + Symbol("x418") + (Integer(-1)*Symbol("y418"))), (Symbol("x419") + Symbol("y419")), (Integer(1) + Symbol("x419") + (Integer(-1)*Symbol("y419"))), (Symbol("x420") + Symbol("y420")), (Integer(1) + Symbol("x420") + (Integer(-1)*Symbol("y420"))), (Symbol("x421") + Symbol("y421")), (Integer(1) + Symbol("x421") + (Integer(-1)*Symbol("y421"))), (Symbol("x422") + Symbol("y422")), (Integer(1) + Symbol("x422") + (Integer(-1)*Symbol("y422"))), (Symbol("x423") + Symbol("y423")), (Integer(1) + Symbol("x423") + (Integer(-1)*Symbol("y423"))), (Symbol("x424") + Symbol("y424")), (Integer(1) + Symbol("x424") + (Integer(-1)*Symbol("y424"))), (Symbol("x425") + Symbol("y425")), (Integer(1) + Symbol("x425") + (Integer(-1)*Symbol("y425"))), (Symbol("x426") + Symbol("y426")), (Integer(1) + Symbol("x426") + (Integer(-1)*Symbol("y426"))), (Symbol("x427") + Symbol("y427")), (Integer(1) + Symbol("x427") + (Integer(-1)*Symbol("y427"))), (Symbol("x428") + Symbol("y428")), (Integer(1) + Symbol("x428") + (Integer(-1)*Symbol("y428"))), (Symbol("x429") + Symbol("y429")), (Integer(1) + Symbol("x429") + (Integer(-1)*Symbol("y429"))), (Symbol("x430") + Symbol("y430")), (Integer(1) + Symbol("x430") + (Integer(-1)*Symbol("y430"))), (Symbol("x431") + Symbol("y431")), (Integer(1) + Symbol("x431") + (Integer(-1)*Symbol("y431"))), (Symbol("x432") + Symbol("y432")), (Integer(1) + Symbol("x432") + (Integer(-1)*Symbol("y432"))), (Symbol("x433") + Symbol("y433")), (Integer(1) + Symbol("x433") + (Integer(-1)*Symbol("y433"))), (Symbol("x434") + Symbol("y434")), (Integer(1) + Symbol("x434") + (Integer(-1)*Symbol("y434"))), (Symbol("x435") + Symbol("y435")), (Integer(1) + Symbol("x435") + (Integer(-1)*Symbol("y435"))), (Symbol("x436") + Symbol("y436")), (Integer(1) + Symbol("x436") + (Integer(-1)*Symbol("y436"))), (Symbol("x437") + Symbol("y437")), (Integer(1) + Symbol("x437") + (Integer(-1)*Symbol("y437"))), (Symbol("x438") + Symbol("y438")), (Integer(1) + Symbol("x438") + (Integer(-1)*Symbol("y438"))), (Symbol("x439") + Symbol("y439")), (Integer(1) + Symbol("x439") + (Integer(-1)*Symbol("y439"))), (Symbol("x440") + Symbol("y440")), (Integer(1) + Symbol("x440") + (Integer(-1)*Symbol("y440"))), (Symbol("x441") + Symbol("y441")), (Integer(1) + Symbol("x441") + (Integer(-1)*Symbol("y441"))), (Symbol("x442") + Symbol("y442")), (Integer(1) + Symbol("x442") + (Integer(-1)*Symbol("y442"))), (Symbol("x443") + Symbol("y443")), (Integer(1) + Symbol("x443") + (Integer(-1)*Symbol("y443"))), (Symbol("x444") + Symbol("y444")), (Integer(1) + Symbol("x444") + (Integer(-1)*Symbol("y444"))), (Symbol("x445") + Symbol("y445")), (Integer(1) + Symbol("x445") + (Integer(-1)*Symbol("y445"))), (Symbol("x446") + Symbol("y446")), (Integer(1) + Symbol("x446") + (Integer(-1)*Symbol("y446"))), (Symbol("x447") + Symbol("y447")), (Integer(1) + Symbol("x447") + (Integer(-1)*Symbol("y447"))), (Symbol("x448") + Symbol("y448")), (Integer(1) + Symbol("x448") + (Integer(-1)*Symbol("y448"))), (Symbol("x449") + Symbol("y449")), (Integer(1) + Symbol("x449") + (Integer(-1)*Symbol("y449"))), (Symbol("x450") + Symbol("y450")), (Integer(1) + Symbol("x450") + (Integer(-1)*Symbol("y450"))), (Symbol("x451") + Symbol("y451")), (Integer(1) + Symbol("x451") + (Integer(-1)*Symbol("y451"))), (Symbol("x452") + Symbol("y452")), (Integer(1) + Symbol("x452") + (Integer(-1)*Symbol("y452"))), (Symbol("x453") + Symbol("y453")), (Integer(1) + Symbol("x453") + (Integer(-1)*Symbol("y453"))), (Symbol("x454") + Symbol("y454")), (Integer(1) + Symbol("x454") + (Integer(-1)*Symbol("y454"))), (Symbol("x455") + Symbol("y455")), (Integer(1) + Symbol("x455") + (Integer(-1)*Symbol("y455"))), (Symbol("x456") + Symbol("y456")), (Integer(1) + Symbol("x456") + (Integer(-1)*Symbol("y456"))), (Symbol("x457") + Symbol("y457")), (Integer(1) + Symbol("x457") + (Integer(-1)*Symbol("y457"))), (Symbol("x458") + Symbol("y458")), (Integer(1) + Symbol("x458") + (Integer(-1)*Symbol("y458"))), (Symbol("x459") + Symbol("y459")), (Integer(1) + Symbol("x459") + (Integer(-1)*Symbol("y459"))), (Symbol("x460") + Symbol("y460")), (Integer(1) + Symbol("x460") + (Integer(-1)*Symbol("y460"))), (Symbol("x461") + Symbol("y461")), (Integer(1) + Symbol("x461") + (Integer(-1)*Symbol("y461"))), (Symbol("x462") + Symbol("y462")), (Integer(1) + Symbol("x462") + (Integer(-1)*Symbol("y462"))), (Symbol("x463") + Symbol("y463")), (Integer(1) + Symbol("x463") + (Integer(-1)*Symbol("y463"))), (Symbol("x464") + Symbol("y464")), (Integer(1) + Symbol("x464") + (Integer(-1)*Symbol("y464"))), (Symbol("x465") + Symbol("y465")), (Integer(1) + Symbol("x465") + (Integer(-1)*Symbol("y465"))), (Symbol("x466") + Symbol("y466")), (Integer(1) + Symbol("x466") + (Integer(-1)*Symbol("y466"))), (Symbol("x467") + Symbol("y467")), (Integer(1) + Symbol("x467") + (Integer(-1)*Symbol("y467"))), (Symbol("x468") + Symbol("y468")), (Integer(1) + Symbol("x468") + (Integer(-1)*Symbol("y468"))), (Symbol("x469") + Symbol("y469")), (Integer(1) + Symbol("x469") + (Integer(-1)*Symbol("y469"))), (Symbol("x470") + Symbol("y470")), (Integer(1) + Symbol("x470") + (Integer(-1)*Symbol("y470"))), (Symbol("x471") + Symbol("y471")), (Integer(1) + Symbol("x471") + (Integer(-1)*Symbol("y471"))), (Symbol("x472") + Symbol("y472")), (Integer(1) + Symbol("x472") + (Integer(-1)*Symbol("y472"))), (Symbol("x473") + Symbol("y473")), (Integer(1) + Symbol("x473") + (Integer(-1)*Symbol("y473"))), (Symbol("x474") + Symbol("y474")), (Integer(1) + Symbol("x474") + (Integer(-1)*Symbol("y474"))), (Symbol("x475") + Symbol("y475")), (Integer(1) + Symbol("x475") + (Integer(-1)*Symbol("y475"))), (Symbol("x476") + Symbol("y476")), (Integer(1) + Symbol("x476") + (Integer(-1)*Symbol("y476"))), (Symbol("x477") + Symbol("y477")), (Integer(1) + Symbol("x477") + (Integer(-1)*Symbol("y477"))), (Symbol("x478") + Symbol("y478")), (Integer(1) + Symbol("x478") + (Integer(-1)*Symbol("y478"))), (Symbol("x479") + Symbol("y479")), (Integer(1) + Symbol("x479") + (Integer(-1)*Symbol("y479"))), (Symbol("x480") + Symbol("y480")), (Integer(1) + Symbol("x480") + (Integer(-1)*Symbol("y480"))), (Symbol("x481") + Symbol("y481")), (Integer(1) + Symbol("x481") + (Integer(-1)*Symbol("y481"))), (Symbol("x482") + Symbol("y482")), (Integer(1) + Symbol("x482") + (Integer(-1)*Symbol("y482"))), (Symbol("x483") + Symbol("y483")), (Integer(1) + Symbol("x483") + (Integer(-1)*Symbol("y483"))), (Symbol("x484") + Symbol("y484")), (Integer(1) + Symbol("x484") + (Integer(-1)*Symbol("y484"))), (Symbol("x485") + Symbol("y485")), (Integer(1) + Symbol("x485") + (Integer(-1)*Symbol("y485"))), (Symbol("x486") + Symbol("y486")), (Integer(1) + Symbol("x486") + (Integer(-1)*Symbol("y486"))), (Symbol("x487") + Symbol("y487")), (Integer(1) + Symbol("x487") + (Integer(-1)*Symbol("y487"))), (Symbol("x488") + Symbol("y488")), (Integer(1) + Symbol("x488") + (Integer(-1)*Symbol("y488"))), (Symbol("x489") + Symbol("y489")), (Integer(1) + Symbol("x489") + (Integer(-1)*Symbol("y489"))), (Symbol("x490") + Symbol("y490")), (Integer(1) + Symbol("x490") + (Integer(-1)*Symbol("y490"))), (Symbol("x491") + Symbol("y491")), (Integer(1) + Symbol("x491") + (Integer(-1)*Symbol("y491"))), (Symbol("x492") + Symbol("y492")), (Integer(1) + Symbol("x492") + (Integer(-1)*Symbol("y492"))), (Symbol("x493") + Symbol("y493")), (Integer(1) + Symbol("x493") + (Integer(-1)*Symbol("y493"))), (Symbol("x494") + Symbol("y494")), (Integer(1) + Symbol("x494") + (Integer(-1)*Symbol("y494"))), (Symbol("x495") + Symbol("y495")), (Integer(1) + Symbol("x495") + (Integer(-1)*Symbol("y495"))), (Symbol("x496") + Symbol("y496")), (Integer(1) + Symbol("x496") + (Integer(-1)*Symbol("y496"))), (Symbol("x497") + Symbol("y497")), (Integer(1) + Symbol("x497") + (Integer(-1)*Symbol("y497"))), (Symbol("x498") + Symbol("y498")), (Integer(1) + Symbol("x498") + (Integer(-1)*Symbol("y498"))), (Symbol("x499") + Symbol("y499")), (Integer(1) + Symbol("x499") + (Integer(-1)*Symbol("y499")))), list(Symbol("x0"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15"), Symbol("x16"), Symbol("x17"), Symbol("x18"), Symbol("x19"), Symbol("x20"), Symbol("x21"), Symbol("x22"), Symbol("x23"), Symbol("x24"), Symbol("x25"), Symbol("x26"), Symbol("x27"), Symbol("x28"), Symbol("x29"), Symbol("x30"), Symbol("x31"), Symbol("x32"), Symbol("x33"), Symbol("x34"), Symbol("x35"), Symbol("x36"), Symbol("x37"), Symbol("x38"), Symbol("x39"), Symbol("x40"), Symbol("x41"), Symbol("x42"), Symbol("x43"), Symbol("x44"), Symbol("x45"), Symbol("x46"), Symbol("x47"), Symbol("x48"), Symbol("x49"), Symbol("x50"), Symbol("x51"), Symbol("x52"), Symbol("x53"), Symbol("x54"), Symbol("x55"), Symbol("x56"), Symbol("x57"), Symbol("x58"), Symbol("x59"), Symbol("x60"), Symbol("x61"), Symbol("x62"), Symbol("x63"), Symbol("x64"), Symbol("x65"), Symbol("x66"), Symbol("x67"), Symbol("x68"), Symbol("x69"), Symbol("x70"), Symbol("x71"), Symbol("x72"), Symbol("x73"), Symbol("x74"), Symbol("x75"), Symbol("x76"), Symbol("x77"), Symbol("x78"), Symbol("x79"), Symbol("x80"), Symbol("x81"), Symbol("x82"), Symbol("x83"), Symbol("x84"), Symbol("x85"), Symbol("x86"), Symbol("x87"), Symbol("x88"), Symbol("x89"), Symbol("x90"), Symbol("x91"), Symbol("x92"), Symbol("x93"), Symbol("x94"), Symbol("x95"), Symbol("x96"), Symbol("x97"), Symbol("x98"), Symbol("x99"), Symbol("x100"), Symbol("x101"), Symbol("x102"), Symbol("x103"), Symbol("x104"), Symbol("x105"), Symbol("x106"), Symbol("x107"), Symbol("x108"), Symbol("x109"), Symbol("x110"), Symbol("x111"), Symbol("x112"), Symbol("x113"), Symbol("x114"), Symbol("x115"), Symbol("x116"), Symbol("x117"), Symbol("x118"), Symbol("x119"), Symbol("x120"), Symbol("x121"), Symbol("x122"), Symbol("x123"), Symbol("x124"), Symbol("x125"), Symbol("x126"), Symbol("x127"), Symbol("x128"), Symbol("x129"), Symbol("x130"), Symbol("x131"), Symbol("x132"), Symbol("x133"), Symbol("x134"), Symbol("x135"), Symbol("x136"), Symbol("x137"), Symbol("x138"), Symbol("x139"), Symbol("x140"), Symbol("x141"), Symbol("x142"), Symbol("x143"), Symbol("x144"), Symbol("x145"), Symbol("x146"), Symbol("x147"), Symbol("x148"), Symbol("x149"), Symbol("x150"), Symbol("x151"), Symbol("x152"), Symbol("x153"), Symbol("x154"), Symbol("x155"), Symbol("x156"), Symbol("x157"), Symbol("x158"), Symbol("x159"), Symbol("x160"), Symbol("x161"), Symbol("x162"), Symbol("x163"), Symbol("x164"), Symbol("x165"), Symbol("x166"), Symbol("x167"), Symbol("x168"), Symbol("x169"), Symbol("x170"), Symbol("x171"), Symbol("x172"), Symbol("x173"), Symbol("x174"), Symbol("x175"), Symbol("x176"), Symbol("x177"), Symbol("x178"), Symbol("x179"), Symbol("x180"), Symbol("x181"), Symbol("x182"), Symbol("x183"), Symbol("x184"), Symbol("x185"), Symbol("x186"), Symbol("x187"), Symbol("x188"), Symbol("x189"), Symbol("x190"), Symbol("x191"), Symbol("x192"), Symbol("x193"), Symbol("x194"), Symbol("x195"), Symbol("x196"), Symbol("x197"), Symbol("x198"), Symbol("x199"), Symbol("x200"), Symbol("x201"), Symbol("x202"), Symbol("x203"), Symbol("x204"), Symbol("x205"), Symbol("x206"), Symbol("x207"), Symbol("x208"), Symbol("x209"), Symbol("x210"), Symbol("x211"), Symbol("x212"), Symbol("x213"), Symbol("x214"), Symbol("x215"), Symbol("x216"), Symbol("x217"), Symbol("x218"), Symbol("x219"), Symbol("x220"), Symbol("x221"), Symbol("x222"), Symbol("x223"), Symbol("x224"), Symbol("x225"), Symbol("x226"), Symbol("x227"), Symbol("x228"), Symbol("x229"), Symbol("x230"), Symbol("x231"), Symbol("x232"), Symbol("x233"), Symbol("x234"), Symbol("x235"), Symbol("x236"), Symbol("x237"), Symbol("x238"), Symbol("x239"), Symbol("x240"), Symbol("x241"), Symbol("x242"), Symbol("x243"), Symbol("x244"), Symbol("x245"), Symbol("x246"), Symbol("x247"), Symbol("x248"), Symbol("x249"), Symbol("x250"), Symbol("x251"), Symbol("x252"), Symbol("x253"), Symbol("x254"), Symbol("x255"), Symbol("x256"), Symbol("x257"), Symbol("x258"), Symbol("x259"), Symbol("x260"), Symbol("x261"), Symbol("x262"), Symbol("x263"), Symbol("x264"), Symbol("x265"), Symbol("x266"), Symbol("x267"), Symbol("x268"), Symbol("x269"), Symbol("x270"), Symbol("x271"), Symbol("x272"), Symbol("x273"), Symbol("x274"), Symbol("x275"), Symbol("x276"), Symbol("x277"), Symbol("x278"), Symbol("x279"), Symbol("x280"), Symbol("x281"), Symbol("x282"), Symbol("x283"), Symbol("x284"), Symbol("x285"), Symbol("x286"), Symbol("x287"), Symbol("x288"), Symbol("x289"), Symbol("x290"), Symbol("x291"), Symbol("x292"), Symbol("x293"), Symbol("x294"), Symbol("x295"), Symbol("x296"), Symbol("x297"), Symbol("x298"), Symbol("x299"), Symbol("x300"), Symbol("x301"), Symbol("x302"), Symbol("x303"), Symbol("x304"), Symbol("x305"), Symbol("x306"), Symbol("x307"), Symbol("x308"), Symbol("x309"), Symbol("x310"), Symbol("x311"), Symbol("x312"), Symbol("x313"), Symbol("x314"), Symbol("x315"), Symbol("x316"), Symbol("x317"), Symbol("x318"), Symbol("x319"), Symbol("x320"), Symbol("x321"), Symbol("x322"), Symbol("x323"), Symbol("x324"), Symbol("x325"), Symbol("x326"), Symbol("x327"), Symbol("x328"), Symbol("x329"), Symbol("x330"), Symbol("x331"), Symbol("x332"), Symbol("x333"), Symbol("x334"), Symbol("x335"), Symbol("x336"), Symbol("x337"), Symbol("x338"), Symbol("x339"), Symbol("x340"), Symbol("x341"), Symbol("x342"), Symbol("x343"), Symbol("x344"), Symbol("x345"), Symbol("x346"), Symbol("x347"), Symbol("x348"), Symbol("x349"), Symbol("x350"), Symbol("x351"), Symbol("x352"), Symbol("x353"), Symbol("x354"), Symbol("x355"), Symbol("x356"), Symbol("x357"), Symbol("x358"), Symbol("x359"), Symbol("x360"), Symbol("x361"), Symbol("x362"), Symbol("x363"), Symbol("x364"), Symbol("x365"), Symbol("x366"), Symbol("x367"), Symbol("x368"), Symbol("x369"), Symbol("x370"), Symbol("x371"), Symbol("x372"), Symbol("x373"), Symbol("x374"), Symbol("x375"), Symbol("x376"), Symbol("x377"), Symbol("x378"), Symbol("x379"), Symbol("x380"), Symbol("x381"), Symbol("x382"), Symbol("x383"), Symbol("x384"), Symbol("x385"), Symbol("x386"), Symbol("x387"), Symbol("x388"), Symbol("x389"), Symbol("x390"), Symbol("x391"), Symbol("x392"), Symbol("x393"), Symbol("x394"), Symbol("x395"), Symbol("x396"), Symbol("x397"), Symbol("x398"), Symbol("x399"), Symbol("x400"), Symbol("x401"), Symbol("x402"), Symbol("x403"), Symbol("x404"), Symbol("x405"), Symbol("x406"), Symbol("x407"), Symbol("x408"), Symbol("x409"), Symbol("x410"), Symbol("x411"), Symbol("x412"), Symbol("x413"), Symbol("x414"), Symbol("x415"), Symbol("x416"), Symbol("x417"), Symbol("x418"), Symbol("x419"), Symbol("x420"), Symbol("x421"), Symbol("x422"), Symbol("x423"), Symbol("x424"), Symbol("x425"), Symbol("x426"), Symbol("x427"), Symbol("x428"), Symbol("x429"), Symbol("x430"), Symbol("x431"), Symbol("x432"), Symbol("x433"), Symbol("x434"), Symbol("x435"), Symbol("x436"), Symbol("x437"), Symbol("x438"), Symbol("x439"), Symbol("x440"), Symbol("x441"), Symbol("x442"), Symbol("x443"), Symbol("x444"), Symbol("x445"), Symbol("x446"), Symbol("x447"), Symbol("x448"), Symbol("x449"), Symbol("x450"), Symbol("x451"), Symbol("x452"), Symbol("x453"), Symbol("x454"), Symbol("x455"), Symbol("x456"), Symbol("x457"), Symbol("x458"), Symbol("x459"), Symbol("x460"), Symbol("x461"), Symbol("x462"), Symbol("x463"), Symbol("x464"), Symbol("x465"), Symbol("x466"), Symbol("x467"), Symbol("x468"), Symbol("x469"), Symbol("x470"), Symbol("x471"), Symbol("x472"), Symbol("x473"), Symbol("x474"), Symbol("x475"), Symbol("x476"), Symbol("x477"), Symbol("x478"), Symbol("x479"), Symbol("x480"), Symbol("x481"), Symbol("x482"), Symbol("x483"), Symbol("x484"), Symbol("x485"), Symbol("x486"), Symbol("x487"), Symbol("x488"), Symbol("x489"), Symbol("x490"), Symbol("x491"), Symbol("x492"), Symbol("x493"), Symbol("x494"), Symbol("x495"), Symbol("x496"), Symbol("x497"), Symbol("x498"), Symbol("x499"), Symbol("y0"), Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("y7"), Symbol("y8"), Symbol("y9"), Symbol("y10"), Symbol("y11"), Symbol("y12"), Symbol("y13"), Symbol("y14"), Symbol("y15"), Symbol("y16"), Symbol("y17"), Symbol("y18"), Symbol("y19"), Symbol("y20"), Symbol("y21"), Symbol("y22"), Symbol("y23"), Symbol("y24"), Symbol("y25"), Symbol("y26"), Symbol("y27"), Symbol("y28"), Symbol("y29"), Symbol("y30"), Symbol("y31"), Symbol("y32"), Symbol("y33"), Symbol("y34"), Symbol("y35"), Symbol("y36"), Symbol("y37"), Symbol("y38"), Symbol("y39"), Symbol("y40"), Symbol("y41"), Symbol("y42"), Symbol("y43"), Symbol("y44"), Symbol("y45"), Symbol("y46"), Symbol("y47"), Symbol("y48"), Symbol("y49"), Symbol("y50"), Symbol("y51"), Symbol("y52"), Symbol("y53"), Symbol("y54"), Symbol("y55"), Symbol("y56"), Symbol("y57"), Symbol("y58"), Symbol("y59"), Symbol("y60"), Symbol("y61"), Symbol("y62"), Symbol("y63"), Symbol("y64"), Symbol("y65"), Symbol("y66"), Symbol("y67"), Symbol("y68"), Symbol("y69"), Symbol("y70"), Symbol("y71"), Symbol("y72"), Symbol("y73"), Symbol("y74"), Symbol("y75"), Symbol("y76"), Symbol("y77"), Symbol("y78"), Symbol("y79"), Symbol("y80"), Symbol("y81"), Symbol("y82"), Symbol("y83"), Symbol("y84"), Symbol("y85"), Symbol("y86"), Symbol("y87"), Symbol("y88"), Symbol("y89"), Symbol("y90"), Symbol("y91"), Symbol("y92"), Symbol("y93"), Symbol("y94"), Symbol("y95"), Symbol("y96"), Symbol("y97"), Symbol("y98"), Symbol("y99"), Symbol("y100"), Symbol("y101"), Symbol("y102"), Symbol("y103"), Symbol("y104"), Symbol("y105"), Symbol("y106"), Symbol("y107"), Symbol("y108"), Symbol("y109"), Symbol("y110"), Symbol("y111"), Symbol("y112"), Symbol("y113"), Symbol("y114"), Symbol("y115"), Symbol("y116"), Symbol("y117"), Symbol("y118"), Symbol("y119"), Symbol("y120"), Symbol("y121"), Symbol("y122"), Symbol("y123"), Symbol("y124"), Symbol("y125"), Symbol("y126"), Symbol("y127"), Symbol("y128"), Symbol("y129"), Symbol("y130"), Symbol("y131"), Symbol("y132"), Symbol("y133"), Symbol("y134"), Symbol("y135"), Symbol("y136"), Symbol("y137"), Symbol("y138"), Symbol("y139"), Symbol("y140"), Symbol("y141"), Symbol("y142"), Symbol("y143"), Symbol("y144"), Symbol("y145"), Symbol("y146"), Symbol("y147"), Symbol("y148"), Symbol("y149"), Symbol("y150"), Symbol("y151"), Symbol("y152"), Symbol("y153"), Symbol("y154"), Symbol("y155"), Symbol("y156"), Symbol("y157"), Symbol("y158"), Symbol("y159"), Symbol("y160"), Symbol("y161"), Symbol("y162"), Symbol("y163"), Symbol("y164"), Symbol("y165"), Symbol("y166"), Symbol("y167"), Symbol("y168"), Symbol("y169"), Symbol("y170"), Symbol("y171"), Symbol("y172"), Symbol("y173"), Symbol("y174"), Symbol("y175"), Symbol("y176"), Symbol("y177"), Symbol("y178"), Symbol("y179"), Symbol("y180"), Symbol("y181"), Symbol("y182"), Symbol("y183"), Symbol("y184"), Symbol("y185"), Symbol("y186"), Symbol("y187"), Symbol("y188"), Symbol("y189"), Symbol("y190"), Symbol("y191"), Symbol("y192"), Symbol("y193"), Symbol("y194"), Symbol("y195"), Symbol("y196"), Symbol("y197"), Symbol("y198"), Symbol("y199"), Symbol("y200"), Symbol("y201"), Symbol("y202"), Symbol("y203"), Symbol("y204"), Symbol("y205"), Symbol("y206"), Symbol("y207"), Symbol("y208"), Symbol("y209"), Symbol("y210"), Symbol("y211"), Symbol("y212"), Symbol("y213"), Symbol("y214"), Symbol("y215"), Symbol("y216"), Symbol("y217"), Symbol("y218"), Symbol("y219"), Symbol("y220"), Symbol("y221"), Symbol("y222"), Symbol("y223"), Symbol("y224"), Symbol("y225"), Symbol("y226"), Symbol("y227"), Symbol("y228"), Symbol("y229"), Symbol("y230"), Symbol("y231"), Symbol("y232"), Symbol("y233"), Symbol("y234"), Symbol("y235"), Symbol("y236"), Symbol("y237"), Symbol("y238"), Symbol("y239"), Symbol("y240"), Symbol("y241"), Symbol("y242"), Symbol("y243"), Symbol("y244"), Symbol("y245"), Symbol("y246"), Symbol("y247"), Symbol("y248"), Symbol("y249"), Symbol("y250"), Symbol("y251"), Symbol("y252"), Symbol("y253"), Symbol("y254"), Symbol("y255"), Symbol("y256"), Symbol("y257"), Symbol("y258"), Symbol("y259"), Symbol("y260"), Symbol("y261"), Symbol("y262"), Symbol("y263"), Symbol("y264"), Symbol("y265"), Symbol("y266"), Symbol("y267"), Symbol("y268"), Symbol("y269"), Symbol("y270"), Symbol("y271"), Symbol("y272"), Symbol("y273"), Symbol("y274"), Symbol("y275"), Symbol("y276"), Symbol("y277"), Symbol("y278"), Symbol("y279"), Symbol("y280"), Symbol("y281"), Symbol("y282"), Symbol("y283"), Symbol("y284"), Symbol("y285"), Symbol("y286"), Symbol("y287"), Symbol("y288"), Symbol("y289"), Symbol("y290"), Symbol("y291"), Symbol("y292"), Symbol("y293"), Symbol("y294"), Symbol("y295"), Symbol("y296"), Symbol("y297"), Symbol("y298"), Symbol("y299"), Symbol("y300"), Symbol("y301"), Symbol("y302"), Symbol("y303"), Symbol("y304"), Symbol("y305"), Symbol("y306"), Symbol("y307"), Symbol("y308"), Symbol("y309"), Symbol("y310"), Symbol("y311"), Symbol("y312"), Symbol("y313"), Symbol("y314"), Symbol("y315"), Symbol("y316"), Symbol("y317"), Symbol("y318"), Symbol("y319"), Symbol("y320"), Symbol("y321"), Symbol("y322"), Symbol("y323"), Symbol("y324"), Symbol("y325"), Symbol("y326"), Symbol("y327"), Symbol("y328"), Symbol("y329"), Symbol("y330"), Symbol("y331"), Symbol("y332"), Symbol("y333"), Symbol("y334"), Symbol("y335"), Symbol("y336"), Symbol("y337"), Symbol("y338"), Symbol("y339"), Symbol("y340"), Symbol("y341"), Symbol("y342"), Symbol("y343"), Symbol("y344"), Symbol("y345"), Symbol("y346"), Symbol("y347"), Symbol("y348"), Symbol("y349"), Symbol("y350"), Symbol("y351"), Symbol("y352"), Symbol("y353"), Symbol("y354"), Symbol("y355"), Symbol("y356"), Symbol("y357"), Symbol("y358"), Symbol("y359"), Symbol("y360"), Symbol("y361"), Symbol("y362"), Symbol("y363"), Symbol("y364"), Symbol("y365"), Symbol("y366"), Symbol("y367"), Symbol("y368"), Symbol("y369"), Symbol("y370"), Symbol("y371"), Symbol("y372"), Symbol("y373"), Symbol("y374"), Symbol("y375"), Symbol("y376"), Symbol("y377"), Symbol("y378"), Symbol("y379"), Symbol("y380"), Symbol("y381"), Symbol("y382"), Symbol("y383"), Symbol("y384"), Symbol("y385"), Symbol("y386"), Symbol("y387"), Symbol("y388"), Symbol("y389"), Symbol("y390"), Symbol("y391"), Symbol("y392"), Symbol("y393"), Symbol("y394"), Symbol("y395"), Symbol("y396"), Symbol("y397"), Symbol("y398"), Symbol("y399"), Symbol("y400"), Symbol("y401"), Symbol("y402"), Symbol("y403"), Symbol("y404"), Symbol("y405"), Symbol("y406"), Symbol("y407"), Symbol("y408"), Symbol("y409"), Symbol("y410"), Symbol("y411"), Symbol("y412"), Symbol("y413"), Symbol("y414"), Symbol("y415"), Symbol("y416"), Symbol("y417"), Symbol("y418"), Symbol("y419"), Symbol("y420"), Symbol("y421"), Symbol("y422"), Symbol("y423"), Symbol("y424"), Symbol("y425"), Symbol("y426"), Symbol("y427"), Symbol("y428"), Symbol("y429"), Symbol("y430"), Symbol("y431"), Symbol("y432"), Symbol("y433"), Symbol("y434"), Symbol("y435"), Symbol("y436"), Symbol("y437"), Symbol("y438"), Symbol("y439"), Symbol("y440"), Symbol("y441"), Symbol("y442"), Symbol("y443"), Symbol("y444"), Symbol("y445"), Symbol("y446"), Symbol("y447"), Symbol("y448"), Symbol("y449"), Symbol("y450"), Symbol("y451"), Symbol("y452"), Symbol("y453"), Symbol("y454"), Symbol("y455"), Symbol("y456"), Symbol("y457"), Symbol("y458"), Symbol("y459"), Symbol("y460"), Symbol("y461"), Symbol("y462"), Symbol("y463"), Symbol("y464"), Symbol("y465"), Symbol("y466"), Symbol("y467"), Symbol("y468"), Symbol("y469"), Symbol("y470"), Symbol("y471"), Symbol("y472"), Symbol("y473"), Symbol("y474"), Symbol("y475"), Symbol("y476"), Symbol("y477"), Symbol("y478"), Symbol("y479"), Symbol("y480"), Symbol("y481"), Symbol("y482"), Symbol("y483"), Symbol("y484"), Symbol("y485"), Symbol("y486"), Symbol("y487"), Symbol("y488"), Symbol("y489"), Symbol("y490"), Symbol("y491"), Symbol("y492"), Symbol("y493"), Symbol("y494"), Symbol("y495"), Symbol("y496"), Symbol("y497"), Symbol("y498"), Symbol("y499")))), "{(-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2)}")
})

test_that("solveset: test_linsolve_immutable matches SymPy", {
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Integer(1), Integer(1), Integer(2)), list(Integer(0), Integer(1), Integer(2)), list(Integer(0), Integer(0), Integer(1)))), Matrix(list(list(Integer(2)), list(Integer(1)), list(Integer(-1))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(1, 3, -1)}")
  expect_identical(.pyrun(linsolve(Matrix(list(list(Integer(1), Integer(1), Integer(7)), list(Integer(1), Integer(-1), Integer(3)))))), "{(5, 2)}")
})

test_that("solveset: test_solve_decomposition matches SymPy", {
  expect_identical(.pyrun(solve_decomposition((Integer(-6) + (Integer(-6)*exp((Integer(2)*Symbol("x")))) + (Integer(11)*exp(Symbol("x"))) + exp((Integer(3)*Symbol("x")))), Symbol("x"), Reals)), "{0, log(2), log(3)}")
  expect_identical(.pyrun(solve_decomposition((Integer(1) + (sin(Symbol("x")))^(Integer(2)) + (Integer(-2)*sin(Symbol("x")))), Symbol("x"), Reals)), "ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)")
  expect_identical(.pyrun(solve_decomposition(((sin(Symbol("x")))^(Integer(2)) + (Integer(-1)*sin(Symbol("x")))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers))")
  expect_identical(.pyrun(solve_decomposition(sin((Integer(1) + Symbol("x"))), Symbol("x"), Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi - 1), Integers), ImageSet(Lambda(_n, 2*_n*pi - 1 + pi), Integers))")
  expect_identical(.pyrun(solve_decomposition((Integer(-1) + exp((Integer(2) + Symbol("x")))), Symbol("x"), Reals)), "{-2}")
  expect_identical(.pyrun(solve_decomposition((log(Symbol("x")))^(Integer(-1)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solve_decomposition((Symbol("x"))^(Integer(-1)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solve_decomposition(Symbol("x"), Symbol("x"), Interval(Integer(1), Integer(2), FALSE, FALSE))), "EmptySet")
})

test_that("solveset: test_nonlinsolve_basic matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(), list())), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Symbol("x")), list(Symbol("x"), Symbol("y")))), "{(0, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Symbol("x"), (Integer(-5) + Symbol("y") + (Integer(-1)*Symbol("x")))), list(Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sin(Symbol("x")))), list(Symbol("x")))), "{(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers),)}")
  expect_identical(.pyrun(nonlinsolve(list(sin(Symbol("x")), (Integer(-1) + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, 2*_n*pi), Integers), 1), (ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), 1)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list(Symbol("x")))), "{(-1,), (1,)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x") + (Integer(-1)*Symbol("y"))), Integer(0)), Symbol("x"), Symbol("y"))), "{(y, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Integer(0), (Symbol("x") + (Integer(-1)*Symbol("y")))), Symbol("x"), Symbol("y"))), "{(y, y)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x") + (Integer(-1)*Symbol("y"))), (Symbol("x") + (Integer(-1)*Symbol("y")))), Symbol("x"), Symbol("y"))), "{(y, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Symbol("x"), Integer(0)), Symbol("x"), Symbol("y"))), "{(0, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Function("f")(Symbol("x")), Integer(0)), Function("f")(Symbol("x")), Symbol("y"))), "{(0, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Function("f")(Symbol("x")), Integer(0)), Function("f")(Symbol("x")), Function("f")(Symbol("y")))), "{(0, f(y))}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list(sin(Symbol("x"))))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), sin(Symbol("x")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), Integer(1))), "{(x^2,)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), (Symbol("x") + Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Integer(1)), list(Symbol("x")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Integer(5)), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Integer(-1)), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Rational(1, 2)), list(Symbol("x")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(2))^(Rational(1, 2))), list(Symbol("x")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(Integer(0)), list(Symbol("x"), Symbol("y")))), "{(x, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Integer(0), Integer(0)), list(Symbol("x"), Symbol("y")))), "{(x, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Integer(0)), list(Symbol("x"), Symbol("y")))), "{(x, y)}")
  expect_identical(.pyrun(nonlinsolve(list(Integer(0)), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(x, y, z)}")
  expect_identical(.pyrun(nonlinsolve(list(Eq(Integer(1), (Symbol("x") + Symbol("y"))), Eq(Integer(1), (Integer(-1) + Symbol("y") + (Integer(-1)*Symbol("x")))), Eq(Integer(1), (Integer(-1) + Symbol("y") + (Integer(-1)*Symbol("x"))))), Symbol("x"), Symbol("y"))), "{(-1/2, 3/2)}")
})

test_that("solveset: test_nonlinsolve_abs matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*Symbol("y")) + Abs(Symbol("x")))), Symbol("x"), Symbol("y"))), "ERR NotImplementedError: nonlinsolve cannot solve equations with Abs in the complex domain")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + Abs(Symbol("x"))), (Symbol("x") + (Integer(-1)*Symbol("y")))), Symbol("x"), Symbol("y"))), "ERR NotImplementedError: nonlinsolve cannot solve equations with Abs in the complex domain")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + Abs(Symbol("x"))), (Integer(-2) + Symbol("y"))), Symbol("x"), Symbol("y"))), "ERR NotImplementedError: nonlinsolve cannot solve equations with Abs in the complex domain")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-2) + Abs(Symbol("x"))), (Symbol("x") + Symbol("y"))), Symbol("x"), Symbol("y"))), "ERR NotImplementedError: nonlinsolve cannot solve equations with Abs in the complex domain")
})

test_that("solveset: test_nonlinsolve_sign matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sign(Symbol("x"))), (Integer(-4) + (Symbol("x")*Symbol("y")))), list(Symbol("x"), Symbol("y")))), "ERR NotImplementedError: nonlinsolve cannot handle solution of type Interval for symbol x. Got Interval.open(0, oo) from equation sign(x) - 1.")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sign(Symbol("x"))), (Symbol("x") + (Integer(-1)*Symbol("y")))), list(Symbol("x"), Symbol("y")))), "ERR NotImplementedError: nonlinsolve cannot handle solution of type Interval for symbol x. Got Interval.open(0, oo) from equation sign(x) - 1.")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sign(Symbol("x")))), list(Symbol("x")))), "{(ConditionSet(x, Eq(sign(x) - 1, 0), Complexes),)}")
})

test_that("solveset: test_raise_exception_nonlinsolve matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list())), "ERR IndexError: \\nSymbols must be given, for which solution of the system is to be\\nfound.")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(2)))))), "ERR ValueError: \\nSymbols must be given, for which solution of the system is to be\\nfound.")
})

test_that("solveset: test_trig_system matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sin(Symbol("x"))), (Integer(-1) + cos(Symbol("x")))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + sin(Symbol("x"))), cos(Symbol("x"))), Symbol("x"))), "{(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers),)}")
})

test_that("solveset: test_trig_system_fail matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x") + Symbol("y") + (Rational(-1, 2)*S_Pi)), (Integer(-1) + sin(Symbol("x")) + sin(Symbol("y")))), list(Symbol("x"), Symbol("y")))), "{(2*_n*pi - asin(sin(y) - 1), y), (2*_n*pi - asin(sin(y) - 1), y), (2*_n*pi + asin(sin(y) - 1) + pi, y), (2*_n*pi + asin(sin(y) - 1) + pi, y)}")
  expect_identical(.pyrun(nonlinsolve(list((Rational(-1, 2) + (Rational(-1, 2)*(Integer(3))^(Rational(1, 2))) + sin(Symbol("x")) + sin(Symbol("y"))), (Rational(1, 2) + (Integer(-1)*sin(Symbol("y"))) + (Rational(-1, 2)*(Integer(3))^(Rational(1, 2))) + sin(Symbol("x")))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, 2*_n*pi + asin(sin(y) - 1/2 + sqrt(3)/2)), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers)), (ImageSet(Lambda(_n, 2*_n*pi + asin(sin(y) - 1/2 + sqrt(3)/2)), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers)), (ImageSet(Lambda(_n, 2*_n*pi - asin(sin(y) - 1/2 + sqrt(3)/2) + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers)), (ImageSet(Lambda(_n, 2*_n*pi - asin(sin(y) - 1/2 + sqrt(3)/2) + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers))}")
})

test_that("solveset: test_nonlinsolve_positive_dimensional matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x", extended_real = TRUE)*Symbol("y", extended_real = TRUE)), ((Integer(-1)*Symbol("x", extended_real = TRUE)) + (Symbol("x", extended_real = TRUE)*Symbol("y", extended_real = TRUE)))), list(Symbol("x", extended_real = TRUE), Symbol("y", extended_real = TRUE)))), "{(0, y)}")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("a", extended_real = TRUE))^(Integer(2)) + (Symbol("a", extended_real = TRUE)*Symbol("c", extended_real = TRUE))), (Symbol("a", extended_real = TRUE) + (Integer(-1)*Symbol("b", extended_real = TRUE)))), list(Symbol("a", extended_real = TRUE), Symbol("b", extended_real = TRUE)))), "{(0, 0), (-c, -c)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("a", extended_real = TRUE) + Symbol("b", extended_real = TRUE) + Symbol("c", extended_real = TRUE) + Symbol("d", extended_real = TRUE)), ((Symbol("a", extended_real = TRUE)*Symbol("b", extended_real = TRUE)) + (Symbol("a", extended_real = TRUE)*Symbol("d", extended_real = TRUE)) + (Symbol("b", extended_real = TRUE)*Symbol("c", extended_real = TRUE)) + (Symbol("c", extended_real = TRUE)*Symbol("d", extended_real = TRUE))), ((Symbol("a", extended_real = TRUE)*Symbol("b", extended_real = TRUE)*Symbol("c", extended_real = TRUE)) + (Symbol("a", extended_real = TRUE)*Symbol("b", extended_real = TRUE)*Symbol("d", extended_real = TRUE)) + (Symbol("a", extended_real = TRUE)*Symbol("c", extended_real = TRUE)*Symbol("d", extended_real = TRUE)) + (Symbol("b", extended_real = TRUE)*Symbol("c", extended_real = TRUE)*Symbol("d", extended_real = TRUE))), (Integer(-1) + (Symbol("a", extended_real = TRUE)*Symbol("b", extended_real = TRUE)*Symbol("c", extended_real = TRUE)*Symbol("d", extended_real = TRUE)))), list(Symbol("a", extended_real = TRUE), Symbol("b", extended_real = TRUE), Symbol("c", extended_real = TRUE), Symbol("d", extended_real = TRUE)))), "FiniteSet((-1/d, -d, 1/d, Complement({d}, {0})), (1/d, -d, -1/d, Complement({d}, {0})))")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x", extended_real = TRUE))^(Integer(4)) + (Integer(-3)*(Symbol("x", extended_real = TRUE))^(Integer(2))) + (Symbol("x", extended_real = TRUE)*Symbol("y", extended_real = TRUE))), (Symbol("x", extended_real = TRUE)*(Symbol("z"))^(Integer(2))), (Integer(-1) + (Symbol("z")*Symbol("y", extended_real = TRUE)))), list(Symbol("x", extended_real = TRUE), Symbol("y", extended_real = TRUE), Symbol("z")))), "{(0, 1/z, z)}")
})

test_that("solveset: test_nonlinsolve_issue_30522 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Integer(-1)*Symbol("x", real = TRUE)) + (((Integer(-1) + Symbol("x", real = TRUE)))^(Integer(-1))*(Integer(-1) + (Symbol("x", real = TRUE))^(Integer(2))))), (Integer(-1) + Symbol("y", real = TRUE))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "FiniteSet((Complement({x}, {1}), 1))")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-2) + (Integer(-1)*Symbol("x", real = TRUE)) + (((Integer(-2) + Symbol("x", real = TRUE)))^(Integer(-1))*(Integer(-4) + (Symbol("x", real = TRUE))^(Integer(2))))), (Integer(-3) + Symbol("y", real = TRUE))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "FiniteSet((Complement({x}, {2}), 3))")
})

test_that("solveset: test_nonlinsolve_issue_30521 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*Symbol("y", real = TRUE)) + (((Integer(-1) + Symbol("x", real = TRUE)))^(Integer(-1))*((Integer(-1)*Symbol("x", real = TRUE)) + (Symbol("x", real = TRUE)*Symbol("y", real = TRUE)))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "FiniteSet((y, Complement({y}, {1})))")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((((Integer(1) + (Symbol("x"))^(Integer(2))))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y"))))), list(Symbol("x"), Symbol("y")))), "FiniteSet((y, Complement({y}, {-I, I})))")
  expect_identical(.pyrun(nonlinsolve(list((((Integer(-3) + Symbol("y")))^(Integer(-1))*((Integer(-2) + Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y"))))), list(Symbol("x"), Symbol("y")))), "FiniteSet((y, Complement({y}, {2, 3})))")
  expect_identical(.pyrun(nonlinsolve(list((((Integer(-1) + Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y")))*(Symbol("x") + (Integer(-2)*Symbol("y"))))), list(Symbol("x"), Symbol("y")))), "FiniteSet((y, Complement({y}, {1})), (2*y, Complement({y}, {1/2})))")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions_near_zero matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((((Integer(-1) + Symbol("x")))^(Integer(-1))*(Integer(-1) + Symbol("y"))*(Rational("-999999999999", "1000000000000") + Symbol("x"))), (Rational("-999999999999", "1000000000000") + Symbol("x")), (Float("-1.00000000000000", precision = 53) + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "{(0.999999999999, 1.0)}")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions_periodic matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x"))^(Integer(-1))*(Integer(-1) + exp(Symbol("x")))), Symbol("y")), list(Symbol("x"), Symbol("y")))), "ConditionSet((x, y), Eq(y, 0) & Ne(x, 0) & Contains(x, ImageSet(Lambda(_n, 2*_n*I*pi), Integers)), ProductSet(Complexes, Complexes))")
})

test_that("solveset: test_nonlinsolve_coupled_denominator_restrictions matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Symbol("x") + Symbol("y")))^(Integer(-1)))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "ConditionSet((x, y, z), Eq(z, 0) & Ne(x + y, 0), ProductSet(Complexes, Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Integer(-1) + (Symbol("x")*Symbol("y"))))^(Integer(-1)))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "ConditionSet((x, y, z), Eq(z, 0) & Ne(x*y - 1, 0), ProductSet(Complexes, Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Symbol("y") + sin(Symbol("x"))))^(Integer(-1)))), list(sin(Symbol("x")), Symbol("y"), Symbol("z")))), "ConditionSet((_X0, y, z), Eq(z, 0) & Ne(_X0 + y, 0), ProductSet(Complexes, Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("zr", real = TRUE)*((Symbol("xr", real = TRUE) + Symbol("yr", real = TRUE)))^(Integer(-1)))), list(Symbol("xr", real = TRUE), Symbol("yr", real = TRUE), Symbol("zr", real = TRUE)))), "ConditionSet((xr, yr, zr), Eq(zr, 0) & Ne(xr + yr, 0), ProductSet(Reals, Reals, Reals))")
  expect_identical(.pyrun(nonlinsolve(list((((Symbol("x") + Symbol("y")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y")))*(Symbol("x") + (Integer(-1)*Symbol("z"))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "Union(ConditionSet((x, y, z), Eq(x, y) & Ne(0, y), ProductSet(Complexes, Complexes, Complexes)), ConditionSet((x, y, z), Eq(x, z) & Ne(y + z, 0), ProductSet(Complexes, Complexes, Complexes)))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*(Symbol("x"))^(Integer(-1))*((Symbol("x") + Symbol("y")))^(Integer(-1)))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "ConditionSet((x, y, z), Eq(z, 0) & Ne(0, x) & Ne(x + y, 0), ProductSet(Complexes, Complexes, Complexes))")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions_recast matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((((Integer(-1) + Function("f")(Symbol("x"))))^(Integer(-1))*((Integer(-1)*Symbol("y")) + Function("f")(Symbol("x"))))), list(Function("f")(Symbol("x")), Symbol("y")))), "FiniteSet((y, Complement({y}, {1})))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Symbol("y") + Function("f")(Symbol("x"))))^(Integer(-1)))), list(Function("f")(Symbol("x")), Symbol("y"), Symbol("z")))), "ConditionSet((_X0, y, z), Eq(z, 0) & Ne(_X0 + y, 0), ProductSet(Complexes, Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + exp(Function("f")(Symbol("x")))), (Integer(-1) + Symbol("y"))), list(Function("f")(Symbol("x")), Symbol("y")))), "{(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), 1)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Symbol("y") + Function("f")(Symbol("x"))))^(Integer(-1))), (Integer(-1) + exp(Function("f")(Symbol("x"))))), list(Function("f")(Symbol("x")), Symbol("y"), Symbol("z")))), "ConditionSet((_X0, y, z), Eq(z, 0) & Ne(_X0 + y, 0) & Contains(_X0, ImageSet(Lambda(_n, 2*_n*I*pi), Integers)), ProductSet(Complexes, Complexes, Complexes))")
})

test_that("solveset: test_substitution_denominator_restrictions_normalization matches SymPy", {
  expect_identical(.pyrun(substitution(list(Integer(0)), list(Symbol("x"), Symbol("y"), Symbol("z")), result = list(list(list(Symbol("x"), Integer(1)), list(Symbol("y"), Integer(1)), list(Symbol("z"), Integer(1))), list(list(Symbol("x"), Symbol("x")), list(Symbol("y"), Integer(0)), list(Symbol("z"), Symbol("z")))), exclude = list((Symbol("x") + Symbol("z"))))), "Union(ConditionSet((x, y, z), Eq(y, 0) & Ne(x + z, 0), ProductSet(Complexes, Complexes, Complexes)), {(1, 1, 1)})")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions_fallback matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Function("f")(Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y"))))), list(Symbol("x"), Symbol("y")))), "ConditionSet((x, y), Eq(x, y) & Ne(f(y), 0), ProductSet(Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((((Symbol("x") + Symbol("y")))^(Integer(-1))*Function("f")(Symbol("x")))), list(Symbol("x"), Symbol("y")))), "ConditionSet((x, y), Eq(f(x), 0) & Ne(x + y, 0), ProductSet(Complexes, Complexes))")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z")*((Symbol("x") + Symbol("y")))^(Integer(-1))), (Integer(-1) + exp(Symbol("x")))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "ConditionSet((x, y, z), Eq(z, 0) & Ne(x + y, 0) & Contains(x, ImageSet(Lambda(_n, 2*_n*I*pi), Integers)), ProductSet(Complexes, Complexes, Complexes))")
})

test_that("solveset: test_nonlinsolve_denominator_restrictions_polynomial matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*Symbol("y")) + (((Integer(-1) + Symbol("x")))^(Integer(-1))*((Integer(-1)*Symbol("x")) + (Symbol("x")*Symbol("y"))))), (Integer(-1) + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((((Integer(-1) + Symbol("x")))^(Integer(-1))*(Symbol("x") + (Integer(-1)*Symbol("y")))), (Integer(-2) + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "{(2, 2)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Integer(-1)*Symbol("x")) + (((Integer(-1) + Symbol("x")))^(Integer(-1))*(Integer(-1) + (Symbol("x"))^(Integer(2)))))), list(Symbol("x"), Symbol("y")))), "FiniteSet((Complement({x}, {1}), y))")
})

test_that("solveset: test_nonlinsolve_polysys matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-2) + Symbol("y", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(2))), (Symbol("y", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(2)))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-4) + ((Symbol("x", real = TRUE) + Symbol("y", real = TRUE)))^(Integer(2))), (Integer(-2) + Symbol("x", real = TRUE) + Symbol("y", real = TRUE))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(2 - y, y)}")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x", real = TRUE))^(Integer(2)) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(-y, y), (y, y), (-Abs(y), y), (Abs(y), y)}")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x", real = TRUE))^(Integer(2)) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))))), list(Symbol("y", real = TRUE), Symbol("x", real = TRUE)))), "{(y, -y), (y, y), (y, -Abs(y)), (y, Abs(y))}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-3) + Symbol("y", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(2))), (Integer(-4) + Symbol("x", real = TRUE) + (Integer(-1)*Symbol("y", real = TRUE)))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(-1/2 + sqrt(29)/2, -9/2 + sqrt(29)/2), (-sqrt(29)/2 - 1/2, -9/2 - sqrt(29)/2)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-3) + Symbol("y", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(2))), (Integer(-4) + Symbol("x", real = TRUE) + (Integer(-1)*Symbol("y", real = TRUE)))), list(Symbol("y", real = TRUE), Symbol("x", real = TRUE)))), "{(-9/2 - sqrt(29)/2, -sqrt(29)/2 - 1/2), (-9/2 + sqrt(29)/2, -1/2 + sqrt(29)/2)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("z", real = TRUE) + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2))) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2)))), (Integer(-2)*Symbol("x", real = TRUE)), (Integer(-2)*Symbol("y", real = TRUE)), Integer(1)), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE), Symbol("z", real = TRUE)))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x", real = TRUE) + Symbol("y", real = TRUE) + Symbol("z", real = TRUE)), Integer(1), Integer(1), Integer(1)), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE), Symbol("z", real = TRUE)))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("y", real = TRUE))^(Integer(4)) + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2))*(Symbol("z", real = TRUE))^(Integer(2))) + (Symbol("x", real = TRUE)*Symbol("y", real = TRUE)*Symbol("z", real = TRUE))), ((Symbol("y", real = TRUE)*Symbol("z", real = TRUE)) + (Integer(-2)*Symbol("x", real = TRUE)*(Symbol("z", real = TRUE))^(Integer(2)))), ((Integer(4)*(Symbol("y", real = TRUE))^(Integer(3))) + (Symbol("x", real = TRUE)*Symbol("z", real = TRUE))), ((Symbol("x", real = TRUE)*Symbol("y", real = TRUE)) + (Integer(-2)*Symbol("z", real = TRUE)*(Symbol("x", real = TRUE))^(Integer(2))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE), Symbol("z", real = TRUE)))), "{(0, 0, z), (x, 0, 0)}")
})

test_that("solveset: test_nonlinsolve_using_substitution matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(2) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))) + (Symbol("n", real = TRUE)*(Symbol("x", real = TRUE) + Symbol("y", real = TRUE))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(-(n*y - y^2 + 2)/n, y)}")
  expect_identical(.pyrun(nonlinsolve(list((((Symbol("x", real = TRUE))^(Integer(2))*(Symbol("z", real = TRUE))^(Integer(2))) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))*(Symbol("z", real = TRUE))^(Integer(2))*exp((Integer(-1)*Symbol("x", real = TRUE)))))), list(Symbol("y", real = TRUE), Symbol("x", real = TRUE), Symbol("z", real = TRUE)))), "{(y, x, 0), (-x*exp(x/2), x, z), (x*exp(x/2), x, z), (-exp(x/2)*Abs(x), x, z), (exp(x/2)*Abs(x), x, z)}")
})

test_that("solveset: test_nonlinsolve_complex matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*sin(Symbol("y"))) + exp(Symbol("x"))), (Integer(-3) + (Symbol("y"))^(Integer(-1)))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, 2*_n*I*pi + log(sin(1/3))), Integers), 1/3)}")
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*sin(Symbol("y"))) + exp(Symbol("x"))), (Integer(-4) + (Symbol("y"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(sin(2))), Integers), -2), (ImageSet(Lambda(_n, 2*_n*I*pi + log(sin(2))), Integers), 2)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-2) + exp(Symbol("x"))), (Integer(-2) + (Symbol("y"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "{(log(2), -sqrt(2)), (log(2), sqrt(2)), (ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers), -sqrt(2)), (ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers), sqrt(2))}")
})

test_that("solveset: test_nonlinsolve_radical matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("y"))^(Rational(1, 2)) + (Integer(-1)*Symbol("x")) + (Integer(-1)*Symbol("z"))), (Integer(-1) + Symbol("y"))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(1 - z, 1, z)}")
})

test_that("solveset: test_nonlinsolve_inexact matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-9) + ((Symbol("x") + Symbol("y")))^(Integer(2))), (Float("-0.750000000000000", precision = 53) + (Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2))))), list(Symbol("x"), Symbol("y")))), "{(-1.625, -1.375), (1.625, 1.375)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-9) + ((Symbol("x") + Symbol("y")))^(Integer(2))), (Float("-0.750000000000000", precision = 53) + ((Symbol("x") + Symbol("y")))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Float("-0.0625000000000000", precision = 53) + (Symbol("y"))^(Integer(2)) + ((Float("-0.500000000000000", precision = 53) + Symbol("x")))^(Integer(2))), (Float("-1.00000000000000", precision = 53) + (Integer(2)*Symbol("x"))), (Integer(2)*Symbol("y"))), list(Symbol("x"), Symbol("y")))), "EmptySet")
})

test_that("solveset: test_solve_nonlinear_trans matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x", real = TRUE))^(Integer(2)) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))*exp((Integer(-1)*Symbol("x", real = TRUE)))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(x, -x*exp(x/2)), (x, x*exp(x/2)), (x, -exp(x/2)*Abs(x)), (x, exp(x/2)*Abs(x))}")
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x", real = TRUE))^(Integer(2)) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))*exp((Integer(-1)*Symbol("x", real = TRUE)))))), list(Symbol("y", real = TRUE), Symbol("x", real = TRUE)))), "{(-x*exp(x/2), x), (x*exp(x/2), x), (-exp(x/2)*Abs(x), x), (exp(x/2)*Abs(x), x)}")
})

test_that("solveset: test_nonlinsolve_issue_29750 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-3)*Symbol("c")) + (Symbol("s")*(Integer(3))^(Rational(1, 2))) + (Integer(-2)*(Integer(3))^(Rational(1, 2))*(Symbol("c"))^(Integer(2))) + (Integer(2)*Symbol("c")*Symbol("s"))), (Integer(-1) + (Symbol("c"))^(Integer(2)) + (Symbol("s"))^(Integer(2)))), list(Symbol("s"), Symbol("c")))), "{(-1/2, -sqrt(3)/2), (1/2, -sqrt(3)/2), (-sqrt(3)/2, -1/2), (sqrt(3)/2, 1/2)}")
})

test_that("solveset: test_issue_14642 matches SymPy", {
  expect_identical(.pyrun(solveset((Rational(1, 2) + sym_I + (Symbol("x"))^(Integer(2)) + (Rational(1, 2)*(Symbol("x"))^(Integer(3)))), Symbol("x"))), "{-((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/3 - 2/3 - 4*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/(3*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)) + I*(-((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/3 + 4*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/(3*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6))), -2/3 - sqrt(3)*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 - 4*re(1/((-1/2 - sqrt(3)*I/2)*(43/2 + 27*I + sqrt(-256 + (43 + 54*I)^2)/2)^(1/3)))/3 + ((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 + I*(-4*im(1/((-1/2 - sqrt(3)*I/2)*(43/2 + 27*I + sqrt(-256 + (43 + 54*I)^2)/2)^(1/3)))/3 + ((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 + sqrt(3)*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6), -2/3 - 4*re(1/((-1/2 + sqrt(3)*I/2)*(43/2 + 27*I + sqrt(-256 + (43 + 54*I)^2)/2)^(1/3)))/3 + sqrt(3)*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 + ((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 + I*(-sqrt(3)*((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*cos(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 + ((3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2)^2 + (27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)^2)^(1/6)*sin(atan((27 + 3*sqrt(3)*31985^(1/4)*cos(atan(172/49)/2)/2)/(3*sqrt(3)*31985^(1/4)*sin(atan(172/49)/2)/2 + 43/2))/3)/6 - 4*im(1/((-1/2 + sqrt(3)*I/2)*(43/2 + 27*I + sqrt(-256 + (43 + 54*I)^2)/2)^(1/3)))/3)}")
})

test_that("solveset: test_issue_13961 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*Symbol("mx")) + (Symbol("ax")*Symbol("q")) + (Integer(-1)*Symbol("lx")*Symbol("q"))), (Symbol("ax") + (Integer(-1)*Symbol("lx")) + (Integer(-1)*Symbol("gx")*Symbol("q"))), ((Integer(-1)*Symbol("nx")) + (Symbol("bx")*(Symbol("q"))^(Integer(2))) + (Symbol("cx")*Symbol("q")) + (Integer(-1)*Symbol("jx")*Symbol("q"))), (Symbol("q")*(Symbol("mx") + (Symbol("lx")*Symbol("q")) + (Integer(-1)*Symbol("ax")*Symbol("q")))), (Symbol("q")*(Symbol("lx") + (Integer(-1)*Symbol("ax")) + (Symbol("gx")*Symbol("q"))))), Symbol("ax"), Symbol("bx"), Symbol("cx"), Symbol("gx"), Symbol("jx"), Symbol("lx"), Symbol("mx"), Symbol("nx"), Symbol("q"))), "FiniteSet((lx + mx/q, (-cx*q + jx*q + nx)/q^2, cx, mx/q^2, jx, lx, mx, nx, Complement({q}, {0})), (lx + mx/q, -(cx*q - jx*q - nx)/q^2, cx, mx/q^2, jx, lx, mx, nx, Complement({q}, {0})))")
})

test_that("solveset: test_issue_14541 matches SymPy", {
  expect_identical(.pyrun(solveset(((Float("-2.00000000000000", precision = 53) + (Integer(-1)*(Symbol("x"))^(Integer(2)))))^(Rational(1, 2)), Symbol("x"))), "{-1.4142135623731*I, 1.4142135623731*I}")
})

test_that("solveset: test_issue_13396 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-2)*Symbol("y")*Abs(Symbol("x"))*exp(((Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-1)*(Symbol("y"))^(Integer(2)))))), Symbol("y"), domain = Reals)), "{0}")
  expect_identical(.pyrun(solveset(((Rational(-1, 2)*S_Pi) + atan(((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2)))))), Symbol("y"), Reals)), "EmptySet")
})

test_that("solveset: test_issue_10876 matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("x"))^(Rational(-1, 2)), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_issue_19050 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x") + Symbol("y")), sin(Symbol("y"))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, -2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers)), (ImageSet(Lambda(_n, -2*_n*pi - pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("x") + Symbol("y")), (cos(Symbol("y")) + sin(Symbol("y")))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, -2*_n*pi - 3*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/4), Integers)), (ImageSet(Lambda(_n, -2*_n*pi - 7*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + 7*pi/4), Integers))}")
})

test_that("solveset: test_issue_16618 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((sin(Symbol("x"))*sin(Symbol("y"))), (Integer(-1) + (cos(Symbol("x"))*cos(Symbol("y"))))), list(Symbol("x"), Symbol("y")))), "{(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers)), (ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))}")
})

test_that("solveset: test_issue_16643 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("x"))^(Integer(2))*sin(Symbol("x"))), Symbol("x"))), "Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))")
})

test_that("solveset: test_issue_5132_1 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((((Symbol("x"))^(Integer(2)) + (Symbol("y"))^(Integer(2))))^(Rational(1, 2)) + (Integer(-1)*(Integer(10))^(Rational(1, 2)))), (Integer(-4) + Symbol("x") + Symbol("y"))), list(Symbol("x"), Symbol("y")))), "{(1, 3), (3, 1)}")
})

test_that("solveset: test_issue_5132_2 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("z"))^(Integer(2)) + (Integer(-1)*sin(Symbol("y", real = TRUE))) + exp((Integer(2)*Symbol("x", real = TRUE))))), list(Symbol("x", real = TRUE), Symbol("z")))), "{(log(-z^2 + sin(y))/2, z), (ImageSet(Lambda(_n, I*(2*_n*pi + arg(-z^2 + sin(y)))/2 + log(Abs(z^2 - sin(y)))/2), Integers), z)}")
  expect_identical(.pyrun(nonlinsolve(list((Symbol("r") + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2))) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2)))), ((Integer(-1)*Symbol("y", real = TRUE)*(Symbol("x", real = TRUE))^(Integer(-1))) + tan(Symbol("t")))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "Union(ConditionSet((x, y), Eq(x, -sqrt(r/(tan(t)^2 + 1))) & Ne(-sqrt(r/(tan(t)^2 + 1)), 0) & Eq(y, -sqrt(r/(tan(t)^2 + 1))*tan(t)), ProductSet(Reals, Reals)), ConditionSet((x, y), Eq(x, sqrt(r/(tan(t)^2 + 1))) & Ne(sqrt(r/(tan(t)^2 + 1)), 0) & Eq(y, sqrt(r/(tan(t)^2 + 1))*tan(t)), ProductSet(Reals, Reals)))")
})

test_that("solveset: test_issue_6752 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("a", real = TRUE) + (Symbol("a", real = TRUE))^(Integer(2))), (Symbol("a", real = TRUE) + (Integer(-1)*Symbol("b", real = TRUE)))), list(Symbol("a", real = TRUE), Symbol("b", real = TRUE)))), "{(-1, -1), (0, 0)}")
})

test_that("solveset: test_issue_2777 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-10) + (((Symbol("x", real = TRUE))^(Integer(2)) + (Symbol("y", real = TRUE))^(Integer(2))))^(Rational(1, 2))), (Integer(-3) + (((Symbol("y", real = TRUE))^(Integer(2)) + ((Integer(10) + (Integer(-1)*Symbol("x", real = TRUE))))^(Integer(2))))^(Rational(1, 2)))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "{(191/20, -3*sqrt(391)/20), (191/20, 3*sqrt(391)/20)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-10) + (((Symbol("x", real = TRUE))^(Integer(2)) + (Symbol("y", real = TRUE))^(Integer(2))))^(Rational(1, 2))), (((Rational(-191, 20) + Symbol("x", real = TRUE)))^(Integer(-1))*(Integer(-3) + (((Symbol("y", real = TRUE))^(Integer(2)) + ((Integer(10) + (Integer(-1)*Symbol("x", real = TRUE))))^(Integer(2))))^(Rational(1, 2))))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-10) + (((Symbol("x", real = TRUE))^(Integer(2)) + (Symbol("y", real = TRUE))^(Integer(2))))^(Rational(1, 2))), (Integer(3) + (((Symbol("y", real = TRUE))^(Integer(2)) + ((Integer(10) + (Integer(-1)*Symbol("x", real = TRUE))))^(Integer(2))))^(Rational(1, 2)))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "EmptySet")
})

test_that("solveset: test_issue_8828 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x"))^(Integer(2)) + ((Integer(620) + Symbol("y")))^(Integer(2)) + (Integer(-1)*((Integer(920) + (Integer(-1)*Symbol("z"))))^(Integer(2)))), (((Integer(126) + (Integer(-1)*Symbol("x"))))^(Integer(2)) + ((Integer(276) + (Integer(-1)*Symbol("y"))))^(Integer(2)) + (Integer(-1)*(Symbol("z"))^(Integer(2)))), (((Integer(-205) + Symbol("y")))^(Integer(2)) + ((Integer(-51) + Symbol("x")))^(Integer(2)) + (Integer(-1)*((Integer(104) + (Integer(-1)*Symbol("z"))))^(Integer(2))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(376099587/4696637 - 23220*sqrt(2228185)/4696637, 1094845540/4696637 - 46580*sqrt(2228185)/4696637, 48545*sqrt(2228185)/4696637 + 296435228/4696637), (23220*sqrt(2228185)/4696637 + 376099587/4696637, 46580*sqrt(2228185)/4696637 + 1094845540/4696637, 296435228/4696637 - 48545*sqrt(2228185)/4696637)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-920) + Symbol("z") + (((Symbol("x"))^(Integer(2)) + ((Integer(620) + Symbol("y")))^(Integer(2))))^(Rational(1, 2))), (((Integer(126) + (Integer(-1)*Symbol("x"))))^(Integer(2)) + ((Integer(276) + (Integer(-1)*Symbol("y"))))^(Integer(2)) + (Integer(-1)*(Symbol("z"))^(Integer(2)))), (Integer(-104) + Symbol("z") + ((((Integer(-205) + Symbol("y")))^(Integer(2)) + ((Integer(-51) + Symbol("x")))^(Integer(2))))^(Rational(1, 2)))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(376099587/4696637 - 23220*sqrt(2228185)/4696637, 1094845540/4696637 - 46580*sqrt(2228185)/4696637, 48545*sqrt(2228185)/4696637 + 296435228/4696637), (23220*sqrt(2228185)/4696637 + 376099587/4696637, 46580*sqrt(2228185)/4696637 + 1094845540/4696637, 296435228/4696637 - 48545*sqrt(2228185)/4696637)}")
})

test_that("solveset: test_nonlinsolve_conditionset matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Rational(-1, 2)*S_Pi) + Function("f")(Symbol("x"))), ((Rational(-3, 2)*S_Pi) + Function("f")(Symbol("y")))), list(Symbol("x"), Symbol("y")))), "ConditionSet((x, y), Eq(2*f(x) - pi, 0) & Eq(2*f(y) - 3*pi, 0), ProductSet(Complexes, Complexes))")
})

test_that("solveset: test_substitution_basic matches SymPy", {
  expect_identical(.pyrun(substitution(list(), list(Symbol("x"), Symbol("y")))), "EmptySet")
  expect_identical(.pyrun(substitution(list(), list())), "EmptySet")
  expect_identical(.pyrun(substitution(list((Integer(-30) + (Integer(2)*(Symbol("x"))^(Integer(2))) + (Integer(3)*(Symbol("y"))^(Integer(2)))), (Integer(-19) + (Integer(-2)*(Symbol("y"))^(Integer(2))) + (Integer(3)*(Symbol("x"))^(Integer(2))))), list(Symbol("x"), Symbol("y")))), "{(-3, -2), (-3, 2), (3, -2), (3, 2)}")
  expect_identical(.pyrun(substitution(list((Symbol("x") + Symbol("y"))), list(Symbol("x")), list(list(list(Symbol("y"), Integer(1)))), list(Symbol("y")), list(), list(Symbol("x"), Symbol("y")))), "{(-1, 1)}")
  expect_identical(.pyrun(substitution(list((Symbol("x") + Symbol("y"))), list(Symbol("x")), list(list(list(Symbol("y"), Integer(1)))), list(Symbol("y")), list((Integer(1) + Symbol("x"))), list(Symbol("y"), Symbol("x")))), "EmptySet")
})

test_that("solveset: test_substitution_incorrect matches SymPy", {
  expect_identical(.pyrun(substitution(list((Integer(-1) + Symbol("h")), (Integer(-1) + Symbol("k")), (Integer(-2) + Symbol("f")), (Integer(-4) + Symbol("f")), (Integer(-2)*Symbol("k"))), list(Symbol("h"), Symbol("k"), Symbol("f")))), "{(1, 1, f)}")
  expect_identical(.pyrun(substitution(list((Symbol("x") + Symbol("y") + Symbol("z")), Integer(1), Integer(1), Integer(1)), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(-y - z, y, z)}")
  expect_identical(.pyrun(substitution(list((Symbol("a") + (Integer(-1)*Symbol("d"))), (Symbol("b") + Symbol("d")), (Symbol("c") + Symbol("d")), (Integer(1) + (Symbol("d"))^(Integer(2)))), list(Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d")))), "{(d, -d, -d, d)}")
  expect_identical(.pyrun(substitution(list((Symbol("a")*(Symbol("a") + (Integer(-1)*log(Symbol("b"))))), (Symbol("a")*(Integer(-2) + Symbol("b")))), list(Symbol("a"), Symbol("b")))), "{(0, b)}")
})

test_that("solveset: test_substitution_redundant matches SymPy", {
  expect_identical(.pyrun(substitution(list(((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2)))), (Integer(-1) + Symbol("z"))), list(Symbol("x"), Symbol("z")))), "{(-y, 1), (y, 1), (-sqrt(y^2), 1), (sqrt(y^2), 1)}")
})

test_that("solveset: test_issue_5132_substitution matches SymPy", {
  expect_identical(.pyrun(substitution(list((Symbol("r", real = TRUE) + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2))) + (Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2)))), ((Integer(-1)*Symbol("y", real = TRUE)*(Symbol("x", real = TRUE))^(Integer(-1))) + tan(Symbol("t", real = TRUE)))), list(Symbol("x", real = TRUE), Symbol("y", real = TRUE)))), "FiniteSet((Complement({-sqrt(r)/sqrt(tan(t)^2 + 1)}, {0}), -sqrt(r)*tan(t)/sqrt(tan(t)^2 + 1)), (Complement({sqrt(r)/sqrt(tan(t)^2 + 1)}, {0}), sqrt(r)*tan(t)/sqrt(tan(t)^2 + 1)))")
  expect_identical(.pyrun(substitution(list(((Symbol("z", real = TRUE))^(Integer(2)) + (Integer(-1)*sin(Symbol("y", real = TRUE))) + exp((Integer(2)*Symbol("x", real = TRUE)))), (Integer(-3) + exp((Integer(-1)*Symbol("y", real = TRUE))))), list(Symbol("y", real = TRUE), Symbol("z", real = TRUE)))), "{(-log(3), -sqrt(-exp(2*x) - sin(log(3)))), (-log(3), sqrt(-exp(2*x) - sin(log(3)))), (ImageSet(Lambda(_n, 2*_n*I*pi - log(3)), Integers), ImageSet(Lambda(_n, -sqrt(-exp(2*x) + sin(2*_n*I*pi - log(3)))), Integers)), (ImageSet(Lambda(_n, 2*_n*I*pi - log(3)), Integers), ImageSet(Lambda(_n, sqrt(-exp(2*x) + sin(2*_n*I*pi - log(3)))), Integers))}")
})

test_that("solveset: test_raises_substitution matches SymPy", {
  expect_identical(.pyrun(substitution(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list())), "ERR ValueError: \\nSymbols must be given, for which solution of the system is to be\\nfound.")
  expect_identical(.pyrun(substitution(list((Integer(-1) + (Symbol("x"))^(Integer(2)))))), "ERR TypeError: substitution() missing 1 required positional argument: 'symbols'")
  expect_identical(.pyrun(substitution(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), list(sin(Symbol("x"))))), "ERR ValueError: \\nIterable of symbols must be given as second argument, not type sin:\\nsin(x)")
  expect_identical(.pyrun(substitution(list((Integer(-1) + (Symbol("x"))^(Integer(2)))), Symbol("x"))), "ERR TypeError: \\nsymbols should be given as a sequence, e.g. a list.Not type <class\\n'sympy.core.symbol.Symbol'>: x")
})

test_that("solveset: test_issue_21022 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-16) + Symbol("k")), (Integer(-8) + Symbol("p")), ((Symbol("y"))^(Integer(2)) + (Symbol("z"))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2)))), (Symbol("d") + Symbol("p") + (Integer(-1)*Symbol("x"))), ((Symbol("d"))^(Integer(2)) + (Symbol("k"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2)))), ((Symbol("z"))^(Integer(2)) + (Integer(-1)*(Symbol("k"))^(Integer(2))) + (Integer(-1)*(Symbol("p"))^(Integer(2)))), (Symbol("abc") + (Integer(-1)*Symbol("efg")))), list(Symbol("abc"), Symbol("d"), Symbol("efg"), Symbol("k"), Symbol("p"), Symbol("x"), Symbol("y"), Symbol("z")))), "{(efg, 32, efg, 16, 8, 40, -16*sqrt(5), -8*sqrt(5)), (efg, 32, efg, 16, 8, 40, -16*sqrt(5), 8*sqrt(5)), (efg, 32, efg, 16, 8, 40, 16*sqrt(5), -8*sqrt(5)), (efg, 32, efg, 16, 8, 40, 16*sqrt(5), 8*sqrt(5))}")
})

test_that("solveset: test_issue_17906 matches SymPy", {
  expect_identical(.pyrun(solveset(((Integer(7))^((Integer(-80) + (Symbol("x"))^(Integer(2)))) + (Integer(-1)*(Integer(49))^(Symbol("x")))), Symbol("x"))), "{-8, 10}")
})

test_that("solveset: test_issue_17933 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x")*sin(Integer(45))) + (Integer(-1)*Symbol("y")*cos(Symbol("q")))), ((Symbol("x")*cos(Integer(45))) + (Integer(-1)*Symbol("y")*sin(Symbol("q")))), ((Symbol("y")*cos(Symbol("q"))) + (Rational(9, 10)*Symbol("x")*sin(Integer(45)))), ((Integer(-1)*Symbol("z")) + (Symbol("y")*sin(Symbol("z"))) + (Rational(9, 10)*Symbol("x")*cos(Integer(45))))), Symbol("x"), Symbol("y"), Symbol("z"), Symbol("q"))), "EmptySet")
})

test_that("solveset: test_issue_17933_bis matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Symbol("x")*sin(Integer(45))) + (Integer(-1)*Symbol("y")*cos(Symbol("zz")))), ((Symbol("x")*cos(Integer(45))) + (Integer(-1)*Symbol("y")*sin(Symbol("zz")))), ((Symbol("y")*cos(Symbol("zz"))) + (Rational(9, 10)*Symbol("x")*sin(Integer(45)))), ((Integer(-1)*Symbol("z")) + (Symbol("y")*sin(Symbol("z"))) + (Rational(9, 10)*Symbol("x")*cos(Integer(45))))), Symbol("x"), Symbol("y"), Symbol("z"), Symbol("zz"))), "{(0, 0, 0, zz)}")
})

test_that("solveset: test_issue_9556 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + Abs(Symbol("x"))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Symbol("b", positive = TRUE) + Abs(Symbol("x"))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(S_false, Symbol("b", positive = TRUE), Reals)), "EmptySet")
})

test_that("solveset: test_issue_9611 matches SymPy", {
  expect_identical(.pyrun(solveset(S_true, Symbol("x"), Reals)), "Reals")
  expect_identical(.pyrun(solveset(S_true, Symbol("y"))), "Complexes")
})

test_that("solveset: test_issue_9557 matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("a") + (Symbol("x"))^(Integer(2))), Symbol("x"), Reals)), "Intersection({-sqrt(-a), sqrt(-a)}, Reals)")
})

test_that("solveset: test_issue_9778 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + (Symbol("x", real = TRUE))^(Integer(3))), Symbol("x", real = TRUE), Reals)), "{-1}")
  expect_identical(.pyrun(solveset((Integer(1) + (Symbol("x", real = TRUE))^(Rational(3, 5))), Symbol("x", real = TRUE), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Symbol("y", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(3))), Symbol("x", real = TRUE), Reals)), "{-Abs(y)^(1/3)*sign(y)}")
})

test_that("solveset: test_issue_10214 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(4) + (Symbol("x"))^(Rational(3, 2))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(4) + (Symbol("x"))^(Rational(-3, 2))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(4) + (Symbol("x"))^(Integer(3))), Symbol("x"), Reals)), "{-2^(2/3)}")
})

test_that("solveset: test_issue_9849 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + Abs(sin(Symbol("x")))), Symbol("x"), Reals)), "EmptySet")
})

test_that("solveset: test_issue_9953 matches SymPy", {
  expect_identical(.pyrun(linsolve(list(), Symbol("x"))), "EmptySet")
})

test_that("solveset: test_issue_9913 matches SymPy", {
  expect_identical(.pyrun(solveset((((Integer(-10) + Symbol("x")))^(Integer(-2)) + (Integer(2)*Symbol("x"))), Symbol("x"), Reals)), "{-(3*sqrt(24081)/4 + 4027/4)^(1/3)/3 - 100/(3*(3*sqrt(24081)/4 + 4027/4)^(1/3)) + 20/3}")
})

test_that("solveset: test_issue_10397 matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("x"))^(Rational(1, 2)), Symbol("x"), Complexes)), "{0}")
})

test_that("solveset: test_issue_14987 matches SymPy", {
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("x"))^(Integer(2))), Symbol("x"))), "ERR NonlinearError: nonlinear term: x^2")
  expect_identical(.pyrun(linear_eq_to_matrix(list(((Integer(-1)*Symbol("a")) + (Integer(2)*Symbol("y")) + (Symbol("x")*(Integer(1) + (Integer(-3)*(Symbol("x"))^(Integer(-1))))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: 1/x")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Integer(-3) + (((Integer(-3) + Symbol("x")))^(Integer(-1))*((Symbol("x"))^(Integer(2)) + (Integer(-3)*Symbol("x")))))), Symbol("x"))), "ERR NonlinearError: nonlinear term: 1/(x - 3)")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Integer(7) + ((Integer(1) + Symbol("x")))^(Integer(3)) + (Integer(-1)*(Symbol("x"))^(Integer(3))) + (Integer(-3)*(Symbol("x"))^(Integer(2))))), Symbol("x"))), "ERR NonlinearError: nonlinear term: (x + 1)^3")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("y") + (Symbol("x")*(Integer(1) + (Symbol("x"))^(Integer(-1)))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: 1/x")
  expect_identical(.pyrun(linear_eq_to_matrix(list((Symbol("y")*(Integer(1) + Symbol("x")))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: \\nnonlinear cross-term: y*(x + 1)")
  expect_identical(.pyrun(linear_eq_to_matrix(list(Eq((Symbol("x"))^(Integer(-1)), (Symbol("y") + (Symbol("x"))^(Integer(-1))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: 1/x")
  expect_identical(.pyrun(linear_eq_to_matrix(list(Eq((Symbol("y")*(Symbol("x"))^(Integer(-1))), (Symbol("y") + (Symbol("y")*(Symbol("x"))^(Integer(-1)))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: nonlinear term: 1/x")
  expect_identical(.pyrun(linear_eq_to_matrix(list(Eq((Symbol("x")*(Integer(1) + Symbol("x"))), (Symbol("y") + (Symbol("x"))^(Integer(2))))), list(Symbol("x"), Symbol("y")))), "ERR NonlinearError: \\nnonlinear cross-term: x*(x + 1)")
})

test_that("solveset: test_simplification matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("x") + ((((Integer(-2)*Symbol("a")) + (Integer(2)*Symbol("b"))))^(Integer(-1))*(Symbol("a") + (Integer(-1)*Symbol("b"))))), Symbol("x"))), "{1/2}")
  expect_identical(.pyrun(solveset((Symbol("x") + ((((Integer(-2)*Symbol("a")) + (Integer(2)*Symbol("b"))))^(Integer(-1))*(Symbol("a") + (Integer(-1)*Symbol("b"))))), Symbol("x"), Reals)), "Intersection({-(a - b)/(-2*a + 2*b)}, Reals)")
  expect_identical(.pyrun(solveset((Symbol("x") + ((((Integer(-2)*Symbol("ap", positive = TRUE)) + (Integer(2)*Symbol("bn", negative = TRUE))))^(Integer(-1))*(Symbol("ap", positive = TRUE) + (Integer(-1)*Symbol("bn", negative = TRUE))))), Symbol("x"))), "{1/2}")
  expect_identical(.pyrun(solveset((Symbol("x") + ((((Integer(-2)*Symbol("ap", positive = TRUE)) + (Integer(2)*Symbol("bn", negative = TRUE))))^(Integer(-1))*(Symbol("ap", positive = TRUE) + (Integer(-1)*Symbol("bn", negative = TRUE))))), Symbol("x"), Reals)), "{1/2}")
})

test_that("solveset: test_integer_domain_relational matches SymPy", {
  expect_identical(.pyrun(solveset(Gt((Integer(3) + (Integer(2)*Symbol("x"))), Integer(0)), Symbol("x"), Integers)), "Range(-1, oo, 1)")
  expect_identical(.pyrun(solveset(Ge((Integer(-2) + (Symbol("x"))^(Integer(2)) + (Integer(3)*Symbol("x"))), Integer(0)), Symbol("x"), Integers)), "Union(Range(-oo, -3, 1), Range(1, oo, 1))")
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + (Symbol("x"))^(Integer(-1))), (Integer(-2) + (Symbol("x"))^(Integer(-1)))), Symbol("x"), Integers)), "Union(Range(-1, 0, 1), Range(1, oo, 1))")
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + ((Integer(-5) + (Symbol("x"))^(Integer(2))))^(Rational(1, 2))), Integer(0)), Symbol("x"), Integers)), "Range(3, oo, 1)")
  expect_identical(.pyrun(solveset(Gt(((log(Symbol("x")))^(Integer(-1)) + log(Symbol("x"))), (Integer(-2) + (log(Symbol("x")))^(Integer(-1)))), Symbol("x"), Integers)), "Range(2, oo, 1)")
  expect_identical(.pyrun(solveset(Le(((Symbol("x"))^(Integer(-1))*log(Symbol("x"))), Integer(0)), Symbol("x"), Integers)), "Range(1, 2, 1)")
  expect_identical(.pyrun(solveset(Lt(((Symbol("x"))^(Integer(-1))*log(Symbol("x"))), Integer(0)), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset(Lt((Symbol("x")*((Integer(-3) + Symbol("x")))^(Integer(-1))), Integer(3)), Symbol("x"), domain = Range(Integer(0), Integer(5), Integer(1)))), "Range(0, 3, 1)")
  expect_identical(.pyrun(solveset(Lt((Symbol("x")*((Integer(-3) + (Symbol("x"))^(Integer(2))))^(Integer(-1))), Integer(3)), Symbol("x"), domain = Range(Integer(0), Integer(5), Integer(1)))), "Union(Range(0, 2, 1), Range(2, 5, 1))")
  expect_identical(.pyrun(solveset(Lt((Integer(2) + Symbol("x")), Integer(0)), Symbol("x"), Integers)), "Range(-oo, -2, 1)")
})

test_that("solveset: test_issue_10555 matches SymPy", {
  expect_identical(.pyrun(solveset(((Rational(-1, 2)*S_Pi) + Function("f")(Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(f(x) - pi/2, 0), Reals)")
  expect_identical(.pyrun(solveset(((Rational(-1, 2)*S_Pi) + Function("f")(Function("g")(Symbol("x")))), Function("g")(Symbol("x")), Reals)), "ConditionSet(g(x), Eq(f(g(x)) - pi/2, 0), Reals)")
})

test_that("solveset: test_issue_8715 matches SymPy", {
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + (Symbol("x"))^(Integer(-1))), (Integer(-2) + (Symbol("x"))^(Integer(-1)))), Symbol("x"), Reals)), "Union(Interval.open(-2, 0), Interval.open(0, oo))")
  expect_identical(.pyrun(solveset(Gt(((log(Symbol("x")))^(Integer(-1)) + log(Symbol("x"))), (Integer(-2) + (log(Symbol("x")))^(Integer(-1)))), Symbol("x"), Reals)), "Union(Interval.open(1, oo), Interval.open(exp(-2), 1))")
})

test_that("solveset: test_issue_11174 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("z"))^(Integer(2)) + (Integer(-1)*sin(Symbol("y"))) + exp((Integer(2)*Symbol("x")))), Symbol("x"), Reals)), "Intersection({log(-z^2 + sin(y))/2}, Reals)")
  expect_identical(.pyrun(solveset(((Symbol("x")*tan(Symbol("t"))) + ((Symbol("r"))^(Rational(1, 2))*((Integer(1) + (tan(Symbol("t")))^(Integer(2))))^(Rational(-1, 2))*Abs(tan(Symbol("t"))))), Symbol("x"), Reals)), "Intersection({-sqrt(r)*Abs(tan(t))/(sqrt(tan(t)^2 + 1)*tan(t))}, Reals)")
})

test_that("solveset: test_issue_11534 matches SymPy", {
  expect_identical(.pyrun(solveset(((Integer(-1)*Symbol("y", real = TRUE)) + (Symbol("x", real = TRUE)*((Integer(1) + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2)))))^(Rational(-1, 2)))), Symbol("x", real = TRUE), Reals)), "Union(Complement(ConditionSet(x, Eq(x/sqrt(1 - x^2) - y, 0), {-y/sqrt(y^2 + 1)}), {-1, 1}), Complement({y/sqrt(y^2 + 1)}, {-1, 1}))")
  expect_identical(.pyrun(solveset(((Integer(-1)*(Symbol("y", real = TRUE))^(Integer(2))) + ((Symbol("x", real = TRUE))^(Integer(2))*((Integer(1) + (Integer(-1)*(Symbol("x", real = TRUE))^(Integer(2)))))^(Integer(-1)))), Symbol("x", real = TRUE), Reals)), "Complement({-y/sqrt(y^2 + 1), y/sqrt(y^2 + 1)}, {-1, 1})")
})

test_that("solveset: test_issue_10477 matches SymPy", {
  expect_identical(.pyrun(solveset(Lt(((Symbol("x"))^(Integer(-1))*(Integer(-3) + (Symbol("x"))^(Integer(2)) + (Integer(4)*Symbol("x")))), Integer(2)), Symbol("x"), Reals)), "Union(Interval.open(-oo, -3), Interval.open(0, 1))")
})

test_that("solveset: test_issue_10671 matches SymPy", {
  expect_identical(.pyrun(solveset(sin(Symbol("y")), Symbol("y"), Interval(Integer(0), S_Pi, FALSE, FALSE))), "{0, pi}")
  expect_identical(.pyrun(solveset(Lt((Integer(-1)*(Symbol("x"))^(Integer(-2))), Integer(0)), Symbol("x"), Interval(Integer(1), Integer(10), FALSE, FALSE))), "Interval(1, 10)")
})

test_that("solveset: test_issue_11064 matches SymPy", {
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + ((Integer(-5) + (Symbol("x"))^(Integer(2))))^(Rational(1, 2))), Integer(0)), Symbol("x"), Reals)), "Interval(sqrt(5), oo)")
  expect_identical(.pyrun(solveset(Lt((Symbol("x") + ((Integer(-5) + (Symbol("x"))^(Integer(2))))^(Rational(1, 2))), Integer(0)), Symbol("x"), Reals)), "Interval(-oo, -sqrt(5))")
  expect_identical(.pyrun(solveset(Gt((Symbol("x") + ((Integer(-5) + (Symbol("x"))^(Integer(2))))^(Rational(1, 2))), (Integer(5))^(Rational(1, 2))), Symbol("x"), Reals)), "Interval.open(sqrt(5), oo)")
})

test_that("solveset: test_sqrt_inequality matches SymPy", {
  expect_identical(.pyrun(solveset(Gt((((Integer(2) + (Integer(4)*Symbol("u"))))^(Integer(-1))*(Integer(1) + Symbol("u"))*(Integer(2) + (Integer(-1)*((Symbol("u")*(Integer(8) + (Integer(17)*Symbol("u")))))^(Rational(1, 2))) + (Integer(5)*Symbol("u")))), Integer(0)), Symbol("u"), Reals)), "Union(Interval.open(-1, -1/2), Interval(0, oo))")
  expect_identical(.pyrun(solveset(Gt((((Integer(2) + (Integer(4)*Symbol("u"))))^(Integer(-1))*(Integer(1) + Symbol("u"))*(Integer(2) + (Integer(-1)*((Symbol("u")*(Integer(8) + (Integer(17)*Symbol("u")))))^(Rational(1, 2))) + (Integer(5)*Symbol("u")))), Integer(1)), Symbol("u"), Reals)), "Union(Interval.open(1 - sqrt(3), -1/2), Interval.open(1 + sqrt(3), oo))")
  expect_identical(.pyrun(solveset(Gt((((Integer(2) + (Integer(4)*Symbol("u"))))^(Integer(-1))*(Integer(1) + Symbol("u"))*(Integer(2) + ((Symbol("u")*(Integer(8) + (Integer(17)*Symbol("u")))))^(Rational(1, 2)) + (Integer(5)*Symbol("u")))), Integer(1)), Symbol("u"), Reals)), "Interval.open(0, oo)")
})

test_that("solveset: test_issue_12478 matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(2) + ((Integer(-2) + Symbol("x")))^(Rational(1, 2))), Symbol("x"))), "EmptySet")
  expect_identical(.pyrun(solveset(Lt((Integer(2) + ((Integer(-2) + Symbol("x")))^(Rational(1, 2))), Integer(0)), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Gt((Integer(2) + ((Integer(-2) + Symbol("x")))^(Rational(1, 2))), Integer(0)), Symbol("x"), Reals)), "Interval(2, oo)")
})

test_that("solveset: test_issue_12429 matches SymPy", {
  expect_identical(.pyrun(solveset(Le(((Symbol("x"))^(Integer(-1))*log(Symbol("x"))), Integer(0)), Symbol("x"), Reals)), "Interval.Lopen(0, 1)")
})

test_that("solveset: test_issue_19506 matches SymPy", {
  expect_identical(.pyrun(solveset(arg((sym_I + Symbol("x"))))), "Intersection(ConditionSet(x, Eq(im(x) + 1, 0), Complexes), ConditionSet(x, re(x) > 0, Complexes))")
})

test_that("solveset: test_solveset_arg matches SymPy", {
  expect_identical(.pyrun(solveset(arg(Symbol("x")), Symbol("x"), Reals)), "Interval.open(0, oo)")
  expect_identical(.pyrun(solveset(arg((Integer(-3) + (Integer(4)*Symbol("x")))), Symbol("x"), Reals)), "Interval.open(3/4, oo)")
})

test_that("solveset: test_issue_13550 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-15) + (Symbol("x"))^(Integer(2)) + (Integer(-2)*Symbol("x"))), symbol = Symbol("x"), domain = Interval((-oo), Integer(0), TRUE, FALSE))), "{-3}")
})

test_that("solveset: test_issue_13849 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-1)*(Integer(2))^(Rational(1, 2))) + (Symbol("t")*((Integer(2))^(Rational(1, 2)) + (Integer(5))^(Rational(1, 2))))), Symbol("t")), Symbol("t"))), "EmptySet")
})

test_that("solveset: test_issue_14223 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-2) + Abs((Symbol("x") + Piecewise(list(Integer(2), Ge(Symbol("x"), Integer(2))), list(Symbol("x"), S_true))))), Symbol("x"), Reals)), "{-1, 1}")
  expect_identical(.pyrun(solveset((Integer(-2) + Abs((Symbol("x") + Piecewise(list(Integer(2), Ge(Symbol("x"), Integer(2))), list(Symbol("x"), S_true))))), Symbol("x"), Interval(Integer(0), Integer(2), FALSE, FALSE))), "{1}")
  expect_identical(.pyrun(solveset(Symbol("x"), Symbol("x"), FiniteSet(Integer(1), Integer(2)))), "EmptySet")
})

test_that("solveset: test_issue_10158 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-10) + (Symbol("x")*Max(Integer(15), Symbol("x")))), Symbol("x"), Reals)), "{2/3}")
  expect_identical(.pyrun(solveset((Integer(-10) + (Symbol("x")*Min(Integer(15), Symbol("x")))), Symbol("x"), Reals)), "{-sqrt(10), sqrt(10)}")
  expect_identical(.pyrun(solveset((Integer(-3) + Max((Integer(2) + Symbol("x")), (Integer(-1) + Abs((Integer(-3) + Symbol("x")))))), Symbol("x"), Reals)), "{-1, 1}")
  expect_identical(.pyrun(solveset(((Integer(-1)*Abs(Symbol("y"))) + Abs((Integer(-1) + Symbol("x")))), Symbol("x"), Reals)), "{1 - Abs(y), Abs(y) + 1}")
  expect_identical(.pyrun(solveset(Abs((Symbol("x") + (Integer(4)*Abs((Integer(1) + Symbol("x")))))), Symbol("x"), Reals)), "{-4/3, -4/5}")
  expect_identical(.pyrun(solveset((Integer(-2) + (Integer(2)*Abs((Symbol("x") + Abs((Symbol("x") + Max(Integer(3), Symbol("x")))))))), Symbol("x"), Reals)), "{-2, -1}")
  expect_identical(.pyrun(solveset((Integer(-10) + (Symbol("x")*Max(Integer(15), Symbol("x")))), Symbol("x"), Complexes)), "ERR ValueError: \\nInequalities in the complex domain are not supported. Try the real\\ndomain by setting domain=S.Reals")
  expect_identical(.pyrun(solveset((Integer(-10) + (Symbol("x")*Min(Integer(15), Symbol("x")))), Symbol("x"), Complexes)), "ERR ValueError: \\nInequalities in the complex domain are not supported. Try the real\\ndomain by setting domain=S.Reals")
  expect_identical(.pyrun(solveset((Integer(-3) + Max((Integer(2) + Symbol("x")), (Integer(-1) + Abs((Integer(-3) + Symbol("x")))))), Symbol("x"), Complexes)), "ERR ValueError: \\nInequalities in the complex domain are not supported. Try the real\\ndomain by setting domain=S.Reals")
  expect_identical(.pyrun(solveset(((Integer(-1)*Abs(Symbol("y"))) + Abs((Integer(-1) + Symbol("x")))), Symbol("x"), Complexes)), "ERR ValueError: \\nAbsolute values cannot be inverted in the complex domain.")
  expect_identical(.pyrun(solveset(Abs((Symbol("x") + (Integer(4)*Abs((Integer(1) + Symbol("x")))))), Symbol("x"), Complexes)), "ERR ValueError: \\nAbsolute values cannot be inverted in the complex domain.")
})

test_that("solveset: test_issue_14300 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(1) + (Integer(-1)*Symbol("y")) + (Integer(-1)*exp((Integer(-18000000)*Symbol("x"))))), Symbol("x"), Reals)), "Intersection({-log(1 - y)/18000000}, Reals)")
  expect_identical(.pyrun(solveset((Integer(1) + (Integer(-1)*Symbol("y")) + (Integer(-1)*exp((Integer(-18000000)*Symbol("x"))))), Symbol("x"))), "ImageSet(Lambda(_n, -I*(2*_n*pi + arg(1 - y))/18000000 - log(Abs(y - 1))/18000000), Integers)")
})

test_that("solveset: test_issue_17882 matches SymPy", {
  expect_identical(.pyrun(solveset(((Rational(4, 3)*((Integer(-1) + (Symbol("x"))^(Integer(2))))^(Rational(-1, 3))) + (Rational(-8, 9)*(Symbol("x"))^(Integer(2))*((Integer(-1) + (Symbol("x"))^(Integer(2))))^(Rational(-4, 3)))), Symbol("x"), Complexes)), "{-sqrt(3), sqrt(3)}")
})

test_that("solveset: test_term_factors matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_term_factors((Integer(-2) + (Integer(3))^(Symbol("x"))))), "[-2, 3^x]")
  expect_identical(.pyrun(symcompR:::.sset_term_factors(((Integer(4))^((Integer(1) + Symbol("x"))) + (Integer(4))^((Integer(-1) + Symbol("x"))) + (Integer(4))^((Integer(2) + Symbol("x"))) + (Integer(-1)*(Integer(3))^((Integer(2) + Symbol("x")))) + (Integer(-1)*(Integer(3))^((Integer(3) + Symbol("x"))))))), "[4^(x + 1), 4^(x - 1), 4^(x + 2), -1, 3^(x + 2), -1, 3^(x + 3)]")
})

test_that("solveset: test_transolve matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_transolve((Integer(3))^(Symbol("x")), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(symcompR:::.sset_transolve(((Integer(3))^(Symbol("x")) + (Integer(-1)*(Integer(9))^((Integer(5) + Symbol("x"))))), Symbol("x"), Reals)), "{-10}")
})

test_that("solveset: test_issue_21276 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("y"))^(Integer(2)) + (Symbol("z"))^(Integer(2)) + ((Symbol("y"))^(Integer(2))*(erf((Symbol("y") + (Integer(-1)*Symbol("z")))))^(Integer(2))) + ((Symbol("z"))^(Integer(2))*(erf((Symbol("y") + (Integer(-1)*Symbol("z")))))^(Integer(2))) + (Integer(-4)*Symbol("x")*(Symbol("y"))^(Integer(2))) + (Integer(-4)*Symbol("x")*(Symbol("z"))^(Integer(2))) + (Integer(-2)*Symbol("y")*Symbol("z")) + (Integer(2)*(Symbol("y"))^(Integer(2))*erf((Symbol("y") + (Integer(-1)*Symbol("z"))))) + (Integer(2)*(Symbol("z"))^(Integer(2))*erf((Symbol("y") + (Integer(-1)*Symbol("z"))))) + (Integer(4)*(Symbol("x"))^(Integer(2))*(Symbol("y"))^(Integer(2))) + (Integer(4)*(Symbol("x"))^(Integer(2))*(Symbol("z"))^(Integer(2))) + (Integer(-8)*Symbol("y")*Symbol("z")*(Symbol("x"))^(Integer(2))) + (Integer(-4)*Symbol("x")*(Symbol("y"))^(Integer(2))*erf((Symbol("y") + (Integer(-1)*Symbol("z"))))) + (Integer(-4)*Symbol("x")*(Symbol("z"))^(Integer(2))*erf((Symbol("y") + (Integer(-1)*Symbol("z"))))) + (Integer(-4)*Symbol("y")*Symbol("z")*erf((Symbol("y") + (Integer(-1)*Symbol("z"))))) + (Integer(-2)*Symbol("y")*Symbol("z")*(erf((Symbol("y") + (Integer(-1)*Symbol("z")))))^(Integer(2))) + (Integer(8)*Symbol("x")*Symbol("y")*Symbol("z")) + (Integer(8)*Symbol("x")*Symbol("y")*Symbol("z")*erf((Symbol("y") + (Integer(-1)*Symbol("z")))))), Symbol("y"))), "{z, z + erfinv(2*x - 1)}")
})

test_that("solveset: test_exponential_real matches SymPy", {
  expect_identical(.pyrun(solveset(((Integer(3))^((Integer(2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Integer(3) + Symbol("x"))))), Symbol("x"), Reals)), "{-3*log(2)/(-2*log(3) + log(2))}")
  expect_identical(.pyrun(solveset(((Integer(4))^((Integer(5) + (Integer(-9)*Symbol("x")))) + (Integer(-1)*(Integer(8))^((Integer(2) + (Integer(-1)*Symbol("x")))))), Symbol("x"), Reals)), "{4/15}")
  expect_identical(.pyrun(solveset(((Integer(2))^(Symbol("x")) + (Integer(4))^(Symbol("x"))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(((Integer(-1)*(Integer(2))^(Symbol("x"))) + exp((Symbol("x")*log(Integer(5))))), Symbol("x"), Reals)), "{0}")
  expect_identical(.pyrun(solveset((Integer(-2) + (exp((Symbol("x")*(Symbol("y"))^(Integer(-1))))*exp((Integer(-1)*Symbol("z")*(Symbol("y"))^(Integer(-1)))))), Symbol("x"), Reals)), "Intersection({y*log(2*exp(z/y))}, Reals)")
  expect_identical(.pyrun(solveset(((Integer(5))^((Rational(1, 2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Rational(1, 3)*Symbol("x"))))), Symbol("x"), Reals)), "{0}")
  expect_identical(.pyrun(solveset(((Integer(4))^((Integer(1) + Symbol("x"))) + (Integer(4))^((Integer(-1) + Symbol("x"))) + (Integer(4))^((Integer(2) + Symbol("x"))) + (Integer(-1)*(Integer(3))^((Integer(2) + Symbol("x")))) + (Integer(-1)*(Integer(3))^((Integer(3) + Symbol("x"))))), Symbol("x"), Reals)), "{2}")
  expect_identical(.pyrun(solveset(((Integer(-9)*exp((Integer(5) + (Integer(-2)*Symbol("x"))))) + (Integer(4)*exp((Integer(1) + (Integer(3)*Symbol("x")))))), Symbol("x"), Reals)), "{-2*log(2)/5 + 2*log(3)/5 + 4/5}")
  expect_identical(.pyrun(solveset((Integer(-84) + (Integer(2))^(Symbol("x")) + (Integer(4))^(Symbol("x")) + (Integer(8))^(Symbol("x"))), Symbol("x"), Reals)), "{2}")
  expect_identical(.pyrun(solveset(((Integer(-123)*(Integer(2726))^(Symbol("x"))) + (Integer(29)*(Integer(2))^((Integer(1) + Symbol("x")))*(Integer(615))^(Symbol("x")))), Symbol("x"), Reals)), "{(-log(29) - log(2) + log(123))/(-log(2726) + log(2) + log(615))}")
  expect_identical(.pyrun(solveset_real(((Integer(2))^((Integer(1) + Symbol("x"))) + (Integer(-9)*exp((Integer(5) + (Integer(-2)*Symbol("x")))))), Symbol("x"))), "{-(-5 - 2*log(3) + log(2))/(log(2) + 2)}")
  expect_identical(.pyrun(solveset_real(((Integer(4))^((Rational(1, 2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Rational(1, 3)*Symbol("x"))))), Symbol("x"))), "{0}")
  expect_identical(.pyrun(solveset_real(((Integer(5))^((Rational(1, 2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Integer(3)*(Symbol("x"))^(Integer(-1)))))), Symbol("x"))), "{-sqrt(6)*sqrt(log(2))/sqrt(log(5)), sqrt(6)*sqrt(log(2))/sqrt(log(5))}")
  expect_identical(.pyrun(solveset_real((Symbol("C1") + (Integer(-1)*exp((Integer(-1)*Function("f")(Symbol("x"))))) + (Symbol("C2")*(Symbol("x"))^(Integer(-2)))), Function("f")(Symbol("x")))), "Intersection({-log(C1 + C2/x^2)}, Reals)")
  expect_identical(.pyrun(solveset_real(((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y", positive = TRUE))^(Integer(2))*exp((Integer(-1)*Symbol("x"))))), Symbol("y", positive = TRUE))), "Intersection({-sqrt(x^2*exp(x)), sqrt(x^2*exp(x))}, Reals)")
  expect_identical(.pyrun(solveset_real(((Integer(1) + (Symbol("p", positive = TRUE))^(Integer(-1))))^((Integer(1) + Symbol("p", positive = TRUE))), Symbol("p", positive = TRUE))), "ConditionSet(_R, Eq((1 + 1/_R)^(_R + 1), 0), Reals)")
  expect_identical(.pyrun(solveset((Integer(12) + (Integer(2))^(Symbol("x")) + (Integer(-1)*(Integer(4))^(Symbol("x")))), Symbol("x"), Reals)), "{2}")
  expect_identical(.pyrun(solveset((Integer(12) + (Integer(2))^(Symbol("x")) + (Integer(-1)*(Integer(2))^((Integer(2)*Symbol("x"))))), Symbol("x"), Reals)), "{2}")
})

test_that("solveset: test_exponential_complex matches SymPy", {
  expect_identical(.pyrun(solveset_complex(((Integer(2))^(Symbol("x")) + (Integer(4))^(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(2^x + 1, 0), Complexes)")
  expect_identical(.pyrun(solveset_complex((Integer(-2) + ((Symbol("x"))^(Symbol("z"))*(Symbol("y"))^(Symbol("z")))), Symbol("z"))), "ConditionSet(z, Eq(x^z*y^z - 2, 0), Complexes)")
  expect_identical(.pyrun(solveset_complex(((Integer(4))^((Rational(1, 2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Rational(1, 3)*Symbol("x"))))), Symbol("x"))), "{0}")
  expect_identical(.pyrun(solveset((Integer(32) + (Integer(2))^(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(2^x + 32, 0), Complexes)")
  expect_identical(.pyrun(solveset_complex((((Integer(15) + (Symbol("x"))^(Integer(2))))^(Integer(-1))*(Integer(2) + (Integer(2))^(exp(((Symbol("x"))^(Integer(-1))*(Symbol("y"))^(Integer(2))))))), Symbol("y"))), "ConditionSet(y, Eq(2^exp(y^2/x) + 2, 0), Complexes)")
  expect_identical(.pyrun(solveset(((Integer(2))^(Symbol("x")) + (Integer(4))^(Symbol("x")) + (Integer(8))^(Symbol("x"))), Symbol("x"))), "ConditionSet(x, Eq(2^(2*x) + 2^x + 1, 0), Complexes)")
  expect_identical(.pyrun(solveset(((Integer(4))^((Integer(1) + Symbol("x"))) + (Integer(4))^((Integer(-1) + Symbol("x"))) + (Integer(4))^((Integer(2) + Symbol("x"))) + (Integer(-1)*(Integer(3))^((Integer(2) + Symbol("x")))) + (Integer(-1)*(Integer(3))^((Integer(3) + Symbol("x"))))), Symbol("x"))), "{2}")
})

test_that("solveset: test_expo_conditionset matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-2) + ((Integer(1) + exp(Symbol("x"))))^(Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, Eq((exp(x) + 1)^x - 2, 0), Reals)")
  expect_identical(.pyrun(solveset((Integer(-3) + (Symbol("x")*((Integer(2) + Symbol("x")))^(Symbol("y")))), Symbol("x"), Reals)), "ConditionSet(x, Eq(x*(x + 2)^y - 3, 0), Reals)")
  expect_identical(.pyrun(solveset((Integer(-3) + (Integer(2))^(Symbol("x")) + (Integer(-1)*exp(Symbol("x")))), Symbol("x"), Reals)), "ConditionSet(x, Eq(2^x - exp(x) - 3, 0), Reals)")
  expect_identical(.pyrun(solveset(((Integer(-1)*exp(Symbol("x"))) + log(Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(-exp(x) + log(x), 0), Reals)")
  expect_identical(.pyrun(solveset(((Integer(2))^(Symbol("x")) + (Integer(3))^(Symbol("x")) + (Integer(-1)*(Integer(5))^(Symbol("x")))), Symbol("x"), Reals)), "ConditionSet(x, Eq(2^x + 3^x - 5^x, 0), Reals)")
})

test_that("solveset: test_exponential_symbols matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("z", positive = TRUE))^(Symbol("x", positive = TRUE)) + (Integer(-1)*Symbol("y", positive = TRUE))), Symbol("x", positive = TRUE), Reals)), "Intersection({log(y)/log(z)}, Reals)")
  expect_identical(.pyrun(solveset(((Integer(-4)*(Symbol("y", positive = TRUE))^(Symbol("w"))) + (Integer(2)*(Symbol("x", positive = TRUE))^(Symbol("w")))), Symbol("w"), Reals)), "Intersection({log(2)/(log(x) - log(y))}, Reals)")
  expect_identical(.pyrun(solveset((Integer(-2) + ((Symbol("x", positive = TRUE)*(Symbol("y", positive = TRUE))^(Integer(-1))))^(Symbol("w"))), Symbol("w"), Reals)), "Intersection({log(2)/log(x/y)}, Reals)")
  expect_identical(.pyrun(solveset((Symbol("x", positive = TRUE))^(Symbol("x", positive = TRUE)), Symbol("x", positive = TRUE), Interval(Integer(0), oo, TRUE, TRUE))), "ConditionSet(_R, Eq(_R^_R, 0), Interval.open(0, oo))")
  expect_identical(.pyrun(solveset((Integer(-1) + (Symbol("x", positive = TRUE))^(Symbol("y", positive = TRUE))), Symbol("y", positive = TRUE), Reals)), "{0}")
  expect_identical(.pyrun(solveset((Integer(-2) + (exp((Symbol("x", positive = TRUE)*(Symbol("y", positive = TRUE))^(Integer(-1))))*exp((Integer(-1)*Symbol("z", positive = TRUE)*(Symbol("y", positive = TRUE))^(Integer(-1)))))), Symbol("y", positive = TRUE), Reals)), "Complement({(x - z)/log(2)}, {0})")
  expect_identical(.pyrun(solveset((Integer(-2) + (exp((Symbol("xr", real = TRUE)*(Symbol("y", positive = TRUE))^(Integer(-1))))*exp((Integer(-1)*Symbol("zr", real = TRUE)*(Symbol("y", positive = TRUE))^(Integer(-1)))))), Symbol("y", positive = TRUE), Reals)), "Complement({(xr - zr)/log(2)}, {0})")
  expect_identical(.pyrun(solveset(((Symbol("a"))^(Symbol("x", positive = TRUE)) + (Integer(-1)*(Symbol("b"))^(Symbol("x", positive = TRUE)))), Symbol("x", positive = TRUE))), "ConditionSet(_C, Ne(a, 0) & Ne(b, 0), {0})")
})

test_that("solveset: test_ignore_assumptions matches SymPy", {
  expect_identical(.pyrun(solveset_complex((Integer(-4) + (Symbol("x", positive = TRUE))^(Integer(2))), Symbol("x", positive = TRUE))), "{-2, 2}")
  expect_identical(.pyrun(solveset_complex((Integer(-4) + (Symbol("x"))^(Integer(2))), Symbol("x"))), "{-2, 2}")
})

test_that("solveset: test_issue_10864 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("x"))^((Symbol("y")*Symbol("z"))) + (Integer(-1)*Symbol("x"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(-x + x^(y*z), 0), Reals)")
})

test_that("solveset: test_solve_only_exp_2 matches SymPy", {
  expect_identical(.pyrun(solveset_real((Integer(-4) + (exp(Symbol("x")))^(Rational(1, 2)) + (exp((Integer(-1)*Symbol("x"))))^(Rational(1, 2))), Symbol("x"))), "ConditionSet(x, Eq(-4*exp(x/2) + exp(x) + 1, 0), Reals)")
})

test_that("solveset: test_is_exponential matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(Symbol("y"), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential((Integer(-2) + (Integer(3))^(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Integer(5))^(Symbol("x")) + (Integer(-1)*(Integer(7))^((Integer(2) + (Integer(-1)*Symbol("x")))))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Integer(-4)*Symbol("x")) + sin((Integer(2))^(Symbol("x")))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^(Symbol("y")) + (Integer(-1)*Symbol("z"))), Symbol("y"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^(Symbol("y")) + (Integer(-1)*Symbol("z"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential((Integer(-1) + (Integer(2))^(Symbol("x")) + (Integer(4))^(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^((Symbol("y")*Symbol("z"))) + (Integer(-1)*Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^((Integer(2)*Symbol("x"))) + (Integer(-1)*(Integer(3))^(Symbol("x")))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^(Symbol("y")) + (Integer(-1)*Symbol("y")*Symbol("z"))), Symbol("y"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((Symbol("x"))^(Symbol("y")) + (Integer(-1)*Symbol("x")*Symbol("z"))), Symbol("y"))), "True")
})

test_that("solveset: test_solve_exponential matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_solve_exponential(((Integer(3))^((Integer(2)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Integer(3) + Symbol("x"))))), Integer(0), Symbol("x"), Reals)), "{-3*log(2)/(-2*log(3) + log(2))}")
  expect_identical(.pyrun(symcompR:::.sset_solve_exponential(((Integer(2))^(Symbol("y")) + (Integer(4))^(Symbol("y"))), Integer(1), Symbol("y"), Reals)), "{log(-1/2 + sqrt(5)/2)/log(2)}")
  expect_identical(.pyrun(symcompR:::.sset_solve_exponential(((Integer(2))^(Symbol("y")) + (Integer(4))^(Symbol("y"))), Integer(0), Symbol("y"), Reals)), "EmptySet")
  expect_identical(.pyrun(symcompR:::.sset_solve_exponential(((Integer(2))^(Symbol("x")) + (Integer(3))^(Symbol("x")) + (Integer(-1)*(Integer(5))^(Symbol("x")))), Integer(0), Symbol("x"), Reals)), "ConditionSet(x, Eq(2^x + 3^x - 5^x, 0), Reals)")
})

test_that("solveset: test_logarithmic matches SymPy", {
  expect_identical(.pyrun(solveset_real((log((Integer(-3) + Symbol("x"))) + log((Integer(3) + Symbol("x")))), Symbol("x"))), "{-sqrt(10), sqrt(10)}")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*log((Integer(-1) + (Integer(2)*Symbol("x"))))) + log((Integer(1) + Symbol("x")))), Symbol("x"))), "{2}")
  expect_identical(.pyrun(solveset_real((Integer(-3) + log((Integer(1) + (Integer(3)*(Symbol("x"))^(Integer(-1))))) + log((Integer(3) + Symbol("x")))), Symbol("x"))), "{-3 + sqrt(-12 + exp(3))*exp(3/2)/2 + exp(3)/2, -sqrt(-12 + exp(3))*exp(3/2)/2 - 3 + exp(3)/2}")
  expect_identical(.pyrun(solveset_real((Symbol("z") + (Integer(-1)*log(Symbol("x"))) + log((Symbol("y")*(Symbol("x"))^(Integer(-1))*((Integer(-1) + ((Symbol("x"))^(Integer(-2))*(Symbol("y"))^(Integer(2)))))^(Integer(-1))))), Symbol("x"))), "Complement(Intersection({-sqrt(y^2 - y*exp(z)), sqrt(y^2 - y*exp(z))}, Reals), Intersection({-sqrt(y^2), sqrt(y^2)}, Reals))")
  expect_identical(.pyrun(solveset_real(((Integer(-1)*log((Integer(1) + (Integer(-1)*Symbol("x"))))) + (Integer(-1)*log((Integer(1) + (Integer(4)*Symbol("x"))))) + log((Integer(3)*Symbol("x")))), Symbol("x"))), "{-1/2, 1/2}")
  expect_identical(.pyrun(solveset(((Integer(-1)*Symbol("y")*log(Symbol("x"))) + log((Symbol("x"))^(Symbol("y")))), Symbol("x"), Reals)), "Reals")
})

test_that("solveset: test_uselogcombine_2 matches SymPy", {
  expect_identical(.pyrun(solveset_real(((Integer(-1)*log(Integer(2))) + log((Integer(1) + (Integer(-1)*tanh(Symbol("x"))))) + log((Integer(1) + exp((Integer(2)*Symbol("x")))))), Symbol("x"))), "ConditionSet(x, Eq((1 - tanh(x))*(exp(2*x) + 1) - 2, 0), Reals)")
  expect_identical(.pyrun(solveset_real((Integer(-2) + (Integer(-1)*log((Integer(1) + (Symbol("x"))^(Rational(1, 2))))) + log((Integer(8)*Symbol("x")))), Symbol("x"))), "{sqrt(exp(2) + 32)*exp(3)/128 + (exp(2) + 16)*exp(2)/128}")
})

test_that("solveset: test_is_logarithmic matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic(Symbol("y"), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic(log(Symbol("x")), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((Integer(-3) + log(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((log(Symbol("x"))*log(Symbol("y"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((log(Symbol("x")))^(Integer(2)), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((log((Integer(-3) + Symbol("x"))) + log((Integer(3) + Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic(((Integer(-1)*Symbol("y")*log(Symbol("x"))) + log((Symbol("x"))^(Symbol("y")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic(sin(log(Symbol("x"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((Symbol("x") + Symbol("y")), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((Integer(4) + (Integer(-1)*log((Integer(1) + (Integer(-1)*Symbol("x"))))) + log((Integer(3)*Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((Symbol("x") + log(Symbol("x")) + log(Symbol("y"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((log((Integer(-3) + Symbol("x"))) + log(log((Integer(-3) + Symbol("x"))))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((log(Symbol("x")) + log((Symbol("x") + log(Integer(3))))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_logarithmic((((Integer(3) + Symbol("y"))*log(Symbol("x"))) + log(Symbol("x"))), Symbol("y"))), "False")
})

test_that("solveset: test_solve_logarithm matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_solve_logarithm(((Integer(-1)*Symbol("y")*log(Symbol("x"))) + log((Symbol("x"))^(Symbol("y")))), Integer(0), Symbol("x"), Reals)), "Reals")
  expect_identical(.pyrun(symcompR:::.sset_solve_logarithm((log(Symbol("x"))*log(Symbol("y", positive = TRUE))), Integer(0), Symbol("x"), Reals)), "{1}")
})

test_that("solveset: test_is_lambert matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Symbol("x"))^(Integer(2)), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Symbol("c") + (Symbol("a"))^((Symbol("x"))^(Integer(2))) + (Symbol("b")*Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert(exp(Integer(2)), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Symbol("x")*exp(Integer(2))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert(((Integer(3)*log(Symbol("x"))) + (Integer(-1)*Symbol("x")*log(Integer(3)))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((log((Integer(-3) + Symbol("x"))) + log(log((Integer(-3) + Symbol("x"))))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Integer(-1) + (Integer(3)*exp((Integer(2) + (Integer(-7)*Symbol("x"))))) + (Integer(5)*Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert(((Rational(1, 4)*exp((Rational(1, 2)*Symbol("x")))) + (Integer(2)*Symbol("a")*(Symbol("x"))^(Integer(-3)))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Integer(3) + ((Integer(-1) + (((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2))))^(Integer(2)) + (Integer(-2)*(((Integer(3)*Symbol("x")) + log(Symbol("x"))))^(Integer(2)))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Integer(-1) + (Symbol("x")*sinh(Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert((Integer(-5) + (Symbol("x")*cos(Symbol("x")))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert(((Integer(-5)*Symbol("x")) + tanh(Symbol("x"))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_lambert(((Integer(-1)*sinh(Symbol("x"))) + cosh(Symbol("x"))), Symbol("x"))), "False")
})

test_that("solveset: test_is_modular matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_modular(Symbol("y"), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-1) + sym_Mod(Symbol("x"), Integer(3))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-1) + sym_Mod((Integer(1) + (Symbol("x"))^(Integer(3)) + (Integer(-1)*Symbol("x")) + (Integer(-3)*(Symbol("x"))^(Integer(2)))), Integer(3))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-2) + sym_Mod(exp((Symbol("x") + Symbol("y"))), Integer(3))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_modular(((Integer(-1)*log(Symbol("x"))) + sym_Mod(exp((Symbol("x") + Symbol("y"))), Integer(3))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-1) + sym_Mod(Symbol("x"), Integer(3))), Symbol("y"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-5) + (sym_Mod(Symbol("x"), Integer(3)))^(Integer(2))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_modular(((sym_Mod(Symbol("x"), Integer(3)))^(Integer(2)) + (Integer(-1)*Symbol("y"))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-1) + exp(sym_Mod(Symbol("x"), Integer(3)))), Symbol("x"))), "False")
  expect_identical(.pyrun(symcompR:::.sset_is_modular((Integer(-1) + sym_Mod(Integer(3), Symbol("y"))), Symbol("y"))), "False")
})

test_that("solveset: test_invert_modular matches SymPy", {
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(Symbol("x"), Integer(12)), Rational(1, 2), .D1, Symbol("x"))})), "(x, EmptySet)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(sin(Symbol("x")), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(Mod(sin(x), 7), 5)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(exp(Symbol("x")), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(Mod(exp(x), 7), 5)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(log(Symbol("x")), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(Mod(log(x), 7), 5)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(Symbol("x"), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(x, ImageSet(Lambda(_n, 7*_n + 5), Integers))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Integer(1) + Symbol("x")), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(x, ImageSet(Lambda(_n, 7*_n + 4), Integers))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Symbol("x") + (Symbol("x"))^(Integer(2))), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(Mod(x^2 + x, 7), 5)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Integer(3)*Symbol("x")), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(x, ImageSet(Lambda(_n, 7*_n + 4), Integers))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod(((Integer(1) + Symbol("x"))*(Integer(2) + Symbol("x"))), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(Mod((x + 1)*(x + 2), 7), 5)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Symbol("x"))^(Integer(4)), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(x, EmptySet)")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Integer(3))^(Symbol("x")), Integer(4)), Integer(3), .D1, Symbol("x"))})), "(x, ImageSet(Lambda(_n, 2*_n + 1), Naturals0))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((Integer(2))^((Integer(1) + Symbol("x") + (Symbol("x"))^(Integer(2)))), Integer(7)), Integer(2), .D1, Symbol("x"))})), "(x^2 + x + 1, ImageSet(Lambda(_n, 3*_n + 1), Naturals0))")
  expect_identical(.pyrun(local({.D1 <- Dummy("n", integer = TRUE); symcompR:::.sset_invert_modular(sym_Mod((sin(Symbol("x")))^(Integer(4)), Integer(7)), Integer(5), .D1, Symbol("x"))})), "(x, EmptySet)")
})

test_that("solveset: test_solve_modular matches SymPy", {
  expect_identical(.pyrun(solveset(((Integer(-1)*Symbol("x")) + sym_Mod(Symbol("x"), Integer(4))), Symbol("x"), Integers)), "ConditionSet(x, Eq(-x + Mod(x, 4), 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod(sin(Symbol("x")), Integer(7)))), Symbol("x"), Integers)), "ConditionSet(x, Eq(Mod(sin(x), 7) - 3, 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod(log(Symbol("x")), Integer(7)))), Symbol("x"), Integers)), "ConditionSet(x, Eq(Mod(log(x), 7) - 3, 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod(exp(Symbol("x")), Integer(7)))), Symbol("x"), Integers)), "ConditionSet(x, Eq(Mod(exp(x), 7) - 3, 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(7) + (Integer(-1)*sym_Mod(Symbol("x"), Integer(5)))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(5) + (Integer(-1)*sym_Mod(Symbol("x"), Integer(5)))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(2) + sym_Mod(Symbol("x"), Integer(-3))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, -3*_n - 2), Integers)")
  expect_identical(.pyrun(solveset((Integer(4) + sym_Mod(Symbol("x"), Integer(-3))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod(Symbol("x"), Integer(5)))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 5*_n + 3), Integers)")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod((Integer(6) + (Integer(5)*Symbol("x"))), Integer(7)))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 7*_n + 5), Integers)")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod((Integer(5)*Symbol("x")), Integer(7)))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 7*_n + 2), Integers)")
  expect_identical(.pyrun(solveset((Integer(-9) + sym_Mod((Symbol("x"))^(Integer(2)), Integer(160))), Symbol("x"), Integers)), "Union(ImageSet(Lambda(_n, 160*_n + 3), Integers), ImageSet(Lambda(_n, 160*_n + 13), Integers), ImageSet(Lambda(_n, 160*_n + 67), Integers), ImageSet(Lambda(_n, 160*_n + 77), Integers), ImageSet(Lambda(_n, 160*_n + 83), Integers), ImageSet(Lambda(_n, 160*_n + 93), Integers), ImageSet(Lambda(_n, 160*_n + 147), Integers), ImageSet(Lambda(_n, 160*_n + 157), Integers))")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod((Symbol("x"))^(Integer(4)), Integer(7)))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(-13) + sym_Mod((Symbol("x"))^(Integer(4)), Integer(17))), Symbol("x"), Integers)), "Union(ImageSet(Lambda(_n, 17*_n + 3), Integers), ImageSet(Lambda(_n, 17*_n + 5), Integers), ImageSet(Lambda(_n, 17*_n + 12), Integers), ImageSet(Lambda(_n, 17*_n + 14), Integers))")
  expect_identical(.pyrun(solveset((Integer(-15) + sym_Mod((Integer(7))^(Symbol("x")), Integer(41))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 40*_n + 3), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(-18) + sym_Mod((Integer(12))^(Symbol("x")), Integer(21))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 6*_n + 2), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(-3) + sym_Mod((Integer(3))^(Symbol("x")), Integer(4))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 2*_n + 1), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(-2) + sym_Mod((Integer(2))^(Symbol("x")), Integer(7))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 3*_n + 1), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(-3) + sym_Mod((Integer(3))^((Integer(3))^(Symbol("x"))), Integer(4))), Symbol("x"), Integers)), "Intersection(ImageSet(Lambda(_n, Intersection({log(2*_n + 1)/log(3)}, Integers)), Naturals0), Integers)")
  expect_identical(.pyrun(solveset((Integer(-2) + sym_Mod((Symbol("x"))^(Integer(3)), Integer(7))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(-1) + sym_Mod((Symbol("x"))^(Integer(3)), Integer(8))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 8*_n + 1), Integers)")
  expect_identical(.pyrun(solveset((Integer(-4) + sym_Mod((Symbol("x"))^(Integer(4)), Integer(9))), Symbol("x"), Integers)), "Union(ImageSet(Lambda(_n, 9*_n + 4), Integers), ImageSet(Lambda(_n, 9*_n + 5), Integers))")
  expect_identical(.pyrun(solveset((Integer(3) + (Integer(-1)*sym_Mod((Integer(6) + (Integer(5)*Symbol("x"))), Integer(7)))), Symbol("x"), Naturals0)), "Intersection(ImageSet(Lambda(_n, 7*_n + 5), Integers), Naturals0)")
  expect_identical(.pyrun(solveset(((Integer(-1)*sym_I) + sym_Mod(Symbol("x"), Integer(3))), Symbol("x"), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(-2) + sym_Mod((sym_I*Symbol("x")), Integer(3))), Symbol("x"), Integers)), "ConditionSet(x, Eq(Mod(I*x, 3) - 2, 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(-2) + sym_Mod((sym_I + Symbol("x")), Integer(3))), Symbol("x"), Integers)), "ConditionSet(x, Eq(Mod(x + I, 3) - 2, 0), Integers)")
  expect_identical(.pyrun(solveset((Integer(-11) + sym_Mod((Symbol("x"))^(Integer(4)), Integer(14))), Symbol("x"), Integers)), "Union(ImageSet(Lambda(_n, 14*_n + 3), Integers), ImageSet(Lambda(_n, 14*_n + 11), Integers))")
  expect_identical(.pyrun(solveset((Integer(-43) + sym_Mod((Symbol("x"))^(Integer(31)), Integer(74))), Symbol("x"), Integers)), "ImageSet(Lambda(_n, 74*_n + 31), Integers)")
  expect_identical(.pyrun(solveset((Integer(20170816) + (Integer(-1)*sym_Mod((Integer(1898888478)*(Integer(742938285))^(Symbol("n", integer = TRUE))), Integer(2147483647)))), Symbol("n", integer = TRUE), Integers)), "ImageSet(Lambda(_n, 2147483646*_n + 100), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(20170816) + (Integer(-1)*sym_Mod((Integer(1898888478)*(Integer(742938285))^(Symbol("n", integer = TRUE))), Integer(2147483647)))), Symbol("n", integer = TRUE), Naturals0)), "Intersection(ImageSet(Lambda(_n, 2147483646*_n + 100), Naturals0), Naturals0)")
  expect_identical(.pyrun(solveset((Integer(20170816) + (Integer(-1)*sym_Mod((Integer(1898888478)*(Integer(742938285))^((Integer(2)*Symbol("n", integer = TRUE)))), Integer(2147483647)))), Symbol("n", integer = TRUE), Integers)), "Intersection(ImageSet(Lambda(_n, 1073741823*_n + 50), Naturals0), Integers)")
  expect_identical(.pyrun(solveset((Integer(20170816) + (Integer(-1)*sym_Mod((Integer(1898888478)*(Integer(742938285))^((Integer(7) + (Integer(2)*Symbol("n", integer = TRUE))))), Integer(2147483647)))), Symbol("n", integer = TRUE), Integers)), "EmptySet")
  expect_identical(.pyrun(solveset((Integer(20170816) + (Integer(-1)*sym_Mod((Integer(1898888478)*(Integer(742938285))^((Integer(-4) + Symbol("n", integer = TRUE)))), Integer(2147483647)))), Symbol("n", integer = TRUE), Integers)), "Intersection(ImageSet(Lambda(_n, 2147483646*_n + 104), Naturals0), Integers)")
})

test_that("solveset: test_issue_10426 matches SymPy", {
  expect_identical(.pyrun(local({.D1 <- Dummy("x"); solveset(((Integer(-1)*sin(.D1)) + sin((Symbol("a") + .D1))), Symbol("a"))})), "Union(ImageSet(Lambda(_n, 2*_n*pi - _x + asin(sin(_x))), Integers), ImageSet(Lambda(_n, 2*_n*pi - _x - asin(sin(_x)) + pi), Integers))")
})

test_that("solveset: test_solveset_conjugate matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-3) + sym_I + conjugate(Symbol("x"))))), "{3 + I}")
})

test_that("solveset: test_issue_18208 matches SymPy", {
  expect_identical(.pyrun(linear_eq_to_matrix(list((Integer(-51) + Symbol("x0") + Symbol("x1") + Symbol("x2") + Symbol("x3")), (Integer(-46) + Symbol("x0") + Symbol("x1") + Symbol("x4") + Symbol("x5")), (Integer(-39) + Symbol("x2") + Symbol("x3") + Symbol("x6") + Symbol("x7")), (Integer(-50) + Symbol("x0") + Symbol("x3") + Symbol("x4") + Symbol("x7")), (Integer(-35) + Symbol("x1") + Symbol("x2") + Symbol("x5") + Symbol("x6")), (Integer(-34) + Symbol("x4") + Symbol("x5") + Symbol("x6") + Symbol("x7")), (Integer(-46) + Symbol("x4") + Symbol("x5") + Symbol("x8") + Symbol("x9")), (Integer(-23) + Symbol("x10") + Symbol("x11") + Symbol("x6") + Symbol("x7")), (Integer(-25) + Symbol("x11") + Symbol("x4") + Symbol("x7") + Symbol("x8")), (Integer(-44) + Symbol("x10") + Symbol("x5") + Symbol("x6") + Symbol("x9")), (Integer(-35) + Symbol("x10") + Symbol("x11") + Symbol("x8") + Symbol("x9")), (Integer(-35) + Symbol("x12") + Symbol("x13") + Symbol("x8") + Symbol("x9")), (Integer(-29) + Symbol("x10") + Symbol("x11") + Symbol("x14") + Symbol("x15")), (Integer(-35) + Symbol("x11") + Symbol("x12") + Symbol("x15") + Symbol("x8")), (Integer(-29) + Symbol("x10") + Symbol("x13") + Symbol("x14") + Symbol("x9")), (Integer(-29) + Symbol("x12") + Symbol("x13") + Symbol("x14") + Symbol("x15")), (Integer(-55) + Symbol("y0") + Symbol("y1") + Symbol("y2") + Symbol("y3")), (Integer(-53) + Symbol("y0") + Symbol("y1") + Symbol("y4") + Symbol("y5")), (Integer(-56) + Symbol("y2") + Symbol("y3") + Symbol("y6") + Symbol("y7")), (Integer(-57) + Symbol("y0") + Symbol("y3") + Symbol("y4") + Symbol("y7")), (Integer(-52) + Symbol("y1") + Symbol("y2") + Symbol("y5") + Symbol("y6")), (Integer(-54) + Symbol("y4") + Symbol("y5") + Symbol("y6") + Symbol("y7")), (Integer(-48) + Symbol("y4") + Symbol("y5") + Symbol("y8") + Symbol("y9")), (Integer(-60) + Symbol("y10") + Symbol("y11") + Symbol("y6") + Symbol("y7")), (Integer(-51) + Symbol("y11") + Symbol("y4") + Symbol("y7") + Symbol("y8")), (Integer(-57) + Symbol("y10") + Symbol("y5") + Symbol("y6") + Symbol("y9")), (Integer(-54) + Symbol("y10") + Symbol("y11") + Symbol("y8") + Symbol("y9")), (Integer(-2) + Symbol("x10")), (Integer(-5) + Symbol("x11")), (Integer(-1) + Symbol("x12")), (Integer(-6) + Symbol("x13")), (Integer(-1) + Symbol("x14")), (Integer(-21) + Symbol("x15")), (Integer(-12) + Symbol("y0")), (Integer(-20) + Symbol("y1"))), list(Symbol("x0"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15"), Symbol("y0"), Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("y7"), Symbol("y8"), Symbol("y9"), Symbol("y10"), Symbol("y11")))), "(Matrix([\\n[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\\n[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), Matrix([\\n[51],\\n[46],\\n[39],\\n[50],\\n[35],\\n[34],\\n[46],\\n[23],\\n[25],\\n[44],\\n[35],\\n[35],\\n[29],\\n[35],\\n[29],\\n[29],\\n[55],\\n[53],\\n[56],\\n[57],\\n[52],\\n[54],\\n[48],\\n[60],\\n[51],\\n[57],\\n[54],\\n[ 2],\\n[ 5],\\n[ 1],\\n[ 6],\\n[ 1],\\n[21],\\n[12],\\n[20]]))")
  expect_identical(.pyrun(linsolve(list((Integer(-51) + Symbol("x0") + Symbol("x1") + Symbol("x2") + Symbol("x3")), (Integer(-46) + Symbol("x0") + Symbol("x1") + Symbol("x4") + Symbol("x5")), (Integer(-39) + Symbol("x2") + Symbol("x3") + Symbol("x6") + Symbol("x7")), (Integer(-50) + Symbol("x0") + Symbol("x3") + Symbol("x4") + Symbol("x7")), (Integer(-35) + Symbol("x1") + Symbol("x2") + Symbol("x5") + Symbol("x6")), (Integer(-34) + Symbol("x4") + Symbol("x5") + Symbol("x6") + Symbol("x7")), (Integer(-46) + Symbol("x4") + Symbol("x5") + Symbol("x8") + Symbol("x9")), (Integer(-23) + Symbol("x10") + Symbol("x11") + Symbol("x6") + Symbol("x7")), (Integer(-25) + Symbol("x11") + Symbol("x4") + Symbol("x7") + Symbol("x8")), (Integer(-44) + Symbol("x10") + Symbol("x5") + Symbol("x6") + Symbol("x9")), (Integer(-35) + Symbol("x10") + Symbol("x11") + Symbol("x8") + Symbol("x9")), (Integer(-35) + Symbol("x12") + Symbol("x13") + Symbol("x8") + Symbol("x9")), (Integer(-29) + Symbol("x10") + Symbol("x11") + Symbol("x14") + Symbol("x15")), (Integer(-35) + Symbol("x11") + Symbol("x12") + Symbol("x15") + Symbol("x8")), (Integer(-29) + Symbol("x10") + Symbol("x13") + Symbol("x14") + Symbol("x9")), (Integer(-29) + Symbol("x12") + Symbol("x13") + Symbol("x14") + Symbol("x15")), (Integer(-55) + Symbol("y0") + Symbol("y1") + Symbol("y2") + Symbol("y3")), (Integer(-53) + Symbol("y0") + Symbol("y1") + Symbol("y4") + Symbol("y5")), (Integer(-56) + Symbol("y2") + Symbol("y3") + Symbol("y6") + Symbol("y7")), (Integer(-57) + Symbol("y0") + Symbol("y3") + Symbol("y4") + Symbol("y7")), (Integer(-52) + Symbol("y1") + Symbol("y2") + Symbol("y5") + Symbol("y6")), (Integer(-54) + Symbol("y4") + Symbol("y5") + Symbol("y6") + Symbol("y7")), (Integer(-48) + Symbol("y4") + Symbol("y5") + Symbol("y8") + Symbol("y9")), (Integer(-60) + Symbol("y10") + Symbol("y11") + Symbol("y6") + Symbol("y7")), (Integer(-51) + Symbol("y11") + Symbol("y4") + Symbol("y7") + Symbol("y8")), (Integer(-57) + Symbol("y10") + Symbol("y5") + Symbol("y6") + Symbol("y9")), (Integer(-54) + Symbol("y10") + Symbol("y11") + Symbol("y8") + Symbol("y9")), (Integer(-2) + Symbol("x10")), (Integer(-5) + Symbol("x11")), (Integer(-1) + Symbol("x12")), (Integer(-6) + Symbol("x13")), (Integer(-1) + Symbol("x14")), (Integer(-21) + Symbol("x15")), (Integer(-12) + Symbol("y0")), (Integer(-20) + Symbol("y1"))), list(Symbol("x0"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15"), Symbol("y0"), Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("y7"), Symbol("y8"), Symbol("y9"), Symbol("y10"), Symbol("y11")))), "{(38 - x3, x3 - 10, 23 - x3, x3, 12 - x7, x7 + 6, 16 - x7, x7, 8, 20, 2, 5, 1, 6, 1, 21, 12, 20, -y11 + y9 + 2, y11 - y9 + 21, -y11 - y7 + y9 + 24, y11 + y7 - y9 - 3, 33 - y7, y7, 27 - y9, y9, 27 - y11, y11)}")
  expect_identical(.pyrun(linsolve(list(Matrix(list(list(Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(1)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)), list(Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)))), Matrix(list(list(Integer(51)), list(Integer(46)), list(Integer(39)), list(Integer(50)), list(Integer(35)), list(Integer(34)), list(Integer(46)), list(Integer(23)), list(Integer(25)), list(Integer(44)), list(Integer(35)), list(Integer(35)), list(Integer(29)), list(Integer(35)), list(Integer(29)), list(Integer(29)), list(Integer(55)), list(Integer(53)), list(Integer(56)), list(Integer(57)), list(Integer(52)), list(Integer(54)), list(Integer(48)), list(Integer(60)), list(Integer(51)), list(Integer(57)), list(Integer(54)), list(Integer(2)), list(Integer(5)), list(Integer(1)), list(Integer(6)), list(Integer(1)), list(Integer(21)), list(Integer(12)), list(Integer(20))))), list(Symbol("x0"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15"), Symbol("y0"), Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("y7"), Symbol("y8"), Symbol("y9"), Symbol("y10"), Symbol("y11")))), "{(38 - x3, x3 - 10, 23 - x3, x3, 12 - x7, x7 + 6, 16 - x7, x7, 8, 20, 2, 5, 1, 6, 1, 21, 12, 20, -y11 + y9 + 2, y11 - y9 + 21, -y11 - y7 + y9 + 24, y11 + y7 - y9 - 3, 33 - y7, y7, 27 - y9, y9, 27 - y11, y11)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-51) + Symbol("x0") + Symbol("x1") + Symbol("x2") + Symbol("x3")), (Integer(-46) + Symbol("x0") + Symbol("x1") + Symbol("x4") + Symbol("x5")), (Integer(-39) + Symbol("x2") + Symbol("x3") + Symbol("x6") + Symbol("x7")), (Integer(-50) + Symbol("x0") + Symbol("x3") + Symbol("x4") + Symbol("x7")), (Integer(-35) + Symbol("x1") + Symbol("x2") + Symbol("x5") + Symbol("x6")), (Integer(-34) + Symbol("x4") + Symbol("x5") + Symbol("x6") + Symbol("x7")), (Integer(-46) + Symbol("x4") + Symbol("x5") + Symbol("x8") + Symbol("x9")), (Integer(-23) + Symbol("x10") + Symbol("x11") + Symbol("x6") + Symbol("x7")), (Integer(-25) + Symbol("x11") + Symbol("x4") + Symbol("x7") + Symbol("x8")), (Integer(-44) + Symbol("x10") + Symbol("x5") + Symbol("x6") + Symbol("x9")), (Integer(-35) + Symbol("x10") + Symbol("x11") + Symbol("x8") + Symbol("x9")), (Integer(-35) + Symbol("x12") + Symbol("x13") + Symbol("x8") + Symbol("x9")), (Integer(-29) + Symbol("x10") + Symbol("x11") + Symbol("x14") + Symbol("x15")), (Integer(-35) + Symbol("x11") + Symbol("x12") + Symbol("x15") + Symbol("x8")), (Integer(-29) + Symbol("x10") + Symbol("x13") + Symbol("x14") + Symbol("x9")), (Integer(-29) + Symbol("x12") + Symbol("x13") + Symbol("x14") + Symbol("x15")), (Integer(-55) + Symbol("y0") + Symbol("y1") + Symbol("y2") + Symbol("y3")), (Integer(-53) + Symbol("y0") + Symbol("y1") + Symbol("y4") + Symbol("y5")), (Integer(-56) + Symbol("y2") + Symbol("y3") + Symbol("y6") + Symbol("y7")), (Integer(-57) + Symbol("y0") + Symbol("y3") + Symbol("y4") + Symbol("y7")), (Integer(-52) + Symbol("y1") + Symbol("y2") + Symbol("y5") + Symbol("y6")), (Integer(-54) + Symbol("y4") + Symbol("y5") + Symbol("y6") + Symbol("y7")), (Integer(-48) + Symbol("y4") + Symbol("y5") + Symbol("y8") + Symbol("y9")), (Integer(-60) + Symbol("y10") + Symbol("y11") + Symbol("y6") + Symbol("y7")), (Integer(-51) + Symbol("y11") + Symbol("y4") + Symbol("y7") + Symbol("y8")), (Integer(-57) + Symbol("y10") + Symbol("y5") + Symbol("y6") + Symbol("y9")), (Integer(-54) + Symbol("y10") + Symbol("y11") + Symbol("y8") + Symbol("y9")), (Integer(-2) + Symbol("x10")), (Integer(-5) + Symbol("x11")), (Integer(-1) + Symbol("x12")), (Integer(-6) + Symbol("x13")), (Integer(-1) + Symbol("x14")), (Integer(-21) + Symbol("x15")), (Integer(-12) + Symbol("y0")), (Integer(-20) + Symbol("y1"))), list(Symbol("x0"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("x10"), Symbol("x11"), Symbol("x12"), Symbol("x13"), Symbol("x14"), Symbol("x15"), Symbol("y0"), Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("y7"), Symbol("y8"), Symbol("y9"), Symbol("y10"), Symbol("y11")))), "{(38 - x3, x3 - 10, 23 - x3, x3, 12 - x7, x7 + 6, 16 - x7, x7, 8, 20, 2, 5, 1, 6, 1, 21, 12, 20, -y5 + y7 - 1, y5 - y7 + 24, 21 - y5, y5, 33 - y7, y7, 27 - y9, y9, -y5 + y7 - y9 + 24, y5 - y7 + y9 + 3)}")
})

test_that("solveset: test_substitution_with_infeasible_solution matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Symbol("c00") + Symbol("c01") + Symbol("m0") + (Integer(-1)*Symbol("c00")*Symbol("l0")) + (Integer(-1)*Symbol("c01")*Symbol("l1"))), (Symbol("m1") + (Integer(-1)*Symbol("c10")*Symbol("l0")) + (Integer(-1)*Symbol("c11")*Symbol("l1"))), (Symbol("c00") + Symbol("c01") + (Integer(-1)*Symbol("c00")*Symbol("l2")) + (Integer(-1)*Symbol("c01")*Symbol("l3"))), (Symbol("m3") + (Integer(-1)*Symbol("c10")*Symbol("l2")) + (Integer(-1)*Symbol("c11")*Symbol("l3"))), (Symbol("p00") + Symbol("p10") + (Integer(-1)*Symbol("l0")*Symbol("p00")) + (Integer(-1)*Symbol("l2")*Symbol("p10"))), (Symbol("p00") + Symbol("p10") + (Integer(-1)*Symbol("l1")*Symbol("p00")) + (Integer(-1)*Symbol("l3")*Symbol("p10"))), ((Integer(-1)*Symbol("l0")*Symbol("p01")) + (Integer(-1)*Symbol("l2")*Symbol("p11"))), ((Integer(-1)*Symbol("l1")*Symbol("p01")) + (Integer(-1)*Symbol("l3")*Symbol("p11"))), ((Integer(-1)*Symbol("a00")) + (Symbol("c00")*Symbol("p00")) + (Symbol("c10")*Symbol("p01"))), ((Integer(-1)*Symbol("a01")) + (Symbol("c01")*Symbol("p00")) + (Symbol("c11")*Symbol("p01"))), ((Integer(-1)*Symbol("a10")) + (Symbol("c00")*Symbol("p10")) + (Symbol("c10")*Symbol("p11"))), ((Integer(-1)*Symbol("a11")) + (Symbol("c01")*Symbol("p10")) + (Symbol("c11")*Symbol("p11"))), (Integer(-1)*Symbol("m0")*Symbol("p00")), (Integer(-1)*Symbol("m1")*Symbol("p01")), (Integer(-1)*Symbol("m2")*Symbol("p10")), (Integer(-1)*Symbol("m3")*Symbol("p11")), (Integer(-1)*Symbol("c00")*Symbol("m4")), (Integer(-1)*Symbol("c01")*Symbol("m5")), (Integer(-1)*Symbol("c10")*Symbol("m6")), (Integer(-1)*Symbol("c11")*Symbol("m7")), Symbol("m2"), Symbol("m4"), Symbol("m5"), Symbol("m6"), Symbol("m7")), list(Symbol("p00"), Symbol("p01"), Symbol("p10"), Symbol("p11"), Symbol("c00"), Symbol("c01"), Symbol("c10"), Symbol("c11"), Symbol("m0"), Symbol("m1"), Symbol("m3"), Symbol("l0"), Symbol("l1"), Symbol("l2"), Symbol("l3")))), "EmptySet")
})

test_that("solveset: test_issue_20097 matches SymPy", {
  expect_identical(.pyrun(solveset((Symbol("x"))^(Rational(-1, 2)))), "EmptySet")
})

test_that("solveset: test_issue_15350 matches SymPy", {
  expect_identical(.pyrun(solveset((((Symbol("x") + (Symbol("x"))^(Integer(-1))))^(Rational(-1, 2))*(Rational(1, 2) + (Rational(-1, 2)*(Symbol("x"))^(Integer(-2))))))), "{-1, 1}")
})

test_that("solveset: test_issue_18359 matches SymPy", {
  expect_identical(.pyrun(solveset((Rational(-1, 2) + Piecewise(list(Integer(0), Lt(Symbol("x"), Integer(0))), list(((Rational(1, 2)*Min(Integer(1), Symbol("x"))) + (Rational(1, 2)*Min(Integer(3), Symbol("x"))) + (Rational(-1, 2)*Min(Integer(2), Symbol("x")))), S_true))), Symbol("x"), Interval(Integer(0), Integer(3), FALSE, FALSE))), "Interval(1, 2)")
  expect_identical(.pyrun(solveset((Rational(-1, 2) + Piecewise(list(Piecewise(list(Integer(0), Lt(Symbol("x"), Integer(0))), list(((Rational(1, 2)*Min(Integer(1), Symbol("x"))) + (Rational(1, 2)*Min(Integer(3), Symbol("x"))) + (Rational(-1, 2)*Min(Integer(2), Symbol("x")))), S_true)), Ge(Symbol("x"), Integer(0))), list(Integer(0), S_true))), Symbol("x"), Interval(Integer(0), Integer(3), FALSE, FALSE))), "Interval(1, 2)")
})

test_that("solveset: test_issue_17604 matches SymPy", {
  expect_identical(.pyrun(symcompR:::.sset_is_exponential(((S_Pi)^((Rational(1, 11)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Rational(3, 11)*Symbol("x")))*exp((Rational(1, 11)*Symbol("x"))))), Symbol("x"))), "True")
  expect_identical(.pyrun(symcompR:::.sset_solve_exponential(((S_Pi)^((Rational(1, 11)*Symbol("x"))) + (Integer(-1)*(Integer(2))^((Rational(3, 11)*Symbol("x")))*exp((Rational(1, 11)*Symbol("x"))))), Integer(0), Symbol("x"), Complexes)), "{0}")
})

test_that("solveset: test_issue_17580 matches SymPy", {
  expect_identical(.pyrun(solveset(((Integer(1) + (Integer(-1)*(Symbol("x"))^(Integer(3)))))^(Integer(-2)), Symbol("x"), Reals)), "EmptySet")
})

test_that("solveset: test_issue_17566_actual matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-3) + (Integer(2))^(Symbol("x")) + (Integer(2))^(Symbol("y"))), (Integer(-5) + (Integer(4))^(Symbol("x")) + (Integer(9))^(Symbol("y")))), Symbol("x"), Symbol("y"))), "{(log(3 - 2^y)/log(2), y)}")
})

test_that("solveset: test_issue_17565 matches SymPy", {
  expect_identical(.pyrun(solveset(Ge(((((Integer(1) + Symbol("x")))^(Rational(2, 3))*(Integer(-4) + (Integer(2)*Symbol("x")))) + (Rational(2, 3)*((Integer(1) + Symbol("x")))^(Rational(-1, 3))*((Integer(-2) + Symbol("x")))^(Integer(2)))), Integer(0)), Symbol("x"), Reals)), "Union(Interval.Lopen(-1, -1/4), Interval(2, oo))")
})

test_that("solveset: test_issue_15024 matches SymPy", {
  expect_identical(.pyrun(solveset(((((Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-10)*Symbol("x"))))^(Rational(-1, 2))*(Integer(5) + Symbol("x"))), Symbol("x"), Reals)), "{-5}")
})

test_that("solveset: test_issue_16877 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + Symbol("x")), sin(Symbol("y"))), Symbol("x"), Symbol("y"))), "{(1, ImageSet(Lambda(_n, 2*_n*pi), Integers)), (1, ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))}")
})

test_that("solveset: test_issue_16876 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(sin(Symbol("x")), ((Integer(-4)*Symbol("y")) + (Integer(2)*Symbol("x")))), Symbol("x"), Symbol("y"))), "{(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, _n*pi), Integers)), (ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), ImageSet(Lambda(_n, _n*pi + pi/2), Integers))}")
})

test_that("solveset: test_issue_21236 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("x"))^(Symbol("y", rational = TRUE)) + (Integer(-1)*Symbol("z"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(x^y - z, 0), Reals)")
  expect_identical(.pyrun(solveset(((Symbol("x"))^((Symbol("e1", even = TRUE)*(Symbol("e2", even = TRUE))^(Integer(-1)))) + (Integer(-1)*Symbol("z"))), Symbol("x"), Reals)), "ConditionSet(x, Eq(x^(e1/e2) - z, 0), Reals)")
})

test_that("solveset: test_issue_21908 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((((Symbol("x"))^(Integer(2)) + (Integer(-1)*(Symbol("y"))^(Integer(2))) + (Integer(2)*Symbol("x")))*exp(Symbol("x"))), (Integer(-2)*Symbol("y")*exp(Symbol("x")))), Symbol("x"), Symbol("y"))), "{(-2, 0), (0, 0)}")
})

test_that("solveset: test_issue_19144 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + Symbol("x") + Symbol("y")), (Integer(1) + (Symbol("y"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "{(1 - I, I), (1 + I, -I)}")
  expect_identical(.pyrun(nonlinsolve(list(Eq((Integer(-1) + Symbol("x") + Symbol("y")), Integer(0)), Eq((Integer(1) + (Symbol("y"))^(Integer(2))), Integer(0))), list(Symbol("x"), Symbol("y")))), "{(1 - I, I), (1 + I, -I)}")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x")*(Symbol("y"))^(Integer(-1)))), (Integer(1) + (Symbol("y"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "{(-I, -I), (I, I)}")
  expect_identical(.pyrun(nonlinsolve(list(Eq((Integer(-1) + (Symbol("x")*(Symbol("y"))^(Integer(-1)))), Integer(0)), Eq((Integer(1) + (Symbol("y"))^(Integer(2))), Integer(0))), list(Symbol("x"), Symbol("y")))), "{(-I, -I), (I, I)}")
  expect_identical(.pyrun(nonlinsolve(list(Eq((Symbol("x") + (Symbol("x"))^(Integer(-1))), (Symbol("x"))^(Integer(-1)))), list(Symbol("x")))), "EmptySet")
})

test_that("solveset: test_issue_22413 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(4)*Symbol("y")*(Integer(1) + (Integer(2)*Symbol("x")) + (Integer(2)*exp(Symbol("y"))))*exp((Integer(2)*Symbol("x")))), (Integer(1) + (Integer(4)*exp((Symbol("y") + (Integer(2)*Symbol("x"))))) + (Integer(4)*Symbol("x")*exp((Integer(2)*Symbol("x")))) + (Integer(4)*Symbol("y")*exp((Symbol("y") + (Integer(2)*Symbol("x"))))))), Symbol("x"), Symbol("y"))), "{(x, 0), (-exp(y) - 1/2, y)}")
})

test_that("solveset: test_issue_23318 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(Eq(Float("53.5780461486929", precision = 53), (Symbol("x")*(Symbol("y"))^(Integer(-1))*log((Integer(1) + (Symbol("y")*((Float("5.00000000000000", precision = 53) + (Integer(-1)*Symbol("y"))))^(Integer(-1))))))), Eq(Symbol("x"), (Float("0.00150000000000000", precision = 53)*Symbol("z"))), Eq(Float("0.00150000000000000", precision = 53), (Float("7845.32000000000", precision = 53)*Symbol("y")*(Symbol("z"))^(Integer(-1))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(266.97755814852, 0.0340301680681629, 177985.03876568)}")
  expect_identical(.pyrun(nonlinsolve(list((Float("53.5780461486929", precision = 53) + (Integer(-1)*Symbol("x")*(Symbol("y"))^(Integer(-1))*log((Integer(1) + (Symbol("y")*((Float("5.00000000000000", precision = 53) + (Integer(-1)*Symbol("y"))))^(Integer(-1))))))), (Symbol("x") + (Float("-0.00150000000000000", precision = 53)*Symbol("z"))), (Float("0.00150000000000000", precision = 53) + (Float("-7845.32000000000", precision = 53)*Symbol("y")*(Symbol("z"))^(Integer(-1))))), list(Symbol("x"), Symbol("y"), Symbol("z")))), "{(266.97755814852, 0.0340301680681629, 177985.03876568)}")
  expect_identical(.pyrun(solveset(((Float("1.02439504345316e-5", precision = 53)*Symbol("z")) + (Float("-0.00150000000000000", precision = 53)*Symbol("z")*log((Integer(1) + (Float("1.91196789933362e-7", precision = 53)*Symbol("z")*((Float("5.00000000000000", precision = 53) + (Float("-1.91196789933362e-7", precision = 53)*Symbol("z"))))^(Integer(-1))))))), Symbol("z"))), "{0, 177985.038765679}")
})

test_that("solveset: test_issue_19814 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(2))^(Symbol("m")) + (Integer(-1)*(Integer(2))^((Integer(2)*Symbol("n"))))), ((Integer(-1)*(Integer(2))^((Integer(4)*Symbol("n")))) + (Integer(4)*(Integer(2))^(Symbol("m"))))), Symbol("m"), Symbol("n"))), "{(log(2^(2*n))/log(2), Complexes)}")
})

test_that("solveset: test_issue_22058 matches SymPy", {
  expect_identical(.pyrun(solveset(((Symbol("t"))^(Rational(1, 2)) + (Integer(2)*Symbol("x")) + (Integer(-1)*(Symbol("t"))^(Rational(1, 2))*(Symbol("x"))^(Integer(2)))), Symbol("x"), Reals)), "Intersection({-sqrt(t + 1)/sqrt(t) + 1/sqrt(t), sqrt(t + 1)/sqrt(t) + 1/sqrt(t)}, Reals)")
})

test_that("solveset: test_issue_11184 matches SymPy", {
  expect_identical(.pyrun(solveset((Integer(-60) + (Integer(20)*(((Symbol("y"))^(Integer(2)) + ((Integer(10) + (((Integer(10) + Symbol("y"))*(Integer(10) + (Integer(-1)*Symbol("y")))))^(Rational(1, 2))))^(Integer(2))))^(Rational(1, 2)))), Symbol("y"), Reals)), "EmptySet")
})

test_that("solveset: test_issue_21890 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-2)*Symbol("y")) + (Integer(4)*(Symbol("x"))^(Integer(3))*(Symbol("y"))^(Integer(4)))), ((Integer(-2)*Symbol("x")) + (Integer(4)*(Symbol("x"))^(Integer(4))*(Symbol("y"))^(Integer(3))))), Symbol("x"), Symbol("y"))), "{(2^(2/3)/(2*y), y), ((-2^(2/3)/4 - 2^(2/3)*sqrt(3)*I/4)/y, y), ((-2^(2/3)/4 + 2^(2/3)*sqrt(3)*I/4)/y, y)}")
  expect_identical(.pyrun(nonlinsolve(list(((Integer(1) + (Integer(-4)*(Symbol("x"))^(Integer(2))))*exp(((Integer(-2)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*(Symbol("y"))^(Integer(2)))))), (Integer(-4)*Symbol("x")*Symbol("y")*exp((Integer(-2)*(Symbol("x"))^(Integer(2))))*exp((Integer(-2)*(Symbol("y"))^(Integer(2)))))), Symbol("x"), Symbol("y"))), "{(-1/2, 0), (1/2, 0)}")
  expect_identical(.pyrun(nonlinsolve(list(((Integer(-2)*Symbol("y", real = TRUE)) + (Integer(4)*(Symbol("x", real = TRUE))^(Integer(3))*(Symbol("y", real = TRUE))^(Integer(4)))), ((Integer(-2)*Symbol("x", real = TRUE)) + (Integer(4)*(Symbol("x", real = TRUE))^(Integer(4))*(Symbol("y", real = TRUE))^(Integer(3))))), Symbol("x", real = TRUE), Symbol("y", real = TRUE))), "{(2^(2/3)/(2*y), y), ((-2^(2/3)/4 - 2^(2/3)*sqrt(3)*I/4)/y, y), ((-2^(2/3)/4 + 2^(2/3)*sqrt(3)*I/4)/y, y)}")
})

test_that("solveset: test_issue_22628 matches SymPy", {
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + Symbol("h")), (Integer(-1) + Symbol("k")), (Integer(-2) + Symbol("f")), (Integer(-4) + Symbol("f")), (Integer(-2)*Symbol("k"))), Symbol("h"), Symbol("k"), Symbol("f"))), "EmptySet")
  expect_identical(.pyrun(nonlinsolve(list((Integer(-1) + (Symbol("x"))^(Integer(3))), (Symbol("x") + Symbol("y")), (Integer(-4) + (Symbol("x"))^(Integer(2)))), list(Symbol("x"), Symbol("y")))), "EmptySet")
})

test_that("solveset: test_issue_30584 matches SymPy", {
  expect_identical(.pyrun(solveset(Gt((((Integer(1) + Symbol("y")))^(Integer(2)) + (Integer(-1)*Symbol("y")*(Integer(2) + Symbol("y")))), Integer(0)), Symbol("x"), Reals)), "Reals")
  expect_identical(.pyrun(solveset(Lt((((Integer(1) + Symbol("y")))^(Integer(2)) + (Integer(-1)*Symbol("y")*(Integer(2) + Symbol("y")))), Integer(0)), Symbol("x"), Reals)), "EmptySet")
})

test_that("solveset: test_solveset_symbol_independent_relational matches SymPy", {
  expect_identical(.pyrun(solveset(Gt(Symbol("y"), Integer(0)), Symbol("x"), Reals)), "ConditionSet(x, y > 0, Reals)")
  expect_identical(.pyrun(solveset(Eq(Symbol("y"), Integer(0)), Symbol("x"), Reals)), "ConditionSet(x, Eq(y, 0), Reals)")
  expect_identical(.pyrun(solveset(Symbol("y"), Symbol("x"), Reals)), "EmptySet")
})

test_that("solveset: test_solveset_symbol_independent_relational_after_expansion matches SymPy", {
  expect_identical(.pyrun(solveset((((Symbol("x") + Symbol("y")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x")*Symbol("y"))), Symbol("x"), Reals)), "EmptySet")
  expect_identical(.pyrun(solveset(Eq((((Symbol("x") + Symbol("y")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x")*Symbol("y"))), Integer(0)), Symbol("x"), Reals)), "ConditionSet(x, Eq(y^2, 0), Reals)")
  expect_identical(.pyrun(solveset(Gt((((Symbol("x") + Symbol("y")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x")*Symbol("y"))), Integer(0)), Symbol("x"), Reals)), "ConditionSet(x, y^2 > 0, Reals)")
  expect_identical(.pyrun(solveset(Ge((((Symbol("x") + Symbol("y")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x")*Symbol("y"))), (Symbol("y"))^(Integer(2))), Symbol("x"), Reals)), "ConditionSet(x, y^2 >= y^2, Reals)")
  expect_identical(.pyrun(solveset(Gt((((Symbol("x") + Symbol("y")))^(Integer(2)) + (Integer(-1)*(Symbol("x"))^(Integer(2))) + (Integer(-2)*Symbol("x")*Symbol("y"))), (Symbol("y"))^(Integer(2))), Symbol("x"), Reals)), "ConditionSet(x, y^2 > y^2, Reals)")
  expect_identical(.pyrun(solveset(Gt((((Integer(1) + Symbol("x")))^(Integer(2)) + (Integer(-1)*Symbol("x")*(Integer(2) + Symbol("x")))), Integer(0)), Symbol("x"), Reals)), "Reals")
  expect_identical(.pyrun(solveset(Lt((((Integer(1) + Symbol("x")))^(Integer(2)) + (Integer(-1)*Symbol("x")*(Integer(2) + Symbol("x")))), Integer(0)), Symbol("x"), Reals)), "EmptySet")
})

test_that("condset: test_CondSet matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(sin(Symbol("x")), Integer(0)), Interval(Integer(0), (Integer(2)*S_Pi), FALSE, TRUE))), "ConditionSet(x, Eq(sin(x), 0), Interval.Ropen(0, 2*pi))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x"))^(Integer(2)), Integer(4)), Reals)), "ConditionSet(x, x^2 > 4, Reals)")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt(Symbol("y"), Integer(5)), Interval(Integer(1), Integer(7), FALSE, FALSE))), "ConditionSet(x, y > 5, Interval(1, 7))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(1)), list(Symbol("x"), Symbol("y")))), "ERR SympifyError: SympifyError: [x, y]")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), (Integer(1) + Symbol("x")), list(Symbol("x"), Symbol("y")))), "ERR SympifyError: SympifyError: [x, y]")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Symbol("x"), Integer(1))), "ERR TypeError: base set should be a Set object, not 1")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(1)))), "ConditionSet(x, x < 1)")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x")), Lt(Symbol("x"), Integer(1)))), "ConditionSet((x,), x < 1)")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), Lt(Symbol("x"), Symbol("y")), ProductSet(Integers, Integers))), "ConditionSet((x, y), x < y, ProductSet(Integers, Integers))")
})

test_that("condset: test_CondSet_intersect matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x"))^(Integer(2)), Integer(4)), Interval(Integer(1), Integer(4), FALSE, FALSE))), "ConditionSet(x, x^2 > 4, Interval(1, 4))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x"))^(Integer(2)), Integer(4)), Interval(Integer(1), Integer(3), FALSE, FALSE))), "ConditionSet(x, x^2 > 4, Interval(1, 3))")
  expect_identical(.pyrun(Intersection(ConditionSet(Symbol("x"), Gt((Symbol("x"))^(Integer(2)), Integer(4)), Interval(Integer(1), Integer(4), FALSE, FALSE)), Interval(Integer(0), Integer(3), FALSE, FALSE))), "ConditionSet(x, x^2 > 4, Interval(1, 3))")
})

test_that("condset: test_issue_9849 matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), S_true, Naturals)), "Naturals")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), S_false, Naturals)), "EmptySet")
})

test_that("condset: test_simplified_FiniteSet_in_CondSet matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), And(Gt(Symbol("x"), Integer(-3)), Lt(Symbol("x"), Integer(1))), FiniteSet(Integer(0), Integer(1), Integer(2)))), "{0}")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(0)), FiniteSet(Integer(0), Integer(1), Integer(2)))), "EmptySet")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(-3)), S_EmptySet)), "EmptySet")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt(Symbol("x"), Integer(0)), FiniteSet(Integer(-1), Integer(0), Integer(1), Symbol("y")))), "Union(ConditionSet(x, x > 0, {y}), {1})")
  expect_identical(.pyrun(Union(FiniteSet(Integer(1)), ConditionSet(Symbol("x"), Gt(Symbol("x"), Integer(0)), FiniteSet(Symbol("y"))))), "Union(ConditionSet(x, x > 0, {y}), {1})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(sym_Mod(Symbol("x"), Integer(3)), Integer(1)), FiniteSet(Integer(1), Integer(2), Integer(4), Symbol("y")))), "Union(ConditionSet(x, Eq(Mod(x, 3), 1), {y}), {1, 4})")
  expect_identical(.pyrun(Union(FiniteSet(Integer(1), Integer(4)), ConditionSet(Symbol("x"), Eq(sym_Mod(Symbol("x"), Integer(3)), Integer(1)), FiniteSet(Symbol("y"))))), "Union(ConditionSet(x, Eq(Mod(x, 3), 1), {y}), {1, 4})")
})

test_that("condset: test_free_symbols matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(Symbol("y"), Integer(0)), FiniteSet(Symbol("z")))), "ConditionSet(x, Eq(y, 0), {z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(Symbol("x"), Integer(0)), FiniteSet(Symbol("z")))), "ConditionSet(x, Eq(x, 0), {z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(Symbol("x"), Integer(0)), FiniteSet(Symbol("x"), Symbol("z")))), "ConditionSet(x, Eq(x, 0), {x, z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(Symbol("x"), Integer(0)), ImageSet(Lambda(Symbol("y"), (Symbol("y"))^(Integer(2))), Integers))), "ConditionSet(x, Eq(x, 0), ImageSet(Lambda(y, y^2), Integers))")
})

test_that("condset: test_bound_symbols matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(Symbol("x"), Integer(0)), FiniteSet(Symbol("x"), Symbol("y")))), "ConditionSet(x, Eq(x, 0), {x, y})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(10)), ImageSet(Lambda(Symbol("y"), (Symbol("y"))^(Integer(2))), Integers))), "ConditionSet(x, x < 10, ImageSet(Lambda(y, y^2), Integers))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(10)), ConditionSet(Symbol("y"), Gt(Symbol("y"), Integer(1)), Integers))), "ConditionSet(x, (x > 1) & (x < 10), Integers)")
})

test_that("condset: test_as_dummy matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(1)), Interval(Symbol("y"), oo, FALSE, TRUE))), "ConditionSet(x, x < 1, Interval(y, oo))")
  expect_identical(.pyrun(ConditionSet(Symbol("_0"), Lt(Symbol("_0"), Integer(1)), Interval(Symbol("y"), oo, FALSE, TRUE))), "ConditionSet(_0, _0 < 1, Interval(y, oo))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(1)), Interval(Symbol("x"), oo, FALSE, TRUE))), "ConditionSet(x, x < 1, Interval(x, oo))")
  expect_identical(.pyrun(ConditionSet(Symbol("_0"), Lt(Symbol("_0"), Integer(1)), Interval(Symbol("x"), oo, FALSE, TRUE))), "ConditionSet(_0, _0 < 1, Interval(x, oo))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(1)), ImageSet(Lambda(Symbol("y"), (Symbol("y"))^(Integer(2))), Integers))), "ConditionSet(x, x < 1, ImageSet(Lambda(y, y^2), Integers))")
  expect_identical(.pyrun(ConditionSet(Symbol("_0"), Lt(Symbol("_0"), Integer(1)), ImageSet(Lambda(Symbol("_1"), (Symbol("_1"))^(Integer(2))), Integers))), "ConditionSet(_0, _0 < 1, ImageSet(Lambda(_1, _1^2), Integers))")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), Le(Symbol("x"), Symbol("y")), ProductSet(Reals, Reals))), "ConditionSet((x, y), x <= y, ProductSet(Reals, Reals))")
  expect_identical(.pyrun(ConditionSet(list(Symbol("_0"), Symbol("_1")), Le(Symbol("_0"), Symbol("_1")), ProductSet(Reals, Reals))), "ConditionSet((_0, _1), _0 <= _1, ProductSet(Reals, Reals))")
  expect_identical(.pyrun(ConditionSet(list(Symbol("y"), Symbol("x")), Le(Symbol("y"), Symbol("x")), ProductSet(Reals, Reals))), "ConditionSet((y, x), y <= x, ProductSet(Reals, Reals))")
})

test_that("condset: test_subs_CondSet matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(2)), FiniteSet(Symbol("y"), Symbol("z")))), "ConditionSet(x, x < 2, {y, z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(2)), FiniteSet(Symbol("y")))), "ConditionSet(x, x < 2, {y})")
  expect_identical(.pyrun(ConditionSet(Symbol("y"), Lt(Symbol("y"), Integer(2)), FiniteSet(Symbol("y"), Symbol("z")))), "ConditionSet(y, y < 2, {y, z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Symbol("y")), FiniteSet(Symbol("y"), Symbol("z")))), "ConditionSet(x, x < y, {y, z})")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Symbol("w")), FiniteSet(Symbol("w"), Symbol("z")))), "ConditionSet(x, x < w, {w, z})")
  expect_identical(.pyrun(ConditionSet(Symbol("n", negative = TRUE), S_false, Integers)), "EmptySet")
  expect_identical(.pyrun(ConditionSet(Symbol("n", negative = TRUE), Lt(Symbol("n", negative = TRUE), Symbol("y")), Integers)), "ConditionSet(n, n < y, Integers)")
  expect_identical(.pyrun(ConditionSet((Integer(1) + Symbol("x")), Lt(Symbol("x"), Integer(1)), Integers)), "ERR ValueError: `x + 1` is not symbol-like")
  expect_identical(.pyrun(ConditionSet(Symbol("p", positive = TRUE), Lt(Symbol("n", negative = TRUE), Symbol("x")), Interval(Integer(-5), Integer(5), FALSE, FALSE))), "ConditionSet(p, n < x, Interval(-5, 5))")
  expect_identical(.pyrun(ConditionSet(Symbol("n", negative = TRUE), Lt(Symbol("n", negative = TRUE), Symbol("x")), Interval((-oo), Integer(0), TRUE, FALSE))), "ConditionSet(n, n < x, Interval(-oo, 0))")
  expect_identical(.pyrun(ConditionSet(Function("f")(Symbol("x")), Lt(Function("f")(Symbol("x")), Integer(1)), list(Symbol("w"), Symbol("z")))), "ERR SympifyError: SympifyError: [w, z]")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("k"), ((Integer(2)*S_Pi*Symbol("k")) + asin(Symbol("y")))), Integers)), "ImageSet(Lambda(k, 2*pi*k + asin(y)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("k"), ((Integer(2)*S_Pi*Symbol("k")) + asin(Rational(1, 3)))), Integers)), "ImageSet(Lambda(k, 2*pi*k + asin(1/3)), Integers)")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), Lt((Symbol("x") + Symbol("y")), Integer(3)), ProductSet(Integers, Integers))), "ConditionSet((x, y), x + y < 3, ProductSet(Integers, Integers))")
  expect_identical(.pyrun(ConditionSet(Symbol("n", negative = TRUE), Lt(Symbol("n", negative = TRUE), Integer(-10)), Interval(Integer(0), Integer(10), FALSE, FALSE))), "ERR TypeError: sym `n` is not in base_set `Interval(0, 10)`")
})

test_that("condset: test_subs_CondSet_tebr matches SymPy", {
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), list((Integer(1) + Symbol("x")), (Symbol("x") + Symbol("y"))), ProductSet(Reals, Reals))), "ERR SympifyError: SympifyError: [x + 1, x + y]")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), And(Eq((Integer(1) + Symbol("x")), Integer(0)), Eq((Symbol("x") + Symbol("y")), Integer(0))), ProductSet(Reals, Reals))), "ConditionSet((x, y), Eq(x + 1, 0) & Eq(x + y, 0), ProductSet(Reals, Reals))")
})

test_that("condset: test_dummy_eq matches SymPy", {
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), And(Eq((Integer(1) + Symbol("x")), Integer(0)), Eq((Symbol("x") + Symbol("y")), Integer(0))), ProductSet(Complexes, Complexes))), "ConditionSet((x, y), Eq(x + 1, 0) & Eq(x + y, 0), ProductSet(Complexes, Complexes))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("m"), (S_Pi*Symbol("m"))), Integers)), "ImageSet(Lambda(m, pi*m), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (S_Pi*Symbol("n"))), Integers)), "ImageSet(Lambda(n, pi*n), Integers)")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Ne(Symbol("a"), Integer(0)), ImageSet(Lambda(Symbol("m"), (S_Pi*Symbol("m"))), Integers))), "ConditionSet(x, Ne(a, 0), ImageSet(Lambda(m, pi*m), Integers))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Ne(Symbol("a"), Integer(0)), ImageSet(Lambda(Symbol("n"), (S_Pi*Symbol("n"))), Integers))), "ConditionSet(x, Ne(a, 0), ImageSet(Lambda(n, pi*n), Integers))")
})

test_that("condset: test_contains matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt(Symbol("x"), Integer(5)), Interval(Integer(1), Integer(7), FALSE, FALSE))), "ConditionSet(x, x > 5, Interval(1, 7))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Ge((Symbol("x"))^(Integer(-1)), Integer(0)), Reals)), "ConditionSet(x, 1/x >= 0, Reals)")
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), Gt((Symbol("x") + Symbol("y")), Integer(1)), ProductSet(Integers, Integers))), "ConditionSet((x, y), x + y > 1, ProductSet(Integers, Integers))")
  expect_identical(.pyrun(ConditionSet(list(Symbol("w"), list(Symbol("x"), Symbol("y"))), Gt((Symbol("w") + Symbol("x") + Symbol("y")), Integer(1)), ProductSet(Integers, ProductSet(Integers, Integers)))), "ConditionSet((w, (x, y)), w + x + y > 1, ProductSet(Integers, ProductSet(Integers, Integers)))")
})

test_that("condset: test_as_relational matches SymPy", {
  expect_identical(.pyrun(ConditionSet(list(Symbol("x"), Symbol("y")), Gt(Symbol("x"), Integer(1)), ProductSet(Integers, Integers))), "ConditionSet((x, y), x > 1, ProductSet(Integers, Integers))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt(Symbol("x"), Integer(1)), Integers)), "ConditionSet(x, x > 1, Integers)")
})

test_that("condset: test_flatten matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x") + sin(Symbol("x"))), Integer(0)))), "ConditionSet(x, x + sin(x) > 0)")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x") + sin(Symbol("x"))), Integer(0)), Reals)), "ConditionSet(x, x + sin(x) > 0, Reals)")
  expect_identical(.pyrun(ConditionSet(Symbol("y"), Gt((Symbol("y") + sin(Symbol("y"))), Integer(0)))), "ConditionSet(y, y + sin(y) > 0)")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Gt((Symbol("x") + sin(Symbol("x"))), Integer(0)), Interval(Integer(-1), Integer(1), FALSE, FALSE))), "ConditionSet(x, x + sin(x) > 0, Interval(-1, 1))")
})

test_that("condset: test_duplicate matches SymPy", {
  expect_identical(.pyrun(ConditionSet(list(Symbol("a"), Symbol("a")), Lt(Symbol("x"), Integer(0)))), "ERR BadSignatureError: Duplicate symbols detected")
})

test_that("condset: test_SetKind_ConditionSet matches SymPy", {
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Eq(sin(Symbol("x")), Integer(0)), Interval(Integer(0), (Integer(2)*S_Pi), FALSE, FALSE))), "ConditionSet(x, Eq(sin(x), 0), Interval(0, 2*pi))")
  expect_identical(.pyrun(ConditionSet(Symbol("x"), Lt(Symbol("x"), Integer(0)))), "ConditionSet(x, x < 0)")
})

test_that("fancysets: test_ImageSet matches SymPy", {
  expect_identical(.pyrun(ImageSet(Symbol("x"), Integers)), "ERR ValueError: First argument must be a Lambda")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), Integer(1)), Integers)), "{1}")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), Symbol("y")), Integers)), "{y}")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), Integer(1)), S_EmptySet)), "EmptySet")
  expect_identical(.pyrun(Intersection(FiniteSet(((S_Pi)^(Integer(-1))*log(Integer(2)))), Integers)), "Intersection({log(2)/pi}, Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Symbol("x"))^(Integer(2))), Naturals)), "ImageSet(Lambda(x, x^2), Naturals)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Symbol("x"))^(Integer(-1))), Naturals)), "ImageSet(Lambda(x, 1/x), Naturals)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Integer(-1)*Symbol("x")), Interval(Integer(0), Integer(1), FALSE, FALSE))), "Interval(-1, 0)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Symbol("x"))^(Integer(2))), Interval(Integer(0), Integer(2), FALSE, FALSE))), "ImageSet(Lambda(x, x^2), Interval(0, 2))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Symbol("x"), Symbol("y")), (Integer(2)*Symbol("x"))), list(Integer(4)), list(Integer(3)))), "ERR SympifyError: SympifyError: [4]")
  expect_identical(.pyrun(ImageSet(Lambda(list(Symbol("x"), Symbol("y")), (Symbol("x") + Symbol("y"))), list(Integer(1), Integer(2), Integer(3)), list(Integer(10), Integer(20), Integer(30)))), "ERR SympifyError: SympifyError: [1, 2, 3]")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), Tuple(Symbol("x"), (Integer(2)*Symbol("y")))), ProductSet(Interval(Integer(1), Integer(3), FALSE, FALSE), Interval(Integer(1), Integer(3), FALSE, FALSE)))), "ImageSet(Lambda(((x, y),), (x, 2*y)), ProductSet(Interval(1, 3), Interval(1, 3)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), Tuple(Symbol("x"), (Symbol("y"))^(Integer(-1)))), ProductSet(Interval(Integer(1), Integer(3), FALSE, FALSE), Interval(Integer(1), Integer(3), FALSE, FALSE)))), "ImageSet(Lambda(((x, y),), (x, 1/y)), ProductSet(Interval(1, 3), Interval(1, 3)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), Tuple(Symbol("x"), (Symbol("y"))^(Integer(2)))), ProductSet(Interval(Integer(1), Integer(3), FALSE, FALSE), Interval(Integer(1), Integer(3), FALSE, FALSE)))), "ImageSet(Lambda(((x, y),), (x, y^2)), ProductSet(Interval(1, 3), Interval(1, 3)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), Tuple(Symbol("x"), Integer(-2))), ProductSet(Interval(Integer(1), Integer(3), FALSE, FALSE), Interval(Integer(1), Integer(3), FALSE, FALSE)))), "ImageSet(Lambda(((x, y),), (x, -2)), ProductSet(Interval(1, 3), Interval(1, 3)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("t"), Symbol("y"), Symbol("x"))), Tuple(Symbol("y"), Symbol("t"), Symbol("x"))), ProductSet(Interval(Integer(3), Integer(7), FALSE, FALSE), Interval(Integer(8), Integer(11), FALSE, FALSE), Interval(Integer(5), Integer(9), FALSE, FALSE)))), "ImageSet(Lambda(((t, y, x),), (y, t, x)), ProductSet(Interval(3, 7), Interval(8, 11), Interval(5, 9)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("t"), Symbol("y"), Symbol("x"))), Tuple((Symbol("y"))^(Integer(-1)), Symbol("t"), Symbol("x"))), ProductSet(Interval(Integer(3), Integer(7), FALSE, FALSE), Interval(Integer(8), Integer(11), FALSE, FALSE), Interval(Integer(5), Integer(9), FALSE, FALSE)))), "ImageSet(Lambda(((t, y, x),), (1/y, t, x)), ProductSet(Interval(3, 7), Interval(8, 11), Interval(5, 9)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), (Integer(2)*(Symbol("x"))^(Integer(-1)))), ProductSet(Interval(Integer(1), Integer(3), FALSE, FALSE), Interval(Integer(1), Integer(3), FALSE, FALSE)))), "ImageSet(Lambda(((x, y),), 2/x), ProductSet(Interval(1, 3), Interval(1, 3)))")
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("x"), Symbol("y"))), (Symbol("x") + Symbol("y"))), list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4))))), "ERR SympifyError: SympifyError: [(1, 2), (3, 4)]")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Symbol("x"))^(Integer(2))), Integer(1))), "ERR TypeError: Set arguments to ImageSet should of type Set")
})

test_that("fancysets: test_image_is_ImageSet matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), (sin(Symbol("x")))^(Rational(1, 2)), Range(Integer(0), Integer(5), Integer(1)))), "ImageSet(Lambda(x, sqrt(sin(x))), Range(0, 5, 1))")
})

test_that("fancysets: test_halfcircle matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(list(Tuple(Symbol("r", real = TRUE), Symbol("theta", real = TRUE))), Tuple((Symbol("r", real = TRUE)*cos(Symbol("theta", real = TRUE))), (Symbol("r", real = TRUE)*sin(Symbol("theta", real = TRUE))))), ProductSet(Interval(Integer(0), Integer(1), FALSE, FALSE), Interval(Integer(0), S_Pi, FALSE, FALSE)))), "ImageSet(Lambda(((r, theta),), (r*cos(theta), r*sin(theta))), ProductSet(Interval(0, 1), Interval(0, pi)))")
})

test_that("fancysets: test_ImageSet_iterator_not_injective matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Symbol("x") + (Integer(-1)*sym_Mod(Symbol("x"), Integer(2))))), Naturals)), "ImageSet(Lambda(x, x - Mod(x, 2)), Naturals)")
})

test_that("fancysets: test_Integers_eval_imageset matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Rational(3, 7) + (Integer(2)*Symbol("x")))), Integers)), "ImageSet(Lambda(x, 2*x + 3/7), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Rational(3, 7) + (Integer(-2)*Symbol("x")))), Integers)), "ImageSet(Lambda(x, 2*x + 3/7), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Rational(-11, 7) + (Integer(-2)*Symbol("x")))), Integers)), "ImageSet(Lambda(x, 2*x + 3/7), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("y") + (Integer(2)*Symbol("x"))), Integers)), "ImageSet(Lambda(x, 2*x + y), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Float("0.0920000000000000", precision = 53) + (Float("0.341000000000000", precision = 53)*Symbol("x"))), Integers)), "ImageSet(Lambda(x, 0.341*x + 0.092), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Float("0.433000000000000", precision = 53) + (Float("0.341000000000000", precision = 53)*Symbol("x"))), Integers)), "ImageSet(Lambda(x, 0.341*x + 0.433), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x", negative = TRUE), (Integer(1) + (Symbol("x", negative = TRUE))^(Integer(2)) + (Integer(-1)*Symbol("x", negative = TRUE))), Integers)), "ImageSet(Lambda(x, x^2 + x + 1), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x", negative = TRUE), (Integer(-1) + (Integer(3)*Symbol("x", negative = TRUE))), Integers)), "ImageSet(Lambda(x, 3*x + 2), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), list(Symbol("x"), (Symbol("x"))^(Integer(-1))), Integers)), "ImageSet(Lambda(x, (x, 1/x)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), Tuple(Symbol("x"), (Symbol("x"))^(Integer(-1)))), Integers)), "ImageSet(Lambda(x, (x, 1/x)), Integers)")
})

test_that("fancysets: test_Range_eval_imageset matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("c") + (Symbol("a")*(Symbol("b") + Symbol("x")))), Range(Integer(0), Integer(3), Integer(1)))), "ImageSet(Lambda(x, a*b + a*x + c), Range(0, 3, 1))")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("c") + (Symbol("a")*Symbol("b")) + (Symbol("a")*Symbol("x"))), Range(Integer(0), Integer(3), Integer(1)))), "ImageSet(Lambda(x, a*b + a*x + c), Range(0, 3, 1))")
  expect_identical(.pyrun(imageset(Symbol("x"), ((Integer(1) + Symbol("x")))^(Integer(2)), Range(Integer(0), Integer(3), Integer(1)))), "ImageSet(Lambda(x, (x + 1)^2), Range(0, 3, 1))")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("c") + (Symbol("a")*(Symbol("b") + Symbol("x")))), Range(Integer(3), Integer(-3), Integer(-2)))), "ImageSet(Lambda(x, a*b - 2*a*x + 3*a + c), Range(0, 3, 1))")
})

test_that("fancysets: test_fun matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), sin((Rational(1, 4)*S_Pi*Symbol("x")))), Range(Integer(-10), Integer(11), Integer(1)))), "ImageSet(Lambda(x, sin(pi*x/4)), Range(-10, 11, 1))")
})

test_that("fancysets: test_infinitely_indexed_set_1 matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(2)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 2*n), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("m"), (Integer(1) + (Integer(2)*Symbol("m")))), Integers)), "ImageSet(Lambda(m, 2*m + 1), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(1) + (Integer(2)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*n + 1), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("m"), (Integer(2)*Symbol("m"))), Integers)), "ImageSet(Lambda(m, 2*m), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(3)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 3*n), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("t"), (Integer(6)*Symbol("t"))), Integers)), "ImageSet(Lambda(t, 6*t), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Rational(1, 3) + (Rational(1, 2)*Symbol("x"))), Integers)), "ImageSet(Lambda(x, x/2 + 1/3), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Rational(1, 2) + (Rational(1, 2)*Symbol("x"))), Integers)), "ImageSet(Lambda(x, x/2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(3) + (Integer(5)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 5*n + 3), Integers)")
})

test_that("fancysets: test_infinitely_indexed_set_2 matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Symbol("n") + Symbol("a", integer = TRUE))), Integers)), "Integers")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (S_Pi + Symbol("n"))), Integers)), "ImageSet(Lambda(n, n - 3 + pi), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (S_Pi + Symbol("n") + Symbol("a", integer = TRUE))), Integers)), "ImageSet(Lambda(n, n - 3 + pi), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Symbol("a", integer = TRUE) + (Integer(-1)*Symbol("n")))), Integers)), "Integers")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(-6)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 6*n), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(6)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 6*n), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (S_Pi + (Integer(2)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*n - 2 + pi), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(-2) + S_Pi + (Integer(2)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*n - 2 + pi), Integers)")
})

test_that("fancysets: test_imageset_intersect_real matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Symbol("n") + (sym_I*(Integer(1) + Symbol("n"))*(Integer(-1) + Symbol("n"))))), Integers)), "ImageSet(Lambda(n, n + I*(n - 1)*(n + 1)), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Symbol("n") + (sym_I*(Rational(1, 2) + Symbol("n"))*(Integer(-1) + Symbol("n"))))), Integers)), "ImageSet(Lambda(n, n + I*(n - 1)*(n + 1/2)), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Symbol("n") + (sym_I*(Integer(1) + Symbol("n"))*(Rational(1, 2) + Symbol("n"))*(Integer(-1) + Symbol("n"))))), Naturals0)), "ImageSet(Lambda(n, n + I*(n - 1)*(n + 1/2)*(n + 1)), Naturals0)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), ((Rational(1, 2)*Symbol("n")) + (sym_I*(Rational(-1, 2) + (Symbol("n"))^(Integer(2)) + (Rational(-1, 2)*Symbol("n")))))), Integers)), "ImageSet(Lambda(n, n/2 + I*(n^2 - n/2 - 1/2)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Rational(1, 2)*Symbol("x"))), ConditionSet(Symbol("n"), Eq((Rational(-1, 2) + (Symbol("n"))^(Integer(2)) + (Rational(-1, 2)*Symbol("n"))), Integer(0)), Integers))), "ImageSet(Lambda(x, x/2), ConditionSet(n, Eq(n^2 - n/2 - 1/2, 0), Integers))")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), ((Symbol("n")*((Integer(-1) + (Symbol("n"))^(Integer(-1))))^(Integer(-1))) + (sym_I*(Integer(1) + Symbol("n"))*(Rational(1, 2) + Symbol("n"))*(Integer(-1) + Symbol("n"))))), Integers)), "ImageSet(Lambda(n, n/(-1 + 1/n) + I*(n - 1)*(n + 1/2)*(n + 1)), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), ((Symbol("n")*((Integer(-6) + Symbol("n")))^(Integer(-1))) + (sym_I*((Integer(2) + (Integer(2)*Symbol("n"))))^(Integer(-1))*(Integer(1) + Symbol("n"))*(Integer(-3) + Symbol("n"))))), Integers)), "ImageSet(Lambda(n, n/(n - 6) + I*(n - 3)*(n + 1)/(2*n + 2)), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), ((Symbol("n")*((Integer(-9) + (Symbol("n"))^(Integer(2))))^(Integer(-1))) + (sym_I*((Integer(2) + (Integer(2)*Symbol("n"))))^(Integer(-1))*(Integer(1) + Symbol("n"))*(Integer(-3) + Symbol("n"))))), Integers)), "ImageSet(Lambda(n, n/(n^2 - 9) + I*(n - 3)*(n + 1)/(2*n + 2)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(-1, 4)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*pi*n - pi/4), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), ((Rational(-1, 4)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*pi*n + 7*pi/4), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(7, 4)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*pi*n + 7*pi/4), Integers)")
})

test_that("fancysets: test_imageset_intersect_interval matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(2)*Symbol("n"))), Interval(Integer(0), S_Pi, FALSE, FALSE))), "ImageSet(Lambda(n, 2*n), Interval(0, pi))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(1, 2)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*pi*n + pi/2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Symbol("n"))^(Integer(2))), Integers)), "ImageSet(Lambda(n, n^2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), Abs(Symbol("n"))), Integers)), "ImageSet(Lambda(n, Abs(n)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), exp(Symbol("n"))), Naturals0)), "ImageSet(Lambda(n, exp(n)), Naturals0)")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("n"), (Symbol("n"))^(Integer(2))), Integers), Interval(Integer(0), Integer(10), FALSE, FALSE))), "Intersection(ImageSet(Lambda(n, n^2), Integers), Interval(0, 10))")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("n"), Abs(Symbol("n"))), Integers), Interval(Integer(0), Integer(2), FALSE, FALSE))), "Intersection(ImageSet(Lambda(n, Abs(n)), Integers), Interval(0, 2))")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("n"), exp(Symbol("n"))), Naturals0), Interval(Integer(1), Integer(2), FALSE, FALSE))), "Intersection(ImageSet(Lambda(n, exp(n)), Naturals0), Interval(1, 2))")
})

test_that("fancysets: test_imageset_intersect_diophantine matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(1) + (Integer(2)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*n + 1), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(1) + (Integer(4)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 4*n + 1), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(2)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 2*n), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Integer(9)*(Symbol("n"))^(Integer(-1))) + (Rational(20, 3)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 20*n/3 + 9/n), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Integer(-2) + Symbol("n")))^(Integer(2))), Integers)), "ImageSet(Lambda(n, (n - 2)^2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(-1)*((Integer(-3) + Symbol("n")))^(Integer(2)))), Integers)), "ImageSet(Lambda(n, -(n - 3)^2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(5) + (Symbol("n"))^(Integer(2)))), Integers)), "ImageSet(Lambda(n, n^2 + 5), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("m"), (Integer(2)*Symbol("m"))), Integers)), "ImageSet(Lambda(m, 2*m), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(6) + (Integer(4)*Symbol("n")) + (Integer(4)*(Symbol("n"))^(Integer(2))))), Integers)), "ImageSet(Lambda(n, 4*n^2 + 4*n + 6), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(-9) + (Symbol("n"))^(Integer(2)))), Integers)), "ImageSet(Lambda(n, n^2 - 9), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("m"), (Integer(-1)*(Symbol("m"))^(Integer(2)))), Integers)), "ImageSet(Lambda(m, -m^2), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("m"), (Integer(40) + (Symbol("m"))^(Integer(2)))), Integers)), "ImageSet(Lambda(m, m^2 + 40), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(41)*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 41*n), Integers)")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("m"), (Integer(40) + (Symbol("m"))^(Integer(2)))), Integers), ImageSet(Lambda(Symbol("n"), (Integer(41)*Symbol("n"))), Integers))), "Intersection(ImageSet(Lambda(m, m^2 + 40), Integers), ImageSet(Lambda(n, 41*n), Integers))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(-16) + (Symbol("n"))^(Integer(4)))), Integers)), "ImageSet(Lambda(n, n^4 - 16), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("m"), (Integer(81) + (Integer(-1)*(Symbol("m"))^(Integer(4))))), Integers)), "ImageSet(Lambda(m, 81 - m^4), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(1, 12)*S_Pi) + (Rational(5, 12)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 5*pi*n/12 + pi/12), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(7, 12)*S_Pi) + (Rational(11, 12)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 11*pi*n/12 + 7*pi/12), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Rational(17, 4)*S_Pi) + (Rational(55, 12)*S_Pi*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 55*pi*n/12 + 17*pi/4), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Symbol("n")*log(Integer(2)))), Integers)), "ImageSet(Lambda(n, n*log(2)), Integers)")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("n"), (Symbol("n")*log(Integer(2)))), Integers), Integers)), "Intersection(ImageSet(Lambda(n, n*log(2)), Integers), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(1) + (Symbol("n"))^(Integer(3)))), Integers)), "ImageSet(Lambda(n, n^3 + 1), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Symbol("n"))^(Integer(3))), Integers)), "ImageSet(Lambda(n, n^3), Integers)")
  expect_identical(.pyrun(Intersection(ImageSet(Lambda(Symbol("n"), (Integer(1) + (Symbol("n"))^(Integer(3)))), Integers), ImageSet(Lambda(Symbol("n"), (Symbol("n"))^(Integer(3))), Integers))), "Intersection(ImageSet(Lambda(n, n^3), Integers), ImageSet(Lambda(n, n^3 + 1), Integers))")
})

test_that("fancysets: test_infinitely_indexed_set_3 matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("m"), (Integer(2)*S_Pi*Symbol("m"))), Integers)), "ImageSet(Lambda(m, 2*pi*m), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(3)*S_Pi*Symbol("n"))), Integers)), "ImageSet(Lambda(n, 3*pi*n), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("t"), (Integer(6)*S_Pi*Symbol("t"))), Integers)), "ImageSet(Lambda(t, 6*pi*t), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(-1) + (Integer(2)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 2*n + 1), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(2) + (Integer(3)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 3*n + 2), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), (Integer(-1) + (Integer(3)*Symbol("n")))), Integers)), "ImageSet(Lambda(n, 3*n + 2), Integers)")
})

test_that("fancysets: test_ImageSet_simplification matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("n"), sin(Symbol("n"))), ImageSet(Lambda(Symbol("m"), tan(Symbol("m"))), Integers))), "ImageSet(Lambda(m, sin(tan(m))), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("m"), sin(tan(Symbol("m")))), Integers)), "ImageSet(Lambda(m, sin(tan(m))), Integers)")
  expect_identical(.pyrun(imageset(Symbol("n"), (Integer(1) + (Integer(2)*Symbol("n"))), Naturals)), "Range(3, oo, 2)")
  expect_identical(.pyrun(imageset(Symbol("n"), (Integer(1) + (Integer(2)*Symbol("n"))), Naturals0)), "Range(1, oo, 2)")
  expect_identical(.pyrun(imageset(Symbol("n"), (Integer(1) + (Integer(-2)*Symbol("n"))), Naturals)), "Range(-1, -oo, -2)")
})

test_that("fancysets: test_ImageSet_contains matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("x") + (Integer(3)*sym_I)), Integers)), "ImageSet(Lambda(x, x + 3*I), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("x") + (sym_I*Symbol("y"))), Integers)), "ImageSet(Lambda(x, x + I*y), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("x") + (sym_I*Symbol("y")*(Symbol("x"))^(Integer(-1)))), Integers)), "ImageSet(Lambda(x, x + I*y/x), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), (Symbol("x") + (sym_I*(Integer(-1) + (cos(Integer(1)))^(Integer(2)) + (sin(Integer(1)))^(Integer(2))))), Integers)), "ImageSet(Lambda(x, x + I*(-1 + cos(1)^2 + sin(1)^2)), Integers)")
})

test_that("fancysets: test_issue_16871 matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Integer(-3) + Symbol("x"))), Integers)), "ImageSet(Lambda(x, x - 3), Integers)")
})

test_that("fancysets: test_issue_18050 matches SymPy", {
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Integer(1) + (sym_I*Symbol("x")))), Integers)), "ImageSet(Lambda(x, I*x + 1), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Integer(1) + (sym_I*Symbol("x")))), Integers)), "ImageSet(Lambda(x, I*x + 1), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Integer(4) + (Integer(8)*sym_I) + (Integer(3)*sym_I*Symbol("x")))), Integers)), "ImageSet(Lambda(x, 3*I*x + 4 + 2*I), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Integer(4) + (Integer(2)*sym_I) + (Integer(3)*sym_I*Symbol("x")))), Integers)), "ImageSet(Lambda(x, 3*I*x + 4 + 2*I), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), ((Integer(2)*Symbol("x")) + (Integer(3)*sym_I))), Integers)), "ImageSet(Lambda(x, 2*x + 3*I), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), ((Integer(2)*Symbol("x")) + (Integer(3)*sym_I))), Integers)), "ImageSet(Lambda(x, 2*x + 3*I), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Integer(10) + (Symbol("x")*Symbol("r", positive = TRUE)))), Integers)), "ImageSet(Lambda(x, r*x + 10), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Integer(10) + (Symbol("x")*Symbol("r", positive = TRUE)))), Integers)), "ImageSet(Lambda(x, r*x + 10), Integers)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), (Integer(8) + (Integer(3)*Symbol("x")) + (Integer(5)*sym_I))), Integers)), "ImageSet(Lambda(x, 3*x + 2 + 5*I), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), (Integer(2) + (Integer(3)*Symbol("x")) + (Integer(5)*sym_I))), Integers)), "ImageSet(Lambda(x, 3*x + 2 + 5*I), Integers)")
})

test_that("fancysets: test_imageset_intersection matches SymPy", {
  expect_identical(.pyrun(local({.D1 <- Dummy("Dummy_4120"); ImageSet(Lambda(.D1, ((Rational(-1, 4)*S_Pi) + (Integer(2)*S_Pi*.D1))), Integers)})), "ImageSet(Lambda(_Dummy_4120, 2*_Dummy_4120*pi - pi/4), Integers)")
  expect_identical(.pyrun(local({.D1 <- Dummy("Dummy_4120"); ImageSet(Lambda(.D1, ((Rational(7, 4)*S_Pi) + (Integer(2)*S_Pi*.D1))), Integers)})), "ImageSet(Lambda(_Dummy_4120, 2*_Dummy_4120*pi + 7*pi/4), Integers)")
})

test_that("sets: test_imageset matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), (Integer(-1) + Symbol("x")), Naturals)), "Naturals0")
  expect_identical(.pyrun(imageset(Symbol("x"), (Integer(1) + Symbol("x")), Naturals0)), "Naturals")
  expect_identical(.pyrun(imageset(Symbol("x"), Abs(Symbol("x")), Naturals0)), "Naturals0")
  expect_identical(.pyrun(imageset(Symbol("x"), Abs(Symbol("x")), Naturals)), "Naturals")
  expect_identical(.pyrun(imageset(Symbol("x"), Abs(Symbol("x")), Integers)), "Naturals0")
  expect_identical(.pyrun(imageset(Symbol("x"), list(Symbol("x"), Symbol("x")), Reals)), "ImageSet(Lambda(x, (x, x)), Reals)")
  expect_identical(.pyrun(imageset(Symbol("x"), (sym_I + Symbol("x")), Reals)), "ImageSet(Lambda(x, x + I), Reals)")
  expect_identical(.pyrun(imageset(Symbol("x"), list(Symbol("x")), Reals)), "ImageSet(Lambda(x, (x,)), Reals)")
  expect_identical(.pyrun(imageset(Symbol("x"), Integers)), "ERR TypeError: \\nexpecting lambda, Lambda, or FunctionClass, not 'Symbol'.")
  expect_identical(.pyrun(imageset(Symbol("x"), Symbol("y"), Symbol("z"), Integers)), "ERR ValueError: arguments after mapping should be sets, not ['Symbol', 'Integers']")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x"), cos(Symbol("x"))), Symbol("y"))), "ERR ValueError: arguments after mapping should be sets, not ['Symbol']")
  expect_identical(.pyrun(imageset(Lambda(list(Symbol("x"), Symbol("y")), Tuple(Symbol("x"), Symbol("y"))), Integers, Integers)), "ImageSet(Lambda((x, y), (x, y)), Integers, Integers)")
  expect_identical(.pyrun(imageset("cos", Integers)), "ImageSet(Lambda(x, cos(x)), Integers)")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), cos(Symbol("x"))), Integers)), "ImageSet(Lambda(x, cos(x)), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), cos(Symbol("x")), Integers)), "ImageSet(Lambda(x, cos(x)), Integers)")
  expect_identical(.pyrun(imageset(Symbol("x"), Integer(1), Integers)), "{1}")
  expect_identical(.pyrun(imageset(Symbol("x"), Symbol("y"), Integers)), "{y}")
  expect_identical(.pyrun(imageset(list(Symbol("x"), Symbol("y")), list(Integer(1), Symbol("z")), Integers, Reals)), "{(1, z)}")
  expect_identical(.pyrun(ImageSet(Lambda(list(Symbol("x1"), Symbol("x2")), (Symbol("x1") + Symbol("x2"))), Interval(Integer(1), Integer(2), FALSE, FALSE), Interval(Integer(2), Integer(3), FALSE, FALSE))), "ImageSet(Lambda((x1, x2), x1 + x2), Interval(1, 2), Interval(2, 3))")
})

test_that("sets: test_is_subset matches SymPy", {
  expect_identical(.pyrun(Union(Interval(Integer(0), Integer(1), FALSE, FALSE), Interval(Integer(2), Integer(5), FALSE, FALSE))), "Union(Interval(0, 1), Interval(2, 5))")
  expect_identical(.pyrun(Union(Interval((-oo), Integer(-2), TRUE, FALSE), Interval(Integer(3), oo, FALSE, TRUE))), "Union(Interval(-oo, -2), Interval(3, oo))")
  expect_identical(.pyrun(imageset(Lambda(Symbol("n", integer = TRUE), (Symbol("n", integer = TRUE))^(Integer(-1))), Integers)), "ImageSet(Lambda(n, 1/n), Integers)")
})

test_that("sets: test_image_interval matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Integer(2)*Symbol("x", real = TRUE)), Interval(Integer(-2), Integer(1), FALSE, FALSE))), "Interval(-4, 2)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Integer(2)*Symbol("x", real = TRUE)), Interval(Integer(-2), Integer(1), TRUE, FALSE))), "Interval.Lopen(-4, 2)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)), Interval(Integer(-2), Integer(1), TRUE, FALSE))), "Interval.Ropen(0, 4)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)), Interval(Integer(-2), Integer(1), FALSE, FALSE))), "Interval(0, 4)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)), Interval(Integer(-2), Integer(1), TRUE, TRUE))), "Interval.Ropen(0, 4)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), ((Integer(-2) + Symbol("x", real = TRUE)))^(Integer(2)), Interval(Integer(1), Integer(3), FALSE, FALSE))), "Interval(0, 1)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), ((Integer(-90)*Symbol("x", real = TRUE)) + (Integer(-26)*(Symbol("x", real = TRUE))^(Integer(3))) + (Integer(3)*(Symbol("x", real = TRUE))^(Integer(4))) + (Integer(78)*(Symbol("x", real = TRUE))^(Integer(2)))), Interval(Integer(0), Integer(4), FALSE, FALSE))), "Interval(-35, 0)")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE) + (Symbol("x", real = TRUE))^(Integer(-1))), Interval((-oo), oo, TRUE, TRUE))), "Union(Interval(-oo, -2), Interval(2, oo))")
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), ((Symbol("x", real = TRUE))^(Integer(-1)) + ((Integer(-1) + Symbol("x", real = TRUE)))^(Integer(-2))), Interval(Integer(0), Integer(2), TRUE, FALSE))), "Interval(3/2, oo)")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x", real = TRUE), (Symbol("a", real = TRUE)*Symbol("x", real = TRUE))), Interval(Integer(0), Integer(1), FALSE, FALSE))), "ImageSet(Lambda(x, a*x), Interval(0, 1))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x", real = TRUE), (Symbol("a", real = TRUE)*Symbol("x", real = TRUE))), Interval(Integer(0), Integer(1), FALSE, FALSE))), "ImageSet(Lambda(x, a*x), Interval(0, 1))")
  expect_identical(.pyrun(imageset(Lambda(Symbol("x", real = TRUE), sin(cos(Symbol("x", real = TRUE)))), Interval(Integer(0), Integer(1), FALSE, FALSE))), "ImageSet(Lambda(x, sin(cos(x))), Interval(0, 1))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x", real = TRUE), sin(cos(Symbol("x", real = TRUE)))), Interval(Integer(0), Integer(1), FALSE, FALSE))), "ImageSet(Lambda(x, sin(cos(x))), Interval(0, 1))")
})

test_that("sets: test_image_piecewise matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), Piecewise(list(Symbol("x"), Le(Symbol("x"), Integer(-1))), list((Symbol("x"))^(Integer(-2)), Le(Symbol("x"), Integer(5))), list((Symbol("x"))^(Integer(3)), S_true)), Interval(Integer(-5), Integer(5), FALSE, FALSE))), "Union(Interval(-5, -1), Interval(1/25, oo))")
  expect_identical(.pyrun(Union(Interval(Integer(-5), Integer(-1), FALSE, FALSE), Interval(Rational(1, 25), oo, FALSE, TRUE))), "Union(Interval(-5, -1), Interval(1/25, oo))")
  expect_identical(.pyrun(imageset(Symbol("x"), Piecewise(list(Integer(0), Le(Symbol("x"), Integer(1))), list(Integer(1), Le(Symbol("x"), Integer(2))), list(Integer(2), S_true)), Interval(Integer(1), Integer(2), FALSE, FALSE))), "{0, 1}")
})

test_that("sets: test_image_Intersection matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)), Intersection(Interval(Integer(-2), Integer(0), FALSE, FALSE), Interval(Symbol("x", real = TRUE), Symbol("y", real = TRUE), FALSE, FALSE)))), "ImageSet(Lambda(x, x^2), Intersection(Interval(-2, 0), Interval(x, y)))")
})

test_that("sets: test_image_FiniteSet matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Integer(2)*Symbol("x", real = TRUE)), FiniteSet(Integer(1), Integer(2), Integer(3)))), "{2, 4, 6}")
})

test_that("sets: test_image_Union matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Symbol("x", real = TRUE))^(Integer(2)), Union(Interval(Integer(-2), Integer(0), FALSE, FALSE), FiniteSet(Integer(1), Integer(2), Integer(3))))), "Union({9}, Interval(0, 4))")
})

test_that("sets: test_image_EmptySet matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x", real = TRUE), (Integer(2)*Symbol("x", real = TRUE)), S_EmptySet)), "EmptySet")
})

test_that("sets: test_issue_10113 matches SymPy", {
  expect_identical(.pyrun(imageset(Symbol("x"), ((Symbol("x"))^(Integer(2))*((Integer(-4) + (Symbol("x"))^(Integer(2))))^(Integer(-1))), Reals)), "Union(Interval(-oo, 0), Interval.open(1, oo))")
  expect_identical(.pyrun(Union(Interval((-oo), Integer(0), TRUE, FALSE), Interval(Integer(1), oo, TRUE, TRUE))), "Union(Interval(-oo, 0), Interval.open(1, oo))")
  expect_identical(.pyrun(imageset(Symbol("x"), ((Symbol("x"))^(Integer(2))*((Integer(-4) + (Symbol("x"))^(Integer(2))))^(Integer(-1))), Interval(Integer(-2), Integer(2), FALSE, FALSE))), "Interval(-oo, 0)")
  expect_identical(.pyrun(imageset(Symbol("x"), ((Symbol("x"))^(Integer(2))*((Integer(-4) + (Symbol("x"))^(Integer(2))))^(Integer(-1))), Interval(Integer(-2), Integer(3), FALSE, FALSE))), "Union(Interval(-oo, 0), Interval(9/5, oo))")
  expect_identical(.pyrun(Union(Interval((-oo), Integer(0), TRUE, FALSE), Interval(Rational(9, 5), oo, FALSE, TRUE))), "Union(Interval(-oo, 0), Interval(9/5, oo))")
})

test_that("sets: test_issue_18505 matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), ((Integer(-1) + (Rational(1, 2)*S_Pi) + (Rational(1, 2)*S_Pi*Symbol("n"))))^(Rational(1, 2))), Integers)), "ImageSet(Lambda(n, sqrt(pi*n/2 - 1 + pi/2)), Integers)")
})

test_that("sets: test_SetKind_ImageSet_Special matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Symbol("n"))^(Integer(2))), Interval(Integer(1), Integer(4), FALSE, FALSE))), "ImageSet(Lambda(n, n^2), Interval(1, 4))")
})

test_that("sets: test_16878 matches SymPy", {
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("x"), Tuple(Symbol("x"), Symbol("x"))), Reals)), "ImageSet(Lambda(x, (x, x)), Reals)")
})

test_that("sets: test_issue_24726 matches SymPy", {
  expect_identical(.pyrun(solveset(sin(Symbol("t")), Symbol("t"), domain = Reals)), "Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))")
  expect_identical(.pyrun(ImageSet(Lambda(Symbol("n"), (Integer(-9) + (Integer(2)*Symbol("n")))), Range(Integer(0), oo, Integer(1)))), "ImageSet(Lambda(n, 2*n - 9), Range(0, oo, 1))")
  expect_identical(.pyrun(Union(ImageSet(Lambda(Symbol("n"), ((Integer(-3)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))), ImageSet(Lambda(Symbol("n"), ((Integer(-2)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))))), "Union(ImageSet(Lambda(n, 2*pi*n - 3*pi), Range(0, oo, 1)), ImageSet(Lambda(n, 2*pi*n - 2*pi), Range(0, oo, 1)))")
  expect_identical(.pyrun(Union(ImageSet(Lambda(Symbol("n"), ((Integer(2)*S_Pi) + (Integer(-2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))), ImageSet(Lambda(Symbol("n"), ((Integer(3)*S_Pi) + (Integer(-2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))))), "Union(ImageSet(Lambda(n, -2*pi*n + 2*pi), Range(0, oo, 1)), ImageSet(Lambda(n, -2*pi*n + 3*pi), Range(0, oo, 1)))")
  expect_identical(.pyrun(Union(ImageSet(Lambda(Symbol("n"), ((Integer(4)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))), ImageSet(Lambda(Symbol("n"), ((Integer(5)*S_Pi) + (Integer(2)*S_Pi*Symbol("n")))), Range(Integer(0), oo, Integer(1))))), "Union(ImageSet(Lambda(n, 2*pi*n + 4*pi), Range(0, oo, 1)), ImageSet(Lambda(n, 2*pi*n + 5*pi), Range(0, oo, 1)))")
})

test_that("docs: mismatches-docs matches SymPy", {
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); solveset(cos(x), x, Reals)})), "Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers))")
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); linsolve(list(x + y + z - 1, x + y + 2*z - 3), list(x, y, z))})), "{(-y - 1, y, 2)}")
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); solve(list(x^2 + y^2 - 1, x - y), list(x, y))})), "[(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)]")
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); nonlinsolve(list(x^2 + y^2 - 1, x - y), list(x, y))})), "{(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)}")
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); solve_poly_system(list(x^2 + y^2 - 1, x - y), x, y)})), "[(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)]")
})

test_that("docs: singular_points matches SymPy", {
  expect_identical(.pyrun(local({x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); singular_points(ImplicitRegion(list(x, y, z), x*y*z + y^4 - x^2*z^2))})), "{(0, 0, z), (x, 0, 0)}")
})
