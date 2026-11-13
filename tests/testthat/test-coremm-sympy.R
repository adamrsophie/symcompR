skip_unless_full_tests()

# Results checked against SymPy (solveset and ODE rows fixed by the core changes)

coremm_run <- function(code) gsub("**", "^", .pyfmt(eval(parse(text = code))), fixed = TRUE)

test_that("solveset matches SymPy: test_solve_trig:963", {
  expect_equal(coremm_run("solveset_real((Integer(1) + (Integer(2)*sin(Symbol(\"x\"))*tan(Symbol(\"x\")))), Symbol(\"x\"))"),
               "Union(ImageSet(Lambda(_n, 2*_n*pi - atan(sqrt(2)*sqrt(-1 + sqrt(17))/(1 - sqrt(17))) + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + atan(sqrt(2)*sqrt(-1 + sqrt(17))/(1 - sqrt(17))) + pi), Integers))")
})

test_that("solveset matches SymPy: test_solve_trig:1010", {
  expect_equal(coremm_run("solveset_real((Integer(-1) + (Integer(2)*cos(Symbol(\"x\"))*cos((Integer(2)*Symbol(\"x\"))))), Symbol(\"x\"))"),
               "Union(ImageSet(Lambda(_n, 2*_n*pi - 2*atan(sqrt(-2*2^(1/3)*(67 + 9*sqrt(57))^(2/3) + 8*2^(2/3) + 11*(67 + 9*sqrt(57))^(1/3))/(3*(67 + 9*sqrt(57))^(1/6))) + 2*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + 2*atan(sqrt(-2*2^(1/3)*(67 + 9*sqrt(57))^(2/3) + 8*2^(2/3) + 11*(67 + 9*sqrt(57))^(1/3))/(3*(67 + 9*sqrt(57))^(1/6)))), Integers))")
})

test_that("solveset matches SymPy: test_old_trig_issues:1150", {
  expect_equal(coremm_run("solveset((Float(\"-0.500000000000000\", precision = 53) + cos((Integer(2)*Symbol(\"x\")))), Symbol(\"x\"), Interval(Integer(0), (Integer(2)*S_Pi), FALSE, FALSE))"),
               "{0.523598775598299, -0.523598775598299 + 1.0*pi, -0.523598775598299 + 2.0*pi, 0.523598775598299 + 1.0*pi}")
})

test_that("solveset matches SymPy: test_issue_9616:1284", {
  expect_equal(coremm_run("solveset((Integer(-1) + (Rational(1, 2)*exp(Symbol(\"x\"))) + (Rational(-1, 2)*exp((Integer(-1)*Symbol(\"x\")))) + (((exp(Symbol(\"x\")) + exp((Integer(-1)*Symbol(\"x\")))))^(Integer(-1))*((Integer(-1)*exp((Integer(-1)*Symbol(\"x\")))) + exp(Symbol(\"x\"))))), Symbol(\"x\"))"),
               "Union(Complement(ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(-sqrt(2)/2 + sqrt(-1/2 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi)/2), Integers)), Complement(ImageSet(Lambda(_n, I*(2*_n*pi - pi + atan(sqrt(2)*sqrt(1/2 + sqrt(2)))) + log(sqrt(1 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi)/2), Integers)), Complement(ImageSet(Lambda(_n, I*(2*_n*pi - atan(sqrt(2)*sqrt(1/2 + sqrt(2))) + pi) + log(sqrt(1 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi)/2), Integers)), Complement(ImageSet(Lambda(_n, 2*_n*I*pi + log(sqrt(2)/2 + sqrt(-1/2 + sqrt(2)))), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + pi)/2), Integers)))")
})

test_that("solveset matches SymPy: test_issue_17566:2419", {
  expect_equal(coremm_run("nonlinsolve(list(((Integer(-1)*(Integer(4))^(Symbol(\"y\"))) + (Integer(32)*(Integer(2))^(Symbol(\"x\"))*(Integer(2))^(Symbol(\"y\")))), ((Integer(-1)*(Integer(3))^((Integer(-1)*Symbol(\"y\")))) + (Integer(27)*(Integer(3))^(Symbol(\"x\"))))), Symbol(\"x\"), Symbol(\"y\"))"),
               "{(-log(81)/log(3), 1)}")
})

test_that("solveset matches SymPy: test_issue_19587:2431", {
  expect_equal(coremm_run("nonlinsolve(list(((Integer(-1)*(Integer(4))^(Symbol(\"n\"))) + (Integer(32)*(Integer(2))^(Symbol(\"m\"))*(Integer(2))^(Symbol(\"n\")))), ((Integer(-1)*(Integer(3))^((Integer(-1)*Symbol(\"n\")))) + (Integer(27)*(Integer(3))^(Symbol(\"m\"))))), Symbol(\"m\"), Symbol(\"n\"))"),
               "{(-log(81)/log(3), 1)}")
})

test_that("dsolve matches SymPy: test_ode_solutions:157", {
  expect_equal(coremm_run("constantsimp(Eq(((Rational(1, 2)*(Symbol(\"x\"))^(Integer(2))*(Function(\"f\")(Symbol(\"x\")))^(Integer(-2))) + log((Symbol(\"x\")*(Integer(2))^(Rational(1, 2))*(Symbol(\"C1\"))^(Integer(-1))*((Symbol(\"x\"))^(Integer(-1)))^(Rational(1, 2))*(Function(\"f\")(Symbol(\"x\")))^(Rational(1, 2))))), Integer(0)), list(Symbol(\"C1\")))"),
               "Eq((x^2/f(x)^2)/2 + log(C1*sqrt(x)*sqrt(f(x))), 0)")
})

test_that("dsolve matches SymPy: test_solve_ics:499", {
  expect_equal(coremm_run("dsolve((Integer(-1) + (Integer(2)*Symbol(\"x\", real = TRUE)*(Function(\"f\")(Symbol(\"x\", real = TRUE)))^(Integer(2))) + (Integer(2)*(Symbol(\"x\", real = TRUE))^(Integer(2))*Derivative(Function(\"f\")(Symbol(\"x\", real = TRUE)), Symbol(\"x\", real = TRUE))*Function(\"f\")(Symbol(\"x\", real = TRUE)))), Function(\"f\")(Symbol(\"x\", real = TRUE)), ics = list(list(Function(\"f\")(Integer(1)), Integer(0))))"),
               "[Eq(f(x), -sqrt(x - 1)/x), Eq(f(x), sqrt(x - 1)/x)]")
})

test_that("dsolve matches SymPy: test_issue_29827:1180", {
  expect_equal(coremm_run("dsolve(Eq(Derivative(Function(\"f\")(Symbol(\"x\", real = TRUE)), Symbol(\"x\", real = TRUE)), log(Function(\"f\")(Symbol(\"x\", real = TRUE)))), func = Function(\"f\")(Symbol(\"x\", real = TRUE)), ics = list(list(Function(\"f\")(Integer(1)), Integer(1))))"),
               "[Eq(li(f(x)), -oo), Eq(f(x), 1)]")
})

test_that("dsolve matches SymPy: test_single.1st_rational_riccati:rational_riccati_09", {
  expect_equal(coremm_run("dsolve(((Integer(-2)*sym_I*(Symbol(\"x\"))^(Integer(-1))*(Integer(1) + (Function(\"f\")(Symbol(\"x\")))^(Integer(2)))) + Derivative(Function(\"f\")(Symbol(\"x\")), Symbol(\"x\"))), Function(\"f\")(Symbol(\"x\")), hint = \"1st_rational_riccati\")"),
               "Eq(f(x), (-I*C1 + I*x^4)/(C1 + x^4))")
})

test_that("dsolve matches SymPy: test_single.separable:separable_24", {
  expect_equal(coremm_run("dsolve((Integer(-1) + (Float(\"51.0500000000000\", precision = 53)*Symbol(\"y\")*Function(\"f\")(Symbol(\"a\"))) + Derivative(Function(\"f\")(Symbol(\"a\")), Symbol(\"a\"))), Function(\"f\")(Symbol(\"a\")), hint = \"separable\")"),
               "Eq(f(a), (0.019588638589618*exp(y*(C1 - 51.05*a)) + 0.019588638589618)/y)")
})

test_that("dsolve matches SymPy: test_single.separable:separable_26", {
  expect_equal(coremm_run("dsolve((Float(\"940.800000000000\", precision = 53) + (Integer(-96)*Derivative(Function(\"v\")(Symbol(\"a\")), Symbol(\"a\"))) + (Float(\"-0.200000000000000\", precision = 53)*(Function(\"v\")(Symbol(\"a\")))^(Integer(2)))), Function(\"v\")(Symbol(\"a\")), hint = \"separable\")"),
               "Eq(v(a), -68.585712797929/tanh(C1 - 0.142886901662352*a))")
})
