# Reference values from SymPy: Mul._eval_power, Integer and Rational
# _eval_power for negative bases, NegativeOne._eval_power, inverse
# trigonometric evaluation, and limits on branch cuts.
syms <- function() list(x = Symbol("x", real = TRUE), y = Symbol("y"), u = Symbol("u", positive = TRUE),
                        v = Symbol("v", positive = TRUE), w = Symbol("w", positive = TRUE), n = Symbol("n", integer = TRUE))

test_that("powers of products split only as SymPy splits them", {
  L <- syms()
  cases <- list(
    c("sqrt(x**2*y)", "sqrt(y)*Abs(x)"),
    c("sqrt(4*y)", "2*sqrt(y)"),
    c("(-2*u)**(1/3)", "(-2)^(1/3)*u^(1/3)"),
    c("(-8*u)**(1/3)", "2*(-1)^(1/3)*u^(1/3)"),
    c("sqrt(-u*v)", "I*sqrt(u)*sqrt(v)"),
    c("sqrt(-u*y)", "sqrt(u)*sqrt(-y)"),
    c("sqrt(-2*u*y)", "sqrt(2)*sqrt(u)*sqrt(-y)"),
    c("(I*u)**(1/2)", "sqrt(I)*sqrt(u)"),
    c("(2*I*u)**(1/2)", "sqrt(2)*sqrt(I)*sqrt(u)"),
    c("(-u*w)**(1/3)", "(-1)^(1/3)*u^(1/3)*w^(1/3)"),
    c("(u*y)**(2/3)", "u^(2/3)*y^(2/3)"),
    c("(u*y)**x", "(u*y)^x"),
    c("(u*v)**n", "(u*v)^n"),
    c("(-u)**u", "(-u)^u"),
    c("(2*u)**pi", "(2*u)^pi"),
    c("sqrt(-x*y)", "sqrt(-x*y)"),
    c("(3*u)**0.5", "1.73205080756888*u^0.5"),
    c("sqrt(2*I)", "1 + I"),
    c("sqrt(-2*I)", "1 - I"),
    c("(u*sqrt(2))**(1/2)", "2^(1/4)*sqrt(u)"),
    c("(2**(1/3)*u)**(3/2)", "sqrt(2)*u^(3/2)"),
    c("sqrt(8*u*y)", "2*sqrt(2)*sqrt(u)*sqrt(y)"),
    c("(x*y)**(1/2)", "sqrt(x*y)"))
  for (cs in cases) expect_equal(sstr(sympify(cs[1], locals = L)), cs[2], label = cs[1])
})

test_that("negative numbers to rational powers follow SymPy", {
  cases <- list(
    c("(-2)**(1/3)", "(-2)^(1/3)"),
    c("(-8)**(1/3)", "2*(-1)^(1/3)"),
    c("(-12)**(1/3)", "(-12)^(1/3)"),
    c("(-16)**(1/3)", "2*(-2)^(1/3)"),
    c("(-2)**(-1/3)", "-(-2)^(2/3)/2"),
    c("(-2)**(2/3)", "(-2)^(2/3)"),
    c("(-4)**(1/4)", "(-1)^(1/4)*sqrt(2)"),
    c("(-1/2)**(1/3)", "(-1)^(1/3)*2^(2/3)/2"),
    c("(-8/27)**(1/3)", "2*(-1)^(1/3)/3"),
    c("(-3/4)**(1/2)", "sqrt(3)*I/2"),
    c("(-2/3)**(2/3)", "(-2)^(2/3)*3^(1/3)/3"),
    c("(-2)**(5/3)", "-2*(-2)^(2/3)"),
    c("(-27)**(-2/3)", "-(-1)^(1/3)/9"),
    c("(-1/8)**(-1/3)", "-2*(-1)^(2/3)"),
    c("(-12)**(1/2)", "2*sqrt(3)*I"))
  for (cs in cases) expect_equal(sstr(sympify(cs[1])), cs[2], label = cs[1])
})

test_that("inverse trigonometric functions evaluate as in SymPy", {
  x <- Symbol("x"); I <- sym_I
  got <- vapply(list(acos(Integer(-2)), acos(-x), acot(I/2), acot(-I/2), asec(Integer(-2)), acsc(Integer(-2)),
                     asec(Integer(2)), acsc(sqrt(Integer(2))), acos(-Rational(1, 2))), sstr, "")
  expect_equal(got, c("acos(-2)", "acos(-x)", "-I*acoth(1/2)", "I*acoth(1/2)", "2*pi/3", "-pi/6", "pi/3", "pi/4", "2*pi/3"))
})

test_that("limits on branch cuts take SymPy's one sided values", {
  x <- Symbol("x")
  cases <- list(
    c("asin(I*x + 2)", "0", "+", "pi - asin(2)"),
    c("asin(I*x - 2)", "0", "-", "-pi + asin(2)"),
    c("acos(I*x + 2)", "0", "+", "-acos(2)"),
    c("acos(I*x - 2)", "0", "+", "acos(-2)"),
    c("acos(I*x - 2)", "0", "-", "2*pi - acos(-2)"),
    c("atan(x + 2*I)", "0", "-", "-pi + I*atanh(2)"),
    c("atan(x - 2*I)", "0", "+", "pi - I*atanh(2)"),
    c("acot(x + I/2)", "0", "+", "pi - I*acoth(1/2)"),
    c("acot(x + I/2)", "0", "-", "-I*acoth(1/2)"),
    c("acot(x - I/2)", "0", "+", "I*acoth(1/2)"),
    c("acot(x - I/2)", "0", "-", "-pi + I*acoth(1/2)"),
    c("asec(I*x + 1/2)", "0", "-", "-asec(1/2)"),
    c("asec(I*x - 1/2)", "0", "+", "2*pi - asec(-1/2)"),
    c("acsc(I*x + 1/2)", "0", "-", "pi - acsc(1/2)"),
    c("acsc(I*x - 1/2)", "0", "+", "-pi + acsc(1/2)"))
  for (cs in cases) expect_equal(sstr(limit(sympify(cs[1], locals = list(x = x)), x, sympify(cs[2]), cs[3])), cs[4], label = cs[1])
  expect_equal(sstr(limit((-x)^x, x, oo)), "zoo")
  expect_equal(sstr(limit((1 - x/2)^(3*x), x, oo)), "zoo")
  expect_equal(sstr(limit((-Rational(1, 4))^x, x, oo)), "0")
})

test_that("powers of Gaussian rationals keep SymPy's Mul forms", {
  x <- Symbol("x")
  got <- vapply(list(sqrt((2 - sym_I)^(-2)), sqrt((1 + sym_I)^(-2)), sqrt((-2 - sym_I)^(-2)),
                     -1 * ((-2 + sym_I) / 5), 5 * ((-2 + sym_I) / 5), 1 / (3 - 4 * sym_I) * 2,
                     2 * ((x + 1) / 3), 3 * (2 * (x + 1)), sqrt(2 * x * (2 - sym_I)^(-2))), sstr, "")
  expect_identical(got, c("(2 + I)/5", "(1 - I)/2", "-(-2 + I)/5", "2/5 - I/5", "-2 + I", "2*(3 + 4*I)/25",
                          "2*x/3 + 2/3", "6*x + 6", "sqrt(2)*sqrt(x/(2 - I)^2)"))
})

test_that("cancel reads I back as a Gaussian coefficient", {
  x <- Symbol("x")
  d <- cancel(-x * exp(-x^2 / (2 * (2 - sym_I)^2)) / (2 - sym_I)^2)
  expect_identical(sstr(d), "-I*x*(4 - 3*I)*exp(-3*x^2/50)*exp(-2*I*x^2/25)/25")
  y <- Symbol("y"); z <- Symbol("z")
  got <- vapply(list(x / ((4 * sym_I - 3) * y * z), (6 * x + 3 * sym_I * y) / (4 * z * sym_I + 2 * y),
                     (2 * x + 2 * sym_I * x) / (x * y)), function(e) sstr(cancel(e)), "")
  expect_identical(got, c("-I*x*(4 - 3*I)/(25*y*z)", "(6*x + 3*I*y)/(2*y + 4*I*z)", "(2 + 2*I)/y"))
  expect_null(heurisch(exp(-x^2 / (2 * (2 - sym_I)^2)), x, hints = list()))
})

test_that("heurisch visits Bessel functions of higher order first", {
  x <- Symbol("x"); nu <- Symbol("nu")
  expect_identical(sstr(heurisch((nu * besselj(nu, x) - x * besselj(nu + 1, x)) / x, x)), "besselj(nu, x)")
  expect_null(heurisch(besselj(nu, x) * x, x))
})

test_that("a rebuilt power stays a nested Mul as in Mul.flatten", {
  got <- vapply(list(sqrt(1 / (2 * (2 - sym_I)^2)), sqrt(3 / (2 - sym_I)^2),
                     (1 / (8 * (1 + sym_I)^2))^Rational(1, 3)), sstr, "")
  expect_identical(got, c("sqrt(2)*((2 + I)/5)/2", "sqrt(3)*((2 + I)/5)", "((1 + I)^(-2))^(1/3)/2"))
})

test_that("Mul rebuilds an unevaluated number power as Mul.flatten does", {
  a <- new_sym("Pow", list(S_Two, Rational(-1, 2)))
  expect_identical(sstr(a), "1/sqrt(2)")
  expect_identical(.k(a * -1), "Mul(Q-1/2,Pow(Z2,Q1/2))")
  expect_identical(.k(a * sqrt(as_sym(3))), "Mul(Q1/2,Pow(Z6,Q1/2))")
  m <- represent(HadamardGate(0), nqubits = 1)
  expect_identical(vapply(m$data, sstr, ""), c("1/sqrt(2)", "1/sqrt(2)", "1/sqrt(2)", "-sqrt(2)/2"))
  s <- Symbol("s")
  expect_identical(sstr(phase_margin(TransferFunction(2, s^3 + 2 * s^2 + 2 * s + 1, s))),
                   "180*(-pi + atan((-2*3^(1/6) + sqrt(3))/(1 - 2*3^(1/3))))/pi + 180")
})
