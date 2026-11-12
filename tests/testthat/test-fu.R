# Expected values are the outputs in the docstrings of sympy/simplify/fu.py
# (powers written with ^).
x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); z <- Symbol("z")
pi <- sym_pi
ta <- tan(a); tb <- tan(b); tc <- tan(c)
I <- function(n) Integer(n)

test_that("the Fu transformation rules and fu() follow SymPy", {
  cases <- list(
   c('TR1(2*csc(x) + sec(x))', '1/cos(x) + 2/sin(x)'), c('TR2(tan(x))', 'sin(x)/cos(x)'), c('TR2(cot(x))', 'cos(x)/sin(x)'),
   c('TR2(tan(tan(x) - sin(x)/cos(x)))', '0'), c('TR2i(sin(x)/cos(x))', 'tan(x)'),
   c('TR2i(sin(x)^2/(cos(x) + 1)^2, half=TRUE)', 'tan(x/2)^2'), c('TR2i(sin(x)^a/(cos(x) + 1)^a)', 'sin(x)^a/(cos(x) + 1)^a'),
   c('TR3(cos(y - x*(y - x)))', 'cos(x*(x - y) + y)'),
   c('TR5(sin(x)^2)', '1 - cos(x)^2'), c('TR5(sin(x)^-2)', 'sin(x)^(-2)'), c('TR5(sin(x)^4)', '(1 - cos(x)^2)^2'),
   c('TR6(cos(x)^2)', '1 - sin(x)^2'), c('TR6(cos(x)^4)', '(1 - sin(x)^2)^2'),
   c('TR5(sin(x)^3)', '(1 - cos(x)^2)*sin(x)'), c('TR5(sin(x)^6, max=6)', '(1 - cos(x)^2)^3'), c('TR5(sin(x)^6, max=6, pow=TRUE)', 'sin(x)^6'), c('TR5(sin(x)^8, max=10, pow=TRUE)', '(1 - cos(x)^2)^4'),
   c('TR7(cos(x)^2)', 'cos(2*x)/2 + 1/2'), c('TR7(cos(x)^2 + 1)', 'cos(2*x)/2 + 3/2'),
   c('TR8(cos(I(2))*cos(I(3)))', 'cos(5)/2 + cos(1)/2'), c('TR8(cos(I(2))*sin(I(3)))', 'sin(5)/2 + sin(1)/2'), c('TR8(sin(I(2))*sin(I(3)))', '-cos(5)/2 + cos(1)/2'),
   c('TR9(cos(I(1)) + cos(I(2)))', '2*cos(1/2)*cos(3/2)'), c('TR9(cos(I(1)) + 2*sin(I(1)) + 2*sin(I(2)))', 'cos(1) + 4*sin(3/2)*cos(1/2)'),
   c('TR9(cos(I(3)) + cos(I(3))*cos(I(2)))', 'cos(3) + cos(2)*cos(3)'),
   c('TR10(cos(a + b))', '-sin(a)*sin(b) + cos(a)*cos(b)'), c('TR10(sin(a + b))', 'sin(a)*cos(b) + sin(b)*cos(a)'),
   c('TR10(sin(a + b + c))', '(-sin(a)*sin(b) + cos(a)*cos(b))*sin(c) + (sin(a)*cos(b) + sin(b)*cos(a))*cos(c)'),
   c('TR10i(cos(I(1))*cos(I(3)) + sin(I(1))*sin(I(3)))', 'cos(2)'), c('TR10i(cos(I(1))*sin(I(3)) + sin(I(1))*cos(I(3)) + cos(I(3)))', 'cos(3) + sin(4)'),
   c('TR10i(sqrt(I(2))*cos(x)*x + sqrt(I(6))*sin(x)*x)', '2*sqrt(2)*x*sin(x + pi/6)'),
   c('TR11(sin(2*x))', '2*sin(x)*cos(x)'), c('TR11(cos(2*x))', '-sin(x)^2 + cos(x)^2'), c('TR11(sin(4*x))', '4*(-sin(x)^2 + cos(x)^2)*sin(x)*cos(x)'),
   c('TR11(sin(4*x/3))', '4*(-sin(x/3)^2 + cos(x/3)^2)*sin(x/3)*cos(x/3)'), c('TR11(cos(I(2)))', 'cos(2)'), c('TR11(cos(I(4)), 2)', '-sin(2)^2 + cos(2)^2'),
   c('TR11(cos(6*pi/7) + cos(3*pi/7), 3*pi/7)', '-sin(3*pi/7)^2 + cos(3*pi/7)^2 + cos(3*pi/7)'),
   c('TR11(sin(x/3)/(cos(x/6)))', 'sin(x/3)/cos(x/6)'),
   c('TR12(tan(x + y))', '(tan(x) + tan(y))/(-tan(x)*tan(y) + 1)'),
   c('TR12i((ta + tb)/(-ta*tb + 1))', 'tan(a + b)'), c('TR12i((ta + tb)/(ta*tb - 1))', '-tan(a + b)'), c('TR12i((-ta - tb)/(ta*tb - 1))', 'tan(a + b)'),
   c('TR12i(expand((ta + tb)/(-ta*tb + 1)^2*(-3*ta - 3*tc)/(2*(ta*tc - 1))))', '-3*tan(a + b)*tan(a + c)/(2*(tan(a) + tan(b) - 1))'),
   c('TR13(tan(I(3))*tan(I(2)))', '-tan(2)/tan(5) - tan(3)/tan(5) + 1'), c('TR13(cot(I(3))*cot(I(2)))', 'cot(2)*cot(5) + 1 + cot(3)*cot(5)'),
   c('TRmorrie(cos(x)*cos(2*x))', 'sin(4*x)/(4*sin(x))'),
   c('TRmorrie(7*cos(I(0))*cos(I(1))*cos(I(2))*cos(I(3))*cos(I(4))*cos(I(5))*cos(I(6))*cos(I(7))*cos(I(8))*cos(I(9)))', '7*sin(12)*sin(16)*cos(5)*cos(7)*cos(9)/(64*sin(1)*sin(3))'),
   c('TRmorrie(cos(pi/7)*cos(2*pi/7)*cos(4*pi/7))', '-sin(3*pi/7)*cos(3*pi/7)/(4*sin(pi/7))'),
   c('TR8(TRmorrie(cos(pi/7)*cos(2*pi/7)*cos(4*pi/7)))', '-1/8'),
   c('TRmorrie(cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9))', '1/16'),
   c('TR8(cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9))', '1/16'),
   c('TR14((cos(x) - 1)*(cos(x) + 1))', '-sin(x)^2'), c('TR14((sin(x) - 1)*(sin(x) + 1))', '-cos(x)^2'),
   c('TR14((cos(x) + 1)*(cos(x) - 1)*(cos(y) - 1)*2*(cos(y) + 1)*(3*(cos(y) - 1))*(3*(cos(y) + 1))*(x - 1))', '-18*(x - 1)*sin(x)^2*sin(y)^4'),
   c('TR15(1 - 1/sin(x)^2)', '-cot(x)^2'), c('TR16(1 - 1/cos(x)^2)', '-tan(x)^2'), c('TR111(1 - 1/tan(x)^2)', '1 - cot(x)^2'),
   c('TR22(1 + tan(x)^2)', 'sec(x)^2'), c('TR22(1 + cot(x)^2)', 'csc(x)^2'),
   c('TRpower(sin(x)^6)', '-15*cos(2*x)/32 + 3*cos(4*x)/16 - cos(6*x)/32 + 5/16'),
   c('TRpower(sin(x)^3*cos(2*x)^4)', '(3*sin(x)/4 - sin(3*x)/4)*(cos(4*x)/2 + cos(8*x)/8 + 3/8)'),
   c('fu(sin(I(50))^2 + cos(I(50))^2 + sin(pi/6))', '3/2'), c('fu(sqrt(I(6))*cos(x) + sqrt(I(2))*sin(x))', '2*sqrt(2)*sin(x + pi/3)'),
   c('fu(sin(x)^4 - cos(y)^2 + sin(y)^2 + 2*cos(x)^2)', 'cos(x)^4 - 2*cos(y)^2 + 2'), c('fu(Rational(1,2) - cos(2*x)/2)', 'sin(x)^2'),
   c('fu(sin(a)*(cos(b) - sin(b)) + cos(a)*(sin(b) + cos(b)))', 'sqrt(2)*sin(a + b + pi/4)'), c('fu(sqrt(I(3))*cos(x)/2 + sin(x)/2)', 'sin(x + pi/3)'),
   c('fu(1 - sin(2*x)^2/4 - sin(y)^2 - cos(x)^4)', '-cos(x)^2 + cos(y)^2'), c('fu(cos(4*pi/9))', 'sin(pi/18)'),
   c('fu(cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9))', '1/16'),
   c('fu(tan(7*pi/18) + tan(5*pi/18) - sqrt(I(3))*tan(5*pi/18)*tan(7*pi/18))', '-sqrt(3)'),
   c('fu(sin(x)/cos(x))', 'tan(x)'), c('fu(sin(x)/cos(x), measure = function(e) -count_ops(e))', 'sin(x)/cos(x)'))
  
  for (cs in cases) expect_equal(sstr(eval(parse(text = cs[1L]))), cs[2L], label = cs[1L])
})

test_that("core changes made for fu follow SymPy", {
  p <- sym_pi
  expect_equal(vapply(list(sin(8*p/7), cos(4*p/7), cos(6*p/7), sin(4*p/7), cos(13*p/7), tan(5*p/7), cot(5*p/7),
                           tan(8*p/7), sin(15*p/7), sin(9*p/7), cos(9*p/7), sec(5*p/7), csc(8*p/7)), sstr, ""),
               c("-sin(pi/7)", "-cos(3*pi/7)", "-cos(pi/7)", "sin(3*pi/7)", "cos(pi/7)", "-tan(2*pi/7)", "-cot(2*pi/7)",
                 "tan(pi/7)", "sin(pi/7)", "-sin(2*pi/7)", "-cos(2*pi/7)", "-sec(2*pi/7)", "-csc(pi/7)"))
  expect_equal(sstr(expand(1/((x + 1)*(y + 1)))), "1/(x*y + x + y + 1)")
  expect_equal(sstr(expand((x + 1)/(2*(y + 1)))), "x/(2*y + 2) + 1/(2*y + 2)")
  ops <- list(list(cos(4*p/9), 3L, "COS + DIV + MUL"), list(sin(p/18), 2L, "DIV + SIN"), list(x^2 + sin(x)/y, 4L, "ADD + DIV + POW + SIN"),
              list(-x - y, 2L, "NEG + SUB"), list(-2*x/3, 3L, "DIV + MUL + NEG"), list(sqrt(x), 2L, "DIV + POW"),
              list(Eq(x, y + 1), 2L, "ADD + EQ"), list(Function("f")(x) + 1, 2L, "ADD + FUNC_F"))
  for (o in ops) { expect_identical(count_ops(o[[1L]]), o[[2L]]); expect_equal(sstr(count_ops(o[[1L]], visual = TRUE)), o[[3L]]) }
  expect_equal(sstr(signsimp(y - x*(y - x))), "x*(x - y) + y")
  expect_equal(sstr(signsimp(x/(-x^2 - x + 1))), "-x/(x^2 + x - 1)")
  expect_equal(sstr(signsimp(x*(1 - x))), "-x*(x - 1)")
  expect_equal(sstr(simplify(x/(1 - x))), "-x/(x - 1)")
  nd <- .as_numer_denom(x/(2*y) + 1/(3*y))
  expect_equal(c(sstr(nd[[1L]]), sstr(nd[[2L]])), c("3*x + 2", "6*y"))
})
