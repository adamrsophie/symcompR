# Fast checks of the special functions evaluated through the mpmath port.
# Expected strings are SymPy's str(N(expr, 15)); the full battery lives in
# test-evalf-mp-sympy.R.

check_n <- function(s, want, n = 15) expect_identical(sstr(N(parse_expr(s), n)), want, label = s)

test_that("Bessel, Hankel and Airy functions match SymPy", {
  check_n("besselj(1/3, 5/2)", "0.198320933418608")
  check_n("bessely(1, 2+I)", "-0.0163154378204725 + 0.599406841766854*I")
  check_n("besseli(1/2, 1+2*I)", "0.171623805371532 + 0.774048495527522*I")
  check_n("besselk(1/3, 1+I)", "0.0752632588864521 - 0.367969362760595*I")
  check_n("hankel1(1/3, 2+I)", "0.17508628430757 + 0.0816831806798347*I")
  check_n("jn(3, 15)", "-0.0644178013257291")
  check_n("yn(1, 1+I)", "-0.472261895218924 + 0.748611845881195*I")
  check_n("airyai(-10)", "0.0402412384864432")
  check_n("airyai(25/2)", "2.39682782607805e-14")
  check_n("airybiprime(1+I)", "0.075662844174966 + 0.783700998785455*I")
})

test_that("error functions and exponential integrals match SymPy", {
  check_n("erf(2-3*I)", "-20.8294614276146 - 8.68731827147016*I")
  check_n("erfc(27)", "5.23704892378926e-319")
  check_n("erfi(1/3+5*I)", "-4.34015867632437e-13 + 1.00000000000166*I")
  check_n("expint(1/2, 1+I)", "-0.0204291294153138 - 0.213024768215182*I")
  check_n("Ei(-2-3*I)", "0.0248262079441994 - 3.16190932850084*I")
  check_n("E1(40)", "1.03677326145166e-19")
  check_n("Si(30)", "1.56675654003035")
  check_n("Ci(-2-I/3)", "0.442265516129874 - 3.21438048781948*I")
  check_n("Chi(1+I)", "0.882172180555936 + 1.28354719327494*I")
  check_n("li(100)", "30.1261415840796")
  check_n("Li(10)", "5.12043572466981")
  check_n("fresnels(20)", "0.484084535925954")
  check_n("fresnelc(2-I/2)", "0.0949915944867347 - 1.78079635294637*I")
})

test_that("incomplete gamma, elliptic integrals and orthogonal polynomials match SymPy", {
  check_n("uppergamma(1+I, 2)", "0.0646320441729588 + 0.112840294049426*I")
  check_n("lowergamma(-3/2, 1+I)", "2.42022659192865 + 0.038271064163532*I")
  check_n("elliptic_k(1+I)", "1.50923695405127 + 0.625146415202697*I")
  check_n("elliptic_e(2, 2)", "0.599070117367796 + 1.0013371715505*I")
  check_n("elliptic_f(5/4, 1+I)", "1.27443708640421 + 0.398695740367898*I")
  check_n("elliptic_pi(-2, 1/3)", "0.974484017926750")
  check_n("legendre(1/3, 2)", "1.18887220307567")
  check_n("chebyshevu(1+I, 1/4)", "1.00100367964513 - 1.56734191050733*I")
  check_n("hermite(5/2, -1/2)", "1.42769027775410")
  check_n("assoc_legendre(5/2, 2, 1/2)", "4.10811266079151")
})

test_that("hypergeometric review examples match SymPy", {
  expect_identical(sstr(N(hyper(list(1 + sym_I), list(2), -100))),
                   "0.00348029442321405 + 0.0186566128060348*I")
  expect_identical(sstr(N(hyper(list(1, 1, 1), list(2, 2), 1))), "1.64493406684823")
  expect_identical(sstr(N(hyper(list(1, 1), list(2), 1))), "oo")
  expect_identical(sstr(N(gamma(sym_I / Integer(2)^200))), "-0.577215664901533 - 1.60693804425899e+60*I")
  check_n("meijerg(((),()),((0,),()),exp_polar(1))", "0.0659880358453125370767901875969", 30)
  check_n("pi", "3.14159265358979323846264338328", 30)
})

test_that("zeta, polylog, dirichlet_eta and LambertW match SymPy", {
  check_n("zeta(3+I)", "1.10721440843141 - 0.148290867178175*I")
  check_n("zeta(1/2, 2)", "-2.46035450880958681288949915252", 30)
  check_n("polylog(2, 1/3)", "0.36621322997706348762", 20)
  check_n("polylog(3, 2+I)", "1.89018684728999 + 1.84738717203211*I")
  check_n("dirichlet_eta(5/3)", "0.785787324658121169625725739737", 30)
  check_n("LambertW(1)", "0.56714329040978387300", 20)
  check_n("LambertW(-1/3, -1)", "-1.51213455165784247389673967807", 30)
})

test_that("special functions are not capped at 15 digits", {
  check_n("besselj(0, 1)", "0.76519768655796655145", 20)
  check_n("airyai(1)", "0.13529241631288141552", 20)
  check_n("Si(1)", "0.94608307036718301494", 20)
  check_n("li(3)", "2.1635885946671919729", 20)
  check_n("hyper((1, 1), (2,), 1/2)", "1.3862943611198906188", 20)
  check_n("airybiprime(1)", "0.932435933392776")
  expect_identical(sstr(N(gamma(1 + sym_I), 15)), "0.498015668118356 - 0.154949828301811*I")
  expect_identical(sstr(N(Matrix(list(list(sqrt(S(2)), pi))))), "Matrix([[1.41421356237310, 3.14159265358979]])")
})
