# Reference values from SymPy.
x <- Symbol("x")

test_that("fourier_series coefficients and truncation", {
  s <- fourier_series(x^2, list(x, -S_Pi, S_Pi))
  expect_equal(sstr(truncate(s)), "-4*cos(x) + cos(2*x) + pi^2/3")
  expect_equal(sstr(truncate(s, 4)), "-4*cos(x) + cos(2*x) - 4*cos(3*x)/9 + pi^2/3")
  s2 <- fourier_series(x, list(x, -S_Pi, S_Pi))
  expect_equal(sstr(truncate(s2, 3)), "2*sin(x) - sin(2*x) + 2*sin(3*x)/3")
  p <- fourier_series(Piecewise(list(0, x < 0), list(1, TRUE)), list(x, -S_Pi, S_Pi))
  expect_equal(sstr(truncate(p, 4)), "2*sin(x)/pi + 2*sin(3*x)/(3*pi) + 2*sin(5*x)/(5*pi) + 1/2")
  expect_equal(sstr(truncate(fourier_series(x, list(x, 0, 1)), 3)), "-sin(2*pi*x)/pi - sin(4*pi*x)/(2*pi) + 1/2")
  expect_equal(sstr(sigma_approximation(s, 4)),
               "-4*cos(x)*sinc(pi/4) + 2*cos(2*x)/pi - 4*cos(3*x)*sinc(3*pi/4)/9 + pi^2/3")
  expect_equal(latex(s), "- 4 \\cos{\\left(x \\right)} + \\cos{\\left(2 x \\right)} + \\frac{\\pi^{2}}{3} + \\ldots")
})

test_that("finite Fourier series and transformations", {
  expect_equal(sstr(fourier_series(sin(x)^2, list(x, -S_Pi, S_Pi))), "FiniteFourierSeries(sin(x)^2, (x, -pi, pi), (1/2, {2: -1/2}, {}))")
  expect_equal(sstr(fourier_series(sin(x) + cos(2*x) + 3)), "FiniteFourierSeries(sin(x) + cos(2*x) + 3, (x, -pi, pi), (3, {2: 1}, {1: 1}))")
  expect_equal(sstr(truncate(fourier_series(cos(x)^3), NULL)), "3*cos(x)/4 + cos(3*x)/4")
  s2 <- fourier_series(x, list(x, -S_Pi, S_Pi))
  expect_equal(sstr(truncate(shift(s2, 1), 2)), "2*sin(x) + 1")
  expect_equal(sstr(truncate(scale(s2, 2), 2)), "4*sin(x) - 2*sin(2*x)")
  expect_equal(sstr(truncate(-s2, 2)), "-2*sin(x) + sin(2*x)")
  expect_equal(sstr(truncate(s2 + fourier_series(x^2, list(x, -S_Pi, S_Pi)), 3)),
               "2*sin(x) - sin(2*x) - 4*cos(x) + cos(2*x) + pi^2/3")
  expect_equal(sstr(sinc(S_Pi/2)), "2/pi")
  expect_equal(sstr(sym_integrate(Piecewise(list(x, x < 1), list(2, TRUE)), list(x, 0, 3))), "9/2")
})

test_that("p12 review fixes for Fourier series, Piecewise integrals and ranges", {
  a <- Symbol("a", real = TRUE); n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(sym_integrate(Piecewise(list(1, a > 0), list(0, TRUE)), list(x, 0, 1))), "Piecewise((1, a > 0), (0, True))")
  expect_identical(sym_integrate(Piecewise(list(1, x < 0)), list(x, -1, 1)), S_NaN)
  expect_equal(fourier_series(sin(1/x))$type, "FourierSeries")
  expect_equal(sstr(truncate(scalex(fourier_series(sin(x) + cos(2*x)), Rational(1, 2)), 2)), "sin(x/2) + cos(x)")
  expect_equal(sstr(fourier_series(sin(Integer("2147483648")*x))),
               "FiniteFourierSeries(sin(2147483648*x), (x, -pi, pi), (0, {}, {2147483648: 1}))")
  expect_equal(sstr(truncate(fourier_series(sin(x), finite = FALSE), 1)), "sin(x)")
  expect_identical(coeff(SeqPer(list(1, 2, 3)), Integer("9007199254740993")), S_One)
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(sqrt(Integer("10000000000000000000000000000000000000001")), Inf))),
               "Range(100000000000000000001, oo, 1)")
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(0, sqrt(Integer("99999999999999999999"))))), "Range(0, 10000000000, 1)")
  expect_equal(sstr(contains(Range(n, n + 6, 2), 2)), "Contains(2, Range(n, n + 6, 2))")
  expect_length(sequence_terms(Range(0, Inf, -1), 3), 0L)
  expect_equal(sstr(Range(-Inf, Inf)), "Range(-oo, oo, 1)")
  expect_length(sequence_terms(Range(1, 5, Inf), 3), 1L)
  expect_length(sequence_terms(Range(n, n + 5, 2), 5), 3L)
  expect_equal(Intersection(Range(0, 10), Interval(a, 5))$type, "Intersection")
})

test_that("p13 h2 fix: resonant coefficients survive shifts and scalings", {
  x <- Symbol("x")
  s <- fourier_series(sin(x), list(x, -sym_pi, sym_pi), finite = FALSE)
  expect_equal(sstr(s$args[[3L]]$args[[3L]]),
               "SeqFormula(Piecewise((-pi, Eq(_n, -1)), (pi, Eq(_n, 1)), (-2*sin(_n*pi)/(_n^2 - 1), True))*sin(_n*x)/pi, (_n, 1, oo))")
  expect_equal(sstr(sigma_approximation(shiftx(s, sym_pi/2), 2)), "2*cos(x)/pi")
  expect_equal(sstr(sigma_approximation(scalex(s, 2), 2)), "2*sin(2*x)/pi")
  expect_equal(sstr(truncate(fourier_series(x*sin(x), list(x, -sym_pi, sym_pi), finite = FALSE), 3)),
               "-cos(x)/2 - 2*cos(2*x)/3 + 1")
})

test_that("p14 h2 fix: parameter-dependent resonances", {
  x <- Symbol("x"); a <- Symbol("a")
  expect_equal(sstr(subs(sigma_approximation(fourier_series(sin(a*x), list(x, -sym_pi, sym_pi), finite = FALSE), 2), a, 1)), "2*sin(x)/pi")
})
