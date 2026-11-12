# Expected values are SymPy's outputs for the same calls.

x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a"); I <- sym_I; oo <- S_Infinity

test_that("limits keep tiny exact coefficients and certify signs", {
  expect_equal(sstr(limit(Float(1e-12) / x, x, 0)), "oo")
  expect_equal(sstr(limit(exp(-Float(1e-12) * x), x, oo)), "0")
  expect_equal(sstr(limit(sin(exp(Integer(-1000))) / x, x, 0)), "oo")
  expect_equal(sstr(limit(I / x, x, 0, "+-")), "zoo")
  expect_true(is_infinite(oo * I))
})

test_that("limits on inverse trigonometric branch cuts follow SymPy", {
  expect_equal(sstr(limit((asin(2 + I * x) - asin(Integer(2))) / x, x, 0)), "oo*sign(pi - 2*asin(2))")
  expect_equal(sstr(limit(asin(2 + I * x) / x, x, 0)), "oo*sign(pi - asin(2))")
  expect_equal(sstr(limit(asin(2 + I * Rational(1, "10000000000000") * x), x, 0)), "pi - asin(2)")
  v <- limit(asin(2 + x + I * x^2), x, 0)
  expect_equal(as.complex(N(v)), as.complex(N(sym_pi - asin(Integer(2)))), tolerance = 1e-12)
})

test_that("log at a pole takes the side of the cut", {
  expect_equal(sstr(limit(log(-1 / x - I) - log(1 / x), x, 0)), "-I*pi")
  expect_equal(sstr(limit(log(-x - I) - log(x), x, oo)), "-I*pi")
  expect_equal(sstr(limit(log(-1 / x + I) - log(1 / x), x, 0)), "I*pi")
  expect_equal(sstr(series(log(-1 / x - I), x, 0, 3)), "-I*pi - log(x) + I*x + x^2/2 + O(x^3)")
})

test_that("gamma and loggamma at poles", {
  expect_equal(sstr(limit(gamma(x), x, -1)), "-oo")
  expect_equal(sstr(limit(gamma(x - 1), x, 0)), "-oo")
  expect_equal(sstr(limit(x * gamma(x - 1), x, 0)), "-1")
  expect_equal(sstr(limit(gamma(x - 2) * x, x, 0)), "1/2")
  expect_equal(sstr(limit(gamma(x), x, -2, "-")), "-oo")
  expect_equal(sstr(gruntz(loggamma(x) / log(x), x, 0)), "-1")
  expect_equal(sstr(limit(loggamma(-x), x, oo)), "zoo")
  expect_equal(sstr(limit(loggamma(x - 1), x, 0)), "oo")
  expect_equal(sstr(series(loggamma(x), x, 0, 2)), "-log(x) - EulerGamma*x + O(x^2)")
  expect_equal(sstr(series(gamma(x - 2), x, 0, 1)), "1/(2*x) + 3/4 - EulerGamma/2 + O(x)")
  expect_error(series(loggamma(x - 1), x, 0, 2))
  expect_equal(sstr(loggamma(Rational(5, 2))), "log(3*sqrt(pi)/4)")
  expect_equal(sstr(loggamma(Integer(-3))), "oo")
  expect_equal(sstr(loggamma(-oo)), "zoo")
})

test_that("limit forms of SymPy's Limit.doit", {
  expect_equal(sstr(limit(S_Infinity / (x^2 - 4), x, oo)), "Limit(oo/(x^2 - 4), x, oo, dir='-')")
  expect_equal(sstr(limit(3 * sec(-4 * sym_pi * x + x / 3), x, 3 * sym_pi / (-2 + 24 * sym_pi))), "-oo")
  expect_equal(sstr(limit((I * x - 1)^Rational(1, 3), x, 0, "-")), "-(-1)^(2/3)")
})

test_that("arbitrary precision keeps tiny and cancelling values", {
  expect_equal(sstr(N(sin(exp(Integer(-1000))), 20)), "5.0759588975494567653e-435")
  expect_equal(sstr(N(sin(Integer(355)), 30)), "-0.0000301443533594884492143302800087")
  expect_equal(sstr(N(log(1 + exp(Integer(-100))), 20)), "3.7200759760208359630e-44")
  expect_equal(sstr(N(S_ComplexInfinity)), "zoo")
  expect_equal(sstr(N(exp(oo * (I + Rational(1, "10000000000000"))))), "zoo")
})

test_that("exp of an infinite multiple decides by the exact real part", {
  expect_equal(sstr(simplify(exp(oo * (-Rational(1, 10^13) + I)))), "0")
  expect_equal(sstr(subs(acot(-x), x, 0)), "pi/2")
  expect_equal(sstr(acot(Integer(-2))), "-acot(2)")
})

test_that("exact rational powers follow Integer._eval_power", {
  cases <- list(
    list(Integer(12), Rational(1, 4), "sqrt(2)*3^(1/4)"),
    list(Integer(504), Rational(1, 6), "sqrt(2)*3^(1/3)*7^(1/6)"),
    list(Integer(18), Rational(5, 6), "3*2^(5/6)*3^(2/3)"),
    list(Integer(6480), Rational(1, 8), "5^(1/8)*sqrt(6)"),
    list(Integer(-4), Rational(1073741825, 2147483647), "2*(-1)^(1073741825/2147483647)*2^(3/2147483647)"),
    list(Integer(-12), Rational(1, 4), "(-3)^(1/4)*sqrt(2)"),
    list(Rational(3, 4), Rational(2, 3), "6^(2/3)/4"),
    list(Integer(4), Rational(-1, 2), "1/2"),
    list(Rational(-1, 2), Rational(-1, 3), "-(-1)^(2/3)*2^(1/3)"))
  for (cs in cases) expect_equal(sstr(cs[[1L]]^cs[[2L]]), cs[[3L]])
  expect_equal(sstr(Integer(-1)^Integer(10000001)), "-1")
})

test_that("RootSum preprocesses symbolic and integer coefficients", {
  expect_equal(sstr(RootSum(x^2 - a, Lambda(x, x^2), x)), "2*a")
  expect_equal(sstr(RootSum(x^2 + a, x = x)), "0")
  expect_equal(sstr(RootSum(x - 2)), "2")
  expect_equal(sstr(RootSum(x)), "0")
  expect_equal(sstr(RootSum(x^5 + 512 * x + 1024, Lambda(x, log(x)), x)),
               "RootSum(w^5 + 2*w + 1, Lambda(w, log(4*w)))")
  expect_equal(sstr(RootSum(x^3 - a^3 * x + a^3, Lambda(x, log(x)), x)),
               "RootSum(-a^3*w + a^3 + w^3, Lambda(w, log(w)))")
  expect_equal(sstr(N(RootSum(x^2 - 2, Lambda(x, 1 / (x + Integer("10000000000000"))), auto = FALSE))),
               "2.00000000000000e-13")
})

test_that("sqf_norm checks every generator", {
  r <- sqf_norm(x * y, sqrt(Integer(2)))
  expect_equal(vapply(r[[1L]], sstr, ""), c("1", "1"))
  expect_equal(sstr(r[[3L]]), "x^2*y^2 - 2*x^2 - 2*y^2 + 4")
})

test_that("symmetrize infers generators and handles constants", {
  r <- symmetrize(sin(x)^2 + sin(y)^3)
  expect_equal(sstr(r[[1L]]), "(sin(x) + sin(y))^2 - 2*sin(x)*sin(y)")
  expect_equal(sstr(r[[2L]]), "sin(y)^3 - sin(y)^2")
  r <- symmetrize(Integer(3), formal = TRUE)
  expect_equal(sstr(r[[1L]]), "3")
  expect_equal(sstr(r[[2L]]), "0")
  expect_length(r[[3L]], 0L)
})

test_that("Dict order, WildFunction and Poly domains", {
  k <- Integer("9007199254740992"); j <- k + 1
  expect_true(Dict(list(k, j), list(0, 0)) == Dict(list(j, k), list(0, 0)))
  F <- WildFunction("F")
  m <- sym_match(F, F)
  expect_equal(sstr(m$F), "F_")
  expect_equal(capture.output(print(Poly(x / 2 + a, x))), "Poly(1/2*x + a, x, domain='QQ[a]')")
})
