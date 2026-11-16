# Regression cases for polynomials, Poly, domains, factorization and polynomial arithmetic.
# Expected values are SymPy's output for the same calls. A case that is slow
# on one core starts with skip_on_cran() or skip_unless_full_tests().

test_that("factor of x^20 - y^20 matches SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(sym_factor(x^20 - y^20)),
               paste0("(x - y)*(x + y)*(x^2 + y^2)*(x^4 - x^3*y + x^2*y^2 - x*y^3 + y^4)*",
                      "(x^4 + x^3*y + x^2*y^2 + x*y^3 + y^4)*(x^8 - x^6*y^2 + x^4*y^4 - x^2*y^6 + y^8)"))
})

test_that("coefficients over a fraction field are a numerator over a denominator", {
  x <- Symbol("x")
  a <- Symbol("a")
  p <- Poly(x^2 * (a + 1) / a + x / (2 * a) - x / a^2, x)
  expect_identical(vapply(all_coeffs(p), sstr, ""), c("(a + 1)/a", "(a - 2)/(2*a^2)", "0"))
})

test_that("exp and power factors read as powers of the generators", {
  x <- Symbol("x"); n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(Poly(exp(Integer(6)) + S_E, S_E)), "Poly(E^6 + E, E, domain='ZZ')")
  expect_equal(sstr(Poly(exp(Integer(6)) * x + 1, x)), "Poly(exp(6)*x + 1, x, domain='ZZ[E]')")
  expect_equal(sstr(as_expr(Poly(exp(2 * x) + exp(x), exp(x)))), "exp(2*x) + exp(x)")
  expect_equal(degree(Poly(2^(2 * n) + 2^n, 2^n)), 2)
})

test_that("Poly reads I as a coefficient in ZZ_I and QQ_I", {
  x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a")
  expect_equal(sstr(Poly(S_I * x + 1, x)), "Poly(I*x + 1, x, domain='ZZ_I')")
  expect_equal(sstr(Poly(S_I * x + y)), "Poly(I*x + y, x, y, domain='ZZ_I')")
  expect_equal(sstr(Poly(S_I * x / 2 + 1)), "Poly(I/2*x + 1, x, domain='QQ_I')")
  expect_equal(sstr(Poly(S_I * a * x + 1, x)), "Poly(I*a*x + 1, x, domain='ZZ_I[a]')")
  expect_equal(sstr(Poly(S_I * sqrt(Integer(3)) * x + 1, x)), "Poly(sqrt(3)*I*x + 1, x, domain='EX')")
})

test_that("terms_gcd follows SymPy's options", {
  x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a"); I <- S_I
  expect_equal(sstr(terms_gcd(2 * x + 2 * I * x)), "x*(2 + 2*I)")
  expect_equal(sstr(terms_gcd(I * x + I * y)), "I*x + I*y")
  expect_equal(sstr(terms_gcd(x^2 * I + x)), "x*(I*x + 1)")
  expect_equal(sstr(terms_gcd((1 + I) * x + (1 + I) * y)), "(1 + I)*(x + y)")
  expect_equal(sstr(terms_gcd(a * x + a * y, x)), "a*(x + y)")
  expect_equal(sstr(terms_gcd(x / 2 + x^2 / 4)), "x*(x + 2)/4")
  expect_equal(sstr(terms_gcd(x * (2 * y + 2 * x), deep = TRUE)), "2*x*(x + y)")
  expect_equal(sstr(terms_gcd(Eq(2 * x + 2 * y, 4 * x * y))), "Eq(2*(x + y), 4*x*y)")
  expect_equal(sstr(terms_gcd(x / 2 + y / 3, clear = FALSE)), "(3*x + 2*y)/6")
  expect_equal(sstr(terms_gcd(x / 2 + 1, clear = FALSE)), "x/2 + 1")
  expect_error(terms_gcd(Lt(x, y)), "Inequalities")
})

test_that("qp_exquo inverts multiplication and rejects non divisors", {
  set.seed(20261115)
  rpoly <- function(n, terms, deg, cmax) {
    e <- matrix(sample(0:deg, terms * n, replace = TRUE), terms, n)
    c <- gmp::as.bigq(floor(stats::runif(terms, 1, cmax + 1)) * sample(c(-1, 1), terms, replace = TRUE))
    symcompR:::qp_norm(e, c)
  }
  slow <- function(f, g) {
    d <- symcompR:::qp_div_multi(f, list(g))
    if (!symcompR:::qp_is_zero(d$r)) NULL else d$q[[1L]]
  }
  same <- function(a, b) identical(a$e, b$e) && all(a$c == b$c)
  for (i in 1:60) {
    n <- sample(1:4, 1L)
    cmax <- if (i %% 3 == 0) 10^12 else 9
    f <- rpoly(n, sample(1:5, 1L), 3L, cmax)
    g <- rpoly(n, sample(1:4, 1L), 2L, cmax)
    if (i %% 4 == 0) g <- symcompR:::qp(g$e, g$c * 3L)
    fg <- symcompR:::qp_mul(f, g)
    expect_true(same(symcompR:::qp_exquo(fg, g), f))
    h <- symcompR:::qp_add(fg, symcompR:::qp_const(1L, n))
    expect_identical(is.null(symcompR:::qp_exquo(h, g)), is.null(slow(h, g)))
    q <- symcompR:::qp_exquo(h, g)
    if (!is.null(q)) expect_true(same(q, slow(h, g)))
  }
})

test_that(".zp_mul agrees with the schoolbook product", {
  set.seed(115)
  school <- function(f, g) {
    out <- gmp::as.bigz(rep(0, length(f) + length(g) - 1L))
    for (i in seq_along(f)) { idx <- i:(i + length(g) - 1L); out[idx] <- out[idx] + f[i] * g }
    out
  }
  rz <- function(n, big) {
    v <- gmp::as.bigz(floor(stats::runif(n, -1000, 1000)))
    if (big) v * gmp::as.bigz(2)^sample(0:90, n, replace = TRUE) else v
  }
  for (i in 1:80) {
    f <- rz(sample(1:20, 1L), i %% 2 == 0); g <- rz(sample(1:20, 1L), i %% 3 == 0)
    if (i %% 5 == 0) f[1L] <- 0
    expect_true(all(symcompR:::.zp_mul(f, g) == school(f, g)))
  }
})

test_that("polynomial extras match SymPy", {
  x <- Symbol("x")
  y <- Symbol("y")
  z <- Symbol("z")
  a <- Symbol("a")
  b <- Symbol("b")
  s <- function(e) sstr(e)
  iv <- function(r) lapply(r, function(o) c(s(o[[1]][[1]]), s(o[[1]][[2]]), o[[2]]))
  expect_equal(vapply(sym_decompose(x^4 + 2*x^3 - x - 1), s, ""), c("x^2 - x - 1", "x^2 + x"))
  expect_equal(vapply(sym_decompose(x^6 - 3*x^4 + 3*x^2 - 1), s, ""), c("x^3 - 3*x^2 + 3*x - 1", "x^2"))
  g <- gff_list(x^5 + 2*x^4 - x^3 - 2*x^2)
  expect_equal(lapply(g, function(p) c(s(p[[1]]), p[[2]])), list(c("x", "1"), c("x + 2", "4")))
  expect_equal(vapply(subresultants(x^3 + x + 1, 3*x^2 + 1, x), s, ""), c("x^3 + x + 1", "3*x^2 + 1", "6*x + 9", "31"))
  expect_equal(s(sym_trunc(2*x^3 + 3*x^2 + 5*x + 7, 3)), "-x^3 - x + 1")
  expect_equal(vapply(symmetrize(x^2 + y^2), s, ""), c("-2*x*y + (x + y)^2", "0"))
  expect_equal(s(symmetrize(x^3 + y^3 + z^3)[[1]]), "3*x*y*z + (x + y + z)^3 - 3*(x + y + z)*(x*y + x*z + y*z)")
  expect_equal(s(rational_interpolate(list(c(1, 1), c(2, 4), c(3, 9), c(4, 16)), 2, x)), "x^2")
  expect_equal(s(swinnerton_dyer_poly(2, x)), "x^4 - 10*x^2 + 1")
  expect_equal(s(nth_power_roots_poly(x^2 - 2, 2)), "x^2 - 4*x + 4")
  expect_equal(iv(intervals(x^2 - 3)), list(c("-2", "-1", "1"), c("1", "2", "1")))
  expect_equal(iv(intervals((x - 1)^2 * (x^2 - 2))), list(c("-2", "-1", "1"), c("1", "1", "2"), c("1", "2", "1")))
  expect_equal(ground_roots(x^6 - 4*x^4 + 4*x^3 - x^2), c("0" = 2, "1" = 2))
  expect_equal(vapply(hurwitz_conditions(x^2 + a*x + b, x), s, ""), c("a", "b"))
  expect_equal(vapply(schur_conditions(x^2 + a*x + b, x), s, ""), c("2*a*b - 2*a - 2*b^2 + 2", "-a^2 + b^2 + 2*b + 1"))
  expect_equal(s(minimal_polynomial(sqrt(S(2)) + sqrt(S(3)), x)), "x^4 - 10*x^2 + 1")
  expect_equal(s(minimal_polynomial(sqrt(S(2)) * sym_I, x)), "x^2 + 2")
  expect_equal(s(minimal_polynomial(cos(sym_pi / 5), x)), "4*x^2 - 2*x - 1")
})

test_that("count_roots takes SymPy's numeric bounds", {
  skip_on_cran()
  x <- Symbol("x"); f <- x^4 - 5 * x^2 + 4
  expect_equal(count_roots(f, 0, 3), 2L)
  expect_equal(count_roots(f), 4L)
  expect_equal(count_roots(f, -3, 0), 2L)
  expect_equal(count_roots(f, 0), 2L)
  expect_equal(count_roots(x^4 - 1, x, -2, 0), 1L)
})

test_that("polytools over GF(p) follow SymPy for any prime and any number of gens", {
  skip_on_cran()
  x <- Symbol("x"); y <- Symbol("y")
  big <- Integer("2305843009213693951")
  s <- function(e) if (is.list(e) && !is_sym(e) && !inherits(e, "Poly")) unname(vapply(e, sstr, "")) else sstr(e)
  expect_equal(s(pdiv(x^2*y + 1, 2*x*y + 3, modulus = 7)), c("2*x*y^2 - 3*y", "-3*y^2 + 2*y"))
  expect_equal(s(gcdex(x^3 + 1, x^2 + 2, modulus = 7)), c("x - 3", "-x^2 + 3*x + 2", "1"))
  expect_equal(s(resultant(x^2 + y, x*y + 1, modulus = 5)), "y^3 + 1")
  expect_equal(s(discriminant(x^2*y + y*x + 3, modulus = 5)), "y^2 - 2*y")
  expect_equal(s(subresultants(x^2 + y, x*y + 1, modulus = 5)), c("x^2 + y", "x*y + 1", "y^3 + 1"))
  expect_equal(s(monic(3*x*y + y + 1, modulus = 5)), "x*y + 2*y + 2")
  expect_equal(s(primitive(3*x*y + 6, modulus = 7)), c("1", "3*x*y - 1"))
  expect_equal(s(terms_gcd(3*x^2*y + 6*x*y^2, modulus = 7)), "x*y*(3*x - y)")
  expect_equal(s(gcd(x^4*y^2 - y^2, x^2*y - y, modulus = 3)), "x^2*y - y")
  expect_equal(s(div(x^2 + y^2, x + y, modulus = 7)), c("x - y", "2*y^2"))
  expect_equal(s(cancel((x^2*y - y)/(x*y + y), modulus = 7)), "x - 1")
  expect_equal(s(compose(x^2 + y, x + y, modulus = 5)), "x^2 + 2*x*y + y^2 + y")
  expect_equal(s(sym_decompose(x^6 + 3*x^4 + 3*x^2 + 2, modulus = 7)), c("x^3 + 3*x^2 + 3*x + 2", "x^2"))
  expect_error(gcdex(x*y + 1, x + y, modulus = 5), "univariate polynomial expected")
  expect_error(exquo(x^2 + 3, x + 1, modulus = 5), "x \\+ 1 does not divide x\\^2 - 2")
  expect_error(resultant(x^2 + y, x*y + 1, y, modulus = 5), "expected an integer, got x\\^2")
  expect_equal(s(gcd(Poly(x^2 - 1, modulus = 5), Poly(x + 1, modulus = 5))), "Poly(x + 1, x, modulus=5)")
  expect_equal(s(resultant(x^3 + x + 1, x^2 + 3*x + 7, modulus = big)), "349")
  expect_equal(s(discriminant(x^3 + x + 1, modulus = big)), "-31")
  expect_equal(s(gcdex(x^3 + 1, x^2 + 2, modulus = big))[3L], "1")
  expect_equal(s(monic(3*x^2 + 1, modulus = big)), "x^2 - 768614336404564650")
  expect_equal(s(subresultants(x^3 + x + 1, x^2 + 3, modulus = big)), c("x^3 + x + 1", "x^2 + 3", "1 - 2*x", "13"))
  # beyond SymPy, which raises NotImplementedError for several generators
  expect_equal(sstr(factor(x^2*y + x, modulus = 5)), "x*(x*y + 1)")
  expect_equal(sstr(factor(x^2 + y^2, modulus = 5)), "(x - 2*y)*(x + 2*y)")
  expect_equal(sstr(factor(x^2 + y^2, modulus = 7)), "x^2 + y^2")
  expect_equal(sstr(factor(expand((x + y)^2*(x - y)), modulus = 2)), "(x + y)^3")
  expect_equal(sstr(factor(expand((x^2 + y^2 + 1)*(x*y + 2)), modulus = big)), "(x*y + 2)*(x^2 + y^2 + 1)")
  fl <- factor_list(expand((x*y + 1)^2*(x + y)^3*(x^2 + y + 3)), modulus = 5)
  expect_equal(vapply(fl$factors, function(f) paste(sstr(f[[1L]]), sstr(f[[2L]])), ""),
               c("x*y + 1 2", "x + y 3", "x^2 + y - 2 1"))
  sl <- sqf_list(expand((x*y + 1)^2*(x + y)^2*(x + 2)), modulus = 5)
  expect_equal(vapply(sl$factors, function(f) paste(sstr(f[[1L]]), sstr(f[[2L]])), ""),
               c("x + 2 1", "x^2*y + x*y^2 + x + y 2"))
  expect_equal(sstr(sqf_part(x^2*y + 2*x*y + y, modulus = 5)), "x*y + y")
})

test_that("polynomial ring conversion accepts Indexed generators", {
  skip_on_cran()
  a0 <- Indexed("A", 0); l <- Symbol("l")
  expect_equal(format(construct_domain(list(a0^2 - 2*a0 + 1))[[1L]]), "ZZ[A[0]]")
  r <- roots(l^2 - 2*a0*l + a0^2, l)
  expect_equal(unname(vapply(r, sstr, "")), "A[0]")
})

test_that("RootSum preprocesses symbolic and integer coefficients", {
  skip_on_cran()
  x <- Symbol("x")
  a <- Symbol("a")
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
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  r <- sqf_norm(x * y, sqrt(Integer(2)))
  expect_equal(vapply(r[[1L]], sstr, ""), c("1", "1"))
  expect_equal(sstr(r[[3L]]), "x^2*y^2 - 2*x^2 - 2*y^2 + 4")
})

test_that("symmetrize infers generators and handles constants", {
  skip_on_cran()
  x <- Symbol("x")
  y <- Symbol("y")
  r <- symmetrize(sin(x)^2 + sin(y)^3)
  expect_equal(sstr(r[[1L]]), "(sin(x) + sin(y))^2 - 2*sin(x)*sin(y)")
  expect_equal(sstr(r[[2L]]), "sin(y)^3 - sin(y)^2")
  r <- symmetrize(Integer(3), formal = TRUE)
  expect_equal(sstr(r[[1L]]), "3")
  expect_equal(sstr(r[[2L]]), "0")
  expect_length(r[[3L]], 0L)
})

test_that("cancel works over the Gaussian integers and rationals", {
  skip_on_cran()
  x <- Symbol("x")
  I <- sym_I
  expect_equal(sstr(cancel((x + I) / (x^2 + 1))), "1/(x - I)")
  expect_equal(sstr(cancel((1 + I) / (2 * x))), "I*(1 - I)/(2*x)")
  expect_equal(sstr(cancel((x^2 + 1) / (x - I))), "x + I")
  expect_equal(sstr(cancel((x^2 + 1) / (2 * x - 2 * I))), "x/2 + I/2")
  expect_equal(sstr(cancel((x / 2 + I / 3) / (x^2 + Rational(4, 9)))), "3/(6*x - 4*I)")
  expect_equal(sstr(cancel(1 / (x - I / 2))), "2/(2*x - I)")
})

test_that("Poly floor division, remainder and equality", {
  x <- Symbol("x"); y <- Symbol("y")
  p <- Poly(x^2 + 2, x); q <- Poly(x, x)
  expect_equal(sstr(p %/% q), "Poly(x, x, domain='ZZ')")
  expect_equal(sstr(p %% q), "Poly(2, x, domain='ZZ')")
  expect_false(Poly(x + y, x) == Poly(x + y, x, y))
  expect_true(Poly(x + y, x) == Poly(x + y, x))
})

test_that("Poly() of a Poly applies requested generators and options", {
  x <- Symbol("x"); y <- Symbol("y")
  p <- Poly(x + 2, x)
  expect_identical(Poly(p, x), p)
  expect_equal(sstr(Poly(p, x, modulus = 2)), "Poly(x, x, modulus=2)")
  expect_equal(sstr(Poly(Poly(x * y + 1, x, y), y, x)), "Poly(y*x + 1, y, x, domain='ZZ')")
})

test_that("arithmetic between Polys unifies their generators", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_equal(sstr(Poly(x, x) + Poly(y, y)), "Poly(x + y, x, y, domain='ZZ')")
  expect_equal(sstr(Poly(y, y) * Poly(x, x)), "Poly(y*x, y, x, domain='ZZ')")
})
