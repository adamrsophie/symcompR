# Expected values are SymPy's str() outputs for the same calls.

pd_str <- function(v) gsub("**", "^", .dmpy_str(v), fixed = TRUE)

test_that("DomainMatrix basics match SymPy", {
  A <- DM(list(list(1, 2), list(3, 4)), ZZ)
  expect_identical(pd_str(A), "DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)")
  expect_identical(pd_str(A$det()), "-2")
  expect_identical(pd_str(A$rref()), "(DomainMatrix([[1, 0], [0, 1]], (2, 2), QQ), (0, 1))")
  expect_identical(pd_str(A$to_field()$inv()), "DomainMatrix([[-2, 1], [3/2, -1/2]], (2, 2), QQ)")
  expect_identical(pd_str(A$charpoly()), "[1, -5, -2]")
  expect_identical(pd_str(A$to_sparse()), "DomainMatrix({0: {0: 1, 1: 2}, 1: {0: 3, 1: 4}}, (2, 2), ZZ)")
  expect_identical(pd_str(A$rref_den()), "(DomainMatrix([[-2, 0], [0, -2]], (2, 2), ZZ), -2, (0, 1))")
  expect_identical(pd_str(DM(list(list(1, 1), list(1, 1)), QQ)$nullspace()), "DomainMatrix([[-1, 1]], (1, 2), QQ)")
  expect_identical(pd_str(DM(list(list(1, 2), list(3, 4)), GF(5))$inv()),
                   "DomainMatrix([[3 mod 5, 1 mod 5], [4 mod 5, 2 mod 5]], (2, 2), GF(5))")
})

test_that("normal forms accept a DomainMatrix", {
  expect_identical(pd_str(hermite_normal_form(DM(list(list(2, 7), list(0, 3)), ZZ))),
                   "DomainMatrix([[2, 1], [0, 3]], (2, 2), ZZ)")
  expect_identical(pd_str(smith_normal_form(DM(list(list(2, 4), list(6, 8)), ZZ))),
                   "DomainMatrix({0: {0: 2}, 1: {1: 4}}, (2, 2), ZZ)")
})

test_that("DomainMatrix errors carry SymPy's classes", {
  A <- DM(list(list(1, 2), list(3, 4)), ZZ)
  expect_error(A + DM(list(list(1, 2)), ZZ), class = "DMShapeError")
  expect_error(A$add(DM(list(list(1, 2), list(3, 4)), QQ)), class = "DMDomainError")
})

test_that("power series rings match SymPy", {
  R <- PythonPowerSeriesRingZZ()
  Q <- PythonPowerSeriesRingQQ(4)
  expect_identical(R$pretty(R$inverse(R(c(1, 1)))), "1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6)")
  expect_identical(Q$pretty(Q$exp(Q$gen)), "1 + x + 1/2*x**2 + 1/6*x**3 + O(x**4)")
  expect_identical(Q$pretty(Q$sqrt(Q(c(4, 4, 1)))), "2 + x")
  expect_identical(pd_str(R$to_list(R$pow_int(R(c(1, 2, 3)), 2))), "[1, 4, 10, 12, 9]")
})

test_that("generators that are not symbols stay independent in a domain", {
  K <- frac_field(QQ, exp(Rational(1, 2)))
  g <- .gen_dummy(exp(Rational(1, 2)))
  expect_identical(sstr(.gd_in(K, sym_E)), sstr(g^2))
  expect_identical(sstr(.gd_out(.gd_in(K, sym_E))), "E")
})

test_that("Poly arithmetic keeps a generator that is not a symbol", {
  x <- Symbol("x")
  p <- Poly(x * exp(Rational(1, 2)) + 1, x)
  expect_identical(paste(utils::capture.output(print(p^2)), collapse = ""),
                   "Poly(E*x^2 + 2*exp(1/2)*x + 1, x, domain='ZZ[exp(1/2)]')")
})

test_that("cancel over fraction fields and GF(p) matches SymPy", {
  x <- Symbol("x"); a <- Symbol("a")
  expect_identical(pd_str(Poly(x / (a + 1) + 1, x)$cancel(Poly(x + a + 1, x))),
                   "(1/(a + 1), Poly(1, x, domain='ZZ(a)'), Poly(1, x, domain='ZZ(a)'))")
  expect_identical(pd_str(Poly(x / (a + 1) + 1, x)$domain), "ZZ(a)")
  expect_identical(sstr(cancel((x^2 - 1) / (x + 1), modulus = 2)), "x + 1")
  F <- field("x,y", ZZ)
  expect_identical(format(Rational(3, 7) + F[[2L]]), "(7*x + 3)/7")
})

test_that("Poly domains, GF(p) gcd, sqf_part and interpolation match SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_identical(pd_str(Poly(x^2 + 1, x, domain = "QQ")), "Poly(x^2 + 1, x, domain='QQ')")
  expect_identical(pd_str(gcd(Poly(x^2 + 1, x, modulus = 2), Poly(x + 1, x, modulus = 2))), "Poly(x + 1, x, modulus=2)")
  expect_identical(sstr(sqf_part(x^3 - 3*x - 2)), "x^2 - x - 2")
  expect_identical(pd_str(sqf_part(Poly(x^3 - 3*x - 2))), "Poly(x^2 - x - 2, x, domain='ZZ')")
  expect_identical(sstr(interpolate(list(1, 8, 27, 64), x)), "x^3")
  expect_identical(sstr(rational_interpolate(list(c(1, 2), c(2, 3), c(3, 5), c(4, 7)), 1, x)),
                   "(57/2 - 9*x/2)/(x^2 - 17*x/2 + 39/2)")
  expect_identical(sstr(interpolating_poly(3, x, list(1, 2, 3), list(1, 4, 9))),
                   "(x - 3)*(x - 2)/2 - 4*(x - 3)*(x - 1) + 9*(x - 2)*(x - 1)/2")
  expect_error(exquo(x^2 + 1, x - 1), class = "ExactQuotientFailed")
})

test_that("PurePoly, Groebner order and all_roots match SymPy", {
  x <- Symbol("x"); y <- Symbol("y")
  expect_true(PurePoly(x^2 + 1) == PurePoly(y^2 + 1))
  expect_identical(pd_str(PurePoly(x^2 + y, x)), "PurePoly(_0^2 + y, _0, domain='ZZ[y]')")
  expect_identical(pd_str(groebner(list(x*y - 2*y, 2*y^2 - x^2), x, y, order = "grlex")),
                   "[y^3 - 2*y, x^2 - 2*y^2, x*y - 2*y]")
  expect_identical(pd_str(all_roots(x^3 - x - 1)),
                   "[CRootOf(x^3 - x - 1, 0), CRootOf(x^3 - x - 1, 1), CRootOf(x^3 - x - 1, 2)]")
  expect_identical(pd_str(Poly(x^2 - 1, x)$factor_list()),
                   "(1, [(Poly(x - 1, x, domain='ZZ'), 1), (Poly(x + 1, x, domain='ZZ'), 1)])")
})
