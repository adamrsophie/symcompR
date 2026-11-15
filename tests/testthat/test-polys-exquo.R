# Exact division of sparse polynomials: the packed fast paths must agree with
# the plain division on random input.

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
