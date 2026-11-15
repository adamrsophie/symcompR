# Dense integer polynomial products by Kronecker substitution must agree with
# the schoolbook product, including negative and very large coefficients.

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
