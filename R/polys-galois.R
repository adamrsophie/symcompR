# Polynomials over GF(p) for small odd primes p < 2^25, stored as double
# vectors (highest degree first), after sympy/polys/galoistools.py.

.gf_trim <- function(a) {
  nz <- which(a != 0)
  if (!length(nz)) return(numeric(0))
  a[nz[1L]:length(a)]
}

.gf_inv <- function(a, p) {
  # modular inverse by extended Euclid
  a <- a %% p
  t0 <- 0
  t1 <- 1
  r0 <- p
  r1 <- a
  while (r1 != 0) {
    q <- r0 %/% r1
    tmp <- r0 - q * r1
    r0 <- r1
    r1 <- tmp
    tmp <- t0 - q * t1
    t0 <- t1
    t1 <- tmp
  }
  if (r0 != 1) stop("not invertible modulo p")
  t0 %% p
}

.gf_mulmod <- function(a, b, p) {
  # a*b mod p for doubles below 2^25 without overflow
  (a * b) %% p
}

.gf_add <- function(a, b, p) {
  la <- length(a)
  lb <- length(b)
  if (la < lb) a <- c(numeric(lb - la), a)
  if (lb < la) b <- c(numeric(la - lb), b)
  .gf_trim((a + b) %% p)
}

.gf_sub <- function(a, b, p) .gf_add(a, (p - b) %% p, p)

.gf_mul <- function(a, b, p) {
  la <- length(a)
  lb <- length(b)
  if (!la || !lb) return(numeric(0))
  out <- numeric(la + lb - 1L)
  for (j in seq_len(lb)) {
    if (b[j] == 0) next
    idx <- j:(j + la - 1L)
    out[idx] <- (out[idx] + (a * b[j]) %% p) %% p
  }
  .gf_trim(out)
}

.gf_monic <- function(a, p) {
  if (!length(a)) return(a)
  if (a[1L] == 1) return(a)
  (a * .gf_inv(a[1L], p)) %% p
}

.gf_divmod <- function(a, b, p) {
  b <- .gf_trim(b)
  a <- .gf_trim(a)
  db <- length(b) - 1L
  if (db < 0L) stop("division by zero polynomial")
  if (length(a) - 1L < db) return(list(q = numeric(0), r = a))
  inv <- .gf_inv(b[1L], p)
  r <- a
  nq <- length(a) - db
  q <- numeric(nq)
  for (i in seq_len(nq)) {
    cf <- (r[i] * inv) %% p
    q[i] <- cf
    if (cf != 0) {
      idx <- i:(i + db)
      r[idx] <- (r[idx] - (b * cf) %% p) %% p
    }
  }
  rr <- if (db > 0L) r[(length(r) - db + 1L):length(r)] else numeric(0)
  list(q = .gf_trim(q), r = .gf_trim(rr))
}

.gf_rem <- function(a, b, p) .gf_divmod(a, b, p)$r

.gf_gcd <- function(a, b, p) {
  a <- .gf_trim(a)
  b <- .gf_trim(b)
  while (length(b)) {
    r <- .gf_rem(a, b, p)
    a <- b
    b <- r
  }
  .gf_monic(a, p)
}

# Extended Euclid: returns s, t, g with s*a + t*b = g (monic).
.gf_gcdex <- function(a, b, p) {
  r0 <- .gf_trim(a)
  r1 <- .gf_trim(b)
  s0 <- 1
  s1 <- numeric(0)
  t0 <- numeric(0)
  t1 <- 1
  while (length(r1)) {
    qr <- .gf_divmod(r0, r1, p)
    r0 <- r1
    r1 <- qr$r
    tmp <- .gf_sub(s0, .gf_mul(qr$q, s1, p), p)
    s0 <- s1
    s1 <- tmp
    tmp <- .gf_sub(t0, .gf_mul(qr$q, t1, p), p)
    t0 <- t1
    t1 <- tmp
  }
  inv <- .gf_inv(r0[1L], p)
  list(s = .gf_trim((s0 * inv) %% p), t = .gf_trim((t0 * inv) %% p),
       g = .gf_trim((r0 * inv) %% p))
}

.gf_diff <- function(a, p) {
  n <- length(a) - 1L
  if (n <= 0L) return(numeric(0))
  .gf_trim(((a[-length(a)] %% p) * ((n:1) %% p)) %% p)
}

# a^e mod f where e is a bigz or number.
.gf_powmod <- function(a, e, f, p) {
  bits <- as.character(gmp::as.bigz(e), b = 2)
  res <- 1
  base <- .gf_rem(a, f, p)
  for (ch in rev(strsplit(bits, "")[[1L]])) {
    if (ch == "1") res <- .gf_rem(.gf_mul(res, base, p), f, p)
    base <- .gf_rem(.gf_mul(base, base, p), f, p)
  }
  .gf_trim(res)
}

.gf_from_bigz <- function(a, p) {
  .gf_trim(as.numeric(gmp::mod.bigz(gmp::as.bigz(a), p)))
}

.gf_sqf_p <- function(f, p) {
  if (length(f) <= 1L) return(TRUE)
  g <- .gf_gcd(f, .gf_diff(f, p), p)
  length(g) == 1L
}

.gf_rng <- new.env()
.gf_rng$state <- 12345

.gf_random <- function(n, p) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    .gf_rng$state <- (.gf_rng$state * 1103515245 + 12345) %% 2147483648
    out[i] <- floor(.gf_rng$state / 2147483648 * p)
  }
  out
}

.gf_ddf <- function(f, p) {
  out <- list()
  fs <- f
  g <- c(1, 0)
  i <- 1L
  x <- c(1, 0)
  while (2L * i <= length(fs) - 1L) {
    g <- .gf_powmod(g, p, fs, p)
    h <- .gf_gcd(fs, .gf_sub(g, x, p), p)
    if (length(h) > 1L) {
      out[[length(out) + 1L]] <- list(h, i)
      fs <- .gf_divmod(fs, h, p)$q
      g <- .gf_rem(g, fs, p)
    }
    i <- i + 1L
  }
  if (length(fs) > 1L) out[[length(out) + 1L]] <- list(fs, length(fs) - 1L)
  out
}

.gf_edf <- function(f, n, p) {
  d <- length(f) - 1L
  if (d <= n) return(list(f))
  e <- (gmp::as.bigz(p)^n - 1) / 2
  repeat {
    r <- .gf_trim(.gf_random(2L * n, p))
    if (length(r) < 2L) next
    h <- .gf_powmod(r, e, f, p)
    g <- .gf_gcd(f, .gf_sub(h, 1, p), p)
    if (length(g) > 1L && length(g) < length(f)) {
      q <- .gf_divmod(f, g, p)$q
      return(c(.gf_edf(g, n, p), .gf_edf(.gf_monic(q, p), n, p)))
    }
  }
}

# Factor a monic square-free polynomial over GF(p) into irreducibles.
.gf_factor_sqf <- function(f, p) {
  f <- .gf_monic(f, p)
  out <- list()
  for (dd in .gf_ddf(f, p)) out <- c(out, .gf_edf(dd[[1L]], dd[[2L]], p))
  out
}

.small_primes <- local({
  n <- 20000L
  is <- rep(TRUE, n)
  is[1L] <- FALSE
  for (i in 2:floor(sqrt(n))) if (is[i]) is[seq(i * i, n, by = i)] <- FALSE
  which(is)
})
