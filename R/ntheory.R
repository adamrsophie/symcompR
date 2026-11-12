# Number theory, after sympy/ntheory and the number theoretic Function classes
# of sympy/functions/combinatorial/numbers.py.
#
# Arguments pass through .nt_int(), which refuses doubles that may already be
# rounded. Internal work uses bigz, and doubles only where every product stays
# below 2^53. Results pass through .nt_out(): plain doubles while exactly
# representable, bigz above that, since bigz vectors walk raw bytes in for()
# and corrupt c() when they come second. SymPy Function classes such as
# totient() return sym integers and stay unevaluated on symbols.

.NT_Z53 <- gmp::as.bigz(2)^53L
.NT_FACTORIZE_MAX <- gmp::as.bigz("10000000000000000000000000000")

.z0 <- function() gmp::as.bigz(integer())

.nt_stop_int <- function(x) {
  stop(if (is_sym(x)) format(x) else paste(format(x), collapse = " "),
       " is not an integer", call. = FALSE)
}

.nt_int <- function(x) {
  if (is_sym(x)) {
    if (x$type != "Integer") .nt_stop_int(x)
    return(gmp::numerator(x$value))
  }
  if (length(x) != 1L) stop("expected a single integer; use vapply() for vectors", call. = FALSE)
  if (inherits(x, "bigz")) {
    if (is.na(x)) stop("integer argument is NA", call. = FALSE)
    attr(x, "mod") <- NULL
    return(x)
  }
  if (inherits(x, "bigq")) {
    if (is.na(x) || gmp::denominator(x) != 1) .nt_stop_int(x)
    return(gmp::numerator(x))
  }
  if (is.character(x)) {
    s <- sub("^\\+", "", trimws(x))
    if (!grepl("^-?[0-9]+$", s)) .nt_stop_int(x)
    return(.str_bigz(s))
  }
  if (is.numeric(x) && is.finite(x) && x == trunc(x)) {
    if (abs(x) >= 2^53)
      stop("doubles at or above 2^53 may be rounded; pass a string or bigz", call. = FALSE)
    return(gmp::as.bigz(x))
  }
  .nt_stop_int(x)
}

.zc <- function(lst) if (length(lst)) gmp::c_bigz(lst) else .z0()

.nt_ints <- function(x) {
  if (inherits(x, "bigz")) {
    attr(x, "mod") <- NULL
    if (anyNA(x)) stop("integer argument is NA", call. = FALSE)
    return(x)
  }
  .zc(lapply(as.list(x), .nt_int))
}

.nt_out <- function(x) {
  if (!length(x)) return(numeric())
  if (inherits(x, "bigz")) {
    if (all(abs(x) < .NT_Z53)) return(as.numeric(x))
    return(x)
  }
  x
}

.bitlen <- function(n) if (n == 0) 0L else gmp::sizeinbase(abs(n), 2L)

# Number of trailing zero bits of a positive integer.
.v2 <- function(n) {
  s <- as.character(n, b = 2L)
  nchar(s) - nchar(sub("0+$", "", s))
}

.isqrt <- function(n) .iroot(n, 2L)$root

.gcd2 <- function(a, b) gmp::gcd.bigz(a, b)

.zinv <- function(a, m) suppressWarnings(gmp::inv.bigz(a %% m, m))

.powm <- function(b, e, m) {
  if (e < 0) {
    b <- .zinv(b, m)
    if (is.na(b)) stop("base is not invertible for the given modulus", call. = FALSE)
    e <- -e
  }
  gmp::powm(b, e, m)
}

# remove(x, f): strip every factor f, as list(x, multiplicity).
.remove <- function(x, f) {
  if (f < 2) stop(.sympy_condition("ValueError", "factor must be > 1"))
  m <- 0
  if (x == 0) return(list(x, 0))
  repeat {
    if (x %% f != 0) break
    x <- x %/% f
    m <- m + 1
  }
  list(x, m)
}

.jacobi <- function(x, y) {
  if (y <= 0 || y %% 2 == 0) stop("y should be an odd positive integer", call. = FALSE)
  x <- x %% y
  if (x == 0) return(if (y == 1) 1L else 0L)
  if (y == 1 || x == 1) return(1L)
  if (.gcd2(x, y) != 1) return(0L)
  j <- 1L
  while (x != 0) {
    while (x > 0 && x %% 2 == 0) {
      x <- x %/% 2
      r <- y %% 8
      if (r == 3 || r == 5) j <- -j
    }
    tmp <- x
    x <- y
    y <- tmp
    if (x %% 4 == 3 && y %% 4 == 3) j <- -j
    x <- x %% y
  }
  j
}

.legendre <- function(x, y) {
  if (y <= 0 || y %% 2 == 0) stop("y should be an odd prime", call. = FALSE)
  x <- x %% y
  if (x == 0) return(0L)
  if (gmp::powm(x, (y - 1) %/% 2, y) == 1) 1L else -1L
}

.kronecker <- function(x, y) {
  if (.gcd2(x, y) != 1) return(0L)
  if (y == 0) return(1L)
  sgn <- if (y < 0 && x < 0) -1L else 1L
  y <- abs(y)
  s <- .v2(y)
  y <- y %/% gmp::as.bigz(2)^s
  r <- x %% 8
  if (s %% 2 == 1 && (r == 3 || r == 5)) sgn <- -sgn
  sgn * .jacobi(x, y)
}

.randint_gen <- function(draw30) {
  function(lo, hi) {
    span <- gmp::as.bigz(hi) - lo + 1
    bits <- .bitlen(span)
    mask <- gmp::as.bigz(2)^bits
    repeat {
      r <- gmp::as.bigz(0)
      for (i in seq_len(ceiling(bits / 30))) r <- r * 2^30 + draw30()
      r <- r %% mask
      if (r < span) return(lo + r)
    }
  }
}

.r_randint <- .randint_gen(function() floor(stats::runif(1) * 2^30))

# Seeded stream private to one call, so .Random.seed stays untouched.
.nt_rng <- function(seed) {
  if (is.null(seed)) return(.r_randint)
  state <- as.numeric(.nt_int(seed) %% 2147483646L) + 1
  .randint_gen(function() {
    state <<- (16807 * state) %% 2147483647
    floor(state / 2)
  })
}

# ---------------------------------------------------------------------------
# Sieve (sympy/ntheory/generate.py)

# Odd primes in (a, b) for even a, from the base primes. Blocks bound memory.
.sieve_segment <- function(a, b, base, interval) {
  if (b %% 2 == 1) b <- b - 1
  out <- list()
  odd <- base[-1L]
  while (a < b) {
    bs <- min(interval, (b - a) %/% 2)
    block <- rep(TRUE, bs)
    lim <- floor(sqrt(a + 2 * bs + 1))
    for (p in odd[odd <= lim]) {
      st <- if (a + 1 < p * p) (p * p - a - 1) / 2 else ((-(a + 1 + p)) %/% 2) %% p
      if (st < bs) block[seq.int(st + 1, bs, by = p)] <- FALSE
    }
    out[[length(out) + 1L]] <- a + 2 * (which(block) - 1) + 1
    a <- a + 2 * bs
  }
  if (length(out)) unlist(out, use.names = FALSE) else numeric()
}

.sieve_extend <- function(s, n) {
  lst <- s$.list
  num <- lst[length(lst)] + 1
  if (n < num) return(lst)
  num2 <- num^2
  while (num2 <= n) {
    lst <- c(lst, .sieve_segment(num, num2, lst, s$.interval))
    num <- num2
    num2 <- num2^2
  }
  lst <- c(lst, .sieve_segment(num, n + 1, lst, s$.interval))
  s$.list <- lst
  lst
}

.sieve_extend_to_no <- function(s, i) {
  lst <- s$.list
  while (length(lst) < i) lst <- .sieve_extend(s, floor(lst[length(lst)] * 1.5))
  lst
}

# Grow in doubling steps so a run of small requests does not rebuild each time.
.sieve_totients <- function(s, b) {
  n <- length(s$.tlist)
  if (b > n) {
    m <- max(b, 2 * n)
    t <- as.numeric(0:(m - 1))
    for (p in .sieve_extend(s, m - 1)) {
      if (p >= m) break
      idx <- seq.int(p, m - 1, by = p) + 1
      t[idx] <- t[idx] - t[idx] %/% p
    }
    s$.tlist <- t
  }
  s$.tlist
}

.sieve_mobius <- function(s, b) {
  n <- length(s$.mlist)
  if (b > n) {
    m <- max(b, 2 * n)
    mu <- rep(1, m)
    mu[1L] <- 0
    for (p in .sieve_extend(s, m - 1)) {
      if (p >= m) break
      idx <- seq.int(p, m - 1, by = p) + 1
      mu[idx] <- -mu[idx]
      if (p * p < m) mu[seq.int(p * p, m - 1, by = p * p) + 1] <- 0
    }
    s$.mlist <- mu
  }
  s$.mlist
}

.sieve_ceiling <- function(x) {
  if (is_sym(x)) x <- fcall("ceiling", list(x))
  if (is.double(x) && length(x) == 1L && is.finite(x)) x <- ceiling(x)
  as.numeric(.nt_int(x))
}

#' Sieve of Eratosthenes
#'
#' `Sieve()` creates a list of primes, totients and Moebius values that grows on
#' demand, and `sieve` is the global instance the number theory functions share.
#' `sieve$extend(n)` grows the primes to at least `n`, and
#' `sieve$extend_to_no(i)` to at least `i` primes. `sieve$primerange(a, b)`,
#' `sieve$totientrange(a, b)` and `sieve$mobiusrange(a, b)` return values for
#' the half open range. `sieve$search(n)` returns the indices of the primes
#' bounding `n`, equal if `n` is prime, and `sieve$contains(n)` tests membership.
#' `sieve[i]` is the `i`-th prime, counting from 1.
#'
#' @param sieve_interval Largest number of odd integers sieved in one block.
#' @return A `Sieve` environment.
#' @examples
#' sieve[1:5]
#' sieve$primerange(10, 30)
#' sieve$totientrange(7, 12)
#' sieve$search(25)
#' @export
Sieve <- function(sieve_interval = 1000000) {
  if (!is.numeric(sieve_interval) || length(sieve_interval) != 1L || sieve_interval <= 0)
    stop("sieve_interval should be a positive integer", call. = FALSE)
  s <- new.env(parent = emptyenv())
  s$.list <- c(2, 3, 5, 7, 11, 13)
  s$.tlist <- c(0, 1, 1, 2, 2, 4)
  s$.mlist <- c(0, 1, -1, -1, 0, -1)
  s$.interval <- sieve_interval
  s$extend <- function(n) {
    .sieve_extend(s, trunc(as.numeric(n)))
    invisible(s)
  }
  s$extend_to_no <- function(i) {
    .sieve_extend_to_no(s, as.numeric(.nt_int(i)))
    invisible(s)
  }
  s$primerange <- function(a, b = NULL) {
    if (is.null(b)) {
      b <- .sieve_ceiling(a)
      a <- 2
    } else {
      a <- max(2, .sieve_ceiling(a))
      b <- .sieve_ceiling(b)
    }
    if (a >= b) return(numeric())
    lst <- .sieve_extend(s, b)
    lst[lst >= a & lst < b]
  }
  s$totientrange <- function(a, b) {
    a <- max(1, .sieve_ceiling(a))
    b <- .sieve_ceiling(b)
    if (a >= b) return(numeric())
    .sieve_totients(s, b)[seq.int(a, b - 1) + 1]
  }
  s$mobiusrange <- function(a, b) {
    a <- max(1, .sieve_ceiling(a))
    b <- .sieve_ceiling(b)
    if (a >= b) return(numeric())
    .sieve_mobius(s, b)[seq.int(a, b - 1) + 1]
  }
  s$search <- function(n) {
    n <- as.numeric(.nt_int(n))
    if (n < 2) stop("n should be >= 2 but got: ", n, call. = FALSE)
    lst <- .sieve_extend(s, n)
    b <- findInterval(n, lst)
    if (lst[b] == n) c(b, b) else c(b, b + 1)
  }
  s$contains <- function(n) {
    if (!is.numeric(n) && !inherits(n, c("bigz", "sym"))) return(FALSE)
    n <- tryCatch(.nt_int(n), error = function(e) NULL)
    if (is.null(n) || n < 2) return(FALSE)
    if (n %% 2 == 0) return(n == 2)
    # ponytail: beyond 2^31 test primality instead of extending the sieve
    if (n > .Machine$integer.max) return(.is_prime(n))
    n <- as.numeric(n)
    r <- s$search(n)
    r[1L] == r[2L]
  }
  s$reset <- function(prime = NULL, totient = NULL, mobius = NULL) {
    if (is.null(prime) && is.null(totient) && is.null(mobius)) prime <- totient <- mobius <- TRUE
    if (isTRUE(prime)) s$.list <- s$.list[1:6]
    if (isTRUE(totient)) s$.tlist <- s$.tlist[1:6]
    if (isTRUE(mobius)) s$.mlist <- s$.mlist[1:6]
    invisible(s)
  }
  class(s) <- "Sieve"
  s
}

#' @export
`[.Sieve` <- function(x, i) {
  i <- as.numeric(i)
  if (any(i < 1)) stop("Sieve indices start at 1.", call. = FALSE)
  .sieve_extend_to_no(x, max(i))[i]
}

#' @export
print.Sieve <- function(x, ...) {
  show <- function(v, name) {
    n <- length(v)
    sprintf("%s sieve (%d): %s, ... %s", name, n,
            paste(format(v[seq_len(min(3L, n))], trim = TRUE), collapse = ", "),
            paste(format(v[max(1L, n - 1L):n], trim = TRUE), collapse = ", "))
  }
  cat(show(x$.list, "prime"), "\n", show(x$.tlist, "totient"), "\n",
      show(x$.mlist, "mobius"), "\n", sep = "")
  invisible(x)
}

#' @rdname Sieve
#' @format `sieve` is the shared `Sieve` instance.
#' @export
sieve <- Sieve()

# ---------------------------------------------------------------------------
# Primality (sympy/ntheory/primetest.py, sympy/external/ntheory.py)

.MERSENNE_EXPS <- c(2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279,
                    2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701,
                    23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
                    1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
                    20996011, 24036583, 25964951, 30402457, 32582657, 37156667,
                    42643801, 43112609, 57885161, 74207281, 77232917, 82589933,
                    136279841)

.is_prime <- function(n) n >= 2 && gmp::isprime(n, 40L) > 0

.check_prp_args <- function(n, a, fn) {
  if (a < 2) stop(fn, "() requires 'a' greater than or equal to 2", call. = FALSE)
  if (n < 1) stop(fn, "() requires 'n' be greater than 0", call. = FALSE)
}

# U, V and Q^k of the Lucas sequence with parameters P, Q at index k, mod n.
.lucas_sequence <- function(n, P, Q, k) {
  if (k == 0) return(list(U = gmp::as.bigz(0), V = gmp::as.bigz(2), Qk = gmp::as.bigz(1)))
  D <- P^2 - 4 * Q
  U <- gmp::as.bigz(1)
  V <- gmp::as.bigz(P)
  Qk <- gmp::as.bigz(Q) %% n
  bits <- strsplit(as.character(gmp::as.bigz(k), b = 2L), "", fixed = TRUE)[[1L]][-1L]
  for (b in bits) {
    U <- (U * V) %% n
    V <- (V * V - 2 * Qk) %% n
    Qk <- (Qk * Qk) %% n
    if (b == "1") {
      Un <- U * P + V
      Vn <- V * P + U * D
      if (Un %% 2 == 1) Un <- Un + n
      if (Vn %% 2 == 1) Vn <- Vn + n
      U <- Un %/% 2
      V <- Vn %/% 2
      Qk <- (Qk * Q) %% n
    }
  }
  list(U = U %% n, V = V %% n, Qk = Qk)
}

.selfridge_params <- function(n) {
  D <- 5
  repeat {
    if (D >= 1000000) stop("appropriate value for D cannot be found", call. = FALSE)
    d <- if (D %% 4 == 3) -D else D
    j <- .jacobi(gmp::as.bigz(d), n)
    if (j == -1L) return(d)
    if (j == 0L && gmp::as.bigz(d) %% n != 0) return(NULL)
    if (D == 13 && .iroot(n, 2L)$exact) return(NULL)
    D <- D + 2
  }
}

.is_selfridge_prp <- function(n) {
  if (n < 1) stop("is_selfridge_prp() requires 'n' be greater than 0", call. = FALSE)
  if (n == 1) return(FALSE)
  if (n %% 2 == 0) return(n == 2)
  d <- .selfridge_params(n)
  if (is.null(d)) return(FALSE)
  .lucas_sequence(n, 1, (1 - d) %/% 4, n + 1)$U == 0
}

.is_strong_selfridge_prp <- function(n) {
  if (n < 1) stop("is_strong_selfridge_prp() requires 'n' be greater than 0", call. = FALSE)
  if (n == 1) return(FALSE)
  if (n %% 2 == 0) return(n == 2)
  d <- .selfridge_params(n)
  if (is.null(d)) return(FALSE)
  s <- .v2(n + 1)
  l <- .lucas_sequence(n, 1, (1 - d) %/% 4, (n + 1) %/% gmp::as.bigz(2)^s)
  if (l$U == 0 || l$V == 0) return(TRUE)
  V <- l$V
  Qk <- l$Qk
  for (i in seq_len(s - 1L)) {
    V <- (V * V - 2 * Qk) %% n
    if (V == 0) return(TRUE)
    Qk <- gmp::powm(Qk, 2L, n)
  }
  FALSE
}

.mr_test <- function(n, base, s, t) {
  b <- gmp::powm(base, t, n)
  if (b == 1 || b == n - 1) return(TRUE)
  for (i in seq_len(s - 1L)) {
    b <- gmp::powm(b, 2L, n)
    if (b == n - 1) return(TRUE)
    if (b == 1) return(FALSE)
  }
  FALSE
}

#' Primality tests
#'
#' `isprime()` tests one integer using GMP (Baillie-PSW plus Miller-Rabin
#' rounds), so its answer is certain below `2^64` and a probable prime above,
#' as in SymPy. `is_gaussian_prime()` tests a Gaussian integer. The remaining
#' functions are the individual pseudoprime tests: Fermat, Euler and
#' Euler-Jacobi to base `a`, Miller-Rabin `mr()` over `bases`, the Selfridge
#' Lucas tests, `proth_test()` for Proth numbers and `is_mersenne_prime()`.
#' [is_square()] tests for a perfect square when given an integer.
#'
#' Integer arguments may be numbers, `gmp` big integers, decimal strings or
#' `sym` integers. Doubles at or above `2^53` are refused, because they may have
#' lost digits already; pass a string or big integer.
#'
#' @param n An integer.
#' @param a Base of the test.
#' @param bases A vector or list of integer bases.
#' @param num A Gaussian integer: an R complex number, a `sym` expression or an
#'   integer.
#' @return `TRUE` or `FALSE`.
#' @examples
#' isprime(97)
#' isprime("170141183460469231731687303715884105727")
#' is_gaussian_prime(complex(real = 2, imaginary = 3))
#' mr(1373651, c(2, 3))
#' @export
isprime <- function(n) .is_prime(.nt_int(n))

#' @rdname isprime
#' @export
is_gaussian_prime <- function(num) {
  num <- as_sym(num)
  a <- .nt_int(Integer(re(num)))
  b <- .nt_int(Integer(im(num)))
  if (a == 0) {
    b <- abs(b)
    return(.is_prime(b) && b %% 4 == 3)
  }
  if (b == 0) {
    a <- abs(a)
    return(.is_prime(a) && a %% 4 == 3)
  }
  .is_prime(a^2 + b^2)
}

#' @rdname isprime
#' @export
is_fermat_pseudoprime <- function(n, a) {
  n <- .nt_int(n)
  a <- .nt_int(a)
  if (a == 1) return(n == 2 || n %% 2 == 1)
  .check_prp_args(n, a, "is_fermat_prp")
  if (n == 1) return(FALSE)
  if (n %% 2 == 0) return(n == 2)
  a <- a %% n
  if (.gcd2(n, a) != 1) stop("is_fermat_prp() requires gcd(n,a) == 1", call. = FALSE)
  gmp::powm(a, n - 1, n) == 1
}

#' @rdname isprime
#' @export
is_euler_pseudoprime <- function(n, a) {
  n <- .nt_int(n)
  a <- .nt_int(a)
  if (a < 1) stop("a should be an integer greater than 0", call. = FALSE)
  if (n < 1) stop("n should be an integer greater than 0", call. = FALSE)
  if (n == 1) return(FALSE)
  if (a == 1) return(n == 2 || n %% 2 == 1)
  if (n %% 2 == 0) return(n == 2)
  if (.gcd2(n, a) != 1) stop(.sympy_condition("ValueError", "The two numbers should be relatively prime"))
  v <- gmp::powm(a, (n - 1) %/% 2, n)
  v == 1 || v == n - 1
}

#' @rdname isprime
#' @export
is_euler_jacobi_pseudoprime <- function(n, a) {
  n <- .nt_int(n)
  a <- .nt_int(a)
  if (a == 1) return(n == 2 || n %% 2 == 1)
  .check_prp_args(n, a, "is_euler_prp")
  if (n == 1) return(FALSE)
  if (n %% 2 == 0) return(n == 2)
  a <- a %% n
  if (.gcd2(n, a) != 1) stop("is_euler_prp() requires gcd(n,a) == 1", call. = FALSE)
  gmp::powm(a, n %/% 2, n) == .jacobi(a, n) %% n
}

.is_square_int <- function(n) {
  n <- .nt_int(n)
  if (n < 0) return(FALSE)
  .iroot(n, 2L)$exact
}

#' @rdname isprime
#' @export
mr <- function(n, bases) {
  n <- .nt_int(n)
  if (n < 2 || (n > 2 && n %% 2 == 0)) return(FALSE)
  s <- .v2(n - 1)
  t <- n %/% gmp::as.bigz(2)^s
  b <- .nt_ints(bases)
  for (i in seq_along(b)) {
    base <- b[i]
    if (base >= n) base <- base %% n
    if (base >= 2 && !.mr_test(n, base, s, t)) return(FALSE)
  }
  TRUE
}

#' @rdname isprime
#' @export
is_lucas_prp <- function(n) {
  n <- .nt_int(n)
  if (n < 2) return(FALSE)
  .is_selfridge_prp(n)
}

#' @rdname isprime
#' @export
is_strong_lucas_prp <- function(n) {
  n <- .nt_int(n)
  if (n < 2) return(FALSE)
  .is_strong_selfridge_prp(n)
}

#' @rdname isprime
#' @export
is_extra_strong_lucas_prp <- function(n) {
  n <- .nt_int(n)
  if (n == 2) return(TRUE)
  if (n < 2 || n %% 2 == 0) return(FALSE)
  if (.iroot(n, 2L)$exact) return(FALSE)
  P <- 3
  repeat {
    D <- P^2 - 4
    j <- .jacobi(gmp::as.bigz(D), n)
    if (j == -1L) break
    if (j == 0L && gmp::as.bigz(D) %% n != 0) return(FALSE)
    P <- P + 1
  }
  s <- .v2(n + 1)
  k <- (n + 1) %/% gmp::as.bigz(2)^s
  l <- .lucas_sequence(n, P, 1, k)
  if (l$U == 0 && (l$V == 2 || l$V == n - 2)) return(TRUE)
  V <- l$V
  for (i in seq_len(s - 1L)) {
    if (V == 0) return(TRUE)
    V <- (V * V - 2) %% n
  }
  FALSE
}

#' @rdname isprime
#' @export
proth_test <- function(n) {
  n <- .nt_int(n)
  if (n < 3) stop(.sympy_condition("ValueError", "n is not Proth number"))
  m <- .v2(n - 1)
  k <- n %/% gmp::as.bigz(2)^m
  if (m < .bitlen(k)) stop(.sympy_condition("ValueError", "n is not Proth number"))
  if (n %% 3 == 0) return(n == 3)
  half <- n %/% 2
  if (k %% 3 != 0) return(gmp::powm(gmp::as.bigz(3), half, n) == n - 1)
  if (.iroot(n, 2L)$exact) return(FALSE)
  a <- gmp::as.bigz(5)
  while (a < n) {
    j <- .jacobi(a, n)
    if (j == -1L) return(gmp::powm(a, half, n) == n - 1)
    if (j == 0L) return(FALSE)
    a <- a + 1
  }
  FALSE
}

#' @rdname isprime
#' @export
is_mersenne_prime <- function(n) {
  n <- .nt_int(n)
  if (n < 1) return(FALSE)
  p <- .bitlen(n)
  if (gmp::as.bigz(2)^p - 1 != n) return(FALSE)
  if (p %in% .MERSENNE_EXPS) return(TRUE)
  if (p < 65000000 || !.is_prime(gmp::as.bigz(p))) return(FALSE)
  # ponytail: GIMPS has verified every exponent below 65 million; an R
  # Lucas-Lehmer run above that is out of reach, so this mirrors the table.
  stop("Lucas-Lehmer test for exponent ", p, " is not practical in R", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Generating primes (sympy/ntheory/generate.py)

# Logarithmic integral, from the convergent series, for bracketing the n-th prime.
.li <- function(x) {
  lx <- log(x)
  s <- 0
  term <- 1
  k <- 1
  repeat {
    term <- term * lx / k
    add <- term / k
    s <- s + add
    if (k > lx && add < s * 1e-17) break
    k <- k + 1
  }
  -digamma(1) + log(lx) + s
}

# Lucy-style prime counting; both updates use the previous pass only, so they
# vectorize directly.
.primepi_count <- function(n) {
  if (n < 2) return(0)
  lst <- sieve$.list
  if (n <= lst[length(lst)]) return(findInterval(n, lst))
  lim <- floor(sqrt(n))
  while (lim * lim > n) lim <- lim - 1
  while ((lim + 1) * (lim + 1) <= n) lim <- lim + 1
  arr1 <- (0:lim + 1) %/% 2
  arr2 <- c(0, (n %/% seq_len(lim) + 1) %/% 2)
  skip <- logical(lim + 1)
  if (lim >= 3) for (i in seq(3, lim, by = 2)) {
    if (skip[i + 1L]) next
    p <- arr1[i]
    skip[seq.int(i, lim, by = i) + 1L] <- TRUE
    up <- min(n %/% (i * i), lim)
    if (up >= 1) {
      js <- seq(1, up, by = 2)
      js <- js[!skip[js + 1L]]
      if (length(js)) {
        st <- i * js
        inner <- st <= lim
        d <- numeric(length(js))
        d[inner] <- arr2[st[inner] + 1L] - p
        d[!inner] <- arr1[n %/% st[!inner] + 1L] - p
        arr2[js + 1L] <- arr2[js + 1L] - d
      }
    }
    if (i * i <= lim) {
      jj <- seq(i * i, lim)
      arr1[jj + 1L] <- arr1[jj + 1L] - (arr1[jj %/% i + 1L] - p)
    }
  }
  arr2[2L]
}

.nextprime <- function(n, ith = 1) {
  if (ith <= 0) stop("ith should be positive", call. = FALSE)
  if (n < 2) n <- gmp::as.bigz(1)
  for (k in seq_len(ith)) n <- gmp::nextprime(n)
  n
}

.prime_nth <- function(n) {
  if (n < 1) stop("nth must be a positive integer; prime(1) == 2", call. = FALSE)
  lst <- sieve$.list
  if (n <= length(lst)) return(lst[n])
  if (n < 1000) return(.sieve_extend(sieve, 8 * n)[n])
  a <- 2
  b <- trunc(n * (log(n) + log(log(n))))
  while (a < b) {
    mid <- (a + b) %/% 2
    if (.li(mid) > n) b <- mid else a <- mid + 1
  }
  .nt_out(.nextprime(gmp::as.bigz(a - 1), n - .primepi_count(a - 1)))
}

#' Generating primes
#'
#' `prime(nth)` is the `nth` prime, `nextprime(n, ith)` the `ith` prime after
#' `n`, `prevprime(n)` the largest prime below `n`, `primerange(a, b)` the
#' primes in `[a, b)` (with one argument, the primes below `a`), and
#' `randprime(a, b)` a random prime in `[a, b)` using the R random number
#' generator, so [set.seed()] applies. `primorial(n, nth)` is the product of the
#' first `n` primes, or with `nth = FALSE` of the primes up to `n`.
#' `composite(nth)` is the `nth` composite and `compositepi(n)` counts the
#' composites up to `n`. `cycle_length(f, x0)` runs Brent's cycle finder on the
#' sequence `x0, f(x0), ...`.
#'
#' Integer results are plain doubles while they are exact (below `2^53`) and
#' `gmp` big integers above that.
#'
#' @param nth,ith Positive index.
#' @param n,a,b Integers. `n`, `a` and `b` of `prevprime()`, `primerange()` and
#'   `randprime()` may be numbers or `gmp` big integers; see [isprime()].
#' @param f A function of one argument.
#' @param x0 Starting value.
#' @param nmax Stop after this many iterations; `NULL` runs to the cycle.
#' @param values If `TRUE`, return the list of visited values instead.
#' @return An integer or vector, as described above. `cycle_length()` returns
#'   `list(lambda, mu)` with `mu` missing when `nmax` stopped the search, or the
#'   list of values. `randprime()` returns `NULL` when `a >= b`.
#' @examples
#' prime(25)
#' nextprime(100, 2)
#' prevprime(100)
#' primerange(2, 30)
#' primorial(5)
#' composite(10)
#' cycle_length(function(x) (x^2 + 1) %% 7, 2)
#' @export
prime <- function(nth) .nt_out(.prime_nth(as.numeric(.nt_int(nth))))

#' @rdname prime
#' @export
nextprime <- function(n, ith = 1) .nt_out(.nextprime(.nt_int(n), as.numeric(.nt_int(ith))))

#' @rdname prime
#' @export
prevprime <- function(n) {
  q <- .nt_rat(n)
  n <- if (!is.null(q)) -((-gmp::numerator(q)) %/% gmp::denominator(q)) else
    if (is.double(n) && length(n) == 1L && is.finite(n)) .nt_int(ceiling(n)) else .nt_int(n)
  if (n < 3) stop(.sympy_condition("ValueError", "no preceding primes"))
  n <- n - 1
  if (n == 2) return(2)
  if (n %% 2 == 0) n <- n - 1
  while (!.is_prime(n)) n <- n - 2
  .nt_out(n)
}

.primerange <- function(a, b) {
  if (a >= b) return(numeric())
  if (b > .NT_Z53) {
    out <- list()
    p <- .nextprime(a - 1)
    while (p < b) {
      out[[length(out) + 1L]] <- p
      p <- .nextprime(p)
    }
    return(.nt_out(.zc(out)))
  }
  a <- max(as.numeric(a), 2)
  b <- as.numeric(b)
  if (b <= 1e7) {
    lst <- .sieve_extend(sieve, b)
    return(lst[lst >= a & lst < b])
  }
  base <- .sieve_extend(sieve, ceiling(sqrt(b)))
  a0 <- if (a %% 2 == 1) a - 1 else a
  c(if (a == 2) 2, .sieve_segment(a0, b, base, sieve$.interval))
}

#' @rdname prime
#' @export
primerange <- function(a, b = NULL) {
  a <- .nt_int(a)
  if (is.null(b)) {
    b <- a
    a <- gmp::as.bigz(2)
  } else {
    b <- .nt_int(b)
  }
  .primerange(a, b)
}

#' @rdname prime
#' @export
randprime <- function(a, b) {
  a <- .nt_int(a)
  b <- .nt_int(b)
  if (a >= b) return(NULL)
  p <- .nextprime(.r_randint(a - 1, b))
  if (p >= b) p <- .nt_int(prevprime(b))
  if (p < a) stop("no primes exist in the specified range", call. = FALSE)
  .nt_out(p)
}

#' @rdname prime
#' @export
primorial <- function(n, nth = TRUE) {
  n <- if (nth) .nt_int(n) else .nt_int(trunc(as.numeric(n)))
  if (n < 1) stop("primorial argument must be >= 1", call. = FALSE)
  ps <- if (nth) .sieve_extend_to_no(sieve, as.numeric(n))[seq_len(as.numeric(n))] else .primerange(gmp::as.bigz(2), n + 1)
  .nt_out(prod(gmp::as.bigz(ps)))
}

#' @rdname prime
#' @export
cycle_length <- function(f, x0, nmax = NULL, values = FALSE) {
  nmax <- if (is.null(nmax)) 0 else as.numeric(nmax)
  power <- 1
  lam <- 1
  tortoise <- x0
  hare <- f(x0)
  i <- 1
  seen <- if (values) list(tortoise) else NULL
  while (tortoise != hare && (nmax == 0 || i < nmax)) {
    i <- i + 1
    if (power == lam) {
      tortoise <- hare
      power <- power * 2
      lam <- 0
    }
    if (values) seen[[length(seen) + 1L]] <- hare
    hare <- f(hare)
    lam <- lam + 1
  }
  if (nmax != 0 && i == nmax)
    return(if (values) seen else list(lambda = nmax, mu = NA_real_))
  if (values) return(seen)
  mu <- 0
  tortoise <- hare <- x0
  for (k in seq_len(lam)) hare <- f(hare)
  while (tortoise != hare) {
    tortoise <- f(tortoise)
    hare <- f(hare)
    mu <- mu + 1
  }
  list(lambda = lam, mu = mu)
}

#' @rdname prime
#' @export
composite <- function(nth) {
  n <- as.numeric(.nt_int(nth))
  if (n < 1) stop("nth must be a positive integer; composite(1) == 4", call. = FALSE)
  arr <- c(4, 6, 8, 9, 10, 12, 14, 15, 16, 18)
  if (n <= 10) return(arr[n])
  lst <- sieve$.list
  a <- 4
  b <- lst[length(lst)]
  if (n <= b - .primepi_count(b) - 1) {
    while (a < b - 1) {
      mid <- (a + b) %/% 2
      if (mid - .primepi_count(mid) - 1 > n) b <- mid else a <- mid
    }
    if (.is_prime(gmp::as.bigz(a))) a <- a - 1
    return(a)
  }
  a <- 4
  b <- trunc(n * (log(n) + log(log(n))))
  while (a < b) {
    mid <- (a + b) %/% 2
    if (mid - .li(mid) - 1 > n) b <- mid else a <- mid + 1
  }
  ncomp <- a - .primepi_count(a) - 1
  while (ncomp > n) {
    if (!.is_prime(gmp::as.bigz(a))) ncomp <- ncomp - 1
    a <- a - 1
  }
  if (.is_prime(gmp::as.bigz(a))) a <- a - 1
  a
}

#' @rdname prime
#' @export
compositepi <- function(n) {
  n <- if (is.double(n) && length(n) == 1L && is.finite(n)) .nt_int(trunc(n)) else .nt_int(n)
  if (n < 4) return(0)
  n <- as.numeric(n)
  n - .primepi_count(n) - 1
}
# ---------------------------------------------------------------------------
# Factoring (sympy/ntheory/factor_.py, sympy/ntheory/ecm.py)

.gcd_num <- function(a, b) {
  while (b) {
    t <- a %% b
    a <- b
    b <- t
  }
  a
}

# Ordering of a bigz vector without passing through doubles.
.zorder <- function(v) {
  if (!length(v)) return(integer())
  s <- as.character(abs(v))
  sg <- ifelse(v < 0, -1, ifelse(v == 0, 0, 1))
  order(sg, nchar(s) * sg, ifelse(sg < 0, chartr("0123456789", "9876543210", s), s),
        method = "radix")
}

.zsort <- function(v) v[.zorder(v)]

# Factor dictionaries are environments keyed by decimal strings while being
# built, and named exponent vectors once returned.
.fac_new <- function() new.env(parent = emptyenv())
.fac_key <- function(p) as.character(gmp::as.bigz(p))
.fac_get <- function(fac, p) get0(.fac_key(p), fac, inherits = FALSE, ifnotfound = 0)
.fac_set <- function(fac, p, e) assign(.fac_key(p), e, envir = fac)
.fac_add <- function(fac, p, e) .fac_set(fac, p, .fac_get(fac, p) + e)
.fac_merge <- function(fac, d, mult = 1) {
  for (k in names(d)) .fac_add(fac, k, d[[k]] * mult)
}
.fac_vec <- function(fac) {
  k <- ls(fac, all.names = TRUE, sorted = FALSE)
  if (!length(k)) return(stats::setNames(numeric(), character()))
  o <- .zorder(gmp::as.bigz(k))
  stats::setNames(vapply(k[o], function(z) get(z, fac), 0, USE.NAMES = FALSE), k[o])
}
.fac_primes <- function(d) gmp::as.bigz(names(d))

.vmsg <- function(verbose, ...) if (verbose) message(sprintf(...))

# ---------------------------------------------------------------------------
# FactorCache

#' Cache of prime factors
#'
#' `factor_cache` remembers one prime factor (the largest seen) for integers
#' that [factorint()] has split, and [factorint()] consults it first.
#' `FactorCache(maxsize)` creates another cache; the oldest entries are dropped
#' beyond `maxsize`. Methods: `get(n)` returns a stored factor or `NULL` (primes
#' within the sieve return themselves), `set(n, factor)` stores a prime factor,
#' `add(n, factors)` stores every prime factor of `n` given in `factors`,
#' `remove(n)` drops an entry, `cache_clear()` empties the cache,
#' `set_maxsize(m)` changes the bound and `size()` counts entries.
#'
#' @param maxsize Maximum number of entries, or `NULL` for no bound.
#' @return A `FactorCache` environment.
#' @examples
#' fc <- FactorCache(10)
#' fc$set(35, 7)
#' fc$get(35)
#' @export
FactorCache <- function(maxsize = NULL) {
  fc <- new.env(parent = emptyenv())
  fc$.store <- new.env(parent = emptyenv())
  fc$.keys <- character()
  fc$.maxsize <- NULL
  trim <- function() {
    m <- fc$.maxsize
    while (!is.null(m) && length(fc$.keys) > m) {
      rm(list = fc$.keys[1L], envir = fc$.store)
      fc$.keys <- fc$.keys[-1L]
    }
  }
  fc$set_maxsize <- function(value) {
    if (!is.null(value) && value <= 0)
      stop("maxsize must be None or a non-negative integer.", call. = FALSE)
    fc$.maxsize <- value
    trim()
    invisible(fc)
  }
  fc$set <- function(n, factor) {
    n <- .nt_int(n)
    factor <- .nt_int(factor)
    if (!(factor > 1 && factor <= n && n %% factor == 0 && .is_prime(factor)))
      stop(factor, " is not a prime factor of ", n, call. = FALSE)
    k <- as.character(n)
    old <- get0(k, fc$.store, inherits = FALSE)
    if (is.null(old)) fc$.keys <- c(fc$.keys, k) else if (old > factor) factor <- old
    assign(k, factor, envir = fc$.store)
    trim()
    invisible(fc)
  }
  fc$get <- function(n, default = NULL) {
    n <- .nt_int(n)
    lst <- sieve$.list
    if (n <= lst[length(lst)] && n >= 2) {
      if (lst[findInterval(as.numeric(n), lst)] == n) return(.nt_out(n))
    }
    k <- as.character(n)
    v <- get0(k, fc$.store, inherits = FALSE)
    if (is.null(v)) return(default)
    fc$.keys <- c(fc$.keys[fc$.keys != k], k)
    .nt_out(v)
  }
  fc$add <- function(n, factors) {
    n <- .nt_int(n)
    fs <- .zsort(.nt_ints(factors))
    for (i in rev(seq_along(fs))) {
      fc$set(n, fs[i])
      n <- .remove(n, fs[i])[[1L]]
    }
    invisible(fc)
  }
  fc$remove <- function(n) {
    k <- as.character(.nt_int(n))
    if (!k %in% fc$.keys) stop(k, " does not exist.", call. = FALSE)
    rm(list = k, envir = fc$.store)
    fc$.keys <- fc$.keys[fc$.keys != k]
    invisible(fc)
  }
  fc$cache_clear <- function() {
    rm(list = ls(fc$.store, all.names = TRUE), envir = fc$.store)
    fc$.keys <- character()
    invisible(fc)
  }
  fc$size <- function() length(fc$.keys)
  fc$set_maxsize(maxsize)
  class(fc) <- "FactorCache"
  fc
}

#' @export
print.FactorCache <- function(x, ...) {
  cat("<FactorCache: ", x$size(), " entries, maxsize ",
      if (is.null(x$.maxsize)) "unbounded" else x$.maxsize, ">\n", sep = "")
  invisible(x)
}

#' @rdname FactorCache
#' @format `factor_cache` is the shared cache, bounded at 1000 entries.
#' @export
factor_cache <- FactorCache(maxsize = 1000)

# Silent internal store: factors found while factoring are true by construction.
.fc_put <- function(n, p) {
  if (n > p && n %% p == 0) tryCatch(factor_cache$set(n, p), error = function(e) NULL)
}

.fc_lookup <- function(n) {
  if (n < 2) return(NULL)
  v <- factor_cache$get(n)
  if (is.null(v)) NULL else gmp::as.bigz(v)
}

# ---------------------------------------------------------------------------
# Small helpers

#' Multiplicities of factors
#'
#' `multiplicity(p, n)` is the largest `m` with `p^m` dividing `n`; rational
#' `p` and `n` (as `sym` rationals or `gmp` big rationals) give the exponent of
#' `p` in `n`, possibly negative, and an unevaluated `sym_factorial(k)` for `n` is
#' handled by `multiplicity_in_factorial(p, k)`, the exponent of `p` in `k!`.
#'
#' @param p Base: an integer greater than 1, or a rational. For an
#'   [ImplicitRegion()], the multiplicity of the point `n` on the region.
#' @param n An integer or rational, or a point of an implicit region.
#' @return An integer.
#' @examples
#' multiplicity(5, 8000)
#' multiplicity(S(3), Rational(1, 9))
#' multiplicity_in_factorial(2, 100)
#' @export
multiplicity <- function(p, n) {
  if (inherits(p, "ImplicitRegion")) return(.ir_multiplicity(p, n))
  pq <- .nt_rat(p)
  nq <- .nt_rat(n)
  if (is.null(pq) && is.null(nq)) {
    if (is_sym(n) && n$type == "Function" && n$name == "factorial" &&
        n$args[[1L]]$type == "Integer" && n$args[[1L]]$value >= 0)
      return(multiplicity_in_factorial(p, n$args[[1L]]))
    p <- .nt_int(p)
    n <- .nt_int(n)
    if (n == 0) stop(.sympy_condition("ValueError", paste0("no such integer exists: multiplicity of ", n, " is not-defined")))
    return(.remove(n, p)[[2L]])
  }
  if (is.null(pq)) pq <- gmp::as.bigq(.nt_int(p))
  if (is.null(nq)) nq <- gmp::as.bigq(.nt_int(n))
  pp <- gmp::numerator(pq)
  pd <- gmp::denominator(pq)
  np <- gmp::numerator(nq)
  nd <- gmp::denominator(nq)
  m <- function(a, b) .remove(b, a)[[2L]]
  if (pd == 1) {
    if (np == 1) return(-m(pp, nd))
    return(m(pp, np) - m(pp, nd))
  }
  if (pp == 1) return(m(pd, nd))
  min(m(pp, np), m(pd, nd)) - min(m(pd, np), m(pp, nd))
}

# bigq of a non-integer rational argument, or NULL.
.nt_rat <- function(x) {
  if (is_sym(x)) return(if (x$type == "Rational") x$value else NULL)
  if (inherits(x, "bigq") && length(x) == 1L && gmp::denominator(x) != 1) return(x)
  if (is.character(x) && length(x) == 1L && grepl("^\\s*[-+]?[0-9]+\\s*/\\s*[0-9]+\\s*$", x)) {
    q <- Rational(gsub("\\s", "", x))
    return(if (q$type == "Rational") q$value else NULL)
  }
  NULL
}

#' @rdname multiplicity
#' @export
multiplicity_in_factorial <- function(p, n) {
  p <- .nt_int(p)
  n <- .nt_int(n)
  if (p <= 0) stop("expecting positive integer got ", p, call. = FALSE)
  if (n < 0) stop("expecting non-negative integer got ", n, call. = FALSE)
  f <- .factorint(p)
  if (!length(f)) stop("expecting an integer greater than 1, got ", p, call. = FALSE)
  best <- list()
  for (k in names(f)) {
    v <- as.character(f[[k]])
    if (is.null(best[[v]]) || gmp::as.bigz(k) > best[[v]]) best[[v]] <- gmp::as.bigz(k)
  }
  res <- NULL
  for (v in names(best)) {
    k <- best[[v]]
    ds <- gmp::as.bigz(0)
    m <- n
    while (m > 0) {
      ds <- ds + m %% k
      m <- m %/% k
    }
    s <- (n - ds) %/% (k - 1) %/% as.numeric(v)
    if (is.null(res) || s < res) res <- s
  }
  .nt_out(res)
}

#' Smoothness of an integer
#'
#' `smoothness(n)` returns `c(B, Bpow)`: the largest prime factor of `n` and the
#' largest prime power dividing it. `smoothness_p(n, m, power, visual)` reports
#' the smoothness of `p + m` for every prime power `p^i` in `n`, sorted by `B`
#' (`power = 0`) or `Bpow` (`power = 1`). It returns `list(m, entries)` where
#' each entry is `c(p, i, B, Bpow)`; with `visual = TRUE` it returns the
#' report as a string, which it also accepts back as input.
#'
#' @param n An integer, a factor vector, a report string or a returned list.
#' @param m Offset added to each prime.
#' @param power 0 or 1, the sort key.
#' @param visual `TRUE` for the string report, `FALSE` for the list.
#' @return As described above.
#' @examples
#' smoothness(2^7 * 3^2)
#' smoothness_p(10431, m = 1)
#' cat(smoothness_p(10431, m = 1, visual = TRUE))
#' @export
smoothness <- function(n) {
  n <- .nt_int(n)
  if (n == 1) return(c(1, 1))
  f <- .factorint(n)
  ps <- .fac_primes(f)
  pw <- .zc(lapply(seq_along(ps), function(i) ps[i]^f[[i]]))
  .nt_out(.zc(list(max(ps), max(pw))))
}

#' @rdname smoothness
#' @export
smoothness_p <- function(n, m = -1, power = 0, visual = NULL) {
  if (identical(visual, 1) || identical(visual, 0)) visual <- as.logical(visual)
  if (!isTRUE(visual) && !isFALSE(visual)) visual <- NULL
  if (is.character(n)) {
    if (isTRUE(visual)) return(n)
    d <- list()
    for (li in strsplit(n, "\n", fixed = TRUE)[[1L]]) {
      kv <- strsplit(strsplit(strsplit(li, "has", fixed = TRUE)[[1L]][1L], "=", fixed = TRUE)[[1L]][2L],
                     "**", fixed = TRUE)[[1L]]
      d[[trimws(kv[1L])]] <- as.numeric(kv[2L])
    }
    d <- unlist(d)
    if (is.null(visual)) return(d)
    return(smoothness_p(d, visual = FALSE))
  }
  is_rv <- is.list(n) && length(n) == 2L && is.null(names(n)) && is.list(n[[2L]])
  if (is_rv) {
    rv <- n
  } else {
    facs <- if (!is.null(names(n))) .fac_refine(n) else .factorint(.nt_int(n))
    k <- if (power) 4L else 3L
    ents <- lapply(names(facs), function(f) {
      sm <- smoothness(gmp::as.bigz(f) + m)
      .nt_out(.zc(list(gmp::as.bigz(f), gmp::as.bigz(facs[[f]]), gmp::as.bigz(sm[1]), gmp::as.bigz(sm[2]))))
    })
    key1 <- .zc(lapply(ents, function(e) gmp::as.bigz(e[k])))
    key2 <- .zc(lapply(ents, function(e) gmp::as.bigz(e[1L])))
    o <- order(match(as.character(key1), unique(as.character(.zsort(key1)))),
               match(as.character(key2), unique(as.character(.zsort(key2)))))
    rv <- list(m, ents[o])
  }
  if (isFALSE(visual) || (!isTRUE(visual) && !is_rv)) return(rv)
  lines <- vapply(rv[[2L]], function(d) {
    d <- as.character(gmp::as.bigz(d))
    mm <- as.numeric(rv[[1L]])
    sprintf("p**i=%s**%s has p%s B=%s, B-pow=%s", d[1L], d[2L],
            if (mm >= 0) paste0("+", mm) else as.character(mm), d[3L], d[4L])
  }, "")
  paste(lines, collapse = "\n")
}

# ---------------------------------------------------------------------------
# Perfect powers

.prod_pow <- function(fac, g) {
  r <- gmp::as.bigz(1)
  for (k in ls(fac, all.names = TRUE)) r <- r * gmp::as.bigz(k)^(get(k, fac) %/% g)
  r
}

# Perfect power detection for integers without factors below next_p,
# returning list(base, exp) with the largest exponent, or FALSE.
.perfect_power_int <- function(n, next_p = 2) {
  if (n <= 3) return(FALSE)
  fac <- .fac_new()
  g <- 0
  multi <- 1
  done <- function(n, g, multi) {
    g <- .gcd_num(g, multi)
    if (g == 1) return(FALSE)
    .fac_set(fac, n, multi)
    list(.prod_pow(fac, g), g)
  }
  if (n <= 1000000) {
    n <- .factorint_small(fac, n, 1000, 1000, next_p)[[1L]]
    if (n > 1) return(FALSE)
    es <- unlist(mget(ls(fac), envir = fac))
    g <- Reduce(.gcd_num, es)
    if (g == 1) return(FALSE)
    return(list(.prod_pow(fac, g), g))
  }
  if (next_p < 3) {
    g <- .v2(n)
    if (g) {
      if (g == 1) return(FALSE)
      n <- n %/% gmp::as.bigz(2)^g
      .fac_set(fac, 2, g)
      if (n == 1) return(list(gmp::as.bigz(2), g))
      r <- .iroot(n, g)
      if (r$exact) return(list(2 * r$root, g))
      if (.is_prime(gmp::as.bigz(g))) return(FALSE)
    }
    next_p <- 3
  }
  while (n %% 8 == 1) {
    r <- .iroot(n, 2L)
    if (!r$exact) break
    n <- r$root
    multi <- multi * 2
  }
  if (n < gmp::as.bigz(next_p)^3) return(done(n, g, multi))
  tf_max <- .bitlen(n) %/% 27 + 24
  if (next_p < tf_max) {
    for (p in .primerange(gmp::as.bigz(next_p), gmp::as.bigz(tf_max))) {
      rm <- .remove(n, p)
      t <- rm[[2L]]
      if (t) {
        n <- rm[[1L]]
        t <- t * multi
        g_ <- .gcd_num(g, t)
        if (g_ == 1) return(FALSE)
        .fac_set(fac, p, t)
        if (n == 1) return(list(.prod_pow(fac, g_), g_))
        if (g == 0 || g_ < g) {
          g <- g_
          r <- .iroot(n^multi, g)
          if (r$exact) return(list(r$root * .prod_pow(fac, g), g))
          if (.is_prime(gmp::as.bigz(g))) return(FALSE)
        }
      }
    }
    next_p <- tf_max
  }
  if (n < gmp::as.bigz(next_p)^3) return(done(n, g, multi))
  if (g) {
    gg <- g / 2^.v2(gmp::as.bigz(g))
    prime_iter <- if (gg > 1) as.numeric(.fac_primes(.factorint(gmp::as.bigz(gg)))) else numeric()
  } else {
    prime_iter <- .primerange(gmp::as.bigz(3), gmp::as.bigz(trunc(log2(n) / log2(next_p)) + 2))
  }
  logn <- log2(n)
  threshold <- logn / 40
  for (p in prime_iter) {
    if (threshold < p) {
      repeat {
        b <- 2^(logn / p)
        rb <- trunc(b + 0.5)
        if (abs(rb - b) < 0.01 && gmp::as.bigz(rb)^p == n) {
          n <- gmp::as.bigz(rb)
          multi <- multi * p
          logn <- log2(n)
        } else break
      }
    } else {
      repeat {
        r <- .iroot(n, p)
        if (!r$exact) break
        n <- r$root
        multi <- multi * p
        logn <- log2(n)
      }
    }
    if (n < gmp::as.bigz(next_p)^(p + 2)) break
  }
  done(n, g, multi)
}

# perfect_power() on bigz or bigq values; list(base, exp) or FALSE.
.perfect_power <- function(n, candidates = NULL, big = TRUE, factor = TRUE) {
  if (n < 0) {
    if (is.null(candidates)) {
      pp <- .perfect_power(-n, NULL, TRUE, factor)
      if (isFALSE(pp)) return(FALSE)
      e <- pp[[2L]]
      e2 <- 2^.v2(gmp::as.bigz(e))
      b <- pp[[1L]]^e2
      e <- e %/% e2
      if (e <= 1) return(FALSE)
      if (big || .is_prime(gmp::as.bigz(e))) return(list(-b, e))
      for (p in .primerange(gmp::as.bigz(3), gmp::as.bigz(e + 1)))
        if (e %% p == 0) return(list(-b^(e %/% p), p))
    }
    odd <- unique(candidates[candidates %% 2 == 1])
    if (!length(odd)) return(FALSE)
    pp <- .perfect_power(-n, odd, big, factor)
    if (isFALSE(pp)) return(FALSE)
    return(list(-pp[[1L]], pp[[2L]]))
  }
  if (inherits(n, "bigq")) {
    if (gmp::denominator(n) == 1) {
      n <- gmp::numerator(n)
    } else {
      p <- gmp::numerator(n)
      q <- gmp::denominator(n)
      if (p == 1) {
        qq <- .perfect_power(q, candidates, big, factor)
        if (isFALSE(qq)) return(FALSE)
        return(list(1 / gmp::as.bigq(qq[[1L]]), qq[[2L]]))
      }
      pp <- .perfect_power(p, NULL, TRUE, factor)
      if (isFALSE(pp)) return(FALSE)
      qq <- .perfect_power(q, NULL, TRUE, factor)
      if (isFALSE(qq)) return(FALSE)
      tup <- function(ex) list(gmp::as.bigq(pp[[1L]]^(pp[[2L]] %/% ex), qq[[1L]]^(qq[[2L]] %/% ex)), ex)
      if (length(candidates)) {
        ok <- candidates[pp[[2L]] %% candidates == 0 & qq[[2L]] %% candidates == 0]
        if (!length(ok)) return(FALSE)
        return(tup(if (big) max(ok) else min(ok)))
      }
      g <- .gcd_num(pp[[2L]], qq[[2L]])
      if (g == 1) return(FALSE)
      if (big) return(tup(g))
      for (p in .primerange(gmp::as.bigz(2), gmp::as.bigz(g + 1))) if (g %% p == 0) return(tup(p))
    }
  }
  if (!is.null(candidates)) candidates <- unique(as.numeric(candidates))
  if (is.null(candidates) && big) return(.perfect_power_int(n))
  if (n <= 3) return(FALSE)
  logn <- log2(n)
  max_possible <- trunc(logn) + 2
  not_square <- as.numeric(n %% 10) %in% c(2, 3, 7, 8)
  min_possible <- 2 + not_square
  if (length(candidates)) {
    cs <- sort(candidates[candidates >= min_possible & candidates < max_possible])
    if (n %% 2 == 0) {
      e <- .v2(n)
      cs <- cs[e %% cs == 0]
    }
    if (big) cs <- rev(cs)
    for (e in cs) {
      r <- .iroot(n, e)
      if (r$exact) return(list(r$root, e))
    }
    return(FALSE)
  }
  cands <- .primerange(gmp::as.bigz(min_possible), gmp::as.bigz(max_possible))
  fac <- gmp::as.bigz(2 + as.numeric(n %% 2))
  for (e in cands) {
    if (factor && n %% fac == 0) {
      e <- .remove(n, fac)[[2L]]
      if (e == 1) return(FALSE)
      r <- .iroot(n, e)
      rt <- r$root
      if (!r$exact) {
        m <- n %/% fac^e
        rE <- .perfect_power(m, .divisors_vec(gmp::as.bigz(e)), TRUE, TRUE)
        if (isFALSE(rE)) return(FALSE)
        E <- rE[[2L]]
        rt <- fac^(e %/% E) * rE[[1L]]
        e <- E
      }
      if (!big) {
        e0 <- as.numeric(.fac_primes(.factorint(gmp::as.bigz(e))))[1L]
        if (e0 != e) {
          rt <- rt^(e %/% e0)
          e <- e0
        }
      }
      return(list(rt, e))
    }
    fac <- gmp::nextprime(fac)
    if (logn / e < 40) {
      b <- 2^(logn / e)
      if (abs(trunc(b + 0.5) - b) > 0.01) next
    }
    r <- .iroot(n, e)
    if (r$exact) {
      rt <- r$root
      if (big) {
        m <- .perfect_power(rt, NULL, big, factor)
        if (!isFALSE(m)) {
          rt <- m[[1L]]
          e <- e * m[[2L]]
        }
      }
      return(list(rt, e))
    }
  }
  FALSE
}

#' Perfect powers
#'
#' `perfect_power(n)` returns `c(b, e)` with `b^e == n` and `e > 1` as large as
#' possible, or `FALSE`. With `big = FALSE` the smallest prime exponent is
#' used. `candidates` restricts the exponents tried; `factor = FALSE` skips the
#' trial division by small primes. A rational `n` returns `list(b, e)` with `b`
#' a `sym` rational.
#'
#' @param n An integer or rational.
#' @param candidates Exponents to try, or `NULL`.
#' @param big Prefer the largest exponent.
#' @param factor Try small prime factors first.
#' @return `c(b, e)`, `list(b, e)` for rationals, or `FALSE`.
#' @examples
#' perfect_power(16)
#' perfect_power(16, big = FALSE)
#' perfect_power(Rational(4, 9))
#' perfect_power(-8)
#' @export
perfect_power <- function(n, candidates = NULL, big = TRUE, factor = TRUE) {
  q <- .nt_rat(n)
  v <- if (is.null(q)) .nt_int(n) else q
  if (!is.null(candidates)) candidates <- as.numeric(.nt_ints(candidates))
  r <- .perfect_power(v, candidates, big, factor)
  if (isFALSE(r)) return(FALSE)
  if (!is.null(q)) return(list(make_rational(gmp::as.bigq(r[[1L]])), r[[2L]]))
  .nt_out(.zc(list(gmp::as.bigz(r[[1L]]), gmp::as.bigz(r[[2L]]))))
}

# ---------------------------------------------------------------------------
# Pollard's rho and p - 1

.pollard_rho <- function(n, s = 2, a = 1, retries = 5, seed = 1234, max_steps = NULL, F = NULL) {
  if (n < 5) stop("pollard_rho should receive n > 4", call. = FALSE)
  randint <- .nt_rng(seed + retries)
  V <- gmp::as.bigz(s)
  a <- gmp::as.bigz(a)
  for (i in seq_len(retries + 1)) {
    U <- V
    f <- if (is.null(F)) function(x) (x * x + a) %% n else function(x) gmp::as.bigz(F(x))
    j <- 0
    repeat {
      if (!is.null(max_steps) && j > max_steps) break
      j <- j + 1
      U <- f(U)
      V <- f(f(V))
      g <- .gcd2(U - V, n)
      if (g == 1) next
      if (g == n) break
      return(g)
    }
    V <- randint(0, n - 1)
    a <- randint(1, n - 3)
    F <- NULL
  }
  NULL
}

#' Pollard's rho and p - 1 factoring
#'
#' `pollard_rho(n)` searches for a nontrivial factor of `n` with the sequence
#' `x -> x^2 + a` started at `s`, or with the function `F`, retrying with random
#' parameters. `pollard_pm1(n, B)` uses Pollard's p - 1 method with
#' smoothness bound `B`. Both return a factor or `NULL`. The random retries
#' use a stream private to the call, seeded by `seed`, so the results are
#' reproducible but differ from SymPy's.
#'
#' @param n An integer greater than 4 (`pollard_rho()`) or 3 (`pollard_pm1()`).
#' @param s Starting value.
#' @param a Constant of the default polynomial, or base for p - 1.
#' @param retries Number of retries with new random parameters.
#' @param seed Seed for the retries.
#' @param max_steps Maximum number of steps per attempt.
#' @param F A function of one big integer to use instead of `x^2 + a`.
#' @param B Smoothness bound, greater than 2.
#' @return A factor, or `NULL`.
#' @examples
#' pollard_rho(8051)
#' pollard_pm1(21477639576571, B = 4000)
#' @export
pollard_rho <- function(n, s = 2, a = 1, retries = 5, seed = 1234, max_steps = NULL, F = NULL) {
  r <- .pollard_rho(.nt_int(n), .nt_int(s), .nt_int(a), retries, seed, max_steps, F)
  if (is.null(r)) NULL else .nt_out(r)
}

.pollard_pm1 <- function(n, B = 10, a = 2, retries = 0, seed = 1234) {
  if (n < 4 || B < 3) stop("pollard_pm1 should receive n > 3 and B > 2", call. = FALSE)
  randint <- .nt_rng(seed + B)
  ps <- .primerange(gmp::as.bigz(2), gmp::as.bigz(B + 1))
  L <- prod(gmp::as.bigz(ps)^trunc(log(B) / log(ps)))
  a <- gmp::as.bigz(a)
  for (i in seq_len(retries + 1)) {
    g <- .gcd2(gmp::powm(a, L, n) - 1, n)
    if (g > 1 && g < n) return(g)
    a <- randint(2, n - 2)
  }
  NULL
}

#' @rdname pollard_rho
#' @export
pollard_pm1 <- function(n, B = 10, a = 2, retries = 0, seed = 1234) {
  r <- .pollard_pm1(.nt_int(n), as.numeric(B), a, retries, seed)
  if (is.null(r)) NULL else .nt_out(r)
}

# ---------------------------------------------------------------------------
# Elliptic curve method

.ec_add <- function(P, Q, diff, n) {
  u <- (P[[1L]] - P[[2L]]) * (Q[[1L]] + Q[[2L]])
  v <- (P[[1L]] + P[[2L]]) * (Q[[1L]] - Q[[2L]])
  ad <- u + v
  sb <- u - v
  list((diff[[2L]] * ad * ad) %% n, (diff[[1L]] * sb * sb) %% n)
}

.ec_double <- function(P, a24, n) {
  s <- P[[1L]] + P[[2L]]
  t <- P[[1L]] - P[[2L]]
  u <- (s * s) %% n
  v <- (t * t) %% n
  d <- u - v
  list((u * v) %% n, (d * (v + a24 * d)) %% n)
}

.ec_ladder <- function(P, k, a24, n) {
  Q <- P
  R <- .ec_double(P, a24, n)
  bits <- strsplit(as.character(gmp::as.bigz(k), b = 2L), "", fixed = TRUE)[[1L]][-1L]
  for (b in bits) {
    if (b == "1") {
      Q <- .ec_add(R, Q, P, n)
      R <- .ec_double(R, a24, n)
    } else {
      R <- .ec_add(Q, R, P, n)
      Q <- .ec_double(Q, a24, n)
    }
  }
  Q
}

# Stage 1 runs the Montgomery ladder on a batch of curves at once (the
# coordinates are bigz vectors, and every curve shares the bits of k); the
# curves are then examined in order, so the result is that of one at a time.
# Batches grow from 4 to 32 curves: small ones answer lucky cases quickly,
# large ones spread the per-call overhead of gmp.
# ponytail: stage 2 still runs one curve at a time.
.ecm_one_factor <- function(n, B1 = 10000, B2 = 100000, max_curve = 200, seed = NULL) {
  randint <- .nt_rng(seed)
  D <- min(as.numeric(.isqrt(gmp::as.bigz(B2))), B1 %/% 2 - 1)
  ps <- .primerange(gmp::as.bigz(2), gmp::as.bigz(B1 + 1))
  k <- prod(gmp::as.bigz(ps)^trunc(log(B1) / log(ps)))
  deltas_list <- list()
  for (r in seq(B1 + 2 * D, B2 + 2 * D - 1, by = 4 * D)) {
    q <- .primerange(gmp::as.bigz(r - 2 * D), gmp::as.bigz(r + 2 * D))
    deltas_list[[length(deltas_list) + 1L]] <- unique(abs(q - r) %/% 2)
  }
  done <- 0
  batch <- 4
  while (done < max_curve) {
    nb <- min(batch, max_curve - done)
    batch <- min(2 * batch, 32)
    done <- done + nb
    sigma <- .zc(lapply(seq_len(nb), function(i) randint(6, n - 1)))
    u <- (sigma * sigma - 5) %% n
    v <- (4 * sigma) %% n
    u3 <- (u * u * u) %% n
    v3 <- (v * v * v) %% n
    den <- .zc(lapply(seq_len(nb), function(i) .zinv(16 * u3[i] * v[i], n)))
    bad <- is.na(den)
    den[bad] <- gmp::as.bigz(1)
    a24 <- (((v - u)^3L %% n) * (3 * u + v) * den) %% n
    Qb <- .ec_ladder(list(u3, v3), k, a24, n)
    for (j in seq_len(nb)) {
      if (bad[j]) {
        g <- .gcd2(2 * u3[j] * v[j], n)
        if (g == n) next
        return(g)
      }
      Q <- list(Qb[[1L]][j], Qb[[2L]][j])
      g <- .gcd2(Q[[2L]], n)
      if (g != 1 && g != n) return(g)
      if (g == n) next
      g <- .ecm_stage2(Q, a24[j], n, B1, D, deltas_list)
      if (g != 1 && g != n) return(g)
    }
  }
  NULL
}

.ecm_stage2 <- function(Q, a24, n, B1, D, deltas_list) {
  S <- vector("list", D)
  S[[1L]] <- Q
  Q2 <- .ec_double(Q, a24, n)
  S[[2L]] <- .ec_add(Q2, Q, Q, n)
  for (d in seq_len(D - 2) + 2) S[[d]] <- .ec_add(S[[d - 1L]], Q2, S[[d - 2L]], n)
  Sx <- .zc(lapply(S, `[[`, 1L))
  Sz <- .zc(lapply(S, `[[`, 2L))
  beta <- (Sx * Sz) %% n
  g <- gmp::as.bigz(1)
  W <- .ec_ladder(Q, 4 * D, a24, n)
  Tp <- .ec_ladder(Q, B1 - 2 * D, a24, n)
  R <- .ec_ladder(Q, B1 + 2 * D, a24, n)
  for (deltas in deltas_list) {
    alpha <- (R[[1L]] * R[[2L]]) %% n
    if (length(deltas)) {
      idx <- deltas + 1
      f <- ((R[[1L]] - Sx[idx]) * (R[[2L]] + Sz[idx]) - alpha + beta[idx]) %% n
      g <- (g * prod(f)) %% n
    }
    nxt <- .ec_add(R, W, Tp, n)
    Tp <- R
    R <- nxt
  }
  .gcd2(n, g)
}

#' Lenstra's elliptic curve factoring
#'
#' `ecm(n)` returns the sorted prime factors of `n`, after trial division by
#' the primes below 100000. `B1` and `B2` are the even stage 1 and stage 2
#' bounds and `max_curve` the number of curves tried per factor. Random curves
#' come from a stream private to the call, seeded by `seed`.
#'
#' @param n An integer.
#' @param B1,B2 Even stage bounds.
#' @param max_curve Maximum number of curves.
#' @param seed Seed for the curves.
#' @return A vector of distinct primes.
#' @examples
#' ecm("25645121643901801")
#' @export
ecm <- function(n, B1 = 10000, B2 = 100000, max_curve = 200, seed = 1234) {
  n <- .nt_int(n)
  if (n < 1) stop("n should be a positive integer", call. = FALSE)
  if (B1 %% 2 != 0 || B2 %% 2 != 0) stop("both bounds must be even", call. = FALSE)
  TF_LIMIT <- 100000
  ps <- .primerange(gmp::as.bigz(2), gmp::as.bigz(TF_LIMIT))
  factors <- list()
  for (p in ps[n %% gmp::as.bigz(ps) == 0]) {
    factors[[length(factors) + 1L]] <- gmp::as.bigz(p)
    n <- .remove(n, p)[[1L]]
  }
  queue <- list()
  check <- function(m) {
    if (m == 1) return(invisible())
    if (.is_prime(m)) {
      factors[[length(factors) + 1L]] <<- m
      return(invisible())
    }
    r <- .perfect_power_int(m, TF_LIMIT)
    if (!isFALSE(r)) return(check(gmp::as.bigz(r[[1L]])))
    queue[[length(queue) + 1L]] <<- m
  }
  check(n)
  while (length(queue)) {
    n <- queue[[length(queue)]]
    queue[[length(queue)]] <- NULL
    f <- .ecm_one_factor(n, B1, B2, max_curve, seed)
    if (is.null(f)) stop("Increase the bounds", call. = FALSE)
    check(f)
    check(n %/% f)
  }
  fs <- .zc(factors)
  .nt_out(.zsort(fs[!duplicated(as.character(fs))]))
}

# ---------------------------------------------------------------------------
# factorint

.TRIAL_CAND <- local({
  k <- 0:5470
  c(rbind(6 * k + 5, 6 * k + 7))
})

# Trial division by 2, 3 and the integers 6k +- 1, stopping after fail_max
# misses in a row or once next_p^2 passes min(n, limit^2). Returns
# list(n, next_p), with next_p = 0 when the remaining n is 1 or prime.
.factorint_small <- function(fac, n, limit, fail_max, next_p = 2) {
  done <- function(n, d) list(n, if (gmp::as.bigz(d)^2L <= n) d else 0)
  limit2 <- gmp::as.bigz(limit)^2L
  threshold2 <- min(n, limit2)
  strip <- function(p) {
    r <- .remove(n, p)
    n <<- r[[1L]]
    .fac_set(fac, p, r[[2L]])
    threshold2 <<- min(n, limit2)
  }
  if (next_p < 3) {
    if (n %% 2 == 0) strip(2)
    next_p <- 3
    if (threshold2 < 9) return(done(n, next_p))
  }
  if (next_p < 5) {
    if (n %% 3 == 0) strip(3)
    next_p <- 5
    if (threshold2 < 25) return(done(n, next_p))
  }
  p6 <- next_p %% 6
  next_p <- next_p + (if (p6 < 2) -1 else 5) - p6
  fails <- 0
  i <- match(next_p, .TRIAL_CAND)
  if (is.na(i)) {
    cand <- function(j) {
      base <- next_p + 6 * ((j - 1) %/% 2)
      base + 2 * ((j - 1) %% 2)
    }
    i <- 1
  } else {
    cand <- function(j) .TRIAL_CAND[j + off]
    off <- i - 1
    i <- 1
  }
  repeat {
    idx <- i:(i + 511)
    cs <- cand(idx)
    hit <- n %% gmp::as.bigz(cs) == 0
    for (j in seq_along(idx)) {
      c0 <- cs[j]
      if (c0 %% 6 == 5 && fails >= fail_max) return(done(n, c0))
      if (hit[j]) {
        strip(c0)
        fails <- 0
      } else {
        fails <- fails + 1
      }
      nxt <- c0 + if (c0 %% 6 == 5) 2 else 4
      if (threshold2 < gmp::as.bigz(nxt)^2L) return(done(n, nxt))
      if (hit[j]) {
        i <- idx[j] + 1
        break
      }
      if (j == length(idx)) i <- idx[j] + 1
    }
  }
}

.trial <- function(fac, n, candidates, verbose = FALSE) {
  found <- FALSE
  if (!length(candidates)) return(list(n, found))
  cz <- if (inherits(candidates, "bigz")) candidates else gmp::as.bigz(candidates)
  hit <- which(n %% cz == 0)
  for (i in hit) {
    d <- cz[i]
    if (n %% d != 0) next
    if (n != d) .fc_put(n, d)
    r <- .remove(n %/% d, d)
    n <- r[[1L]]
    .fac_set(fac, d, r[[2L]] + 1)
    .vmsg(verbose, "\t%s ** %s", as.character(d), r[[2L]] + 1)
    found <- TRUE
  }
  list(n, found)
}

.check_termination <- function(fac, n, limit, opts, next_p) {
  .vmsg(opts$verbose, "Check for termination")
  if (n == 1) return(TRUE)
  if (n < gmp::as.bigz(next_p)^2L || .is_prime(n)) {
    .fc_put(n, n)
    .fac_set(fac, n, 1)
    return(TRUE)
  }
  if (is.null(limit) && opts$all_methods && n < .NT_FACTORIZE_MAX) {
    # GMP's compiled Pollard rho beats the R pipeline below 10^28; without a
    # limit the factorization is unique, so the result is the same.
    .fac_merge(fac, .gmp_factor(n))
    return(TRUE)
  }
  p <- .perfect_power_int(n, next_p)
  if (isFALSE(p)) return(FALSE)
  base <- gmp::as.bigz(p[[1L]])
  ex <- p[[2L]]
  if (base < gmp::as.bigz(next_p)^2L || .is_prime(base)) {
    .fc_put(n, base)
    .fac_set(fac, base, ex)
  } else {
    facs <- .factorint_pos(base, limit, opts)
    for (b in names(facs)) .fac_set(fac, b, ex * facs[[b]])
  }
  TRUE
}

.gmp_factor <- function(n) {
  f <- gmp::factorize(n)
  ks <- as.character(f)
  tab <- table(ks)
  d <- stats::setNames(as.numeric(tab), names(tab))
  d[.zorder(gmp::as.bigz(names(d)))]
}

.factorint_pos <- function(n, limit, opts) {
  if (n < 10) {
    return(switch(as.character(n), "0" = c("0" = 1), "1" = stats::setNames(numeric(), character()),
                  "2" = c("2" = 1), "3" = c("3" = 1), "4" = c("2" = 2), "5" = c("5" = 1),
                  "6" = c("2" = 1, "3" = 1), "7" = c("7" = 1), "8" = c("2" = 3), "9" = c("3" = 2)))
  }
  fac <- .fac_new()
  verbose <- opts$verbose
  if (verbose) {
    sn <- as.character(n)
    if (nchar(sn) > 50) message(sprintf("Factoring %s..(%d other digits)..%s",
                                        substr(sn, 1, 5), nchar(sn) - 10, substring(sn, nchar(sn) - 4)))
    else message("Factoring ", sn)
  }
  small <- min(2^15, if (is.null(limit)) 2^15 else limit)
  .vmsg(verbose, "Trial division with ints [%d ... %d] and fail_max=%d", 2, small, 600)
  r <- .factorint_small(fac, n, small, 600)
  n <- r[[1L]]
  next_p <- r[[2L]]
  if (next_p == 0) {
    if (n > 1) .fac_set(fac, n, 1)
    .vmsg(verbose, "Factorization is complete.")
    return(.fac_vec(fac))
  }
  repeat {
    p <- .fc_lookup(n)
    if (is.null(p)) break
    rm <- .remove(n, p)
    n <- rm[[1L]]
    .fac_set(fac, p, rm[[2L]])
  }
  if (!is.null(limit) && next_p > limit) {
    .vmsg(verbose, "Exceeded limit: %s", format(limit))
    if (.check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
    if (n > 1) .fac_set(fac, n, 1)
    return(.fac_vec(fac))
  }
  if (.check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
  sqrt_n <- .isqrt(n)
  a <- sqrt_n + 1
  if ((n %% 4 == 1) != (a %% 2 == 1)) a <- a + 1
  b2 <- a^2L - n
  for (i in 1:3) {
    rb <- .iroot(b2, 2L)
    if (rb$exact) {
      .vmsg(verbose, "Close factors satisfying Fermat condition found.")
      for (rr in list(a - rb$root, a + rb$root)) .fac_merge(fac, .factorint_pos(rr, limit, opts))
      return(.fac_vec(fac))
    }
    b2 <- b2 + 4 * (a + 1)
    a <- a + 2
  }
  low <- next_p
  high <- 2 * next_p
  lim <- (if (is.null(limit)) sqrt_n else gmp::as.bigz(limit)) + 1
  iteration <- 0
  sub_factor <- function(cand) {
    if (cand < gmp::as.bigz(next_p)^2L || .is_prime(cand)) cand else .fac_primes(.factorint_pos(cand, limit, opts))
  }
  repeat {
    high_ <- min(gmp::as.bigz(high), lim)
    found_trial <- FALSE
    if (opts$use_trial) {
      .vmsg(verbose, "Trial division with primes [%s ... %s]", format(low), as.character(high_))
      tr <- .trial(fac, n, .primerange(gmp::as.bigz(low), high_), verbose)
      n <- tr[[1L]]
      found_trial <- tr[[2L]]
      next_p <- as.numeric(high_)
      if (found_trial && .check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
    }
    if (high > lim) {
      .vmsg(verbose, "Exceeded limit: %s", as.character(lim))
      if (n > 1) .fac_set(fac, n, 1)
      return(.fac_vec(fac))
    }
    if (!found_trial) {
      if (opts$use_pm1) {
        .vmsg(verbose, "Pollard's p-1 with smoothness bound %s and seed %s", format(low), as.character(high_))
        cand <- .pollard_pm1(n, B = low, seed = as.numeric(high_))
        if (!is.null(cand)) {
          n <- .trial(fac, n, sub_factor(cand))[[1L]]
          if (.check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
        }
      }
      if (opts$use_rho) {
        .vmsg(verbose, "Pollard's rho with retries %d, max_steps %s and seed %s", 1, format(low), as.character(high_))
        cand <- .pollard_rho(n, retries = 1, max_steps = low, seed = as.numeric(high_))
        if (!is.null(cand)) {
          n <- .trial(fac, n, sub_factor(cand))[[1L]]
          if (.check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
        }
      }
    }
    iteration <- iteration + 1
    if (opts$use_ecm && iteration >= 3 && nchar(as.character(n)) >= 24) break
    low <- high
    high <- high * 2
  }
  B1 <- 10000
  num_curves <- 50
  repeat {
    B2 <- 100 * B1
    .vmsg(verbose, "Elliptic Curve with B1 bound %d, B2 bound %d, num_curves %d", B1, B2, num_curves)
    f <- .ecm_one_factor(n, B1, B2, num_curves, seed = B1)
    if (!is.null(f)) {
      n <- .trial(fac, n, sub_factor(f))[[1L]]
      if (.check_termination(fac, n, limit, opts, next_p)) return(.fac_vec(fac))
    }
    B1 <- B1 * 5
    num_curves <- num_curves * 4
  }
}

# Integer factorization for internal callers: named exponent vector.
.factorint <- function(n, limit = NULL) {
  n <- gmp::as.bigz(n)
  if (n < 0) {
    f <- .factorint(-n, limit)
    return(c(c("-1" = 1), f))
  }
  if (is.null(limit) && n < .NT_FACTORIZE_MAX && n > 1) return(.gmp_factor(n))
  .factorint_pos(n, limit, list(use_trial = TRUE, use_rho = TRUE, use_pm1 = TRUE,
                                use_ecm = is.null(limit), verbose = FALSE, all_methods = TRUE))
}

# Refine a factor dictionary so every key is prime, merging exponents.
.fac_refine <- function(d, factor_fn = function(k) .factorint(k)) {
  d <- unlist(d)
  fac <- .fac_new()
  for (k in names(d)) {
    kz <- .nt_int(k)
    if (.is_prime(kz)) .fac_add(fac, kz, d[[k]])
    else .fac_merge(fac, factor_fn(kz), d[[k]])
  }
  .fac_vec(fac)
}

.fac_visual <- function(d) {
  if (!length(d)) return(S_One)
  args <- list()
  if ("-1" %in% names(d)) {
    args <- list(S_NegativeOne)
    d <- d[names(d) != "-1"]
  }
  for (k in names(d)) {
    e <- d[[k]]
    # p^-e is kept as 1/(p^e) so that the printer shows the factored power
    pw <- new_sym("Pow", list(Integer(k), Integer(abs(e))))
    args[[length(args) + 1L]] <- if (e < 0) new_sym("Pow", list(pw, S_NegativeOne)) else pw
  }
  if (length(args) == 1L) return(args[[1L]])
  new_sym("Mul", args)
}

.sym_powers <- function(e) {
  terms <- if (e$type == "Mul") e$args else list(e)
  fac <- .fac_new()
  for (t in terms) {
    ex <- gmp::as.bigz(1)
    while (t$type == "Pow") {
      ex <- ex * .nt_int(t$args[[2L]])
      t <- t$args[[1L]]
    }
    if (abs(ex) >= .NT_Z53) stop("exponents at or above 2^53 are not supported", call. = FALSE)
    ex <- as.numeric(ex)
    .fac_add(fac, .nt_int(t), ex)
  }
  .fac_vec(fac)
}

.fac_multiple <- function(d) {
  if (any(d < 0)) {
    out <- list()
    for (k in names(d)) {
      e <- d[[k]]
      item <- if (e > 0) Integer(k) else Rational(1, k)
      out <- c(out, rep(list(item), abs(e)))
    }
    return(out)
  }
  if (!length(d)) return(numeric())
  .nt_out(gmp::as.bigz(rep(names(d), times = d)))
}

#' Integer factorization
#'
#' `factorint(n)` returns the prime factorization of `n` as a named vector of
#' exponents whose names are the primes, in increasing order; `-1` appears for
#' negative `n` and `0` for zero. The method follows SymPy: trial division,
#' a perfect power test, Fermat's method, then Pollard p - 1 and rho
#' alongside trial division, and finally Lenstra's elliptic curve method. With
#' no `limit`, numbers below `10^28` go to the compiled Pollard rho of GMP,
#' which gives the same unique result faster.
#'
#' `limit` stops trial division at that bound (and disables the elliptic curve
#' method), so the last factor may be composite. The `use_*` switches choose
#' methods. `multiple = TRUE` lists the factors with repetition. A named vector
#' of factors (keys need not be prime) or an unevaluated product is refactored,
#' and returned as an unevaluated `sym` product unless `visual = FALSE`;
#' `visual = TRUE` gives that product for any input. An unevaluated
#' `sym_factorial(k)` is factored with Legendre's formula.
#'
#' `factorrat(r)` factors a rational, with negative exponents for the
#' denominator. `primefactors(n)` lists the distinct primes.
#'
#' @param n An integer, a named factor vector, a `sym` product, or a
#'   `sym_factorial()` expression.
#' @param rat A rational number.
#' @param limit Trial division bound, or `NULL`.
#' @param use_trial,use_rho,use_pm1,use_ecm Method switches.
#' @param verbose Report progress with [message()].
#' @param visual `TRUE` for a `sym` product, `FALSE` for the vector, `NULL`
#'   for the default.
#' @param multiple List factors with repetition.
#' @param ... Passed to `factorint()`.
#' @return A named numeric vector, a vector of factors, or a `sym` product.
#' @examples
#' factorint(2000)
#' factorint(-60)
#' factorint(2^20 + 1, multiple = TRUE)
#' factorint(c("4" = 2, "12" = 3))
#' factorint(100, visual = TRUE)
#' factorrat(Rational(8, 9))
#' primefactors(6936)
#' @export
factorint <- function(n, limit = NULL, use_trial = TRUE, use_rho = TRUE, use_pm1 = TRUE,
                      use_ecm = TRUE, verbose = FALSE, visual = NULL, multiple = FALSE) {
  if (multiple) {
    d <- factorint(n, limit, use_trial, use_rho, use_pm1, use_ecm, verbose, visual = FALSE)
    return(.fac_multiple(d))
  }
  is_dict <- !is_sym(n) && !is.null(names(n))
  is_mul <- is_sym(n) && n$type %in% c("Mul", "Pow")
  factordict <- NULL
  if (is_dict) factordict <- .fac_refine(n, function(k)
    factorint(k, limit, use_trial, use_rho, use_pm1, use_ecm, verbose, visual = FALSE))
  if (is_mul) factordict <- .fac_refine(.sym_powers(n), function(k)
    factorint(k, limit, use_trial, use_rho, use_pm1, use_ecm, verbose, visual = FALSE))
  if (isTRUE(visual) && !is_dict && !is_mul)
    factordict <- factorint(n, limit, use_trial, use_rho, use_pm1, use_ecm, verbose, visual = FALSE)
  if (isTRUE(visual) || (is_dict && is.null(visual))) return(.fac_visual(factordict))
  if (is_dict || is_mul) return(factordict)
  if (!(use_trial || use_rho || use_pm1 || use_ecm)) stop("at least one method must be enabled", call. = FALSE)
  if (is_sym(n) && n$type == "Function" && n$name == "factorial") {
    x <- .nt_int(n$args[[1L]])
    if (x >= 20) {
      ps <- .primerange(gmp::as.bigz(2), x + 1)
      ex <- vapply(ps, function(p) {
        m <- 0
        q <- x %/% p
        while (q != 0) {
          m <- m + as.numeric(q)
          q <- q %/% p
        }
        m
      }, 0)
      .vmsg(verbose, "Factorization is complete.")
      return(stats::setNames(ex, as.character(gmp::as.bigz(ps))))
    }
    n <- gmp::factorialZ(as.numeric(x))
  }
  n <- .nt_int(n)
  if (!is.null(limit)) {
    limit <- as.numeric(.nt_int(trunc(as.numeric(limit))))
    if (limit == 0) limit <- NULL else use_ecm <- FALSE
  }
  opts <- list(use_trial = use_trial, use_rho = use_rho, use_pm1 = use_pm1,
               use_ecm = use_ecm, verbose = verbose,
               all_methods = use_trial && use_rho && use_pm1 && use_ecm)
  if (n < 0) {
    f <- factorint(-n, limit, use_trial, use_rho, use_pm1, use_ecm, verbose, visual = FALSE)
    return(c(c("-1" = 1), f))
  }
  if (!is.null(limit) && limit < 2) {
    if (n == 1) return(stats::setNames(numeric(), character()))
    return(stats::setNames(1, as.character(n)))
  }
  if (is.null(limit) && opts$all_methods && n >= 10 && n < .NT_FACTORIZE_MAX && !verbose) return(.gmp_factor(n))
  .factorint_pos(n, limit, opts)
}

#' @rdname factorint
#' @export
factorrat <- function(rat, limit = NULL, use_trial = TRUE, use_rho = TRUE, use_pm1 = TRUE,
                      verbose = FALSE, visual = NULL, multiple = FALSE) {
  q <- .nt_rat(rat)
  q <- if (is.null(q)) gmp::as.bigq(.nt_int(rat)) else q
  f1 <- factorint(gmp::numerator(q), limit, use_trial, use_rho, use_pm1, verbose = verbose, visual = FALSE)
  f2 <- factorint(gmp::denominator(q), limit, use_trial, use_rho, use_pm1, verbose = verbose, visual = FALSE)
  fac <- .fac_new()
  .fac_merge(fac, f1)
  .fac_merge(fac, f2, -1)
  d <- .fac_vec(fac)
  if (length(d) > 1 && "1" %in% names(d)) d <- d[names(d) != "1"]
  if (multiple) {
    # sorted by value: keys below 1, then reciprocals 1/p, then the primes
    keys <- names(d)
    kz <- gmp::as.bigz(keys)
    o <- c(which(d > 0 & kz < 1), rev(which(d < 0)), which(d > 0 & kz >= 1))
    item <- function(i) if (d[[i]] > 0) Integer(keys[i]) else Rational(1, keys[i])
    return(unlist(lapply(o, function(i) rep(list(item(i)), abs(d[[i]]))), recursive = FALSE))
  }
  if (isTRUE(visual)) return(.fac_visual(d))
  d
}

#' @rdname factorint
#' @export
primefactors <- function(n, limit = NULL, verbose = FALSE, ...) {
  n <- .nt_int(n)
  fs <- .fac_primes(factorint(n, limit = limit, verbose = verbose, visual = NULL, multiple = FALSE, ...))
  if (!length(fs)) return(numeric())
  k <- length(fs)
  s <- fs[seq_len(k - 1L)]
  s <- s[!(as.character(s) %in% c("-1", "0", "1"))]
  if (.is_prime(fs[k])) s <- c(s, fs[k])
  if (!length(s)) return(numeric())
  .nt_out(s)
}

# ---------------------------------------------------------------------------
# Divisors

.divisors_from <- function(f) {
  ps <- .fac_primes(f)
  ps <- ps[ps > 1]
  ex <- f[as.character(ps)]
  tot <- gmp::as.bigz(1)
  for (i in seq_along(ps)) tot <- tot * ps[i]^ex[[i]]
  if (tot < .NT_Z53) {
    d <- 1
    for (i in seq_along(ps)) d <- as.vector(outer(d, as.numeric(ps[i])^(0:ex[[i]])))
    return(sort(d))
  }
  d <- gmp::as.bigz(1)
  for (i in seq_along(ps)) {
    parts <- lapply(0:ex[[i]], function(k) d * ps[i]^k)
    d <- do.call(c, parts)
  }
  .zsort(d)
}

.divisors_vec <- function(n) {
  n <- abs(n)
  if (n == 0) return(numeric())
  if (n == 1) return(1)
  .divisors_from(.factorint(n))
}

#' Divisors
#'
#' `divisors(n)` lists the positive divisors of `n` in increasing order, and
#' `proper_divisors(n)` leaves out `n`. `divisor_count(n, modulus)` counts the
#' divisors that are multiples of `modulus`. `udivisors()` lists the unitary
#' divisors `d` with `gcd(d, n / d) = 1`, and `antidivisors()` the integers
#' that divide `n` as unevenly as possible. The `*_count()` functions count
#' them.
#'
#' @param n An integer.
#' @param proper Leave out `n`.
#' @param modulus Count only multiples of this.
#' @return A vector of divisors, or a count.
#' @examples
#' divisors(24)
#' proper_divisors(24)
#' divisor_count(24, 2)
#' udivisors(60)
#' antidivisors(24)
#' @export
divisors <- function(n, proper = FALSE) {
  n <- abs(.nt_int(n))
  d <- .divisors_vec(n)
  if (proper) d <- d[d != n]
  .nt_out(if (inherits(d, "bigz")) d else gmp::as.bigz(d))
}

#' @rdname divisors
#' @export
proper_divisors <- function(n) divisors(n, proper = TRUE)

#' @rdname divisors
#' @export
divisor_count <- function(n, modulus = 1, proper = FALSE) {
  n <- .nt_int(n)
  modulus <- .nt_int(modulus)
  if (modulus == 0) return(0)
  if (modulus != 1) {
    if (n %% modulus != 0) return(0)
    n <- n %/% modulus
  }
  if (n == 0) return(0)
  f <- .factorint(n)
  f <- f[!(names(f) %in% c("-1", "0", "1"))]
  r <- prod(gmp::as.bigz(f + 1))
  if (proper) r <- r - 1
  .nt_out(r)
}

#' @rdname divisors
#' @export
proper_divisor_count <- function(n, modulus = 1) divisor_count(n, modulus, proper = TRUE)

#' @rdname divisors
#' @export
udivisors <- function(n) {
  n <- abs(.nt_int(n))
  if (n == 0) return(numeric())
  if (n == 1) return(1)
  f <- .factorint(n)
  pw <- .zc(lapply(names(f), function(k) gmp::as.bigz(k)^f[[k]]))
  d <- gmp::as.bigz(1)
  for (i in seq_along(pw)) d <- c(d, d * pw[i])
  .nt_out(.zsort(d))
}

#' @rdname divisors
#' @export
udivisor_count <- function(n) {
  n <- .nt_int(n)
  if (n == 0) return(0)
  f <- .factorint(n)
  .nt_out(gmp::as.bigz(2)^sum(!(names(f) %in% c("-1", "0", "1"))))
}

#' @rdname divisors
#' @export
antidivisors <- function(n) {
  n <- abs(.nt_int(n))
  if (n <= 2) return(numeric())
  d1 <- gmp::as.bigz(.divisors_vec(n)) * 2
  d1 <- d1[n > d1 & n %% d1 != 0]
  odd <- function(m) {
    d <- gmp::as.bigz(.divisors_vec(m))
    d[n > d & d >= 2 & n %% d != 0]
  }
  out <- c(d1, odd(2 * n - 1), odd(2 * n + 1))
  .nt_out(.zsort(out))
}

#' @rdname divisors
#' @export
antidivisor_count <- function(n) {
  n <- abs(.nt_int(n))
  if (n <= 2) return(0)
  divisor_count(2 * n - 1) + divisor_count(2 * n + 1) + divisor_count(n) -
    divisor_count(n, 2) - 5
}

.divisor_sigma <- function(n, k = 1) {
  f <- .factorint(n)
  if (k == 0) return(prod(gmp::as.bigz(f + 1)))
  r <- gmp::as.bigz(1)
  for (p in names(f)) {
    pz <- gmp::as.bigz(p)
    r <- r * ((pz^(k * (f[[p]] + 1)) - 1) %/% (pz^k - 1))
  }
  r
}

#' Square free and other cores, and digital roots
#'
#' `core(n, t)` is the `t`-free part of `n`: the product of `p^(e %% t)` over
#' the prime powers `p^e` of `n`. `dra(n, b)` is the additive digital root of
#' `n` in base `b` and `drm(n, b)` the multiplicative one.
#'
#' @param n A positive integer (any integer for `dra()` and `drm()`).
#' @param t Power, at least 2.
#' @param b Base, at least 2.
#' @return An integer.
#' @examples
#' core(24)
#' core(5^6 * 7^5, 3)
#' dra(19, 12)
#' drm(277777788888899, 10)
#' @export
core <- function(n, t = 2) {
  n <- .nt_int(n)
  t <- as.numeric(.nt_int(t))
  if (n <= 0) stop(.sympy_condition("ValueError", "n must be a positive integer"))
  if (t <= 1) stop(.sympy_condition("ValueError", "t must be >= 2"))
  f <- .factorint(n)
  y <- gmp::as.bigz(1)
  for (p in names(f)) y <- y * gmp::as.bigz(p)^(f[[p]] %% t)
  .nt_out(y)
}

#' @rdname core
#' @export
dra <- function(n, b) {
  num <- abs(.nt_int(n))
  b <- .nt_int(b)
  if (b <= 1) stop("Base should be an integer greater than 1", call. = FALSE)
  if (num == 0) return(0)
  .nt_out(1 + (num - 1) %% (b - 1))
}

#' @rdname core
#' @export
drm <- function(n, b) {
  n <- abs(.nt_int(n))
  b <- .nt_int(b)
  if (b <= 1) stop("Base should be an integer greater than 1", call. = FALSE)
  while (n > b) {
    mul <- gmp::as.bigz(1)
    while (n > 1) {
      r <- n %% b
      n <- n %/% b
      if (r == 0) return(0)
      mul <- mul * r
    }
    n <- mul
  }
  .nt_out(n)
}

# ---------------------------------------------------------------------------
# Perfect, abundant, amicable and Carmichael numbers

#' Perfect, abundant, deficient, amicable and Carmichael numbers
#'
#' `abundance(n)` is `sigma(n) - 2n`. `is_perfect()`, `is_abundant()` and
#' `is_deficient()` test its sign; `is_amicable(m, n)` tests an amicable pair
#' and `is_carmichael(n)` a Carmichael number. `find_carmichael_numbers_in_range(x, y)`
#' lists the Carmichael numbers in `[x, y)` and `find_first_n_carmichaels(n)` the
#' first `n`. `mersenne_prime_exponent(nth)` is the exponent of the `nth`
#' Mersenne prime, from the list of the 52 known ones.
#'
#' @param n,m Integers.
#' @param x,y Range bounds, `0 <= x <= y`.
#' @param nth Index from 1 to 52.
#' @return A logical, an integer or a vector.
#' @examples
#' is_perfect(28)
#' abundance(12)
#' is_amicable(220, 284)
#' find_first_n_carmichaels(3)
#' mersenne_prime_exponent(10)
#' @export
mersenne_prime_exponent <- function(nth) {
  n <- as.numeric(.nt_int(nth))
  if (n < 1) stop("nth must be a positive integer; mersenne_prime_exponent(1) == 2", call. = FALSE)
  if (n > length(.MERSENNE_EXPS))
    stop("There are only ", length(.MERSENNE_EXPS), " perfect numbers; nth must be less than or equal to ",
         length(.MERSENNE_EXPS), call. = FALSE)
  .MERSENNE_EXPS[n]
}

#' @rdname mersenne_prime_exponent
#' @export
is_perfect <- function(n) {
  n <- .nt_int(n)
  if (n < 1) return(FALSE)
  if (n %% 2 == 0) {
    m <- (.bitlen(n) + 1) %/% 2
    two <- gmp::as.bigz(2)
    if (two^(m - 1) * (two^m - 1) != n) return(FALSE)
    return(m %in% .MERSENNE_EXPS || is_mersenne_prime(two^m - 1))
  }
  if (n < gmp::as.bigz(10)^2000L) return(FALSE)
  if (n %% 105 == 0) return(FALSE)
  if (all(c(n %% 12 != 1, n %% 468 != 117, n %% 324 != 81))) return(FALSE)
  result <- .divisor_sigma(n) == 2 * n
  if (result) stop("an odd perfect number was found: ", paste(names(.factorint(n)), collapse = ", "),
                   call. = FALSE)
  result
}

#' @rdname mersenne_prime_exponent
#' @export
abundance <- function(n) {
  n <- .nt_int(n)
  .nt_out(.divisor_sigma(n) - 2 * n)
}

#' @rdname mersenne_prime_exponent
#' @export
is_abundant <- function(n) {
  n <- .nt_int(n)
  if (is_perfect(n)) return(FALSE)
  n %% 6 == 0 || .divisor_sigma(n) - 2 * n > 0
}

#' @rdname mersenne_prime_exponent
#' @export
is_deficient <- function(n) {
  n <- .nt_int(n)
  if (is_perfect(n)) return(FALSE)
  .divisor_sigma(n) - 2 * n < 0
}

#' @rdname mersenne_prime_exponent
#' @export
is_amicable <- function(m, n) {
  m <- .nt_int(m)
  n <- .nt_int(n)
  if (m == n) return(FALSE)
  s <- m + n
  s == .divisor_sigma(m) && s == .divisor_sigma(n)
}

.is_carmichael <- function(n) {
  if (n < 561 || n %% 2 == 0 || .is_prime(n)) return(FALSE)
  f <- .factorint(n)
  all(f == 1) && all((n - 1) %% (gmp::as.bigz(names(f)) - 1) == 0)
}

#' @rdname mersenne_prime_exponent
#' @export
is_carmichael <- function(n) .is_carmichael(.nt_int(n))

#' @rdname mersenne_prime_exponent
#' @export
find_carmichael_numbers_in_range <- function(x, y) {
  x <- .nt_int(x)
  y <- .nt_int(y)
  if (!(0 <= x && x <= y))
    stop("The provided range is not valid. x and y must be non-negative integers and x <= y", call. = FALSE)
  i <- if (x %% 2 == 0) x + 1 else x
  out <- list()
  while (i < y) {
    if (.is_carmichael(i)) out[[length(out) + 1L]] <- i
    i <- i + 2
  }
  .nt_out(.zc(out))
}

#' @rdname mersenne_prime_exponent
#' @export
find_first_n_carmichaels <- function(n) {
  n <- as.numeric(.nt_int(n))
  out <- numeric()
  i <- 561
  while (length(out) < n) {
    if (.is_carmichael(gmp::as.bigz(i))) out <- c(out, i)
    i <- i + 2
  }
  out
}
# ---------------------------------------------------------------------------
# Chinese remainders (sympy/ntheory/modular.py, galoistools gf_crt*)

.symres <- function(a, m) if (a <= m %/% 2) a else a - m

.crt_raw <- function(m, v) {
  p <- prod(m)
  r <- gmp::as.bigz(0)
  for (i in seq_along(m)) {
    e <- p %/% m[i]
    s <- gmp::gcdex(e, m[i])[2L]
    r <- r + e * ((v[i] * s) %% m[i])
  }
  r %% p
}

#' Chinese remainder theorem
#'
#' `crt(m, v)` returns `c(x, M)` with `x %% m[i] == v[i] %% m[i]` for every `i`
#' and `M` the product of the moduli, or `NULL` when the congruences are
#' inconsistent; `symmetric = TRUE` puts `x` within half of `M`. Moduli that are
#' not coprime are handled by [solve_congruence()], whose arguments are pairs
#' `c(remainder, modulus)`. `crt1(m)` precomputes for repeated use with
#' `crt2(m, v, mm, e, s)`. `symmetric_residue(a, m)` maps `a` in `[0, m)` to the
#' residue nearest zero.
#'
#' @param m Moduli.
#' @param v Remainders.
#' @param symmetric Return the residue nearest zero.
#' @param check Check the result, falling back to [solve_congruence()].
#' @param mm,e,s Values from `crt1()`.
#' @param a A residue.
#' @param ... Pairs `c(remainder, modulus)`.
#' @return `c(x, M)` or `NULL`; `crt1()` returns `list(mm, e, s)`.
#' @examples
#' crt(c(99, 97, 95), c(49, 76, 65))
#' solve_congruence(c(2, 3), c(3, 5), c(2, 7))
#' symmetric_residue(3, 6)
#' @export
crt <- function(m, v, symmetric = FALSE, check = TRUE) {
  m <- .nt_ints(m)
  v <- .nt_ints(v)
  res <- .crt(m, v, symmetric, check)
  if (is.null(res)) NULL else .nt_out(res)
}

.crt <- function(m, v, symmetric = FALSE, check = TRUE) {
  result <- .crt_raw(m, v)
  mm <- prod(m)
  if (check && !all(v %% m == result %% m)) {
    pairs <- lapply(seq_along(m), function(i) .zc(list(v[i], m[i])))
    r <- .solve_congruence(pairs, symmetric, FALSE)
    return(r)
  }
  if (symmetric) result <- .symres(result, mm)
  .zc(list(result, mm))
}

#' @rdname crt
#' @export
crt1 <- function(m) {
  m <- .nt_ints(m)
  p <- prod(m)
  E <- p %/% m
  S <- .zc(lapply(seq_along(m), function(i) gmp::gcdex(E[i], m[i])[2L] %% m[i]))
  list(mm = .nt_out(p), e = .nt_out(E), s = .nt_out(S))
}

#' @rdname crt
#' @export
crt2 <- function(m, v, mm, e, s, symmetric = FALSE) {
  m <- .nt_ints(m)
  v <- .nt_ints(v)
  mm <- .nt_int(mm)
  e <- .nt_ints(e)
  s <- .nt_ints(s)
  r <- sum(e * ((v * s) %% m)) %% mm
  if (symmetric) r <- .symres(r, mm)
  .nt_out(.zc(list(r, mm)))
}

#' @rdname crt
#' @export
symmetric_residue <- function(a, m) .nt_out(.symres(.nt_int(a), .nt_int(m)))

.solve_congruence <- function(rm, symmetric = FALSE, check = TRUE) {
  if (check) {
    keys <- character()
    vals <- list()
    for (pr in rm) {
      r <- pr[1L] %% pr[2L]
      k <- as.character(pr[2L])
      if (k %in% keys) {
        if (vals[[k]] != r) return(NULL)
        next
      }
      keys <- c(keys, k)
      vals[[k]] <- r
    }
    rm <- lapply(keys, function(k) .zc(list(vals[[k]], gmp::as.bigz(k))))
    if (all(vapply(rm, function(pr) .is_prime(pr[2L]), TRUE))) {
      return(.crt(.zc(lapply(rm, `[`, 2L)), .zc(lapply(rm, `[`, 1L)), symmetric, FALSE))
    }
  }
  n <- gmp::as.bigz(0)
  m <- gmp::as.bigz(1)
  for (pr in rm) {
    a1 <- n
    m1 <- m
    a <- m1
    b <- pr[1L] - a1
    mc <- pr[2L]
    g <- .gcd2(.gcd2(a, b), mc)
    a <- a %/% g
    b <- b %/% g
    mc <- mc %/% g
    if (a != 1) {
      gx <- gmp::gcdex(a, mc)
      if (gx[1L] != 1) return(NULL)
      b <- b * gx[2L]
    }
    n <- a1 + m1 * b
    m <- m1 * mc
    n <- n %% m
  }
  if (symmetric) n <- .symres(n, m)
  .zc(list(n, m))
}

#' @rdname crt
#' @export
solve_congruence <- function(..., symmetric = FALSE, check = TRUE) {
  rm <- lapply(list(...), .nt_ints)
  r <- .solve_congruence(rm, symmetric, check)
  if (is.null(r)) NULL else .nt_out(r)
}

# ---------------------------------------------------------------------------
# Residues (sympy/ntheory/residue_ntheory.py)

.n_order <- function(a, n) {
  if (n <= 1) stop(.sympy_condition("ValueError", "n should be an integer greater than 1"))
  a <- a %% n
  if (a == 1) return(gmp::as.bigz(1))
  if (.gcd2(a, n) != 1) stop(.sympy_condition("ValueError", "The two numbers should be relatively prime"))
  a_order <- gmp::as.bigz(1)
  f <- .factorint(n)
  for (pk in names(f)) {
    p <- gmp::as.bigz(pk)
    e <- f[[pk]]
    pe <- p^e
    pe_order <- (p - 1) * p^(e - 1)
    fs <- .factorint(p - 1)
    if (e > 1) fs[[pk]] <- e - 1
    order <- gmp::as.bigz(1)
    for (qk in names(fs)) {
      px <- gmp::as.bigz(qk)
      x <- gmp::powm(a, pe_order %/% px^fs[[qk]], pe)
      while (x != 1) {
        x <- gmp::powm(x, px, pe)
        order <- order * px
      }
    }
    a_order <- a_order * order %/% .gcd2(a_order, order)
  }
  a_order
}

#' Orders and primitive roots
#'
#' `n_order(a, n)` is the multiplicative order of `a` modulo `n`.
#' `is_primitive_root(a, p)` tests whether `a` generates the units modulo `p`,
#' and `primitive_root(p)` returns the smallest primitive root, or with
#' `smallest = FALSE` one found quickly; it returns `NULL` when none exists.
#'
#' @param a An integer coprime to the modulus.
#' @param n,p Modulus, greater than 1.
#' @param smallest Return the smallest primitive root.
#' @return An integer, logical or `NULL`.
#' @examples
#' n_order(3, 7)
#' is_primitive_root(3, 10)
#' primitive_root(19)
#' primitive_root(50)
#' @export
n_order <- function(a, n) .nt_out(.n_order(.nt_int(a), .nt_int(n)))

# Primitive roots of an odd prime p, in SymPy's order: the first one, then its
# powers with exponents coprime to p - 1. `take` bounds the count.
.prim_roots_prime <- function(p, take = Inf) {
  if (p == 3) return(list(gmp::as.bigz(2)))
  g_min <- if (as.numeric(p %% 8) %in% c(1, 7)) 3 else 2
  if (p < 41) {
    g <- gmp::as.bigz(if (p == 23) 5 else g_min)
  } else {
    v <- (p - 1) %/% .fac_primes(.factorint(p - 1))
    g <- gmp::as.bigz(g_min)
    while (g < p) {
      if (all(vapply(seq_along(v), function(i) gmp::powm(g, v[i], p) != 1, TRUE))) break
      g <- g + 1
    }
  }
  out <- list(g)
  k <- gmp::as.bigz(3)
  while (length(out) < take && k < p) {
    if (.gcd2(p - 1, k) == 1) out[[length(out) + 1L]] <- gmp::powm(g, k, p)
    k <- k + 2
  }
  out
}

.prim_root_prime <- function(p) .prim_roots_prime(p, 1)[[1L]]

.is_primitive_root <- function(a, p) {
  if (p <= 1) stop(.sympy_condition("ValueError", "p should be an integer greater than 1"))
  a <- a %% p
  if (.gcd2(a, p) != 1) stop(.sympy_condition("ValueError", "The two numbers should be relatively prime"))
  if (p <= 4) return(a == p - 1)
  if (p %% 2 == 1) q <- p else if (p %% 4 != 0) q <- p %/% 2 else return(FALSE)
  if (.is_prime(q)) {
    group_order <- q - 1
    fs <- .fac_primes(.factorint(q - 1))
  } else {
    m <- .perfect_power_int(q, 3)
    if (isFALSE(m)) return(FALSE)
    q <- gmp::as.bigz(m[[1L]])
    e <- m[[2L]]
    if (!.is_prime(q)) return(FALSE)
    group_order <- q^(e - 1) * (q - 1)
    fs <- .fac_primes(.factorint(q - 1))
    if (!(as.character(q) %in% as.character(fs))) fs <- c(fs, q)
  }
  all(vapply(seq_along(fs), function(i) gmp::powm(a, group_order %/% fs[i], p) != 1, TRUE))
}

#' @rdname n_order
#' @export
is_primitive_root <- function(a, p) .is_primitive_root(.nt_int(a), .nt_int(p))

.primitive_root <- function(p, smallest = TRUE) {
  if (p <= 1) stop(.sympy_condition("ValueError", "p should be an integer greater than 1"))
  if (p <= 4) return(p - 1)
  p_even <- p %% 2 == 0
  if (!p_even) q <- p else if (p %% 4 != 0) q <- p %/% 2 else return(NULL)
  if (.is_prime(q)) {
    e <- 1
  } else {
    m <- .perfect_power_int(q, 3)
    if (isFALSE(m)) return(NULL)
    q <- gmp::as.bigz(m[[1L]])
    e <- m[[2L]]
    if (!.is_prime(q)) return(NULL)
  }
  if (!smallest) {
    g <- .prim_root_prime(q)
    if (e > 1) {
      p2 <- q^2L
      t <- (g - .powm(g, 2 - q, p2)) %% p2
      g <- if (t == 0) g + q else g
    }
    if (p_even && g %% 2 == 0) g <- g + q^e
    return(g)
  }
  if (p_even) {
    i <- gmp::as.bigz(3)
    while (i < p) {
      if (i %% q != 0 && .is_primitive_root(i, p)) return(i)
      i <- i + 2
    }
  }
  g <- .prim_root_prime(q)
  if (e == 1 || gmp::powm(g, q - 1, q^2L) != 1) return(g)
  i <- g + 1
  while (i < p) {
    if (i %% q != 0 && .is_primitive_root(i, p)) return(i)
    i <- i + 1
  }
  NULL
}

#' @rdname n_order
#' @export
primitive_root <- function(p, smallest = TRUE) {
  r <- .primitive_root(.nt_int(p), smallest)
  if (is.null(r)) NULL else .nt_out(r)
}

.sqrt_mod_tonelli_shanks <- function(a, p) {
  s <- .v2(p - 1)
  t <- p %/% gmp::as.bigz(2)^s
  if (p %% 3 == 2) {
    d <- gmp::as.bigz(3)
  } else if (as.numeric(p %% 5) %in% c(2, 3)) {
    d <- gmp::as.bigz(5)
  } else {
    d <- gmp::as.bigz(7)
    while (.jacobi(p, d) != -1L) d <- gmp::nextprime(d)
  }
  A <- gmp::powm(a, t, p)
  D <- gmp::powm(d, t, p)
  m <- gmp::as.bigz(0)
  two <- gmp::as.bigz(2)
  for (i in 0:(s - 1)) {
    adm <- (A * gmp::powm(D, m, p)) %% p
    adm <- gmp::powm(adm, two^(s - 1 - i), p)
    if (adm == p - 1) m <- m + two^i
  }
  (gmp::powm(a, (t + 1) %/% 2, p) * gmp::powm(D, m %/% 2, p)) %% p
}

# Square roots of a (coprime to p) modulo p^k, sorted, or NULL.
.sqrt_mod_prime_power <- function(a, p, k) {
  pk <- p^k
  a <- a %% pk
  if (p == 2) {
    if (a %% 8 != 1) return(NULL)
    if (k <= 3) return(gmp::as.bigz(seq(1, as.numeric(pk) - 1, by = 2)))
    r <- gmp::as.bigz(1)
    two <- gmp::as.bigz(2)
    for (nx in 3:(k - 1)) if (((r^2L - a) %/% two^nx) %% 2 == 1) r <- r + two^(nx - 1)
    h <- two^(k - 1)
    return(.zsort(.zc(list(r, pk - r, (r + h) %% pk, (-(r + h)) %% pk))))
  }
  if (.jacobi(a, p) != 1L) return(NULL)
  if (p %% 4 == 3) {
    res <- gmp::powm(a, (p + 1) %/% 4, p)
  } else if (p %% 8 == 5) {
    res <- gmp::powm(a, (p + 3) %/% 8, p)
    if (gmp::powm(res, 2L, p) != a %% p) res <- (res * gmp::powm(gmp::as.bigz(2), (p - 1) %/% 4, p)) %% p
  } else {
    res <- .sqrt_mod_tonelli_shanks(a, p)
  }
  if (k > 1) {
    px <- p
    for (i in seq_len(.bitlen(gmp::as.bigz(k)) - 1L)) {
      px <- px^2L
      res <- (res - (res^2L - a) * .zinv(2 * res, px)) %% px
    }
    if (bitwAnd(k, k - 1) != 0) res <- (res - (res^2L - a) * .zinv(2 * res, pk)) %% pk
  }
  .zsort(.zc(list(res, pk - res)))
}

# Root sets are lists of arithmetic progressions list(start, step, count),
# so that a single root can be taken without listing p^(n/2) of them.
.prog <- function(start, step, count) list(start = gmp::as.bigz(start), step = gmp::as.bigz(step),
                                           count = gmp::as.bigz(count))

.prog_set <- function(v) lapply(seq_along(v), function(i) .prog(v[i], 0, 1))

.set_len <- function(set) sum(.zc(lapply(set, `[[`, "count")))

.set_get <- function(set, j) {
  for (pr in set) {
    if (j < pr$count) return(pr$start + pr$step * j)
    j <- j - pr$count
  }
  stop("index out of range", call. = FALSE)
}

.set_expand <- function(set) {
  if (.set_len(set) > 1e7) stop("too many roots to list", call. = FALSE)
  .zc(lapply(set, function(pr) pr$start + pr$step * gmp::as.bigz(seq_len(as.numeric(pr$count)) - 1)))
}

.sqrt_mod1 <- function(a, p, n) {
  pn <- p^n
  a <- a %% pn
  if (a == 0) {
    st <- p^((n + 1) %/% 2)
    return(list(.prog(0, st, pn %/% st)))
  }
  rr <- .remove(a, p)
  a <- rr[[1L]]
  r <- rr[[2L]]
  if (r %% 2 == 1) return(NULL)
  res <- .sqrt_mod_prime_power(a, p, n - r)
  if (is.null(res)) return(NULL)
  m <- r %/% 2
  lapply(seq_along(res), function(i) .prog(res[i] * p^m, p^(n - m), p^m))
}

# Root sets and CRT data for x^2 = a (mod p), or NULL when there is no root.
.sqrt_mod_sets <- function(a, p) {
  p <- abs(p)
  f <- .factorint(p)
  sets <- list()
  pv <- list()
  for (pk in names(f)) {
    px <- gmp::as.bigz(pk)
    ex <- f[[pk]]
    rx <- if (a %% px != 0) {
      r <- .sqrt_mod_prime_power(a, px, ex)
      if (is.null(r)) NULL else .prog_set(r)
    } else .sqrt_mod1(a, px, ex)
    if (is.null(rx) || !length(rx)) return(NULL)
    sets[[length(sets) + 1L]] <- rx
    pv[[length(pv) + 1L]] <- px^ex
  }
  pv <- .zc(pv)
  mm <- prod(pv)
  E <- mm %/% pv
  S <- .zc(lapply(seq_along(pv), function(i) gmp::gcdex(E[i], pv[i])[2L] %% pv[i]))
  list(sets = sets, pv = pv, mm = mm, E = E, S = S)
}

# All square roots of a modulo p in SymPy's iteration order (last factor fastest).
.sqrt_mod_all <- function(a, p) {
  d <- .sqrt_mod_sets(a, p)
  if (is.null(d)) return(.z0())
  if (length(d$sets) == 1L) return(.set_expand(d$sets[[1L]]))
  acc <- gmp::as.bigz(0)
  for (i in seq_along(d$sets)) {
    term <- d$E[i] * ((.set_expand(d$sets[[i]]) * d$S[i]) %% d$pv[i])
    acc <- do.call(c, lapply(seq_along(acc), function(j) acc[j] + term))
  }
  acc %% d$mm
}

#' Square roots modulo n
#'
#' `sqrt_mod(a, p)` returns one root of `x^2 = a (mod p)`, at most `p / 2`, or
#' `NULL`; `all_roots = TRUE` returns every root in increasing order.
#' `sqrt_mod_iter(a, p)` returns every root in the order SymPy's iterator
#' produces them. `quadratic_residues(p)` lists the squares modulo `p`, and
#' `is_quad_residue(a, p)` tests one.
#'
#' @param a An integer.
#' @param p Modulus.
#' @param all_roots Return all roots.
#' @param domain Function applied to each root, for example [Integer()].
#' @return An integer, a vector, `NULL` or a logical.
#' @examples
#' sqrt_mod(11, 43)
#' sqrt_mod(17, 32, all_roots = TRUE)
#' quadratic_residues(13)
#' is_quad_residue(21, 100)
#' @export
sqrt_mod <- function(a, p, all_roots = FALSE) {
  a <- .nt_int(a)
  p <- abs(.nt_int(p))
  if (all_roots) return(.nt_out(.zsort(.sqrt_mod_all(a, p))))
  d <- .sqrt_mod_sets(a, p)
  if (is.null(d)) return(NULL)
  k <- length(d$sets)
  lens <- .zc(lapply(d$sets, .set_len))
  idx <- rep(list(gmp::as.bigz(0)), k)
  halfp <- p %/% 2
  x <- NULL
  repeat {
    r <- if (k == 1L) .set_get(d$sets[[1L]], idx[[1L]]) else
      sum(.zc(lapply(seq_len(k), function(i) d$E[i] * ((.set_get(d$sets[[i]], idx[[i]]) * d$S[i]) %% d$pv[i])))) %% d$mm
    if (r < halfp) return(.nt_out(r))
    if (r > halfp) return(.nt_out(p - r))
    x <- r
    i <- k
    while (i >= 1L) {
      idx[[i]] <- idx[[i]] + 1
      if (idx[[i]] < lens[i]) break
      idx[[i]] <- gmp::as.bigz(0)
      i <- i - 1L
    }
    if (i < 1L) break
  }
  .nt_out(x)
}

#' @rdname sqrt_mod
#' @export
sqrt_mod_iter <- function(a, p, domain = NULL) {
  roots <- .nt_out(.sqrt_mod_all(.nt_int(a), .nt_int(p)))
  if (is.null(domain)) return(roots)
  lapply(seq_along(roots), function(i) domain(roots[i]))
}

#' @rdname sqrt_mod
#' @export
quadratic_residues <- function(p) {
  p <- as.numeric(.nt_int(p))
  i <- gmp::as.bigz(0:(p %/% 2))
  .nt_out(.zsort(unique(i^2L %% p)))
}

.is_quad_residue <- function(a, p) {
  if (p < 1) stop("p must be > 0", call. = FALSE)
  a <- a %% p
  if (a < 2 || p < 3) return(TRUE)
  t <- .v2(p)
  if (t) {
    a_ <- a %% gmp::as.bigz(2)^t
    if (a_ != 0) {
      r <- .v2(a_)
      if (r %% 2 == 1 || (a_ %/% gmp::as.bigz(2)^r) %% 8 %in% c(3, 5, 7)) return(FALSE)
    }
    p <- p %/% gmp::as.bigz(2)^t
    a <- a %% p
    if (a < 2 || p < 3) return(TRUE)
  }
  j <- .jacobi(a, p)
  if (j == -1L || .is_prime(p)) return(j == 1L)
  f <- .factorint(p)
  for (pk in names(f)) {
    px <- gmp::as.bigz(pk)
    if (a %% px != 0) {
      if (.jacobi(a, px) != 1L) return(FALSE)
    } else {
      a_ <- a %% px^f[[pk]]
      if (a_ == 0) next
      rr <- .remove(a_, px)
      if (rr[[2L]] %% 2 == 1 || .jacobi(rr[[1L]], px) != 1L) return(FALSE)
    }
  }
  TRUE
}

#' @rdname sqrt_mod
#' @export
is_quad_residue <- function(a, p) .is_quad_residue(.nt_int(a), .nt_int(p))

.is_nthpow_pp <- function(a, n, p, k) {
  while (a %% p == 0) {
    a <- a %% p^k
    if (a == 0) return(TRUE)
    rr <- .remove(a, p)
    a <- rr[[1L]]
    mu <- rr[[2L]]
    if (mu %% n != 0) return(FALSE)
    k <- k - mu
  }
  if (p != 2) {
    f <- p^(k - 1) * (p - 1)
    return(gmp::powm(a, f %/% .gcd2(f, n), p^k) == 1)
  }
  if (n %% 2 == 1) return(TRUE)
  c <- min(.v2(n) + 2, k)
  a %% gmp::as.bigz(2)^c == 1
}

#' Power residues and n-th roots modulo m
#'
#' `is_nthpow_residue(a, n, m)` tests whether `x^n = a (mod m)` has a solution,
#' and `nthroot_mod(a, n, p)` returns the smallest solution, every solution
#' with `all_roots = TRUE`, or `NULL` (an empty vector) when there is none.
#'
#' @param a An integer.
#' @param n Exponent.
#' @param m,p Modulus.
#' @param all_roots Return all roots.
#' @return A logical, an integer, a vector or `NULL`.
#' @examples
#' is_nthpow_residue(2, 3, 7)
#' nthroot_mod(11, 4, 19)
#' nthroot_mod(68, 3, 109, all_roots = TRUE)
#' @export
is_nthpow_residue <- function(a, n, m) {
  a <- .nt_int(a)
  n <- .nt_int(n)
  m <- .nt_int(m)
  if (m <= 0) stop("m must be > 0", call. = FALSE)
  if (n < 0) stop("n must be >= 0", call. = FALSE)
  a <- a %% m
  if (n == 0) return(if (m == 1) FALSE else a == 1)
  if (a == 0 || n == 1) return(TRUE)
  if (n == 2) return(.is_quad_residue(a, m))
  f <- .factorint(m)
  all(vapply(names(f), function(pk) .is_nthpow_pp(a, n, gmp::as.bigz(pk), f[[pk]]), TRUE))
}

.nthroot_mod1 <- function(s, q, p) {
  g <- .prim_root_prime(p)
  r <- s
  f <- .factorint(q)
  for (qk in names(f)) {
    qx <- gmp::as.bigz(qk)
    ex <- f[[qk]]
    fx <- (p - 1) %/% qx^ex
    while (fx %% qx == 0) fx <- fx %/% qx
    z <- fx * .zinv(-fx, qx)
    x <- (1 + z) %/% qx
    t <- .discrete_log(p, gmp::powm(r, fx, p), gmp::powm(g, fx * qx, p))
    for (i in seq_len(ex)) {
      r <- (gmp::powm(r, x, p) * gmp::powm(g, (-z * t) %% (p - 1), p)) %% p
      t <- t %/% qx
    }
  }
  h <- gmp::powm(g, (p - 1) %/% q, p)
  res <- list(r)
  hx <- r
  for (i in seq_len(as.numeric(q) - 1)) {
    hx <- (hx * h) %% p
    res[[length(res) + 1L]] <- hx
  }
  .zsort(.zc(res))
}

.nthroot_mod_pp <- function(a, n, p, k) {
  if (!.is_nthpow_pp(a, n, p, k)) return(.z0())
  amp <- a %% p
  if (amp == 0) {
    base <- gmp::as.bigz(0)
  } else if ((p - 1) %% n == 0) {
    base <- .nthroot_mod1(amp, n, p)
  } else {
    pa <- n
    pb <- p - 1
    b <- gmp::as.bigz(1)
    if (pa < pb) {
      tmp <- amp; amp <- b; b <- tmp
      tmp <- pa; pa <- pb; pb <- tmp
    }
    while (pb != 0) {
      q <- pa %/% pb
      pc <- pa %% pb
      cc <- (.powm(b, -q, p) * amp) %% p
      pa <- pb
      pb <- pc
      amp <- b
      b <- cc
    }
    base <- if (pa == 1) amp else if (pa == 2) .zsort(.sqrt_mod_all(amp, p)) else .nthroot_mod1(amp, pa, p)
  }
  if (k == 1) return(base)
  a <- a %% p^k
  tot <- list()
  for (i in seq_along(base)) {
    root <- base[i]
    d <- (gmp::powm(root, n - 1, p) * n) %% p
    new_base <- p
    if (d != 0) {
      m_inv <- .zinv(d, p)
      for (j in seq_len(k - 1)) {
        new_base <- new_base * p
        tmp <- (gmp::powm(root, n, new_base) - a) * m_inv
        root <- (root - tmp) %% new_base
      }
      tot[[length(tot) + 1L]] <- root
    } else {
      rib <- root
      for (j in seq_len(k - 1)) {
        new_base <- new_base * p
        nr <- list()
        for (t in seq_along(rib)) {
          k_ <- rib[t]
          if (gmp::powm(k_, n, new_base) != a %% new_base) next
          step <- new_base %/% p
          nr <- c(nr, lapply(0:(as.numeric(p) - 1), function(s) (k_ + s * step) %% new_base))
        }
        rib <- .zc(nr)
        if (length(rib)) rib <- rib[!duplicated(as.character(rib))]
      }
      tot <- c(tot, if (length(rib)) lapply(seq_along(rib), function(t) rib[t]))
    }
  }
  r <- .zc(tot)
  if (!length(r)) return(r)
  .zsort(r[!duplicated(as.character(r))])
}

#' @rdname is_nthpow_residue
#' @export
nthroot_mod <- function(a, n, p, all_roots = FALSE) {
  a <- .nt_int(a)
  n <- .nt_int(n)
  p <- .nt_int(p)
  if (n < 1) stop(.sympy_condition("ValueError", "n should be positive"))
  if (p < 1) stop("p should be positive", call. = FALSE)
  a <- a %% p
  if (n == 1) return(.nt_out(a))
  if (n == 2) return(sqrt_mod(a, p, all_roots))
  f <- .factorint(p)
  base <- list()
  pp <- list()
  for (qk in names(f)) {
    r <- .nthroot_mod_pp(a, n, gmp::as.bigz(qk), f[[qk]])
    if (!length(r)) return(if (all_roots) numeric() else NULL)
    base[[length(base) + 1L]] <- r
    pp[[length(pp) + 1L]] <- gmp::as.bigz(qk)^f[[qk]]
  }
  pp <- .zc(pp)
  mm <- prod(pp)
  E <- mm %/% pp
  S <- .zc(lapply(seq_along(pp), function(i) gmp::gcdex(E[i], pp[i])[2L] %% pp[i]))
  acc <- gmp::as.bigz(0)
  for (i in seq_along(base)) {
    term <- E[i] * ((base[[i]] * S[i]) %% pp[i])
    acc <- do.call(c, lapply(seq_along(acc), function(j) acc[j] + term))
  }
  acc <- acc %% mm
  ret <- .zsort(acc[!duplicated(as.character(acc))])
  if (all_roots) return(.nt_out(ret))
  .nt_out(ret[1L])
}

# ---------------------------------------------------------------------------
# Discrete logarithms

.dlog_trial <- function(n, a, b, order = NULL) {
  a <- a %% n
  b <- b %% n
  if (is.null(order)) order <- n
  x <- gmp::as.bigz(1)
  i <- 0
  while (i < order) {
    if (x == a) return(gmp::as.bigz(i))
    x <- (x * b) %% n
    i <- i + 1
  }
  stop(.sympy_condition("ValueError", "Log does not exist"))
}

.dlog_shanks <- function(n, a, b, order = NULL) {
  a <- a %% n
  b <- b %% n
  if (is.null(order)) order <- .n_order(b, n)
  m <- as.numeric(.isqrt(order)) + 1
  T <- new.env(hash = TRUE, parent = emptyenv())
  x <- gmp::as.bigz(1)
  for (i in seq_len(m) - 1) {
    k <- as.character(x)
    if (is.null(T[[k]])) assign(k, i, envir = T)
    x <- (x * b) %% n
  }
  z <- .powm(b, -m, n)
  x <- a
  for (i in seq_len(m) - 1) {
    hit <- T[[as.character(x)]]
    if (!is.null(hit)) return(gmp::as.bigz(i) * m + hit)
    x <- (x * z) %% n
  }
  stop(.sympy_condition("ValueError", "Log does not exist"))
}

.dlog_rho <- function(n, a, b, order = NULL, retries = 10, rseed = NULL) {
  a <- a %% n
  b <- b %% n
  if (is.null(order)) order <- .n_order(b, n)
  randint <- .nt_rng(rseed)
  step <- function(x, e1, e2) {
    c <- as.numeric(x %% 3)
    if (c == 0) list((a * x) %% n, e1, (e2 + 1) %% order)
    else if (c == 1) list((x * x) %% n, (e1 + e1) %% order, (e2 + e2) %% order)
    else list((b * x) %% n, (e1 + 1) %% order, e2)
  }
  for (i in seq_len(retries)) {
    aa <- randint(1, order - 1)
    ba <- randint(1, order - 1)
    xa <- (gmp::powm(b, aa, n) * gmp::powm(a, ba, n)) %% n
    s <- step(xa, aa, ba)
    xb <- s[[1L]]
    ab <- s[[2L]]
    bb <- s[[3L]]
    j <- gmp::as.bigz(0)
    while (j < order) {
      s <- step(xa, aa, ba)
      xa <- s[[1L]]; aa <- s[[2L]]; ba <- s[[3L]]
      s <- step(xb, ab, bb)
      s <- step(s[[1L]], s[[2L]], s[[3L]])
      xb <- s[[1L]]; ab <- s[[2L]]; bb <- s[[3L]]
      if (xa == xb) {
        r <- (ba - bb) %% order
        ri <- .zinv(r, order)
        if (!is.na(ri)) {
          e <- (ri * (ab - aa)) %% order
          if ((gmp::powm(b, e, n) - a) %% n == 0) return(e)
        }
        break
      }
      j <- j + 1
    }
  }
  stop("Pollard's Rho failed to find logarithm", call. = FALSE)
}

.dlog_smooth <- function(n, fb) {
  if (n == 0) return(NULL)
  out <- numeric(length(fb))
  for (i in seq_along(fb)) {
    while (n %% fb[i] == 0) {
      out[i] <- out[i] + 1
      n <- n %/% fb[i]
    }
  }
  if (n != 1) NULL else out
}

.dlog_index_calculus <- function(n, a, b, order, rseed = NULL) {
  randint <- .nt_rng(rseed)
  a <- a %% n
  b <- b %% n
  ln <- log(n)
  B <- trunc(exp(0.5 * sqrt(ln * log(ln)) * (1 + 1 / log(ln))))
  mx <- 5 * B * B
  fb <- .primerange(gmp::as.bigz(2), gmp::as.bigz(B))
  lf <- length(fb)
  abx <- a
  rela <- NULL
  x <- gmp::as.bigz(0)
  while (x < order) {
    if (abx == 1) return((order - x) %% order)
    r <- .dlog_smooth(abx, fb)
    if (!is.null(r)) {
      rela <- c(lapply(r, function(v) gmp::as.bigz(v) %% order), list(x))
      break
    }
    abx <- (abx * b) %% n
    x <- x + 1
  }
  if (is.null(rela)) stop("Index Calculus failed", call. = FALSE)
  rela <- .zc(rela)
  rels <- vector("list", lf)
  k <- 1
  kk <- 0
  while (k < 3 * lf && kk < mx) {
    x <- randint(1, order - 1)
    rel <- .dlog_smooth(gmp::powm(b, x, n), fb)
    if (is.null(rel)) {
      kk <- kk + 1
      next
    }
    k <- k + 1
    kk <- 0
    rel <- c(gmp::as.bigz(rel), x)
    index <- lf + 1
    for (i in seq_len(lf)) {
      ri <- rel[i] %% order
      if (ri > 0 && !is.null(rels[[i]])) rel <- (rel - ri * rels[[i]]) %% order
      else rel[i] <- ri
      if (rel[i] > 0 && index == lf + 1) index <- i
    }
    if (index == lf + 1 || !is.null(rels[[index]])) next
    rinv <- .zinv(rel[index], order)
    idx <- index:(lf + 1)
    rel[idx] <- (rinv * rel[idx]) %% order
    rels[[index]] <- rel
    complete <- TRUE
    for (i in seq_len(lf)) {
      if (rela[i] > 0 && !is.null(rels[[i]])) rela <- (rela - rela[i] * rels[[i]]) %% order
      if (rela[i] > 0) {
        complete <- FALSE
        break
      }
    }
    if (complete) {
      x <- (order - rela[lf + 1]) %% order
      if (gmp::powm(b, x, n) == a) return(x)
      stop("Index Calculus failed", call. = FALSE)
    }
  }
  stop("Index Calculus failed", call. = FALSE)
}

.dlog_pohlig_hellman <- function(n, a, b, order = NULL, order_factors = NULL) {
  a <- a %% n
  b <- b %% n
  if (is.null(order)) order <- .n_order(b, n)
  if (is.null(order_factors)) order_factors <- .factorint(order)
  ps <- .fac_primes(order_factors)
  l <- rep(list(gmp::as.bigz(0)), length(ps))
  for (i in seq_along(ps)) {
    pi <- ps[i]
    for (j in seq_len(order_factors[[i]]) - 1) {
      aj <- gmp::powm((a * .powm(b, -l[[i]], n)) %% n, order %/% pi^(j + 1), n)
      bj <- gmp::powm(b, order %/% pi, n)
      cj <- .discrete_log(n, aj, bj, pi, TRUE)
      l[[i]] <- l[[i]] + cj * pi^j
    }
  }
  mods <- .zc(lapply(seq_along(ps), function(i) ps[i]^order_factors[[i]]))
  .crt_raw(mods, .zc(l))
}

.discrete_log <- function(n, a, b, order = NULL, prime_order = NULL) {
  if (n < 1) stop(.sympy_condition("ValueError", "n should be positive"))
  if (n == 1) return(gmp::as.bigz(0))
  # powers of a unit are units
  if (.gcd2(b, n) == 1 && .gcd2(a, n) != 1) stop(.sympy_condition("ValueError", "Log does not exist"))
  order_factors <- NULL
  if (is.null(order)) {
    fac <- .fac_new()
    fn <- .factorint(n)
    for (pk in names(fn)) {
      if (fn[[pk]] > 1) .fac_add(fac, pk, fn[[pk]] - 1)
      .fac_merge(fac, .factorint(gmp::as.bigz(pk) - 1))
    }
    factors <- .fac_vec(fac)
    order <- prod(gmp::as.bigz(names(factors))^factors)
    if (!length(factors)) order <- gmp::as.bigz(1)
    of <- c()
    for (pk in names(factors)) {
      p <- gmp::as.bigz(pk)
      i <- 0
      for (r in seq_len(factors[[pk]])) {
        if (gmp::powm(b, order %/% p, n) == 1) {
          order <- order %/% p
          i <- i + 1
        } else break
      }
      if (i < factors[[pk]]) of[pk] <- factors[[pk]] - i
    }
    order_factors <- of
  }
  if (is.null(prime_order)) prime_order <- .is_prime(order)
  if (order < 1000) return(.dlog_trial(n, a, b, order))
  if (prime_order) {
    if (4 * sqrt(log(n) * log(log(n))) < log(order) - 10) return(.dlog_index_calculus(n, a, b, order))
    if (order < 1e12) return(.dlog_shanks(n, a, b, order))
    return(.dlog_rho(n, a, b, order))
  }
  .dlog_pohlig_hellman(n, a, b, order, order_factors)
}

#' Discrete logarithm
#'
#' `discrete_log(n, a, b)` finds `x` with `b^x = a (mod n)`, choosing among
#' trial multiplication, baby step giant step, Pollard's rho, index calculus
#' and Pohlig-Hellman as SymPy does. `order` is the order of `b` if known, and
#' `prime_order` states whether it is prime.
#'
#' @param n Modulus.
#' @param a Target.
#' @param b Base.
#' @param order Order of `b`, or `NULL`.
#' @param prime_order Whether `order` is prime, or `NULL`.
#' @return An integer.
#' @examples
#' discrete_log(41, 15, 7)
#' discrete_log(587, 2^9, 2)
#' @export
discrete_log <- function(n, a, b, order = NULL, prime_order = NULL) {
  if (!is.null(order)) order <- .nt_int(order)
  .nt_out(.discrete_log(.nt_int(n), .nt_int(a), .nt_int(b), order, prime_order))
}

# ---------------------------------------------------------------------------
# Congruences

.linear_congruence <- function(a, b, m) {
  if (a %% m == 0) {
    if (b %% m == 0) return(gmp::as.bigz(seq(0, as.numeric(m) - 1)))
    return(.z0())
  }
  gx <- gmp::gcdex(a, m)
  g <- gx[1L]
  r <- gx[2L]
  if (b %% g != 0) return(.z0())
  t <- gmp::as.bigz(seq(0, as.numeric(g) - 1))
  .zsort((r * b %/% g + t * m %/% g) %% m)
}

.quadratic_congruence <- function(a, b, c, n) {
  if (n <= 1) stop(.sympy_condition("ValueError", "n should be an integer greater than 1"))
  a <- a %% n
  b <- b %% n
  c <- c %% n
  if (a == 0) return(.linear_congruence(b, -c, n))
  if (n == 2) {
    roots <- list()
    if (c == 0) roots[[1L]] <- gmp::as.bigz(0)
    if ((b + c) %% 2 == 1) roots[[length(roots) + 1L]] <- gmp::as.bigz(1)
    return(.zc(roots))
  }
  if (.gcd2(2 * a, n) == 1) {
    inv_a <- .zinv(a, n)
    b <- b * inv_a
    c <- c * inv_a
    if (b %% 2 == 1) b <- b + n
    b <- b %/% 2
    r <- .sqrt_mod_all(b^2L - c, n)
    if (!length(r)) return(.z0())
    return(.zsort((r - b) %% n))
  }
  r <- .sqrt_mod_all(b^2L - 4 * a * c, 4 * a * n)
  if (!length(r)) return(.z0())
  num <- r - b
  ok <- num %% (2 * a) == 0
  res <- (num[ok] %/% (2 * a)) %% n
  if (!length(res)) return(.z0())
  .zsort(res[!duplicated(as.character(res))])
}

#' Polynomial congruences
#'
#' `quadratic_congruence(a, b, c, n)` solves `a x^2 + b x + c = 0 (mod n)`.
#' `polynomial_congruence(expr, m)` solves a univariate polynomial with integer
#' coefficients modulo `m`. `binomial_mod(n, m, k)` is `choose(n, m) %% k`,
#' by Granville's extension of Lucas' theorem.
#'
#' @param a,b,c Coefficients.
#' @param n,m,k Integers; `n` and `m` of `quadratic_congruence()` and
#'   `polynomial_congruence()` are moduli.
#' @param expr A `sym` polynomial or a string.
#' @return A vector of solutions in increasing order, or an integer.
#' @examples
#' quadratic_congruence(2, 5, 3, 7)
#' polynomial_congruence("x^6 - 2*x^5 - 35", 6125)
#' binomial_mod(10, 2, 6)
#' @export
quadratic_congruence <- function(a, b, c, n) {
  .nt_out(.quadratic_congruence(.nt_int(a), .nt_int(b), .nt_int(c), .nt_int(n)))
}

# Coefficients of a polynomial (highest first) over GF(p), evaluated at x.
.gf_eval <- function(f, x, p) {
  r <- gmp::as.bigz(0)
  for (i in seq_along(f)) r <- (r * x + f[i]) %% p
  r
}

.gf_trim_z <- function(f) {
  nz <- which(f != 0)
  if (!length(nz)) return(.z0())
  f[nz[1L]:length(f)]
}

.gf_rem_z <- function(f, g, p) {
  f <- .gf_trim_z(f %% p)
  g <- .gf_trim_z(g %% p)
  dg <- length(g)
  inv <- .zinv(g[1L], p)
  while (length(f) >= dg) {
    q <- (f[1L] * inv) %% p
    idx <- seq_len(dg)
    f[idx] <- (f[idx] - q * g) %% p
    f <- .gf_trim_z(f)
  }
  f
}

.gf_mul_z <- function(f, g, p) {
  if (!length(f) || !length(g)) return(.z0())
  out <- rep(gmp::as.bigz(0), length(f) + length(g) - 1L)
  for (i in seq_along(f)) {
    idx <- seq_along(g) + i - 1L
    out[idx] <- (out[idx] + f[i] * g) %% p
  }
  .gf_trim_z(out)
}

.gf_gcd_z <- function(f, g, p) {
  f <- .gf_trim_z(f %% p)
  g <- .gf_trim_z(g %% p)
  while (length(g)) {
    r <- .gf_rem_z(f, g, p)
    f <- g
    g <- r
  }
  if (!length(f)) return(f)
  (f * .zinv(f[1L], p)) %% p
}

.gf_div_z <- function(f, g, p) {
  f <- .gf_trim_z(f %% p)
  dq <- length(f) - length(g)
  if (dq < 0) return(.z0())
  inv <- .zinv(g[1L], p)
  q <- rep(gmp::as.bigz(0), dq + 1L)
  for (i in seq_len(dq + 1L)) {
    c0 <- (f[i] * inv) %% p
    q[i] <- c0
    idx <- seq_along(g) + i - 1L
    f[idx] <- (f[idx] - c0 * g) %% p
  }
  .gf_trim_z(q)
}

.gf_pow_mod_z <- function(g, e, f, p) {
  r <- gmp::as.bigz(1)
  bits <- strsplit(as.character(gmp::as.bigz(e), b = 2L), "", fixed = TRUE)[[1L]]
  for (b in bits) {
    r <- .gf_rem_z(.gf_mul_z(r, r, p), f, p)
    if (b == "1") r <- .gf_rem_z(.gf_mul_z(r, g, p), f, p)
  }
  r
}

.gf_sub1 <- function(g, p) {
  if (!length(g)) return(gmp::as.bigz(p - 1))
  g[length(g)] <- (g[length(g)] - 1) %% p
  .gf_trim_z(g)
}

# Roots of f modulo prime p with f(0) != 0 and degree below p (Las Vegas).
.csolve_las_vegas <- function(f, p, seed = NULL) {
  randint <- .nt_rng(seed)
  root <- list()
  g <- .gf_sub1(.gf_pow_mod_z(gmp::as.bigz(c(1, 0)), p - 1, f, p), p)
  factors <- list(.gf_gcd_z(f, g, p))
  while (length(factors)) {
    f <- factors[[length(factors)]]
    factors[[length(factors)]] <- NULL
    lf <- length(f)
    if (lf <= 1) next
    if (lf == 2) {
      root[[length(root) + 1L]] <- (-.zinv(f[1L], p) * f[2L]) %% p
      next
    }
    if (lf == 3) {
      inv <- .zinv(f[1L], p)
      b <- (f[2L] * inv) %% p
      b <- (b + p * (b %% 2)) %/% 2
      r <- .sqrt_mod_all(b^2L - f[3L] * inv, p)
      if (length(r)) root <- c(root, lapply(seq_along(r), function(i) (r[i] - b) %% p))
      next
    }
    repeat {
      a <- randint(0, p - 1)
      g <- .gf_sub1(.gf_pow_mod_z(.zc(list(gmp::as.bigz(1), a)), (p - 1) %/% 2, f, p), p)
      g <- .gf_gcd_z(f, g, p)
      if (length(g) > 1 && length(g) < length(f)) {
        factors[[length(factors) + 1L]] <- g
        factors[[length(factors) + 1L]] <- .gf_div_z(f, g, p)
        break
      }
    }
  }
  r <- .zc(root)
  if (!length(r)) return(r)
  .zsort(r[!duplicated(as.character(r))])
}

.gf_diff_z <- function(f, p) {
  n <- length(f)
  if (n <= 1) return(.z0())
  .gf_trim_z((f[-n] * gmp::as.bigz((n - 1):1)) %% p)
}

.csolve_prime <- function(f, p, e = 1) {
  g <- f
  if (length(g) > p) for (i in seq_len(length(g) - as.numeric(p))) {
    g[i + as.numeric(p) - 1] <- g[i + as.numeric(p) - 1] + g[i]
    g[i] <- 0
  }
  g <- .gf_trim_z(g %% p)
  k <- 0
  while (k < length(g) && g[length(g) - k] == 0) k <- k + 1
  root_zero <- .z0()
  if (k) {
    g <- g[seq_len(length(g) - k)]
    root_zero <- gmp::as.bigz(0)
  }
  if (!length(g)) {
    X1 <- gmp::as.bigz(seq(0, as.numeric(p) - 1))
  } else if (length(g)^2 < p) {
    X1 <- c(root_zero, .csolve_las_vegas(g, p))
  } else {
    xs <- seq(0, as.numeric(p) - 1)
    keep <- vapply(xs, function(x) .gf_eval(g, gmp::as.bigz(x), p) == 0, TRUE)
    X1 <- c(root_zero, gmp::as.bigz(xs[keep]))
  }
  if (!length(X1)) return(.z0())
  if (e == 1) return(X1)
  fd <- .gf_diff_z(f %% p, p)
  X <- list()
  S <- lapply(seq_along(X1), function(i) list(X1[i], 1))
  while (length(S)) {
    top <- S[[length(S)]]
    S[[length(S)]] <- NULL
    x <- top[[1L]]
    s <- top[[2L]]
    if (s == e) {
      X[[length(X) + 1L]] <- x
      next
    }
    ps <- p^s
    alpha <- .gf_eval(fd, x, ps * p)
    fx <- gmp::as.bigz(0)
    for (i in seq_along(f)) fx <- fx * x + f[i]
    beta <- -fx %/% ps
    vs <- .linear_congruence(alpha, beta, p)
    for (i in seq_along(vs)) S[[length(S) + 1L]] <- list(x + vs[i] * ps, s + 1)
  }
  .zsort(.zc(X))
}

.gf_csolve <- function(f, n) {
  P <- .factorint(n)
  X <- lapply(names(P), function(pk) .csolve_prime(f, gmp::as.bigz(pk), P[[pk]]))
  if (any(lengths(X) == 0)) return(.z0())
  mods <- .zc(lapply(names(P), function(pk) gmp::as.bigz(pk)^P[[pk]]))
  combos <- list(.z0())
  for (pool in X) {
    nc <- list()
    for (cb in combos) for (i in seq_along(pool)) nc[[length(nc) + 1L]] <- c(cb, pool[i])
    combos <- nc
  }
  .zsort(.zc(lapply(combos, function(cb) .crt_raw(mods, cb))))
}

#' @rdname quadratic_congruence
#' @export
polynomial_congruence <- function(expr, m) {
  expr <- as_sym(expr)
  m <- .nt_int(m)
  fs <- free_symbols(expr)
  if (length(fs) != 1L) stop("The expression should be univariate", call. = FALSE)
  if (!is_polynomial(expr, fs[[1L]])) stop("The expression should be a polynomial", call. = FALSE)
  cf <- all_coeffs(expr, fs[[1L]])
  if (!all(vapply(cf, function(z) z$type == "Integer", TRUE)))
    stop("The expression should should have integer coefficients", call. = FALSE)
  co <- .zc(lapply(cf, .nt_int)) %% m
  rank <- length(co)
  if (rank == 3) return(.nt_out(.quadratic_congruence(co[1L], co[2L], co[3L], m)))
  if (rank == 2) return(.nt_out(.quadratic_congruence(gmp::as.bigz(0), co[1L], co[2L], m)))
  if (co[1L] == 1 && 1 + co[rank] == sum(co))
    return(nthroot_mod(-co[rank], rank - 1, m, all_roots = TRUE))
  .nt_out(.gf_csolve(co, m))
}

.binomial_mod_pp <- function(n, m, p, q) {
  modulo <- p^q
  fcache <- list(gmp::as.bigz(1))
  factorial_mod <- function(v) {
    v <- as.numeric(v)
    while (length(fcache) <= v) {
      k <- length(fcache)
      fcache[[k + 1L]] <<- (k * fcache[[k]]) %% modulo
    }
    fcache[[v + 1L]]
  }
  bj <- function(u, j, r) {
    prod <- u
    for (i in seq_len(r)) if (i != j) prod <- prod * (u * u - i * i)
    for (i in seq_len(r)) if (i != j) prod <- prod %/% (j * j - i * i)
    prod %/% j
  }
  up_factorial <- function(u) {
    r <- q %/% 2
    fac <- gmp::as.bigz(1)
    prd <- gmp::as.bigz(1)
    if ((r == 1 && p == 2) || (2 * r + 1) == p || (2 * r + 1) == p * p) {
      if (q %% 2 == 1) r <- r + 1
      md <- p^(2 * r)
      div <- p^(2 * r - q)
    } else {
      md <- p^(2 * r + 1)
      div <- p^((2 * r + 1) - q)
    }
    for (j in seq_len(r)) {
      lo <- (j - 1) * p + 1
      hi <- j * p - 1
      if (hi >= lo) for (mul in seq(as.numeric(lo), as.numeric(hi))) fac <- (fac * mul) %% md
      e <- bj(u, j, r)
      prd <- (prd * .powm(fac, e, md)) %% md
    }
    if (p == 2) {
      sm <- u %/% 2
      for (j in seq_len(r)) sm <- sm + (j %/% 2) * bj(u, j, r)
      if (sm %% 2 == 1) prd <- -prd
    }
    prd <- prd %% (md %/% div)
    prd %% md
  }
  up_plus_v_binom <- function(u, v) {
    prd <- gmp::as.bigz(1)
    div <- .zinv(factorial_mod(v), modulo)
    for (j in seq_len(q - 1)) {
      b <- div
      if (v >= 1) for (v_ in seq(as.numeric(j * p + 1), as.numeric(j * p + v))) b <- (b * v_) %% modulo
      aj <- u
      for (i in seq_len(q - 1)) if (i != j) aj <- aj * (u - i)
      for (i in seq_len(q - 1)) if (i != j) aj <- aj %/% (j - i)
      aj <- aj %/% j
      prd <- (prd * .powm(b, aj, modulo)) %% modulo
    }
    prd
  }
  factorial_p <- function(N) {
    u <- N %/% p
    v <- N %% p
    (factorial_mod(v) * up_factorial(u) * up_plus_v_binom(u, v)) %% modulo
  }
  prd <- gmp::as.bigz(1)
  Nj <- n
  Mj <- m
  Rj <- n - m
  e0 <- 0
  carry <- 0
  eq_1 <- 0
  j <- 0
  while (Nj != 0) {
    numerator <- factorial_p(Nj %% modulo)
    denominator <- (factorial_p(Mj %% modulo) * factorial_p(Rj %% modulo)) %% modulo
    Nj <- Nj %/% p
    mj <- Mj %% p
    Mj <- Mj %/% p
    rj <- Rj %% p
    Rj <- Rj %/% p
    carry <- as.numeric((mj + rj + carry) %/% p)
    e0 <- e0 + carry
    if (j >= q - 1) eq_1 <- eq_1 + carry
    prd <- (prd * numerator * .zinv(denominator, modulo)) %% modulo
    j <- j + 1
  }
  mul <- gmp::powm(gmp::as.bigz(if (p == 2 && q >= 3) 1 else -1) %% modulo, eq_1, modulo)
  (gmp::powm(p, e0, modulo) * mul * prd) %% modulo
}

#' @rdname quadratic_congruence
#' @export
binomial_mod <- function(n, m, k) {
  n <- .nt_int(n)
  m <- .nt_int(m)
  k <- .nt_int(k)
  if (k < 1) stop("k is required to be positive", call. = FALSE)
  if (n < 0 || m < 0 || m > n) return(0)
  f <- .factorint(k)
  if (!length(f)) return(0)
  ps <- .fac_primes(f)
  res <- .zc(lapply(seq_along(ps), function(i) .binomial_mod_pp(n, m, ps[i], f[[i]])))
  .nt_out(.crt_raw(.zc(lapply(seq_along(ps), function(i) ps[i]^f[[i]])), res))
}
# ---------------------------------------------------------------------------
# Multinomial coefficients (sympy/ntheory/multinomial.py)

# Named vector of exact doubles, or a named list of big integers when a value
# reaches 2^53 (bigz vectors do not keep names).
.nt_named <- function(keys, vals) {
  if (!length(vals)) return(stats::setNames(numeric(), character()))
  if (length(keys) > 1L && all(nzchar(keys))) {
    m <- do.call(rbind, lapply(strsplit(keys, ",", fixed = TRUE), as.numeric))
    o <- do.call(order, lapply(seq_len(ncol(m)), function(j) m[, j]))
    keys <- keys[o]
    vals <- vals[o]
  }
  if (all(abs(vals) < .NT_Z53)) return(stats::setNames(as.numeric(vals), keys))
  stats::setNames(lapply(seq_along(vals), function(i) vals[i]), keys)
}

.tkey <- function(t) paste(format(t, scientific = FALSE, trim = TRUE), collapse = ",")

#' Binomial and multinomial coefficients
#'
#' `binomial_coefficients(n)` returns the coefficients of `(x + y)^n` named by
#' the exponent pair, as in `"1,2"`, and `binomial_coefficients_list(n)` the
#' row of Pascal's triangle. `multinomial_coefficients(m, n)` returns the
#' coefficients of `(x1 + ... + xm)^n` named by exponent tuples.
#' `multinomial_coefficients_iterator(m, n)` gives the same coefficients in the
#' order of SymPy's iterator, which is cheaper when `m` is much larger than
#' `n`. Values are doubles while exact, otherwise a named list of big integers.
#'
#' @param n Power.
#' @param m Number of terms.
#' @return A named vector (or list), or a vector for `binomial_coefficients_list()`.
#' @examples
#' binomial_coefficients(4)
#' binomial_coefficients_list(5)
#' multinomial_coefficients(3, 2)
#' @export
binomial_coefficients <- function(n) {
  n <- as.numeric(.nt_int(n))
  keys <- c(.tkey(c(0, n)), .tkey(c(n, 0)))
  vals <- list(gmp::as.bigz(1), gmp::as.bigz(1))
  a <- gmp::as.bigz(1)
  for (k in seq_len(n %/% 2)) {
    a <- (a * (n - k + 1)) %/% k
    keys <- c(keys, .tkey(c(k, n - k)), .tkey(c(n - k, k)))
    vals <- c(vals, list(a, a))
  }
  keep <- !duplicated(keys)
  .nt_named(keys[keep], .zc(vals[keep]))
}

#' @rdname binomial_coefficients
#' @export
binomial_coefficients_list <- function(n) {
  n <- as.numeric(.nt_int(n))
  d <- rep(list(gmp::as.bigz(1)), n + 1)
  a <- gmp::as.bigz(1)
  for (k in seq_len(n %/% 2)) {
    a <- (a * (n - k + 1)) %/% k
    d[[k + 1]] <- a
    d[[n - k + 1]] <- a
  }
  .nt_out(.zc(d))
}

#' @rdname binomial_coefficients
#' @export
multinomial_coefficients <- function(m, n) {
  m <- as.numeric(.nt_int(m))
  n <- .nt_int(n)
  if (m == 1 && n >= 0) return(stats::setNames(1, as.character(n)))
  n <- as.numeric(n)
  if (!m) {
    if (n) return(stats::setNames(numeric(), character()))
    return(stats::setNames(1, ""))
  }
  if (m == 2) return(binomial_coefficients(n))
  if (m >= 2 * n && n > 1) return(multinomial_coefficients_iterator(m, n))
  t <- c(n, rep(0, m - 1))
  r <- new.env(hash = TRUE, parent = emptyenv())
  keys <- .tkey(t)
  assign(keys, gmp::as.bigz(1), envir = r)
  j <- if (n) 0 else m
  while (j < m - 1) {
    tj <- t[j + 1]
    if (j) {
      t[j + 1] <- 0
      t[1] <- tj
    }
    if (tj > 1) {
      t[j + 2] <- t[j + 2] + 1
      j <- 0
      start <- 1
      v <- gmp::as.bigz(0)
    } else {
      j <- j + 1
      start <- j + 1
      v <- get(.tkey(t), envir = r)
      t[j + 1] <- t[j + 1] + 1
    }
    for (k in seq_len(m - start) + start - 1) {
      if (t[k + 1]) {
        t[k + 1] <- t[k + 1] - 1
        v <- v + get(.tkey(t), envir = r)
        t[k + 1] <- t[k + 1] + 1
      }
    }
    t[1] <- t[1] - 1
    key <- .tkey(t)
    assign(key, (v * tj) %/% (n - t[1]), envir = r)
    keys <- c(keys, key)
  }
  .nt_named(keys, .zc(lapply(keys, get, envir = r)))
}

#' @rdname binomial_coefficients
#' @export
multinomial_coefficients_iterator <- function(m, n) {
  m <- as.numeric(.nt_int(m))
  n <- as.numeric(.nt_int(n))
  if (m < 2 * n || n == 1) return(multinomial_coefficients(m, n))
  if (n == 0) return(.nt_named(.tkey(rep(0, m)), gmp::as.bigz(1)))
  mc <- multinomial_coefficients(n, n)
  small <- list()
  for (k in names(mc)) {
    nz <- as.numeric(strsplit(k, ",", fixed = TRUE)[[1L]])
    small[[.tkey(nz[nz != 0])]] <- gmp::as.bigz(mc[[k]])
  }
  t <- c(n, rep(0, m - 1))
  keys <- .tkey(t)
  vals <- list(small[[.tkey(t[t != 0])]])
  j <- if (n) 0 else m
  while (j < m - 1) {
    tj <- t[j + 1]
    if (j) {
      t[j + 1] <- 0
      t[1] <- tj
    }
    if (tj > 1) {
      t[j + 2] <- t[j + 2] + 1
      j <- 0
    } else {
      j <- j + 1
      t[j + 1] <- t[j + 1] + 1
    }
    t[1] <- t[1] - 1
    keys <- c(keys, .tkey(t))
    vals[[length(vals) + 1L]] <- small[[.tkey(t[t != 0])]]
  }
  .nt_named(keys, .zc(vals))
}

# ---------------------------------------------------------------------------
# Partitions

.part_cache <- new.env(parent = emptyenv())
.part_cache$p <- list(gmp::as.bigz(1), gmp::as.bigz(1))

# p(n) by Euler's pentagonal recurrence, cached across calls; also backs
# partition().
# ponytail: O(n^1.5) big additions; SymPy's Rademacher series needs
# arbitrary precision floats, which base R lacks. Fine to about n = 10^5.
.npartitions <- function(n) {
  if (n < 0) return(gmp::as.bigz(0))
  p <- .part_cache$p
  have <- length(p) - 1
  if (n > have) {
    p <- c(p, vector("list", n - have))
    for (m in (have + 1):n) {
      k <- seq_len(floor((1 + sqrt(1 + 24 * m)) / 6))
      g1 <- k * (3 * k - 1) / 2
      g2 <- k * (3 * k + 1) / 2
      sg <- ifelse(k %% 2 == 1, 1, -1)
      i1 <- m - g1
      i2 <- m - g2
      pos <- c(i1[sg > 0 & i1 >= 0], i2[sg > 0 & i2 >= 0])
      neg <- c(i1[sg < 0 & i1 >= 0], i2[sg < 0 & i2 >= 0])
      s <- sum(.zc(p[pos + 1]))
      if (length(neg)) s <- s - sum(.zc(p[neg + 1]))
      p[[m + 1]] <- s
    }
    .part_cache$p <- p
  }
  p[[n + 1]]
}

#' Number of partitions
#'
#' `npartitions(n)` is the number of ways of writing `n` as a sum of positive
#' integers, as [partition()] gives symbolically.
#'
#' @param n A nonnegative integer.
#' @return An integer (a `gmp` big integer from `n = 300` on).
#' @examples
#' npartitions(25)
#' npartitions(1000)
#' @export
npartitions <- function(n) .nt_out(.npartitions(as.numeric(.nt_int(n))))

# ---------------------------------------------------------------------------
# Continued fractions (sympy/ntheory/continued_fraction.py)

# floor((p + sgn * sqrt(D)) / q) for nonsquare D, with t = isqrt(D).
.floor_quad <- function(p, sgn, t, q) {
  if (q < 0) {
    p <- -p
    sgn <- -sgn
    q <- -q
  }
  if (sgn > 0) (p + t) %/% q else (p - t - 1) %/% q
}

.cf_rational <- function(q) {
  out <- list()
  repeat {
    i <- gmp::numerator(q) %/% gmp::denominator(q)
    out[[length(out) + 1L]] <- i
    q <- q - i
    if (q == 0) break
    q <- 1 / q
  }
  .nt_out(.zc(out))
}

#' Continued fractions
#'
#' `continued_fraction(a)` returns the terms of a rational or quadratic
#' irrational number. A rational gives a vector of terms; a quadratic
#' irrational gives a list whose last element is the repeating block, so
#' `continued_fraction(sqrt(S(7)))` is `list(2, c(1, 1, 1, 4))`. Wherever
#' terms are accepted, a list means that its last element repeats.
#' `continued_fraction_periodic(p, q, d, s)` does the same for
#' `(p + s * sqrt(d)) / q`. `continued_fraction_reduce(cf)` turns terms back
#' into the number. `continued_fraction_iterator(x, n)` returns up to `n` terms
#' of any real `x`, exactly for rationals and in double precision otherwise
#' (reliable for roughly 10 terms). `continued_fraction_convergents(cf, n)`
#' returns the convergents of the terms, the first `n` of them when the
#' fraction repeats.
#'
#' @param a A `sym` rational or quadratic irrational, or something [as_sym()]
#'   accepts.
#' @param p,q,d,s Integers with `d >= 0` and `q != 0`.
#' @param cf Terms, as returned above.
#' @param x A real number.
#' @param n Maximum number of terms or convergents.
#' @return A vector or list of terms, or a `sym` number, or a list of `sym`
#'   convergents.
#' @examples
#' continued_fraction(Rational(415, 93))
#' continued_fraction(sqrt(S(7)))
#' continued_fraction_periodic(3, 2, 7)
#' continued_fraction_reduce(list(1, 2))
#' continued_fraction_reduce(list(1, c(2)))
#' continued_fraction_iterator(sym_pi, 6)
#' continued_fraction_convergents(c(0, 2, 1, 2))
#' @export
continued_fraction <- function(a) {
  e <- as_sym(a)
  bad <- function() stop("expecting a rational or quadratic irrational, not ", sstr(e), call. = FALSE)
  if (!all(vapply(atoms(e), is_exact_number, TRUE))) bad()
  if (e$type == "Integer") return(continued_fraction_periodic(e, 1, 0))
  if (e$type == "Rational") return(continued_fraction_periodic(num_numer(e), num_denom(e), 0))
  is_sqrt_int <- function(z) z$type == "Pow" && z$args[[1L]]$type == "Integer" &&
    z$args[[2L]]$type == "Rational" && z$args[[2L]]$value == gmp::as.bigq(1, 2)
  if (is_sqrt_int(e)) return(continued_fraction_periodic(0, 1, e$args[[1L]]))
  if (e$type == "Mul" && length(e$args) == 2L && is_exact_number(e$args[[1L]]) && is_sqrt_int(e$args[[2L]])) {
    co <- e$args[[1L]]
    return(continued_fraction_periodic(0, num_denom(co), e$args[[2L]]$args[[1L]], num_numer(co)))
  }
  # (a + b*sqrt(c))/d with integer a, b, d: clear the rational coefficients
  e <- expand(e)
  den <- gmp::as.bigz(1)
  for (z in atoms(e)) if (z$type == "Rational") den <- den * gmp::denominator(z$value) %/% .gcd2(den, gmp::denominator(z$value))
  d <- Integer(den)
  p <- expand(.mul(list(d, e)))
  if (d$type == "Integer") {
    if (is_exact_number(p)) return(continued_fraction_periodic(p, d))
    if (p$type == "Add" && length(p$args) == 2L) {
      a0 <- p$args[[1L]]
      bc <- p$args[[2L]]
    } else {
      a0 <- S_Zero
      bc <- p
    }
    if (a0$type == "Integer") {
      b <- NULL
      if (bc$type == "Mul" && length(bc$args) == 2L) {
        b <- bc$args[[1L]]
        c <- bc$args[[2L]]
      } else if (bc$type == "Pow") {
        b <- S_One
        c <- bc
      }
      if (!is.null(b) && b$type == "Integer" && is_sqrt_int(c))
        return(continued_fraction_periodic(a0, d, c$args[[1L]], b))
    }
  }
  bad()
}

#' @rdname continued_fraction
#' @export
continued_fraction_periodic <- function(p, q, d = 0, s = 1) {
  p <- .nt_int(p)
  q <- .nt_int(q)
  d <- .nt_int(d)
  s <- .nt_int(s)
  if (d < 0) stop("expected non-negative for `d` but got ", d, call. = FALSE)
  if (q == 0) stop("The denominator cannot be 0.", call. = FALSE)
  if (s == 0) d <- gmp::as.bigz(0)
  sq <- .iroot(d, 2L)
  if (sq$exact) return(.cf_rational(gmp::as.bigq(p + s * sq$root, q)))
  if (q < 0) {
    p <- -p
    q <- -q
    s <- -s
  }
  D <- s^2L * d
  sgn <- if (s > 0) 1 else -1
  fl <- .floor_quad(p, sgn, .iroot(D, 2L)$root, q)
  if (fl < 0) {
    # n = (p + s*sqrt(d))/q < 0: expand 1 - frac(-n), then shift the first term
    w <- -fl - 1
    one_f <- continued_fraction_periodic(q * (w + 1) + p, q, d, s)
    one_f[[1L]] <- one_f[[1L]] - .nt_out(w + 1)
    return(one_f)
  }
  if ((D - p^2L) %% q != 0) {
    D <- D * q^2L
    p <- p * q
    q <- q * q
  }
  tD <- .iroot(D, 2L)$root
  terms <- list()
  seen <- new.env(hash = TRUE, parent = emptyenv())
  repeat {
    key <- paste(as.character(p), as.character(q))
    if (!is.null(seen[[key]])) break
    assign(key, length(terms), envir = seen)
    a <- .floor_quad(p, sgn, tD, q)
    terms[[length(terms) + 1L]] <- a
    p <- a * q - p
    q <- (D - p^2L) %/% q
  }
  i <- seen[[paste(as.character(p), as.character(q))]]
  pre <- terms[seq_len(i)]
  per <- .nt_out(.zc(terms[(i + 1):length(terms)]))
  c(lapply(pre, .nt_out), list(per))
}

# A partial quotient as an exact bigq: integers, rationals and exact doubles.
.cf_term <- function(z) {
  if (is_sym(z) && is_exact_number(z)) return(z$value)
  if (is.double(z) && length(z) == 1L && is.finite(z) && z != trunc(z)) return(gmp::as.bigq(z))
  q <- .nt_rat(z)
  if (is.null(q)) gmp::as.bigq(.nt_int(z)) else q
}

#' @rdname continued_fraction
#' @export
continued_fraction_reduce <- function(cf) {
  if (!length(cf)) return(S_Zero)
  is_period <- is.list(cf) && length(cf) > 0
  last <- if (is_period) unlist(cf[[length(cf)]]) else NULL
  conv <- function(terms) {
    p2 <- S_Zero; q2 <- S_One; p1 <- S_One; q1 <- S_Zero
    for (a in terms) {
      a <- as_sym(a)
      p <- .add(list(.mul(list(a, p1)), p2))
      q <- .add(list(.mul(list(a, q1)), q2))
      p2 <- p1; q2 <- q1; p1 <- p; q1 <- q
    }
    list(p1, q1, p2, q2)
  }
  if (!is_period) {
    r <- conv(lapply(as.list(cf), function(z) make_rational(.cf_term(z))))
    return(.div(r[[1L]], r[[2L]]))
  }
  pre <- lapply(cf[-length(cf)], function(z) make_rational(gmp::as.bigq(.nt_int(z))))
  per <- .nt_ints(last)
  per <- lapply(seq_along(per), function(i) make_rational(gmp::as.bigq(per[i])))
  # purely periodic y = [per..., y] solves Q y^2 + (Q1 - P) y - P1 = 0
  r <- conv(per)
  P <- r[[1L]]$value; Q <- r[[2L]]$value; P1 <- r[[3L]]$value; Q1 <- r[[4L]]$value
  A <- Q
  B <- Q1 - P
  C <- -P1
  disc <- B^2 - 4 * A * C
  y <- .add(list(make_rational(-B / (2 * A)),
                 .mul(list(make_rational(1 / (2 * A)), .sqrt(make_rational(disc))))))
  if (!length(pre)) {
    rv <- y
  } else {
    r2 <- conv(pre)
    rv <- expand(radsimp(.div(.add(list(.mul(list(r2[[1L]], y)), r2[[3L]])),
                              .add(list(.mul(list(r2[[2L]], y)), r2[[4L]])))))
  }
  if (rv$type == "Add") rv <- factor_terms(rv)
  rv
}

#' @rdname continued_fraction
#' @export
continued_fraction_iterator <- function(x, n = NULL) {
  x <- as_sym(x)
  if (is_exact_number(x)) {
    r <- .cf_rational(x$value)
    return(if (is.null(n)) r else utils::head(r, n))
  }
  if (x$type == "Float") {
    # a double is an exact binary rational
    r <- .cf_rational(gmp::as.bigq(x$value))
    return(if (is.null(n)) r else utils::head(r, n))
  }
  if (is.null(n)) stop("n is needed for an irrational number", call. = FALSE)
  cf <- tryCatch(continued_fraction(x), error = function(e) NULL)
  if (!is.null(cf)) {
    pre <- lapply(cf[-length(cf)], gmp::as.bigz)
    per <- gmp::as.bigz(cf[[length(cf)]])
    need <- max(0, n - length(pre))
    terms <- .zc(c(pre, if (need) list(rep(per, length.out = need)) else list()))
    return(.nt_out(terms[seq_len(min(n, length(terms)))]))
  }
  # Other irrationals: run the expansion on exact rational bounds around the
  # double value and keep only terms both bounds agree on.
  b <- .nt_bounds(x)
  lo <- b[[1L]]; hi <- b[[2L]]
  out <- list()
  while (length(out) < n) {
    a <- .bigq_floor(lo)
    if (a != .bigq_floor(hi) || lo == a) break
    out[[length(out) + 1L]] <- a
    nlo <- 1 / (hi - a)
    hi <- 1 / (lo - a)
    lo <- nlo
  }
  if (length(out) < n)
    warning("only ", length(out), " terms are certain at double precision", call. = FALSE)
  .nt_out(.zc(lapply(out, gmp::as.bigz)))
}

# Rough amplification of the relative rounding error in evaluating x:
# cancellation in sums and large arguments of exp() and powers widen it.
.nt_cond <- function(x) {
  if (!length(x$args) || is_number_type(x)) return(1)
  k <- prod(vapply(x$args, .nt_cond, 0))
  mag <- function(z) abs(tryCatch(.neval(z, list()), error = function(e) NaN))
  if (x$type == "Add") {
    tot <- sum(vapply(x$args, mag, 0)); v <- mag(x)
    k <- k * if (is.finite(tot) && v > 0) max(1, tot / v) else Inf
  } else if (x$type == "Function" && length(x$args) == 1L && x$name %in% c("floor", "ceiling", "frac", "sign", "Heaviside")) {
    # discontinuous: exact only when the argument stays clear of a jump
    b <- tryCatch(.nt_bounds(x$args[[1L]]), error = function(e) NULL)
    if (is.null(b)) return(Inf)
    ok <- if (x$name %in% c("sign", "Heaviside")) b[[1L]] > 0 || b[[2L]] < 0
          else .bigq_floor(b[[1L]]) == .bigq_floor(b[[2L]]) && b[[1L]] != .bigq_floor(b[[1L]])
    if (!ok) return(Inf)
  } else if (x$type == "Function" && length(x$args) == 1L) {
    # condition number |a f'(a) / f(a)| of the function at its argument; the
    # linear estimate needs a small absolute error in the argument
    a <- x$args[[1L]]
    if (abs(tryCatch(.neval(a, list()), error = function(e) NaN)) * 64 * .Machine$double.eps * k > 1e-6) return(Inf)
    u <- Dummy("u")
    d <- tryCatch(.neval(subs(.diff1(fcall(x$name, list(u)), u), u, a), list()), error = function(e) NaN)
    v <- mag(x); av <- mag(a)
    cn <- if (is.finite(abs(d)) && v > 0) abs(av * d) / v else Inf
    k <- k * max(1, cn)
  } else if (x$type == "Pow") {
    k <- k * max(1, mag(x$args[[2L]]) * max(1, abs(log(mag(x$args[[1L]])))))
  }
  k
}

# Exact rational bounds around the double value of a real constant.
.nt_bounds <- function(x) {
  v <- .neval(x, list())
  if (is.complex(v)) v <- Re(v)
  eps <- 64 * .Machine$double.eps * max(1, abs(v)) * .nt_cond(as_sym(x))
  if (!is.finite(eps)) stop("cannot bound ", format(x), " at double precision", call. = FALSE)
  list(gmp::as.bigq(v - eps), gmp::as.bigq(v + eps))
}

.nt_certified_floor <- function(x) {
  b <- .nt_bounds(x)
  a <- .bigq_floor(b[[1L]])
  if (a != .bigq_floor(b[[2L]])) stop("cannot decide floor(", format(x), ") at double precision", call. = FALSE)
  gmp::as.bigz(a)
}

#' @rdname continued_fraction
#' @export
continued_fraction_convergents <- function(cf, n = NULL) {
  if (is.list(cf)) {
    last <- cf[[length(cf)]]
    pre <- lapply(cf[-length(cf)], .cf_term)
    per <- lapply(as.list(unlist(last)), .cf_term)
    if (is.null(n)) stop("n is needed for a periodic continued fraction", call. = FALSE)
    terms <- c(pre, rep_len(per, max(0, n - length(pre))))
  } else {
    terms <- lapply(as.list(cf), .cf_term)
  }
  if (!is.null(n)) terms <- utils::head(terms, n)
  p2 <- gmp::as.bigq(0); q2 <- gmp::as.bigq(1); p1 <- gmp::as.bigq(1); q1 <- gmp::as.bigq(0)
  out <- vector("list", length(terms))
  for (i in seq_along(terms)) {
    a <- terms[[i]]
    p <- a * p1 + p2
    q <- a * q1 + q2
    p2 <- p1; q2 <- q1; p1 <- p; q1 <- q
    out[[i]] <- make_rational(p / q)
  }
  out
}

# ---------------------------------------------------------------------------
# Digits (sympy/ntheory/digits.py)

.digits_raw <- function(n, b) {
  x <- abs(n)
  out <- list()
  while (x >= b) {
    out[[length(out) + 1L]] <- x %% b
    x <- x %/% b
  }
  out[[length(out) + 1L]] <- x
  rev(as.numeric(.zc(out)))
}

#' Digits of an integer
#'
#' `digits(n, b)` returns the base `b` digits of `n`, most significant first,
#' preceded by `b` (or `-b` for negative `n`) as in SymPy; `digits` pads with
#' zeros to that many digits. `count_digits(n, b)` counts each digit value and
#' `is_palindromic(n, b)` tests whether the digits read the same backwards.
#'
#' @param n An integer.
#' @param b Base, at least 2.
#' @param digits Minimum number of digits, or `NULL`.
#' @return A numeric vector, a named count vector, or a logical.
#' @examples
#' digits(35)
#' digits(27, b = 2)
#' digits(-3958, 27)
#' count_digits(1111339)
#' is_palindromic(88, 3)
#' @export
digits <- function(n, b = 10, digits = NULL) {
  b <- .nt_int(b)
  n <- .nt_int(n)
  if (b < 2) stop("b must be greater than 1", call. = FALSE)
  if (b >= .NT_Z53) stop("bases at or above 2^53 are not supported", call. = FALSE)
  y <- .digits_raw(n, b)
  if (!is.null(digits)) {
    if (length(y) > digits) stop("For ", n, ", at least ", length(y), " digits are needed.", call. = FALSE)
    y <- c(rep(0, digits - length(y)), y)
  }
  c(if (n < 0) -as.numeric(b) else as.numeric(b), y)
}

#' @rdname digits
#' @export
count_digits <- function(n, b = 10) {
  d <- digits(n, b)[-1L]
  tab <- table(factor(d, levels = sort(unique(d))))
  stats::setNames(as.numeric(tab), names(tab))
}

#' @rdname digits
#' @export
is_palindromic <- function(n, b = 10) {
  d <- digits(n, b)[-1L]
  identical(d, rev(d))
}

# ---------------------------------------------------------------------------
# Egyptian fractions (sympy/ntheory/egyptian_fraction.py)

.egypt_greedy <- function(x, y) {
  out <- list()
  while (x != 1) {
    a <- (-y) %% x
    b <- y * (y %/% x + 1)
    out[[length(out) + 1L]] <- y %/% x + 1
    g <- .gcd2(a, b)
    x <- a %/% g
    y <- b %/% g
  }
  c(.zc(out), y)
}

.egypt_dups <- function(l, fix) {
  repeat {
    l <- .zsort(l)
    s <- as.character(l)
    i <- which(duplicated(s))
    if (!length(i)) return(l)
    l <- fix(l, i[1L] - 1L)
  }
}

# x copies of y, SymPy's `[y] * x` (sympy/ntheory/egyptian_fraction.py:164 and
# 194): CPython raises OverflowError past sys.maxsize and MemoryError below it;
# rep() would turn x into NA.
.egypt_copies <- function(x, y) {
  if (x > gmp::as.bigz("9223372036854775807")) stop(.sympy_condition("OverflowError", "cannot fit 'int' into an index-sized integer"))
  if (x > .Machine$integer.max) stop(.sympy_condition("MemoryError", ""))
  rep(y, as.numeric(x))
}

.egypt_graham_jewett <- function(x, y) {
  l <- .egypt_copies(x, y)
  .egypt_dups(l, function(l, i) {
    v <- l[i]
    l[i + 1L] <- v + 1
    c(l, v * (v + 1))
  })
}

.egypt_takenouchi <- function(x, y) {
  if (x == 3) {
    if (y %% 2 == 0) return(.zc(list(y %/% 2, y)))
    i <- (y - 1) %/% 2
    j <- i + 1
    k <- j + i
    return(.zc(list(j, k, j * k)))
  }
  l <- .egypt_copies(x, y)
  .egypt_dups(l, function(l, i) {
    k <- l[i]
    if (k %% 2 == 0) {
      l[i] <- k %/% 2
      return(l[-(i + 1L)])
    }
    l[i] <- (k + 1) %/% 2
    l[i + 1L] <- k * (k + 1) %/% 2
    l
  })
}

.egypt_golomb <- function(x, y) {
  if (x == 1) return(y)
  xp <- .zinv(x, y)
  .zsort(c(.zc(list(xp * y)), .egypt_golomb((x * xp - 1) %/% y, xp)))
}

#' Egyptian fractions
#'
#' `egyptian_fraction(r)` writes a positive rational as a sum of distinct unit
#' fractions and returns their denominators. A rational above 1 first takes
#' the harmonic series `1 + 1/2 + 1/3 + ...` as far as it fits. The
#' `algorithm` for the rest is `"Greedy"`, `"Graham Jewett"`, `"Takenouchi"`
#' or `"Golomb"`.
#'
#' @param r A positive rational, or `c(p, q)`.
#' @param algorithm Algorithm name.
#' @return A vector of denominators.
#' @examples
#' egyptian_fraction(Rational(3, 7))
#' egyptian_fraction(c(3, 7), "Graham Jewett")
#' egyptian_fraction(Rational(7, 3))
#' @export
egyptian_fraction <- function(r, algorithm = "Greedy") {
  q <- if ((is.numeric(r) || inherits(r, "bigz")) && length(r) == 2L) {
    gmp::as.bigq(.nt_int(r[1L]), .nt_int(r[2L]))
  } else if (is_sym(r) && is_exact_number(r)) {
    r$value
  } else if ((is.numeric(r) || inherits(r, "bigz") || is.character(r)) && length(r) == 1L) {
    q1 <- .nt_rat(r)
    if (is.null(q1)) gmp::as.bigq(.nt_int(r)) else q1
  } else {
    stop("Value must be a Rational or tuple of ints", call. = FALSE)
  }
  if (q <= 0) stop(.sympy_condition("ValueError", "Value must be positive"))
  x <- gmp::numerator(q)
  y <- gmp::denominator(q)
  if (y == 1 && x == 2) return(c(1, 2, 3, 6))
  if (x == y + 1) return(.nt_out(.zc(list(gmp::as.bigz(1), y))))
  prefix <- list()
  acc <- gmp::as.bigq(0)
  d <- gmp::as.bigz(1)
  while (acc + 1 / gmp::as.bigq(d) <= q) {
    acc <- acc + 1 / gmp::as.bigq(d)
    prefix[[length(prefix) + 1L]] <- d
    d <- d + 1
  }
  rem <- q - acc
  if (rem == 0) return(.nt_out(.zc(prefix)))
  x <- gmp::numerator(rem)
  y <- gmp::denominator(rem)
  post <- switch(algorithm,
                 "Greedy" = .egypt_greedy(x, y),
                 "Graham Jewett" = .egypt_graham_jewett(x, y),
                 "Takenouchi" = .egypt_takenouchi(x, y),
                 "Golomb" = .egypt_golomb(x, y),
                 stop("Entered invalid algorithm", call. = FALSE))
  .nt_out(c(.zc(prefix), post))
}

# ---------------------------------------------------------------------------
# Hexadecimal digits of pi (sympy/ntheory/bbp_pi.py)

.bbp_dn <- function(n, prec) (.bitlen(gmp::as.bigz(n + 1 + prec)) - 1) %/% 4 + prec + 3

.bbp_series <- function(j, n, prec) {
  D <- .bbp_dn(n, prec)
  D4 <- 4 * D
  two_D4 <- gmp::as.bigz(2)^D4
  s <- gmp::as.bigz(0)
  if (n >= 0) {
    k <- 0:n
    d <- gmp::as.bigz(j + 8 * k)
    pw <- .zc(lapply(seq_along(k), function(i) gmp::powm(gmp::as.bigz(16), n - k[i], d[i])))
    s <- sum((pw * two_D4) %/% d)
  }
  t <- gmp::as.bigz(0)
  e <- D4 - 4
  d <- 8 * (n + 1) + j
  while (e >= 0) {
    dt <- gmp::as.bigz(2)^e %/% d
    if (dt == 0) break
    t <- t + dt
    e <- e - 4
    d <- d + 8
  }
  s + t
}

#' Hexadecimal digits of pi
#'
#' `pi_hex_digits(n, prec)` returns `prec` hexadecimal digits of pi starting
#' after the `n`-th, by the Bailey-Borwein-Plouffe formula, so
#' `pi_hex_digits(0)` begins with the leading 3.
#'
#' @param n Starting position, 0 or more.
#' @param prec Number of digits.
#' @return A string.
#' @examples
#' pi_hex_digits(0)
#' pi_hex_digits(1, 10)
#' @export
pi_hex_digits <- function(n, prec = 14) {
  n <- as.numeric(.nt_int(n))
  prec <- as.numeric(.nt_int(prec))
  if (n < 0) stop("n cannot be negative", call. = FALSE)
  if (prec < 0) stop("prec cannot be negative", call. = FALSE)
  if (prec == 0) return("")
  n <- n - 1
  D <- .bbp_dn(n, prec)
  x <- (4 * .bbp_series(1, n, prec) - 2 * .bbp_series(4, n, prec) -
          .bbp_series(5, n, prec) - .bbp_series(6, n, prec)) %% gmp::as.bigz(16)^D
  h <- as.character(x %/% gmp::as.bigz(16)^(D - prec), b = 16L)
  paste0(strrep("0", max(0, prec - nchar(h))), h)
}
# ---------------------------------------------------------------------------
# Self-initializing quadratic sieve (sympy/ntheory/qs.py)

.QS_SHIFT <- 2^10
.QS_LOG2 <- 710

# Modular inverses of x modulo primes p below 2^25, vectorized over doubles.
.inv_vec <- function(x, p) {
  r0 <- p
  r1 <- x %% p
  t0 <- numeric(length(p))
  t1 <- rep(1, length(p))
  while (any(r1 != 0)) {
    act <- r1 != 0
    q <- r0 %/% ifelse(act, r1, 1)
    nr <- ifelse(act, r0 - q * r1, r1)
    nt <- ifelse(act, t0 - q * t1, t1)
    r0 <- ifelse(act, r1, r0)
    t0 <- ifelse(act, t1, t0)
    r1 <- ifelse(act, nr, 0)
    t1 <- nt
  }
  t0 %% p
}

.qs_factor_base <- function(prime_bound, n) {
  ps <- .primerange(gmp::as.bigz(3), gmp::as.bigz(prime_bound))
  nm <- as.numeric(n %% gmp::as.bigz(ps))
  keep <- vapply(seq_along(ps), function(i) .jacobi(gmp::as.bigz(nm[i]), gmp::as.bigz(ps[i])) == 1L, TRUE)
  ps <- ps[keep]
  nm <- nm[keep]
  tm <- vapply(seq_along(ps), function(i)
    as.numeric(.sqrt_mod_prime_power(gmp::as.bigz(nm[i]), gmp::as.bigz(ps[i]), 1)[1L]), 0)
  prime <- c(2, ps)
  # 0-based index of the last factor base prime below 1000 (and 5000)
  idx_1000 <- if (any(ps > 1000)) which(ps > 1000)[1L] - 1 else NA
  idx_5000 <- if (any(ps > 5000)) which(ps > 5000)[1L] - 1 else NA
  list(prime = prime, tmem = c(1, tm), logp = c(.QS_LOG2, round(log(ps) * .QS_SHIFT)),
       idx_1000 = idx_1000, idx_5000 = idx_5000)
}

# Polynomials (a x + b)^2 - N, the a values chosen close to sqrt(2N)/M and the
# b values stepped through a Gray code. Returns a closure yielding the next
# polynomial; solutions of the sieve congruences live in fb$soln1/soln2
# (0-based residues, NA where the prime divides a).
.qs_polys <- function(N, M, fb, randint) {
  approx_val <- log(2 * N) / 2 - log(M)
  nfb <- length(fb$prime)
  start <- if (is.na(fb$idx_1000) || fb$idx_1000 == 0) 0 else fb$idx_1000
  end <- if (is.na(fb$idx_5000) || fb$idx_5000 == 0) nfb - 1 else fb$idx_5000
  st <- new.env(parent = emptyenv())
  st$i <- 0
  st$count <- 0
  st$fb <- fb
  pz <- gmp::as.bigz(fb$prime[-1L])
  new_a <- function() {
    best <- NULL
    for (rep in 1:50) {
      a <- gmp::as.bigz(1)
      q <- integer()
      while (log(a) < approx_val) {
        rp <- 0
        while (rp == 0 || rp %in% q) rp <- as.numeric(randint(start, end))
        a <- a * fb$prime[rp + 1]
        q <- c(q, rp)
      }
      ratio <- exp(log(a) - approx_val)
      if (is.null(best) || abs(ratio - 1) < abs(best$ratio - 1)) best <- list(a = a, q = q, ratio = ratio)
    }
    a <- best$a
    B <- list()
    for (val in best$q) {
      ql <- fb$prime[val + 1]
      g <- (fb$tmem[val + 1] * as.numeric(.zinv(a %/% ql, gmp::as.bigz(ql)))) %% ql
      if (2 * g > ql) g <- ql - g
      B[[length(B) + 1L]] <- a %/% ql * g
    }
    B <- .zc(B)
    b <- sum(B)
    p <- fb$prime[-1L]
    am <- as.numeric(a %% pz)
    div <- am == 0
    ainv <- .inv_vec(ifelse(div, 1, am), p)
    bm <- as.numeric(b %% pz)
    s1 <- ((ainv * ((fb$tmem[-1L] - bm) %% p)) %% p)
    s2 <- ((ainv * ((-fb$tmem[-1L] - bm) %% p)) %% p)
    s1[div] <- NA
    s2[div] <- NA
    st$bainv <- vapply(seq_along(B), function(j) ((2 * as.numeric(B[j] %% pz)) %% p * ainv) %% p,
                       numeric(length(p)))
    if (!is.matrix(st$bainv)) st$bainv <- matrix(st$bainv, ncol = length(B))
    st$s1 <- s1
    st$s2 <- s2
    st$a <- a
    st$b <- b
    st$B <- B
    st$i <- 0
  }
  function() {
    if (st$count == 0 || st$i >= 2^(length(st$B) - 1) - 1) {
      new_a()
    } else {
      st$i <- st$i + 1
      i <- st$i
      v <- .v2(gmp::as.bigz(i))
      neg <- 2 * ((i %/% 2^(v + 1)) %% 2) - 1
      st$b <- st$b + 2 * neg * st$B[v + 1]
      p <- fb$prime[-1L]
      st$s1 <- (st$s1 - neg * st$bainv[, v + 1]) %% p
      st$s2 <- (st$s2 - neg * st$bainv[, v + 1]) %% p
    }
    st$count <- st$count + 1
    a <- st$a
    b <- st$b
    list(a = a, b = b, a2 = a^2L, ab = 2 * a * b, b2 = b^2L - N, s1 = st$s1, s2 = st$s2)
  }
}

.qs_sieve <- function(M, fb, g) {
  if (g$a %% 2 == 0) {
    arr <- numeric(2 * M)
  } else if ((M - g$b) %% 2 == 1) {
    arr <- rep(c(.QS_LOG2, 0), M)
  } else {
    arr <- rep(c(0, .QS_LOG2), M)
  }
  p <- fb$prime[-1L]
  lp <- fb$logp[-1L]
  for (i in seq_along(p)) {
    if (is.na(g$s1[i])) next
    for (s in c(g$s1[i], g$s2[i])) {
      st <- (M + s) %% p[i]
      if (st < 2 * M) {
        idx <- seq(st, 2 * M - 1, by = p[i]) + 1
        arr[idx] <- arr[idx] + lp[i]
      }
    }
  }
  arr
}

.qs_smooth <- function(num, fbz, nfb) {
  vec <- logical(nfb + 1L)
  if (num < 0) {
    num <- -num
    vec[1L] <- TRUE
  }
  for (i in which(num %% fbz == 0)) {
    e <- 0
    while (num %% fbz[i] == 0) {
      num <- num %/% fbz[i]
      e <- e + 1
    }
    if (e %% 2 == 1) vec[i + 1L] <- TRUE
  }
  list(vec, num)
}

.qs_trial_stage <- function(N, M, fb, fbz, arr, g, partial, ERROR_TERM) {
  acc <- (log(M) + log(N) / 2 - ERROR_TERM) * .QS_SHIFT
  smooth <- list()
  proper <- list()
  bound <- 128 * fb$prime[length(fb$prime)]
  for (pos in which(arr >= acc)) {
    x <- pos - 1 - M
    v <- (g$a2 * x + g$ab) * x + g$b2
    sm <- .qs_smooth(v, fbz, length(fb$prime))
    vec <- sm[[1L]]
    num <- sm[[2L]]
    if (num == 1) {
      smooth[[length(smooth) + 1L]] <- list(g$a * x + g$b, v, vec)
    } else if (num < bound && .is_prime(num)) {
      if (N %% num == 0) {
        proper[[length(proper) + 1L]] <- num
        next
      }
      u <- g$a * x + g$b
      k <- as.character(num)
      prev <- partial[[k]]
      if (!is.null(prev)) {
        rm(list = k, envir = partial)
        smooth[[length(smooth) + 1L]] <- list((u * prev[[1L]] * .zinv(num, N)) %% N,
                                              (v * prev[[2L]]) %/% num^2L, xor(vec, prev[[3L]]))
      } else {
        assign(k, list(u, v, vec), envir = partial)
      }
    }
  }
  list(smooth, proper)
}

# Gaussian elimination over GF(2) on the parity vectors, then a square
# congruence from each dependent relation; returns the factors found, in order.
.qs_find_factor <- function(N, rels, col) {
  row <- length(rels)
  mat <- do.call(rbind, lapply(rels, `[[`, 3L))
  mark <- logical(row)
  for (pos in seq_len(col)) {
    hits <- which(mat[, pos])
    if (!length(hits)) next
    i <- hits[1L]
    add_col <- mat[i, ]
    add_col[pos] <- FALSE
    mat[i, ] <- FALSE
    mat[i, pos] <- TRUE
    mark[i] <- TRUE
    later <- hits[hits > i]
    if (length(later)) mat[later, ] <- xor(mat[later, , drop = FALSE], matrix(add_col, length(later), col, byrow = TRUE))
  }
  pivot <- rep(NA_integer_, row)
  pivot[mark] <- apply(mat[mark, , drop = FALSE], 1L, function(r) which(r)[1L])
  out <- list()
  for (i in which(!mark)) {
    u <- rels[[i]][[1L]]
    v <- rels[[i]][[2L]]
    for (k in which(mark)) {
      if (mat[i, pivot[k]]) {
        u <- u * rels[[k]][[1L]]
        v <- v * rels[[k]][[2L]]
      }
    }
    v <- .isqrt(v)
    g <- .gcd2(u - v, N)
    if (g > 1 && g < N) out[[length(out) + 1L]] <- g
  }
  out
}

#' Quadratic sieve
#'
#' `qs_factor(N, prime_bound, M)` factors `N` with the self-initializing
#' quadratic sieve, using the primes below `prime_bound` as the factor base and
#' a sieve interval of `2 * M`; it returns a named vector of exponents whose
#' names are the factors found (not always prime). `qs()` returns the factors
#' alone. Random choices come from a stream private to the call, seeded by
#' `seed`.
#'
#' @param N An integer greater than 1.
#' @param prime_bound Bound for the factor base.
#' @param M Half the sieve interval.
#' @param ERROR_TERM Slack in the sieve threshold.
#' @param seed Seed.
#' @return A named vector, or a vector of factors.
#' @examples
#' qs(10009202107, 100, 2000)
#' @export
qs_factor <- function(N, prime_bound, M, ERROR_TERM = 25, seed = 1234) {
  N <- .nt_int(N)
  prime_bound <- as.numeric(prime_bound)
  M <- as.numeric(M)
  if (N < 2) stop("N should be greater than 1", call. = FALSE)
  fac <- .fac_new()
  if (N %% 2 == 0) {
    r <- .remove(N, 2)
    N <- r[[1L]]
    .fac_set(fac, 2, r[[2L]])
    if (N == 1) return(.fac_vec(fac))
  }
  if (.is_prime(N)) {
    .fac_set(fac, N, 1)
    return(.fac_vec(fac))
  }
  pp <- .perfect_power_int(N, 3)
  if (!isFALSE(pp)) {
    .fac_set(fac, pp[[1L]], pp[[2L]])
    return(.fac_vec(fac))
  }
  N_copy <- N
  randint <- .nt_rng(seed)
  fb <- .qs_factor_base(prime_bound, N)
  fbz <- gmp::as.bigz(fb$prime)
  threshold <- (length(fb$prime) * 105) %/% 100
  partial <- new.env(hash = TRUE, parent = emptyenv())
  smooth <- list()
  nextpoly <- .qs_polys(N, M, fb, randint)
  strip <- function(p) {
    if (N_copy %% p != 0) return(FALSE)
    r <- .remove(N_copy, p)
    N_copy <<- r[[1L]]
    .fac_set(fac, p, r[[2L]])
    TRUE
  }
  repeat {
    g <- nextpoly()
    arr <- .qs_sieve(M, fb, g)
    tr <- .qs_trial_stage(N, M, fb, fbz, arr, g, partial, ERROR_TERM)
    smooth <- c(smooth, tr[[1L]])
    for (p in tr[[2L]]) strip(p)
    if (threshold <= length(smooth)) break
  }
  for (f in .qs_find_factor(N, smooth, length(fb$prime) + 1L)) {
    if (strip(f) && (N_copy == 1 || .is_prime(N_copy))) break
  }
  if (N_copy != 1) .fac_set(fac, N_copy, 1)
  .fac_vec(fac)
}

#' @rdname qs_factor
#' @export
qs <- function(N, prime_bound, M, ERROR_TERM = 25, seed = 1234) {
  .nt_out(gmp::as.bigz(names(qs_factor(N, prime_bound, M, ERROR_TERM, seed))))
}

# ---------------------------------------------------------------------------
# Elliptic curves (sympy/ntheory/elliptic_curve.py)

.ec_domain <- function(modulus) {
  if (modulus == 0) {
    from <- function(v) {
      if (is_sym(v)) return(gmp::as.bigq(.to_bigq(v)))
      gmp::as.bigq(.to_bigq(v))
    }
    return(list(from = from, add = function(a, b) a + b, sub = function(a, b) a - b,
                mul = function(a, b) a * b, div = function(a, b) a / b, neg = function(a) -a,
                zero = function(a) a == 0, sym = function(a) make_rational(a),
                int = function(a) gmp::denominator(a) == 1))
  }
  p <- gmp::as.bigz(modulus)
  if (!.is_prime(p)) stop("modulus must be a prime", call. = FALSE)
  from <- function(v) {
    q <- gmp::as.bigq(.to_bigq(v))
    d <- .zinv(gmp::denominator(q), p)
    if (is.na(d)) stop("denominator not invertible modulo ", p, call. = FALSE)
    (gmp::numerator(q) * d) %% p
  }
  list(from = from, add = function(a, b) (a + b) %% p, sub = function(a, b) (a - b) %% p,
       mul = function(a, b) (a * b) %% p,
       div = function(a, b) {
         i <- .zinv(b, p)
         if (is.na(i)) stop("division by zero modulo ", p, call. = FALSE)
         (a * i) %% p
       },
       neg = function(a) (-a) %% p, zero = function(a) a %% p == 0,
       sym = function(a) Integer(.symres(a %% p, p)), int = function(a) TRUE)
}

#' Elliptic curves
#'
#' `EllipticCurve(a4, a6, a1, a2, a3, modulus)` is the curve
#' `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6` over the rationals
#' (`modulus = 0`) or over the integers modulo a prime. The result holds the
#' fields `characteristic`, `discriminant`, `is_singular` and `j_invariant`, and
#' the methods `point(x, y, z = 1)` (also `EllipticCurvePoint()`),
#' `point_at_infinity()`, `points()` (a matrix of the affine points over a
#' finite field), `points_x(x)`, `torsion_points()`, `minimal()`, `order()` and
#' `rank()`. Points support `+`, `-`, `*` by an integer and `==`, and
#' [point_order()] gives the order of a point (`Inf` for infinite order).
#'
#' Several SymPy defects are fixed: `order()` counts the point at infinity, so
#' it is the order of the group; a point `P` with `P == -P` has order 2 even
#' when `a1` or `a3` is nonzero; [point_order()] over the rationals relies on
#' Mazur's bound of 12 rather than integral coordinates; `torsion_points()`
#' works through an integral short model, so it also finds the torsion of
#' curves with rational coefficients; and the discriminant of a curve with
#' rational coefficients is exact rather than truncated. Over a finite field,
#' `points()` lists the affine points as a numeric matrix, so the modulus must
#' be below `2^53`, and `points_x(x)` returns a list of `list(x, y)`.
#'
#' @param a4,a6,a1,a2,a3 Coefficients.
#' @param modulus 0, or a prime.
#' @param x,y,z Projective coordinates.
#' @param curve An `EllipticCurve`.
#' @param p A point.
#' @return An `EllipticCurve` environment, or an `EllipticCurvePoint`.
#' @examples
#' e <- EllipticCurve(-1, 9)
#' p <- e$point(0, 3)
#' q <- e$point(-1, 3)
#' p + q
#' 2 * p + q
#' e$discriminant
#' EllipticCurve(0, 1)$torsion_points()
#' f <- EllipticCurve(1, 0, modulus = 19)
#' f$order()
#' point_order(f$point(3, 7))
#' @export
EllipticCurve <- function(a4, a6, a1 = 0, a2 = 0, a3 = 0, modulus = 0) {
  modulus <- .nt_int(modulus)
  D <- .ec_domain(modulus)
  e <- new.env(parent = emptyenv())
  e$domain <- D
  e$modulus <- .nt_out(modulus)
  e$characteristic <- .nt_out(modulus)
  a <- lapply(list(a1, a2, a3, a4, a6), D$from)
  names(a) <- c("a1", "a2", "a3", "a4", "a6")
  for (nm in names(a)) assign(nm, a[[nm]], envir = e)
  m <- D$mul
  ad <- D$add
  sb <- D$sub
  k <- function(v) D$from(v)
  b2 <- ad(m(a$a1, a$a1), m(k(4), a$a2))
  b4 <- ad(m(k(2), a$a4), m(a$a1, a$a3))
  b6 <- ad(m(a$a3, a$a3), m(k(4), a$a6))
  b8 <- sb(ad(sb(ad(m(m(a$a1, a$a1), a$a6), m(m(k(4), a$a2), a$a6)), m(m(a$a1, a$a3), a$a4)),
              m(a$a2, m(a$a3, a$a3))), m(a$a4, a$a4))
  disc <- ad(sb(sb(D$neg(m(m(b2, b2), b8)), m(k(8), m(b4, m(b4, b4)))), m(k(27), m(b6, b6))),
             m(k(9), m(b2, m(b4, b6))))
  e$b2 <- b2
  e$b4 <- b4
  e$b6 <- b6
  e$b8 <- b8
  e$.disc <- disc
  e$discriminant <- if (modulus == 0) {
    if (gmp::denominator(disc) == 1) .nt_out(gmp::numerator(disc)) else make_rational(disc)
  } else .nt_out(disc)
  e$is_singular <- D$zero(disc)
  c4 <- sb(m(b2, b2), m(k(24), b4))
  e$j_invariant <- if (e$is_singular) NULL else D$sym(D$div(m(c4, m(c4, c4)), disc))
  e$point <- function(x, y, z = 1) EllipticCurvePoint(x, y, z, e)
  e$point_at_infinity <- function() .ec_point(D$from(0), D$from(1), D$from(0), e)
  e$minimal <- function() {
    if (modulus == 2) return(e)
    if (modulus == 3) {
      return(EllipticCurve(.ec_num(D$div(b4, k(2)), e), .ec_num(D$div(b6, k(4)), e),
                           a2 = .ec_num(D$div(b2, k(4)), e), modulus = modulus))
    }
    c6 <- D$neg(ad(ad(m(b2, m(b2, b2)), D$neg(m(k(36), m(b2, b4)))), m(k(216), b6)))
    EllipticCurve(.ec_num(D$neg(m(k(27), c4)), e), .ec_num(D$neg(m(k(54), c6)), e), modulus = modulus)
  }
  e$points <- function() {
    if (modulus == 0) stop("Infinitely many points", call. = FALSE)
    rows <- list()
    for (xi in seq(0, as.numeric(modulus) - 1)) {
      ys <- .ec_y_mod(e, gmp::as.bigz(xi))
      for (j in seq_along(ys)) rows[[length(rows) + 1L]] <- c(xi, as.numeric(ys[j]))
    }
    if (!length(rows)) return(matrix(numeric(), 0, 2, dimnames = list(NULL, c("x", "y"))))
    out <- do.call(rbind, rows)
    colnames(out) <- c("x", "y")
    out
  }
  e$points_x <- function(x) {
    if (modulus != 0) {
      xv <- D$from(x)
      ys <- .ec_y_mod(e, xv)
      return(lapply(seq_along(ys), function(i) list(x = .nt_out(xv), y = .nt_out(ys[i]))))
    }
    lapply(.ec_y_rational(e, as_sym(x)), function(yv) list(x = as_sym(x), y = yv))
  }
  e$torsion_points <- function() {
    if (modulus != 0) stop("No torsion point for Finite Field.", call. = FALSE)
    .ec_torsion(e)
  }
  e$order <- function() {
    if (modulus == 0) stop("Still not implemented", call. = FALSE)
    nrow(e$points()) + 1
  }
  e$rank <- function() {
    if (modulus != 0) return(0)
    stop("Still not implemented", call. = FALSE)
  }
  class(e) <- "EllipticCurve"
  e
}

.ec_num <- function(v, e) if (e$modulus == 0) make_rational(v) else v

# y with (x, y) on a curve over GF(p).
.ec_y_mod <- function(e, x) {
  p <- gmp::as.bigz(e$modulus)
  D <- e$domain
  lin <- D$add(D$mul(e$a1, x), e$a3)
  rhs <- D$add(D$add(D$add(D$mul(x, D$mul(x, x)), D$mul(e$a2, D$mul(x, x))), D$mul(e$a4, x)), e$a6)
  .quadratic_congruence(gmp::as.bigz(1), lin, -rhs, p)
}

.ec_sym <- function(v) make_rational(gmp::as.bigq(v))

# Torsion points over Q: scale to integral coefficients (x = X/u^2,
# y = Y/u^3), move to the short model Y^2 = X^3 - 27 c4 X - 54 c6 by
# (36x + 3 b2, 108 (2y + a1 x + a3)), apply Nagell-Lutz there (Y = 0, or Y^2
# dividing 4A^3 + 27B^2), keep the points of finite order and map them back.
# SymPy searches the given model directly, which misses points of curves with
# rational coefficients or with a1, a3 nonzero.
.ec_torsion <- function(e) {
  a <- list(e$a1, e$a2, e$a3, e$a4, e$a6)
  u <- gmp::as.bigz(1)
  for (v in a) {
    d <- gmp::denominator(v)
    u <- u * d %/% .gcd2(u, d)
  }
  uq <- gmp::as.bigq(u)
  sc <- Map(function(v, k) v * uq^k, a, c(1, 2, 3, 4, 6))
  s1 <- sc[[1L]]; s2 <- sc[[2L]]; s3 <- sc[[3L]]; s4 <- sc[[4L]]; s6 <- sc[[5L]]
  b2 <- s1^2 + 4 * s2
  b4 <- 2 * s4 + s1 * s3
  b6 <- s3^2 + 4 * s6
  c4 <- b2^2 - 24 * b4
  c6 <- -b2^3 + 36 * b2 * b4 - 216 * b6
  A <- -27 * c4
  B <- -54 * c6
  short <- EllipticCurve(make_rational(A), make_rational(B))
  X <- Dummy("X")
  xroots <- function(Y2) {
    f <- .add(list(.pow(X, S(3L)), .mul(list(make_rational(A), X)), make_rational(B - Y2)))
    Filter(is_exact_number, solve(f, X))
  }
  back <- function(Xr, Yr) {
    xs <- (Xr - 3 * b2) / 36
    ys <- (Yr / 108 - s1 * xs - s3) / 2
    e$point(make_rational(xs / uq^2), make_rational(ys / uq^3))
  }
  out <- list(e$point_at_infinity())
  for (xr in xroots(gmp::as.bigq(0))) out[[length(out) + 1L]] <- back(xr$value, gmp::as.bigq(0))
  D <- abs(gmp::numerator(4 * A^3 + 27 * B^2))
  dv <- .nt_ints(divisors(D))
  for (ii in seq_along(dv)) {
    r <- .iroot(dv[ii], 2L)
    if (!r$exact) next
    Y <- gmp::as.bigq(r$root)
    for (xr in xroots(Y^2)) {
      if (is.finite(point_order(short$point(xr, make_rational(Y))))) {
        pt <- back(xr$value, Y)
        out <- c(out, list(pt, -pt))
      }
    }
  }
  out
}

.ec_y_rational <- function(e, x0) {
  y <- Dummy("y")
  f <- .sub(.add(list(.pow(y, S_Two), .mul(list(.ec_sym(e$a1), x0, y)), .mul(list(.ec_sym(e$a3), y)))),
            .add(list(.pow(x0, S(3L)), .mul(list(.ec_sym(e$a2), .pow(x0, S_Two))),
                      .mul(list(.ec_sym(e$a4), x0)), .ec_sym(e$a6))))
  solve(f, y)
}

.ec_contains <- function(e, x, y, z) {
  D <- e$domain
  if (e$modulus == 0 && D$zero(z)) return(TRUE)
  m <- D$mul
  ad <- D$add
  lhs <- ad(ad(m(m(y, y), z), m(m(e$a1, x), m(y, z))), m(m(e$a3, y), m(z, z)))
  rhs <- ad(ad(ad(m(x, m(x, x)), m(m(e$a2, m(x, x)), z)), m(m(e$a4, x), m(z, z))), m(e$a6, m(z, m(z, z))))
  D$zero(D$sub(lhs, rhs))
}

.ec_point <- function(x, y, z, curve) structure(list(x = x, y = y, z = z, curve = curve),
                                                class = "EllipticCurvePoint")

#' @rdname EllipticCurve
#' @export
EllipticCurvePoint <- function(x, y, z, curve) {
  D <- curve$domain
  x <- D$from(x)
  y <- D$from(y)
  z <- D$from(z)
  if (!.ec_contains(curve, x, y, z)) stop("The curve does not contain this point", call. = FALSE)
  .ec_point(x, y, z, curve)
}

.ec_is_inf <- function(p) p$curve$domain$zero(p$z)

.ec_neg <- function(p) {
  e <- p$curve
  D <- e$domain
  if (.ec_is_inf(p)) return(e$point_at_infinity())
  .ec_point(p$x, D$sub(D$neg(p$y), D$add(D$mul(e$a1, p$x), D$mul(e$a3, p$z))), p$z, e)
}

.ec_affine <- function(p) {
  D <- p$curve$domain
  list(D$div(p$x, p$z), D$div(p$y, p$z))
}

.ec_eq <- function(p, q) {
  if (!identical(p$curve, q$curve)) return(FALSE)
  if (.ec_is_inf(p) || .ec_is_inf(q)) return(.ec_is_inf(p) && .ec_is_inf(q))
  a <- .ec_affine(p)
  b <- .ec_affine(q)
  D <- p$curve$domain
  D$zero(D$sub(a[[1L]], b[[1L]])) && D$zero(D$sub(a[[2L]], b[[2L]]))
}

.ec_addp <- function(p, q) {
  if (!identical(p$curve, q$curve)) stop("points lie on different curves", call. = FALSE)
  if (.ec_is_inf(p)) return(q)
  if (.ec_is_inf(q)) return(p)
  e <- p$curve
  D <- e$domain
  m <- D$mul
  ad <- D$add
  sb <- D$sub
  k <- D$from
  P1 <- .ec_affine(p)
  P2 <- .ec_affine(q)
  x1 <- P1[[1L]]; y1 <- P1[[2L]]; x2 <- P2[[1L]]; y2 <- P2[[2L]]
  if (!D$zero(sb(x1, x2))) {
    slope <- D$div(sb(y1, y2), sb(x1, x2))
    yint <- D$div(sb(m(y1, x2), m(y2, x1)), sb(x2, x1))
  } else {
    if (D$zero(ad(ad(ad(y1, y2), m(e$a1, x1)), e$a3))) return(e$point_at_infinity())
    den <- ad(ad(m(e$a1, x1), e$a3), m(k(2), y1))
    slope <- D$div(sb(ad(ad(m(k(3), m(x1, x1)), m(m(k(2), e$a2), x1)), e$a4), m(e$a1, y1)), den)
    yint <- D$div(sb(ad(ad(D$neg(m(x1, m(x1, x1))), m(e$a4, x1)), m(k(2), e$a6)), m(e$a3, y1)), den)
  }
  x3 <- sb(sb(sb(ad(m(slope, slope), m(e$a1, slope)), e$a2), x1), x2)
  y3 <- sb(sb(D$neg(m(ad(slope, e$a1), x3)), yint), e$a3)
  .ec_point(x3, y3, k(1), e)
}

.ec_mulp <- function(p, n) {
  n <- .nt_int(n)
  r <- p$curve$point_at_infinity()
  if (n == 0) return(r)
  if (n < 0) return(.ec_neg(.ec_mulp(p, -n)))
  bits <- rev(strsplit(as.character(n, b = 2L), "", fixed = TRUE)[[1L]])
  for (b in bits) {
    if (b == "1") r <- .ec_addp(r, p)
    p <- .ec_addp(p, p)
  }
  r
}

#' @export
Ops.EllipticCurvePoint <- function(e1, e2) {
  pt <- function(z) inherits(z, "EllipticCurvePoint")
  if (missing(e2)) {
    if (.Generic == "-") return(.ec_neg(e1))
    if (.Generic == "+") return(e1)
  }
  switch(.Generic,
         "+" = .ec_addp(e1, e2),
         "-" = .ec_addp(e1, .ec_neg(e2)),
         "*" = if (pt(e1)) .ec_mulp(e1, e2) else .ec_mulp(e2, e1),
         "==" = .ec_eq(e1, e2),
         "!=" = !.ec_eq(e1, e2),
         stop("operator ", .Generic, " is not defined for points", call. = FALSE))
}

#' @rdname EllipticCurve
#' @export
point_order <- function(p) {
  if (.ec_is_inf(p)) return(1)
  if (.ec_eq(p, .ec_neg(p))) return(2)
  q <- .ec_mulp(p, 2)
  if (.ec_eq(q, .ec_neg(p))) return(3)
  # over Q a point of finite order has order at most 12 (Mazur)
  i <- 2
  finite <- p$curve$modulus != 0
  repeat {
    q <- .ec_addp(p, q)
    i <- i + 1
    if (.ec_is_inf(q)) return(i)
    if (!finite && i >= 12) return(Inf)
  }
}

#' @export
chooseOpsMethod.EllipticCurvePoint <- function(x, y, mx, my, cl, reverse) TRUE

.ec_str <- function(v, e) {
  if (e$modulus == 0) return(sstr(make_rational(v)))
  sstr(e$domain$sym(v))
}

#' @export
format.EllipticCurvePoint <- function(x, ...) {
  if (.ec_is_inf(x)) return("O")
  a <- .ec_affine(x)
  paste0("(", .ec_str(a[[1L]], x$curve), ", ", .ec_str(a[[2L]], x$curve), ")")
}

#' @export
print.EllipticCurvePoint <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @export
print.EllipticCurve <- function(x, ...) {
  X <- Symbol("x")
  Y <- Symbol("y")
  cf <- function(v) if (x$modulus == 0) make_rational(v) else Integer(v)
  lhs <- .add(list(.pow(Y, S_Two), .mul(list(cf(x$a1), X, Y)), .mul(list(cf(x$a3), Y))))
  rhs <- .add(list(.pow(X, S(3L)), .mul(list(cf(x$a2), .pow(X, S_Two))), .mul(list(cf(x$a4), X)), cf(x$a6)))
  cat("EllipticCurve: ", sstr(lhs), " = ", sstr(rhs),
      if (x$modulus != 0) paste0(" over GF(", format(x$modulus), ")") else " over QQ", "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Number theoretic functions (sympy/functions/combinatorial/numbers.py)

.isF <- function(e, fact) isFALSE(ask_fact(e, fact))
.isT <- function(e, fact) isTRUE(ask_fact(e, fact))
.zsym <- function(z) make_rational(gmp::as.bigq(z))

.nt_check <- function(e, name, positive = FALSE, nonneg = FALSE) {
  # TypeError and ValueError, as sympy/functions/combinatorial/numbers.py:1764 to 1770
  if (.isF(e, "integer") || e$type %in% c("Rational", "Float"))
    stop(.sympy_condition("TypeError", paste0(name, " should be an integer")))
  if (positive && .isF(e, "positive")) stop(.sympy_condition("ValueError", paste0(name, " should be a positive integer")))
  if (nonneg && .isF(e, "nonnegative")) stop(.sympy_condition("ValueError", paste0(name, " should be a nonnegative integer")))
}

.nt_facts <- function(...) {
  extra <- c(...)
  function(args, afs) c(commutative = TRUE, extra)
}

.nt_latex1 <- function(sym) function(a) paste0(sym, "\\left(", a[1], "\\right)")

def_function("totient", eval = function(args) {
  n <- args[[1L]]
  .nt_check(n, "n", positive = TRUE)
  if (is_one_num(n)) return(S_One)
  if (.isT(n, "prime")) return(.sub(n, S_One))
  if (n$type == "Integer") {
    f <- .factorint(.nt_int(n))
    ps <- .fac_primes(f)
    return(.zsym(prod(ps^(f - 1) * (ps - 1))))
  }
  NULL
}, facts = .nt_facts(integer = TRUE, positive = TRUE), latex = .nt_latex1("\\phi"))

.reduced_totient <- function(n) {
  rr <- .remove(n, 2)
  n <- rr[[1L]]
  t <- rr[[2L]]
  t <- if (t == 0) gmp::as.bigz(1) else if (t > 2) gmp::as.bigz(2)^(t - 2) else gmp::as.bigz(t)
  f <- .factorint(n)
  for (p in names(f)) {
    pz <- gmp::as.bigz(p)
    v <- (pz - 1) * pz^(f[[p]] - 1)
    t <- t * v %/% .gcd2(t, v)
  }
  t
}

def_function("reduced_totient", eval = function(args) {
  n <- args[[1L]]
  .nt_check(n, "n", positive = TRUE)
  if (is_one_num(n)) return(S_One)
  if (.isT(n, "prime")) return(.sub(n, S_One))
  if (n$type == "Integer") return(.zsym(.reduced_totient(.nt_int(n))))
  NULL
}, facts = .nt_facts(integer = TRUE, positive = TRUE), latex = .nt_latex1("\\lambda"))

def_function("divisor_sigma", nargs = c(1L, 2L), eval = function(args) {
  n <- args[[1L]]
  k <- if (length(args) > 1L) args[[2L]] else S_One
  .nt_check(n, "n", positive = TRUE)
  .nt_check(k, "k", nonneg = TRUE)
  if (.isT(n, "prime")) return(.add(list(S_One, .pow(n, k))))
  if (is_one_num(n)) return(S_One)
  if (n$type == "Integer") {
    f <- .factorint(.nt_int(n))
    if (.isT(k, "zero")) return(.zsym(prod(gmp::as.bigz(f + 1))))
    if (k$type == "Integer") return(.zsym(.divisor_sigma(.nt_int(n), as.numeric(.nt_int(k)))))
    if (.isF(k, "zero")) {
      terms <- lapply(names(f), function(p) {
        pz <- Integer(p)
        cancel(.div(.sub(.pow(pz, .mul(list(k, Integer(f[[p]] + 1)))), S_One), .sub(.pow(pz, k), S_One)))
      })
      return(.mul(terms))
    }
  }
  NULL
}, facts = .nt_facts(integer = TRUE, positive = TRUE),
latex = function(a) if (length(a) > 1L) paste0("\\sigma_{", a[2], "}\\left(", a[1], "\\right)")
                    else paste0("\\sigma\\left(", a[1], "\\right)"))

def_function("udivisor_sigma", nargs = c(1L, 2L), eval = function(args) {
  n <- args[[1L]]
  k <- if (length(args) > 1L) args[[2L]] else S_One
  .nt_check(n, "n", positive = TRUE)
  .nt_check(k, "k", nonneg = TRUE)
  if (.isT(n, "prime")) return(.add(list(S_One, .pow(n, k))))
  if (n$type == "Integer") {
    f <- .factorint(.nt_int(n))
    return(.mul(lapply(names(f), function(p) .add(list(S_One, .pow(Integer(p), .mul(list(k, Integer(f[[p]])))))))))
  }
  NULL
}, facts = .nt_facts(integer = TRUE, positive = TRUE),
latex = function(a) if (length(a) > 1L) paste0("\\sigma^*_{", a[2], "}\\left(", a[1], "\\right)")
                    else paste0("\\sigma^*\\left(", a[1], "\\right)"))

def_function("primenu", eval = function(args) {
  n <- args[[1L]]
  .nt_check(n, "n", positive = TRUE)
  if (.isT(n, "prime")) return(S_One)
  if (is_one_num(n)) return(S_Zero)
  if (n$type == "Integer") return(Integer(length(.factorint(.nt_int(n)))))
  NULL
}, facts = .nt_facts(integer = TRUE, nonnegative = TRUE), latex = .nt_latex1("\\nu"))

def_function("primeomega", eval = function(args) {
  n <- args[[1L]]
  .nt_check(n, "n", positive = TRUE)
  if (.isT(n, "prime")) return(S_One)
  if (is_one_num(n)) return(S_Zero)
  if (n$type == "Integer") return(Integer(sum(.factorint(.nt_int(n)))))
  NULL
}, facts = .nt_facts(integer = TRUE, nonnegative = TRUE), latex = .nt_latex1("\\Omega"))

def_function("mobius", eval = function(args) {
  n <- args[[1L]]
  .nt_check(n, "n", positive = TRUE)
  if (.isT(n, "prime")) return(S_NegativeOne)
  if (is_one_num(n)) return(S_One)
  terms <- if (n$type == "Mul") n$args else list(n)
  result <- 1
  unresolved <- FALSE
  for (t in terms) {
    be <- .base_exp(t)
    m <- be[[1L]]
    e <- be[[2L]]
    if (!(.isT(m, "integer") && .isT(m, "positive") && .isT(e, "integer") && .isT(e, "positive"))) return(NULL)
    lt <- ask_fact(.sub(e, S_One), "positive")
    if (isTRUE(lt)) return(S_Zero)
    if (m$type == "Integer") {
      f <- .factorint(.nt_int(m))
      if (any(f > 1)) return(S_Zero)
      if (isFALSE(lt)) result <- result * (-1)^length(f) else unresolved <- TRUE
    } else {
      # SymPy returns a value here, ignoring the unknown factor
      unresolved <- TRUE
    }
  }
  if (unresolved) NULL else Integer(result)
}, facts = .nt_facts(integer = TRUE, prime = FALSE), latex = .nt_latex1("\\mu"))

def_function("legendre_symbol", nargs = 2L, eval = function(args) {
  a <- args[[1L]]
  p <- args[[2L]]
  .nt_check(a, "a")
  .nt_check(p, "p")
  if (.isF(p, "prime") || .isF(p, "odd")) stop(.sympy_condition("ValueError", "p should be an odd prime integer"))
  if (a$type == "Integer" && p$type == "Integer") {
    if (.nt_int(a) %% .nt_int(p) == 0) return(S_Zero)
    return(Integer(.legendre(.nt_int(a), .nt_int(p))))
  }
  if (is_one_num(a)) return(S_One)
  NULL
}, facts = .nt_facts(integer = TRUE, prime = FALSE),
latex = function(a) paste0("\\left(\\frac{", a[1], "}{", a[2], "}\\right)"))

def_function("jacobi_symbol", nargs = 2L, eval = function(args) {
  m <- args[[1L]]
  n <- args[[2L]]
  .nt_check(m, "m")
  .nt_check(n, "n")
  if (.isF(n, "positive") || .isF(n, "odd")) stop(.sympy_condition("ValueError", "n should be an odd positive integer"))
  if (is_one_num(m) || is_one_num(n)) return(S_One)
  if (m$type == "Integer" && n$type == "Integer") return(Integer(.jacobi(.nt_int(m), .nt_int(n))))
  NULL
}, facts = .nt_facts(integer = TRUE, prime = FALSE),
latex = function(a) paste0("\\left(\\frac{", a[1], "}{", a[2], "}\\right)"))

def_function("kronecker_symbol", nargs = 2L, eval = function(args) {
  a <- args[[1L]]
  n <- args[[2L]]
  .nt_check(a, "a")
  .nt_check(n, "n")
  if (is_one_num(a) || is_one_num(n)) return(S_One)
  if (a$type == "Integer" && n$type == "Integer") return(Integer(.kronecker(.nt_int(a), .nt_int(n))))
  NULL
}, facts = .nt_facts(integer = TRUE, prime = FALSE),
latex = function(a) paste0("\\left(\\frac{", a[1], "}{", a[2], "}\\right)"))

def_function("primepi", eval = function(args) {
  n <- args[[1L]]
  if (is_const(n, "oo")) return(S_Infinity)
  if (is_const(n, "-oo")) return(S_Zero)
  if (.isF(n, "real")) stop("n should be a real", call. = FALSE)
  if (isTRUE(ask_fact(.sub(n, S_Two), "negative"))) return(S_Zero)
  if (length(free_symbols(n))) return(NULL)
  v <- if (is_exact_number(n)) gmp::numerator(n$value) %/% gmp::denominator(n$value) else
    if (n$type == "Float") gmp::as.bigz(floor(n$value)) else .nt_certified_floor(n)
  v <- gmp::as.bigz(v)
  if (v < 2) return(S_Zero)
  if (v >= .NT_Z53) stop("primepi() is limited to arguments below 2^53", call. = FALSE)
  Integer(.primepi_count(as.numeric(v)))
}, facts = .nt_facts(integer = TRUE, nonnegative = TRUE), latex = .nt_latex1("\\pi"))

.nt_arg <- function(x) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x == trunc(x) && abs(x) >= 2^53) .nt_int(x)
  if (is.character(x) && length(x) == 1L && grepl("^\\s*[-+]?[0-9]+\\s*$", x)) return(Integer(.nt_int(x)))
  as_sym(x)
}

#' Number theoretic functions
#'
#' Symbolic functions that evaluate on integers and stay unevaluated on
#' symbols, using their assumptions where possible (`totient(p)` is `p - 1`
#' for a prime symbol `p`): Euler's `totient(n)`, Carmichael's
#' `reduced_totient(n)`, `divisor_sigma(n, k)` (the sum of the `k`-th powers of
#' the divisors) and its unitary analogue `udivisor_sigma(n, k)`, `primenu(n)`
#' and `primeomega(n)` (distinct and repeated prime factors), the Moebius
#' function `mobius(n)`, the Legendre, Jacobi and Kronecker symbols, and the
#' prime counting function `primepi(n)`. `totient()` and `reduced_totient()`
#' also accept a factorization as returned by [factorint()].
#'
#' @param n,a,m,p Arguments.
#' @param k Power; `NULL` (the default) means 1.
#' @return A `sym` expression.
#' @examples
#' totient(25)
#' totient(Symbol("p", prime = TRUE))
#' divisor_sigma(18, 2)
#' mobius(30)
#' legendre_symbol(3, 7)
#' jacobi_symbol(45, 77)
#' primepi(100)
#' @export
totient <- function(n) {
  if (!is_sym(n) && !is.null(names(n))) {
    if (!length(n)) return(S_One)
    ps <- gmp::as.bigz(names(n))
    ex <- as.numeric(unlist(n))
    return(.zsym(prod(ps^(ex - 1) * (ps - 1))))
  }
  fcall("totient", list(.nt_arg(n)))
}

#' @rdname totient
#' @export
reduced_totient <- function(n) {
  if (!is_sym(n) && !is.null(names(n))) {
    if (!length(n)) return(S_One)
    ps <- gmp::as.bigz(names(n))
    ex <- as.numeric(unlist(n))
    t <- gmp::as.bigz(1)
    for (i in seq_along(ex)) {
      v <- if (ps[i] == 2 && ex[i] > 2) gmp::as.bigz(2)^(ex[i] - 2) else (ps[i] - 1) * ps[i]^(ex[i] - 1)
      t <- t * v %/% .gcd2(t, v)
    }
    return(.zsym(t))
  }
  fcall("reduced_totient", list(.nt_arg(n)))
}

#' @rdname totient
#' @export
divisor_sigma <- function(n, k = NULL) fcall("divisor_sigma", c(list(.nt_arg(n)), if (!is.null(k)) list(.nt_arg(k))))

#' @rdname totient
#' @export
udivisor_sigma <- function(n, k = NULL) fcall("udivisor_sigma", c(list(.nt_arg(n)), if (!is.null(k)) list(.nt_arg(k))))

#' @rdname totient
#' @export
primenu <- function(n) fcall("primenu", list(.nt_arg(n)))

#' @rdname totient
#' @export
primeomega <- function(n) fcall("primeomega", list(.nt_arg(n)))

#' @rdname totient
#' @export
mobius <- function(n) fcall("mobius", list(.nt_arg(n)))

#' @rdname totient
#' @export
legendre_symbol <- function(a, p) fcall("legendre_symbol", list(.nt_arg(a), .nt_arg(p)))

#' @rdname totient
#' @export
jacobi_symbol <- function(m, n) fcall("jacobi_symbol", list(.nt_arg(m), .nt_arg(n)))

#' @rdname totient
#' @export
kronecker_symbol <- function(a, n) fcall("kronecker_symbol", list(.nt_arg(a), .nt_arg(n)))

#' @rdname totient
#' @export
primepi <- function(n) fcall("primepi", list(.nt_arg(n)))
