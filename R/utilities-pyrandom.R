# Python's random module (MT19937 with CPython's seeding, getrandbits and
# randrange), so that functions taking SymPy's `seed` arguments (an integer
# seeds sympy.core.random's generator) reproduce SymPy's values exactly.
# Words are 32 bit unsigned integers held in doubles; products are split in
# 16 bit halves to stay exact.

.U32 <- 4294967296
.lo16 <- function(a) a %% 65536
.hi16 <- function(a) a %/% 65536
.xor32 <- function(a, b) bitwXor(.hi16(a), .hi16(b)) * 65536 + bitwXor(.lo16(a), .lo16(b))
.and32 <- function(a, b) bitwAnd(.hi16(a), .hi16(b)) * 65536 + bitwAnd(.lo16(a), .lo16(b))
.or32 <- function(a, b) bitwOr(.hi16(a), .hi16(b)) * 65536 + bitwOr(.lo16(a), .lo16(b))
.mul32 <- function(a, b) {
  ah <- .hi16(a); al <- .lo16(a); bh <- .hi16(b); bl <- .lo16(b)
  (((ah * bl + al * bh) %% 65536) * 65536 + al * bl) %% .U32
}
.shr32 <- function(a, k) a %/% 2^k
.shl32 <- function(a, k) (a * 2^k) %% .U32

.py_rng <- new.env()

.mt_init_genrand <- function(s) {
  mt <- numeric(624L)
  mt[1L] <- s %% .U32
  for (i in 2:624) {
    p <- mt[i - 1L]
    mt[i] <- (.mul32(1812433253, .xor32(p, .shr32(p, 30))) + (i - 1)) %% .U32
  }
  mt
}

# init_by_array of CPython's _randommodule.c
.mt_init_by_array <- function(key) {
  mt <- .mt_init_genrand(19650218)
  N <- 624L; i <- 2L; j <- 1L; kl <- length(key)
  for (k in seq_len(max(N, kl))) {
    p <- mt[i - 1L]
    mt[i] <- (.xor32(mt[i], .mul32(.xor32(p, .shr32(p, 30)), 1664525)) + key[j] + (j - 1)) %% .U32
    i <- i + 1L; j <- j + 1L
    if (i > N) { mt[1L] <- mt[N]; i <- 2L }
    if (j > kl) j <- 1L
  }
  for (k in seq_len(N - 1L)) {
    p <- mt[i - 1L]
    mt[i] <- (.xor32(mt[i], .mul32(.xor32(p, .shr32(p, 30)), 1566083941)) - (i - 1)) %% .U32
    i <- i + 1L
    if (i > N) { mt[1L] <- mt[N]; i <- 2L }
  }
  mt[1L] <- 2147483648
  mt
}

# random.seed(a) for an integer a: the 32 bit words of abs(a), least
# significant first
.py_seed <- function(a) {
  a <- abs(gmp::as.bigz(a))
  key <- numeric()
  repeat {
    key <- c(key, as.numeric(a %% .U32))
    a <- a %/% .U32
    if (a == 0) break
  }
  .py_rng$mt <- .mt_init_by_array(key)
  .py_rng$mti <- 625L
  invisible(NULL)
}

.mt_twist <- function(mt) {
  N <- 624L; M <- 397L
  for (kk in seq_len(N)) {
    y <- .or32(.and32(mt[kk], 2147483648), .and32(mt[(kk %% N) + 1L], 2147483647))
    v <- .xor32(mt[((kk - 1L + M) %% N) + 1L], .shr32(y, 1))
    if (y %% 2 == 1) v <- .xor32(v, 2567483615)
    mt[kk] <- v
  }
  mt
}

.py_genrand <- function() {
  if (is.null(.py_rng$mt)) .py_seed(as.numeric(Sys.time()) * 1e6)
  if (.py_rng$mti > 624L) { .py_rng$mt <- .mt_twist(.py_rng$mt); .py_rng$mti <- 1L }
  y <- .py_rng$mt[.py_rng$mti]
  .py_rng$mti <- .py_rng$mti + 1L
  y <- .xor32(y, .shr32(y, 11))
  y <- .xor32(y, .and32(.shl32(y, 7), 2636928640))
  y <- .xor32(y, .and32(.shl32(y, 15), 4022730752))
  .xor32(y, .shr32(y, 18))
}

# random.getrandbits(k), as a bigz
.py_getrandbits <- function(k) {
  if (k <= 32) return(gmp::as.bigz(.shr32(.py_genrand(), 32 - k)))
  out <- gmp::as.bigz(0); sh <- gmp::as.bigz(1)
  while (k > 0) {
    r <- .py_genrand()
    if (k < 32) r <- .shr32(r, 32 - k)
    out <- out + gmp::as.bigz(r) * sh
    sh <- sh * .U32
    k <- k - 32
  }
  out
}

.bit_length <- function(n) {
  n <- gmp::as.bigz(n)
  if (n == 0) return(0L)
  nchar(as.character(n, b = 2))
}

# random.randrange(start, stop)
.py_randrange <- function(a, b = NULL) {
  if (is.null(b)) { b <- a; a <- 0 }
  a <- gmp::as.bigz(a); w <- gmp::as.bigz(b) - a
  if (w <= 0) stop("empty range for randrange()", call. = FALSE)
  k <- .bit_length(w)
  r <- .py_getrandbits(k)
  while (r >= w) r <- .py_getrandbits(k)
  a + r
}

.as_out <- function(v) if (abs(v) < 2^52) as.numeric(v) else v

# sympy.core.random._randrange(seed): NULL draws from the generator as it
# stands, an integer reseeds it, a vector supplies the values in order
.sympy_randrange <- function(seed = NULL) {
  if (is.null(seed)) return(function(a, b = NULL) .as_out(.py_randrange(a, b)))
  if (length(seed) == 1L && !is.list(seed) && !isTRUE(attr(seed, "sequence"))) {
    .py_seed(seed)
    return(function(a, b = NULL) .as_out(.py_randrange(a, b)))
  }
  seq <- as.list(seed)
  function(a, b = NULL) {
    if (is.null(b)) { b <- a; a <- 0 }
    if (b - a < 1) stop("_randrange got empty range", call. = FALSE)
    repeat {
      if (!length(seq)) stop("_randrange sequence was too short", call. = FALSE)
      x <- seq[[1L]]; seq <<- seq[-1L]
      if (a <= x && x < b) return(x)
    }
  }
}

# sympy.core.random._randint(seed)
.sympy_randint <- function(seed = NULL) {
  rr <- .sympy_randrange(seed)
  if (is.null(seed) || (length(seed) == 1L && !is.list(seed) && !isTRUE(attr(seed, "sequence"))))
    return(function(a, b) rr(a, b + 1))
  seq <- as.list(seed)
  function(a, b) {
    if (b - a < 0) stop("_randint got empty range", call. = FALSE)
    repeat {
      if (!length(seq)) stop("_randint sequence was too short", call. = FALSE)
      x <- seq[[1L]]; seq <<- seq[-1L]
      if (a <= x && x <= b) return(x)
    }
  }
}
