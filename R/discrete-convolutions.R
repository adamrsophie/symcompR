# Convolutions, after sympy/discrete/convolutions.py, with the transforms of
# sympy/discrete/transforms.py they run on.

.dc_not_subscriptable <- function(x) {
  stop(.sympy_condition("TypeError", sprintf("'%s' object is not subscriptable", .sympy_class(x))))
}

# a Python sequence as a list: numbers stay R numbers (Python int or float)
.dc_seq <- function(x) if (is_sym(x)) .dc_not_subscriptable(x) else as.list(x)

.dcv_sym <- function(v) {
  if (is_sym(v)) return(v)
  if (is.numeric(v) && v == round(v) && abs(v) < 2^53) return(.int(v))
  if (is.complex(v)) return(as_sym(v))
  as_sym(v)
}

.dc_pow2 <- function(n) {
  if (n > 0 && bitwAnd(n, n - 1L)) 2L^(floor(log2(n)) + 1L) else n
}

.dc_bitrev <- function(a, b) {
  n <- length(a)
  for (i in seq_len(n - 1L)) {
    j <- sum(as.integer(intToBits(i))[seq_len(b)] * 2^(rev(seq_len(b)) - 1))
    if (i < j) { t <- a[[i + 1L]]; a[[i + 1L]] <- a[[j + 1L]]; a[[j + 1L]] <- t }
  }
  a
}

# sympy.discrete.transforms._fourier_transform
.dc_fourier <- function(seq, dps = NULL, inverse = FALSE) {
  if (is_sym(seq))
    stop(.sympy_condition("TypeError", "Expected a sequence of numeric coefficients for Fourier Transform"))
  a <- lapply(as.list(seq), .dcv_sym)
  if (any(vapply(a, function(x) length(free_symbols(x)) > 0L, TRUE)))
    stop(.sympy_condition("ValueError", "Expected non-symbolic coefficients"))
  n <- length(a)
  if (n < 2L) return(a)
  b <- floor(log2(n))
  if (bitwAnd(n, n - 1L)) { b <- b + 1L; n <- 2L^b }
  a <- c(a, rep(list(S_Zero), n - length(a)))
  a <- .dc_bitrev(a, b)
  ang <- if (inverse) -2 * S_Pi / .int(n) else 2 * S_Pi / .int(n)
  if (!is.null(dps)) ang <- N(ang, dps + 2)
  w <- lapply(seq_len(n %/% 2L) - 1L, function(i) {
    t <- ang * .int(i)
    fcall("cos", list(t)) + S_I * fcall("sin", list(t))
  })
  h <- 2L
  while (h <= n) {
    hf <- h %/% 2L; ut <- n %/% h
    for (i in seq(0L, n - 1L, by = h)) for (j in seq_len(hf) - 1L) {
      u <- a[[i + j + 1L]]
      v <- expand_mul(a[[i + j + hf + 1L]] * w[[ut * j + 1L]])
      a[[i + j + 1L]] <- u + v
      a[[i + j + hf + 1L]] <- u - v
    }
    h <- h * 2L
  }
  if (inverse) a <- if (!is.null(dps)) lapply(a, function(x) N(x / .int(n), dps)) else lapply(a, function(x) x / .int(n))
  a
}

# sympy.discrete.transforms._number_theoretic_transform, on gmp integers
.dc_ntt <- function(seq, prime, inverse = FALSE) {
  if (is_sym(seq))
    stop(.sympy_condition("TypeError", "Expected a sequence of integer coefficients for Number Theoretic Transform"))
  p <- .as_int(prime)
  if (!isprime(p)) stop(.sympy_condition("ValueError", "Expected prime modulus for Number Theoretic Transform"))
  if (!length(seq)) return(NULL)
  a <- .zc(lapply(as.list(seq), function(x) .as_int(x) %% p))
  n <- length(a)
  b <- floor(log2(n))
  if (bitwAnd(n, n - 1L)) { b <- b + 1L; n <- 2L^b }
  if ((p - 1) %% n != 0) stop(.sympy_condition("ValueError", "Expected prime modulus of the form (m*2**k + 1)"))
  a <- c(a, gmp::as.bigz(rep(0, n - length(a))))
  for (i in seq_len(n - 1L)) {
    j <- sum(as.integer(intToBits(i))[seq_len(b)] * 2^(rev(seq_len(b)) - 1))
    if (i < j) { t <- a[i + 1L]; a[i + 1L] <- a[j + 1L]; a[j + 1L] <- t }
  }
  pr <- .primitive_root(p, TRUE)
  if (is.null(pr)) stop(.sympy_condition("ValueError", "Expected prime modulus to have a primitive root"))
  rt <- gmp::powm(pr, (p - 1) %/% n, p)
  if (inverse) rt <- gmp::powm(rt, p - 2, p)
  w <- gmp::as.bigz(rep(1, n %/% 2L))
  if (n %/% 2L >= 2L) for (i in 2:(n %/% 2L)) w[i] <- (w[i - 1L] * rt) %% p
  h <- 2L
  while (h <= n) {
    hf <- h %/% 2L; ut <- n %/% h
    for (i in seq(0L, n - 1L, by = h)) for (j in seq_len(hf) - 1L) {
      u <- a[i + j + 1L]
      v <- a[i + j + hf + 1L] * w[ut * j + 1L]
      a[i + j + 1L] <- (u + v) %% p
      a[i + j + hf + 1L] <- (u - v) %% p
    }
    h <- h * 2L
  }
  if (inverse) {
    rv <- gmp::powm(gmp::as.bigz(n), p - 2, p)
    a <- (a * rv) %% p
  }
  a
}

# sympy.discrete.transforms._walsh_hadamard_transform
.dc_wht <- function(seq, inverse = FALSE) {
  if (is_sym(seq))
    stop(.sympy_condition("TypeError", "Expected a sequence of coefficients for Walsh Hadamard Transform"))
  a <- lapply(as.list(seq), .dcv_sym)
  n <- length(a)
  if (n < 2L) return(a)
  n <- .dc_pow2(n)
  a <- c(a, rep(list(S_Zero), n - length(a)))
  h <- 2L
  while (h <= n) {
    hf <- h %/% 2L
    for (i in seq(0L, n - 1L, by = h)) for (j in seq_len(hf) - 1L) {
      u <- a[[i + j + 1L]]; v <- a[[i + j + hf + 1L]]
      a[[i + j + 1L]] <- u + v
      a[[i + j + hf + 1L]] <- u - v
    }
    h <- h * 2L
  }
  if (inverse) a <- lapply(a, function(x) x / .int(n))
  a
}

# sympy.discrete.transforms._mobius_transform
.dc_mobius <- function(seq, sgn, subset) {
  if (is_sym(seq)) stop(.sympy_condition("TypeError", "Expected a sequence of coefficients"))
  a <- lapply(as.list(seq), .dcv_sym)
  n <- length(a)
  if (n < 2L) return(a)
  n <- .dc_pow2(n)
  a <- c(a, rep(list(S_Zero), n - length(a)))
  i <- 1L
  while (i < n) {
    for (j in seq_len(n) - 1L) {
      if (xor(bitwAnd(j, i) != 0L, !subset)) a[[j + 1L]] <- a[[j + 1L]] + .int(sgn) * a[[bitwXor(j, i) + 1L]]
    }
    i <- i * 2L
  }
  a
}

#' Convolutions
#'
#' `convolution(a, b)` is SymPy's `convolution`: the linear convolution of
#' two sequences, cyclic of length `cycle` when `cycle > 0`. Integer and
#' rational sequences use `convolution_int()` (Kronecker substitution); other
#' sequences use `convolution_fft()`, with `dps` decimal digits when given.
#' `prime` asks for `convolution_ntt()` modulo a prime of the form
#' `m 2^k + 1`, `dyadic = TRUE` for the dyadic (XOR) convolution
#' `convolution_fwht()` and `subset = TRUE` for the subset convolution
#' `convolution_subset()`. At most one of these may be given.
#'
#' @param a,b Sequences (vectors or lists).
#' @param cycle Length of a cyclic convolution (0 for linear).
#' @param dps Decimal digits for a floating point transform.
#' @param prime A prime modulus for the number theoretic transform.
#' @param dyadic Use the dyadic convolution.
#' @param subset Use the subset convolution.
#' @return A list of expressions.
#' @examples
#' convolution(c(1, 2, 3), c(4, 5, 6))
#' convolution(c(1, 2, 3), c(4, 5, 6), cycle = 3)
#' convolution(c(1, 2, 3), c(4, 5, 6), prime = 19 * 2^10 + 1)
#' convolution_fwht(c(1, 2, 3), c(4, 5, 6))
#' convolution_int(c(10, -5, 1, 3), c(-5, 6, 7))
#' @rdname convolution
#' @export
convolution <- function(a, b, cycle = 0, dps = NULL, prime = NULL, dyadic = NULL, subset = NULL) {
  c <- as.integer(.as_int(cycle))
  if (c < 0L) stop(.sympy_condition("ValueError", "The length for cyclic convolution must be non-negative"))
  dyadic <- if (isTRUE(dyadic)) TRUE else NULL
  subset <- if (isTRUE(subset)) TRUE else NULL
  if (sum(!vapply(list(prime, dps, dyadic, subset), is.null, TRUE)) > 1L)
    stop(.sympy_condition("TypeError", "Ambiguity in determining the type of convolution"))
  if (!is.null(prime)) {
    ls <- convolution_ntt(a, b, prime = prime)
    if (!c) return(ls)
    p <- .as_int(prime)
    return(lapply(seq_len(c), function(i) {
      s <- if (i <= length(ls)) sum(.zc(lapply(ls[seq(i, length(ls), by = c)], function(z) gmp::as.bigz(z$value)))) else gmp::as.bigz(0)
      .int(s %% p)
    }))
  }
  if (isTRUE(dyadic)) ls <- convolution_fwht(a, b)
  else if (isTRUE(subset)) ls <- convolution_subset(a, b)
  else {
    loop <- function(a) {
      dens <- list()
      for (i in a) {
        if (is_sym(i) && i$type == "Rational") dens[[length(dens) + 1L]] <- gmp::denominator(i$value)
        else if (!((is.numeric(i) && i == round(i)) || (is_sym(i) && i$type == "Integer"))) return(NULL)
      }
      ints <- lapply(a, function(i) if (is_sym(i)) gmp::numerator(i$value) else gmp::as.bigz(i))
      if (length(dens)) {
        l <- Reduce(gmp::lcm.bigz, dens)
        return(list(lapply(a, function(i) if (is_sym(i)) gmp::numerator(i$value) * (l %/% gmp::denominator(i$value)) else gmp::as.bigz(i) * l), l))
      }
      list(ints, gmp::as.bigz(1))
    }
    al <- .dc_seq(a); bl <- .dc_seq(b)
    ls <- NULL
    da <- loop(al)
    if (!is.null(da)) {
      db <- loop(bl)
      if (!is.null(db)) {
        den <- da[[2L]] * db[[2L]]
        li <- .dc_conv_int(da[[1L]], db[[1L]])
        li <- c(li, rep(list(gmp::as.bigz(0)), length(al) + length(bl) - 1L - length(li)))
        ls <- lapply(li, function(z) make_rational(gmp::as.bigq(z, den)))
      }
    }
    if (is.null(ls)) ls <- convolution_fft(a, b, dps)
  }
  if (!c) return(ls)
  lapply(seq_len(c), function(i) {
    if (i > length(ls)) return(S_Zero)
    Reduce(`+`, ls[seq(i, length(ls), by = c)], S_Zero)
  })
}

#' @rdname convolution
#' @export
convolution_fft <- function(a, b, dps = NULL) {
  a <- .dc_seq(a); b <- .dc_seq(b)
  n <- m <- length(a) + length(b) - 1L
  n <- .dc_pow2(n)
  if (n > 0L) {
    a <- c(a, rep(list(S_Zero), n - length(a)))
    b <- c(b, rep(list(S_Zero), n - length(b)))
  }
  a <- .dc_fourier(a, dps); b <- .dc_fourier(b, dps)
  a <- Map(function(x, y) expand_mul(x * y), a, b)
  a <- .dc_fourier(a, dps, inverse = TRUE)
  a[seq_len(max(m, 0L))]
}

#' @rdname convolution
#' @export
convolution_ntt <- function(a, b, prime) {
  a <- .dc_seq(a); b <- .dc_seq(b)
  p <- .as_int(prime)
  n <- m <- length(a) + length(b) - 1L
  n <- .dc_pow2(n)
  if (n > 0L) {
    a <- c(a, rep(list(0), n - length(a)))
    b <- c(b, rep(list(0), n - length(b)))
  }
  a <- .dc_ntt(a, p); b <- .dc_ntt(b, p)
  if (m <= 0L) return(list())
  a <- (a * b) %% p
  a <- .dc_ntt(a, p, inverse = TRUE)
  lapply(seq_len(max(m, 0L)), function(i) .int(a[i]))
}

#' @rdname convolution
#' @export
convolution_fwht <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  a <- .dc_seq(a); b <- .dc_seq(b)
  n <- .dc_pow2(max(length(a), length(b)))
  a <- c(a, rep(list(S_Zero), n - length(a)))
  b <- c(b, rep(list(S_Zero), n - length(b)))
  a <- .dc_wht(a); b <- .dc_wht(b)
  a <- Map(function(x, y) expand_mul(x * y), a, b)
  .dc_wht(a, inverse = TRUE)
}

#' @rdname convolution
#' @export
convolution_subset <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  if (is_sym(a) || is_sym(b)) stop(.sympy_condition("TypeError", "Expected a sequence of coefficients for convolution"))
  a <- lapply(as.list(a), .dcv_sym); b <- lapply(as.list(b), .dcv_sym)
  n <- .dc_pow2(max(length(a), length(b)))
  a <- c(a, rep(list(S_Zero), n - length(a)))
  b <- c(b, rep(list(S_Zero), n - length(b)))
  out <- rep(list(S_Zero), n)
  for (mask in seq_len(n) - 1L) {
    smask <- mask
    while (smask > 0L) {
      out[[mask + 1L]] <- out[[mask + 1L]] + expand_mul(a[[smask + 1L]] * b[[bitwXor(mask, smask) + 1L]])
      smask <- bitwAnd(smask - 1L, mask)
    }
    out[[mask + 1L]] <- out[[mask + 1L]] + expand_mul(a[[smask + 1L]] * b[[bitwXor(mask, smask) + 1L]])
  }
  out
}

#' @rdname convolution
#' @export
convolution_int <- function(a, b) {
  r <- .dc_conv_int(lapply(as.list(a), .as_int), lapply(as.list(b), .as_int))
  lapply(r, .int)
}

# Kronecker substitution on gmp integers
.dc_conv_int <- function(a, b) {
  a <- .zc(a); b <- .zc(b)
  B <- max(abs(a)) * max(abs(b)) * (1 + min(length(a) - 1L, length(b) - 1L))
  x <- gmp::as.bigz(1); power <- 0L
  while (x <= 2 * B) { x <- x * 2; power <- power + 1L }
  sh <- gmp::as.bigz(2)^power
  to_integer <- function(poly) {
    n <- gmp::as.bigz(0); mul <- 0L
    for (c in rev(as.list(poly))) {
      if (c != 0 && !mul) mul <- if (c < 0) -1L else 1L
      n <- n * sh + mul * c
    }
    list(mul, n)
  }
  pa <- to_integer(a); pb <- to_integer(b)
  result <- pa[[2L]] * pb[[2L]]
  mul <- pa[[1L]] * pb[[1L]]
  half <- x %/% 2
  borrow <- 0L
  poly <- list()
  while (result != 0 || borrow) {
    coeff <- (result %% x) + borrow
    result <- result %/% sh
    borrow <- as.integer(coeff >= half)
    poly[[length(poly) + 1L]] <- mul * (if (coeff < half) coeff else coeff - x)
  }
  if (!length(poly)) list(gmp::as.bigz(0)) else poly
}
