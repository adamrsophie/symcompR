# nsimplify() and _real_to_rational() of sympy/simplify/simplify.py, with
# the number identification of mpmath (mpmath/identification.py: pslq,
# findpoly and identify) in exact integer fixed-point arithmetic.

# ---------------------------------------------------------------------------
# mpmath helpers on mpf values list(m, e)

.id_mpf <- function(x, prec) {
  if (is.list(x) && !is_sym(x)) return(.mp_round(x$m, x$e, prec))
  if (inherits(x, "bigz") || inherits(x, "bigq")) return(.mp_from_q(gmp::as.bigq(x), prec))
  .mp_round(.mp_from_dbl(as.numeric(x))$m, .mp_from_dbl(as.numeric(x))$e, prec)
}
.id_abs <- function(a) .mpf(abs(a$m), a$e)
.id_neg <- function(a) .mpf(-a$m, a$e)
.id_lt <- function(a, b) .mp_sign(.mp_sub(a, b, 8L)) < 0
.id_zero <- function(a) a$m == 0
.id_one <- .mpf(1, 0)
# to_fixed: floor(x * 2^prec)
.id_fixed <- function(a, prec) {
  off <- a$e + prec
  if (off >= 0) a$m * .pow2(off) else a$m %/% .pow2(-off)
}
.ID_P2 <- new.env()
.id_p2 <- function(k) {
  key <- as.character(k)
  v <- .ID_P2[[key]]
  if (is.null(v)) { v <- .pow2(k); assign(key, v, envir = .ID_P2) }
  v
}
.id_shr <- function(a, k) a %/% .id_p2(k)
.id_shl <- function(a, k) a * .id_p2(k)
.id_isqrt <- function(n) if (n <= 0) gmp::as.bigz(0) else .iroot(gmp::as.bigz(n), 2L)$root
.id_sqrt_fixed <- function(x, prec) .id_isqrt(.id_shl(x, prec))
.id_round_fixed <- function(x, prec) .id_shl(.id_shr(x + .id_p2(prec - 1), prec), prec)
.id_pow <- function(a, k, prec) { r <- .id_one; for (i in seq_len(k)) r <- .mp_mul(r, a, prec); .mp_round(r$m, r$e, prec) }

# mpmath pslq: an integer relation among the numbers x (mpf values) at
# working precision prec, or NULL
.id_pslq <- function(x, prec, tol = NULL, maxcoeff = 1000, maxsteps = 100) {
  n <- length(x)
  if (n < 2L) stop("ValueError: n cannot be less than 2", call. = FALSE)
  if (prec < 53) stop("ValueError: prec cannot be less than 53", call. = FALSE)
  target <- as.integer(prec * 0.75)
  tol <- if (is.null(tol)) .mpf(1, -target) else .id_mpf(tol, prec)
  prec <- prec + 60L
  tol <- .id_fixed(tol, prec)
  x <- lapply(x, function(xk) .id_fixed(.id_mpf(xk, prec - 60L), prec))
  minx <- Reduce(function(a, b) if (b < a) b else a, lapply(x, abs))
  if (minx == 0) stop("ValueError: PSLQ requires a vector of nonzero numbers", call. = FALSE)
  if (minx < tol %/% 100) return(NULL)
  g <- .id_sqrt_fixed(.id_shl(gmp::as.bigz(4), prec) %/% 3, prec)
  one <- .pow2(prec)
  zero <- gmp::as.bigz(0)
  A <- B <- H <- matrix(list(zero), n, n)
  for (i in seq_len(n)) { A[[i, i]] <- one; B[[i, i]] <- one }
  s <- vector("list", n)
  for (k in seq_len(n)) {
    t <- zero
    for (j in k:n) t <- t + .id_shr(x[[j]]^2, prec)
    s[[k]] <- .id_sqrt_fixed(t, prec)
  }
  t <- s[[1L]]
  y <- x
  for (k in seq_len(n)) {
    y[[k]] <- .id_shl(x[[k]], prec) %/% t
    s[[k]] <- .id_shl(s[[k]], prec) %/% t
  }
  for (i in seq_len(n)) {
    if (i <= n - 1L) H[[i, i]] <- if (s[[i]] != 0) .id_shl(s[[i + 1L]], prec) %/% s[[i]] else zero
    if (i > 1L) for (j in 1:(i - 1L)) {
      sjj1 <- s[[j]] * s[[j + 1L]]
      H[[i, j]] <- if (sjj1 != 0) .id_shl(-y[[i]] * y[[j]], prec) %/% sjj1 else zero
    }
  }
  reduce_row <- function(i, j) {
    t <- .id_round_fixed(.id_shl(H[[i, j]], prec) %/% H[[j, j]], prec)
    y[[j]] <<- y[[j]] + .id_shr(t * y[[i]], prec)
    for (k in seq_len(j)) H[[i, k]] <<- H[[i, k]] - .id_shr(t * H[[j, k]], prec)
    for (k in seq_len(n)) {
      A[[i, k]] <<- A[[i, k]] - .id_shr(t * A[[j, k]], prec)
      B[[k, j]] <<- B[[k, j]] + .id_shr(t * B[[k, i]], prec)
    }
  }
  if (n >= 2L) for (i in 2:n) for (j in (i - 1L):1) if (H[[j, j]] != 0) reduce_row(i, j)
  for (REP in seq_len(maxsteps)) {
    m <- -1L; szmax <- gmp::as.bigz(-1)
    for (i in seq_len(n - 1L)) {
      sz <- .id_shr(g^i * abs(H[[i, i]]), prec * (i - 1L))
      if (sz > szmax) { m <- i; szmax <- sz }
    }
    tmp <- y[[m]]; y[[m]] <- y[[m + 1L]]; y[[m + 1L]] <- tmp
    tmp <- H[m, ]; H[m, ] <- H[m + 1L, ]; H[m + 1L, ] <- tmp
    tmp <- A[m, ]; A[m, ] <- A[m + 1L, ]; A[m + 1L, ] <- tmp
    tmp <- B[, m]; B[, m] <- B[, m + 1L]; B[, m + 1L] <- tmp
    if (m <= n - 2L) {
      t0 <- .id_sqrt_fixed(.id_shr(H[[m, m]]^2 + H[[m, m + 1L]]^2, prec), prec)
      if (t0 == 0) break
      t1 <- .id_shl(H[[m, m]], prec) %/% t0
      t2 <- .id_shl(H[[m, m + 1L]], prec) %/% t0
      for (i in m:n) {
        t3 <- H[[i, m]]; t4 <- H[[i, m + 1L]]
        H[[i, m]] <- .id_shr(t1 * t3 + t2 * t4, prec)
        H[[i, m + 1L]] <- .id_shr(-t2 * t3 + t1 * t4, prec)
      }
    }
    if (m + 1L <= n) for (i in (m + 1L):n) {
      for (j in min(i - 1L, m + 1L):1) {
        if (H[[j, j]] == 0) break
        reduce_row(i, j)
      }
    }
    for (i in seq_len(n)) {
      if (abs(y[[i]]) < tol) {
        vec <- lapply(seq_len(n), function(j) .id_shr(.id_round_fixed(B[[j, i]], prec), prec))
        if (all(vapply(vec, function(v) abs(v) < maxcoeff, NA))) return(vec)
      }
    }
    recnorm <- Reduce(function(a, b) if (b > a) b else a, lapply(H, abs))
    if (recnorm != 0) {
      nrm <- .id_shr(.pow2(2L * prec) %/% recnorm, prec) %/% 100
      if (nrm >= maxcoeff) break
    }
  }
  NULL
}

.id_findpoly <- function(x, n, prec) {
  if (.id_zero(x)) return(list(gmp::as.bigz(1), gmp::as.bigz(0)))
  xs <- list(.id_one)
  for (i in seq_len(n)) {
    xs[[i + 1L]] <- .id_pow(x, i, prec)
    a <- .id_pslq(xs, prec)
    if (!is.null(a)) return(a)
  }
  NULL
}

.id_fracgcd <- function(p, q) {
  g <- gmp::gcd.bigz(p, q)
  if (g != 1) { p <- p %/% g; q <- q %/% g }
  if (q == 1) list(p) else list(p, q)
}
.id_pslqstring <- function(r, constants) {
  q <- r[[1L]]; r <- r[-1L]
  s <- character()
  for (i in seq_along(r)) {
    p <- r[[i]]
    if (p != 0) {
      z <- .id_fracgcd(-p, q)
      cs <- constants[[i]]$name
      cs <- if (cs == "1") "" else paste0("*", cs)
      term <- if (length(z) == 1L) { if (z[[1L]] > 0) paste0(as.character(z[[1L]]), cs) else paste0("(", as.character(z[[1L]]), ")", cs) }
              else paste0("(", as.character(z[[1L]]), "/", as.character(z[[2L]]), ")", cs)
      s <- c(s, term)
    }
  }
  s <- paste(s, collapse = " + ")
  if (grepl("+", s, fixed = TRUE) || grepl("*", s, fixed = TRUE)) s <- paste0("(", s, ")")
  if (nzchar(s)) s else "0"
}
.id_prodstring <- function(r, constants) {
  q <- r[[1L]]; r <- r[-1L]
  num <- character(); den <- character()
  for (i in seq_along(r)) {
    p <- r[[i]]
    if (p != 0) {
      z <- .id_fracgcd(-p, q)
      cs <- constants[[i]]$name
      if (length(z) == 1L) {
        t <- if (abs(z[[1L]]) == 1) cs else paste0(cs, "**", as.character(abs(z[[1L]])))
        if (z[[1L]] < 0) den <- c(den, t) else num <- c(num, t)
      } else {
        t <- paste0(cs, "**(", as.character(abs(z[[1L]])), "/", as.character(z[[2L]]), ")")
        if (z[[1L]] < 0) den <- c(den, t) else num <- c(num, t)
      }
    }
  }
  num <- paste(num, collapse = "*"); den <- paste(den, collapse = "*")
  if (nzchar(num) && nzchar(den)) return(paste0("(", num, ")/(", den, ")"))
  if (nzchar(num)) return(num)
  if (nzchar(den)) return(paste0("1/(", den, ")"))
  NULL
}
.id_quadraticstring <- function(t, a, b, c, prec) {
  if (c < 0) { a <- -a; b <- -b; c <- -c }
  D <- b^2 - 4 * a * c
  sq <- .mp_sqrt(.mp_from_q(gmp::as.bigq(D), prec), prec)
  u1 <- .mp_div(.mp_add(.mp_from_q(gmp::as.bigq(-b), prec), sq, prec), .mp_from_q(gmp::as.bigq(2 * c), prec), prec)
  u2 <- .mp_div(.mp_sub(.mp_from_q(gmp::as.bigq(-b), prec), sq, prec), .mp_from_q(gmp::as.bigq(2 * c), prec), prec)
  ch <- as.character
  if (.id_lt(.id_abs(.mp_sub(u1, t, prec)), .id_abs(.mp_sub(u2, t, prec)))) {
    if (b != 0) sprintf("((%s+sqrt(%s))/%s)", ch(-b), ch(D), ch(2 * c)) else sprintf("(sqrt(%s)/%s)", ch(-4 * a * c), ch(2 * c))
  } else {
    if (b != 0) sprintf("((%s-sqrt(%s))/%s)", ch(-b), ch(D), ch(2 * c)) else sprintf("(-sqrt(%s)/%s)", ch(-4 * a * c), ch(2 * c))
  }
}

.ID_TRANSFORMS <- list(
  list(function(x, c, p) .mp_mul(x, c, p), "$y/$c", FALSE),
  list(function(x, c, p) .mp_div(x, c, p), "$c*$y", TRUE),
  list(function(x, c, p) .mp_div(c, x, p), "$c/$y", FALSE),
  list(function(x, c, p) .id_pow(.mp_mul(x, c, p), 2L, p), "sqrt($y)/$c", FALSE),
  list(function(x, c, p) .id_pow(.mp_div(x, c, p), 2L, p), "$c*sqrt($y)", TRUE),
  list(function(x, c, p) .id_pow(.mp_div(c, x, p), 2L, p), "$c/sqrt($y)", FALSE),
  list(function(x, c, p) .mp_mul(c, .id_pow(x, 2L, p), p), "sqrt($y)/sqrt($c)", TRUE),
  list(function(x, c, p) .mp_div(.id_pow(x, 2L, p), c, p), "sqrt($c)*sqrt($y)", TRUE),
  list(function(x, c, p) .mp_div(c, .id_pow(x, 2L, p), p), "sqrt($c)/sqrt($y)", TRUE),
  list(function(x, c, p) .mp_sqrt(.mp_mul(x, c, p), p), "$y**2/$c", FALSE),
  list(function(x, c, p) .mp_sqrt(.mp_div(x, c, p), p), "$c*$y**2", TRUE),
  list(function(x, c, p) .mp_sqrt(.mp_div(c, x, p), p), "$c/$y**2", FALSE),
  list(function(x, c, p) .mp_mul(c, .mp_sqrt(x, p), p), "$y**2/$c**2", TRUE),
  list(function(x, c, p) .mp_div(.mp_sqrt(x, p), c, p), "$c**2*$y**2", TRUE),
  list(function(x, c, p) .mp_div(c, .mp_sqrt(x, p), p), "$c**2/$y**2", TRUE),
  list(function(x, c, p) .mp_exp(.mp_mul(x, c, p), p), "log($y)/$c", FALSE),
  list(function(x, c, p) .mp_exp(.mp_div(x, c, p), p), "$c*log($y)", TRUE),
  list(function(x, c, p) .mp_exp(.mp_div(c, x, p), p), "$c/log($y)", FALSE),
  list(function(x, c, p) .mp_mul(c, .mp_exp(x, p), p), "log($y/$c)", TRUE),
  list(function(x, c, p) .mp_div(.mp_exp(x, p), c, p), "log($c*$y)", TRUE),
  list(function(x, c, p) .mp_div(c, .mp_exp(x, p), p), "log($c/$y)", FALSE),
  list(function(x, c, p) .id_log(.mp_mul(x, c, p), p), "exp($y)/$c", FALSE),
  list(function(x, c, p) .id_log(.mp_div(x, c, p), p), "$c*exp($y)", TRUE),
  list(function(x, c, p) .id_log(.mp_div(c, x, p), p), "$c/exp($y)", FALSE),
  list(function(x, c, p) .mp_mul(c, .id_log(x, p), p), "exp($y/$c)", TRUE),
  list(function(x, c, p) .mp_div(.id_log(x, p), c, p), "exp($c*$y)", TRUE),
  list(function(x, c, p) .mp_div(c, .id_log(x, p), p), "exp($c/$y)", TRUE))
# ln of a positive mpf; NULL (skipped) otherwise, where mpmath gives a
# complex value whose comparisons fail
.id_log <- function(a, p) if (a$m > 0) .mp_log(a, p) else NULL

# mpmath identify: a string formula for x, or NULL
.id_identify <- function(x, constants, tol, maxcoeff, full, prec) {
  x <- .id_mpf(x, prec)
  if (.id_zero(x)) return(if (full) "0" else "0")
  if (x$m < 0) {
    sol <- .id_identify(.id_neg(x), constants, tol, maxcoeff, full, prec)
    if (is.null(sol)) return(NULL)
    return(paste0("-(", sol, ")"))
  }
  tol <- if (!is.null(tol)) .id_mpf(tol, prec) else .mpf(1, as.integer(floor(0.7 * (1 - prec))))
  M <- maxcoeff
  constants <- lapply(constants, function(cc) list(value = .id_mpf(cc$value, prec), name = cc$name))
  if (!any(vapply(constants, function(cc) cc$value$m == 1 && cc$value$e == 0, NA)))
    constants <- c(list(list(value = .id_one, name = "1")), constants)
  M2 <- .mp_from_q(gmp::as.bigq(M)^2, prec)
  sols <- character()
  for (tr in .ID_TRANSFORMS) {
    for (cc in constants) {
      if (tr[[3L]] && cc$name == "1") next
      t <- tr[[1L]](x, cc$value, prec)
      if (is.null(t)) next
      if (.id_lt(M2, .id_abs(t)) || .id_lt(.id_abs(t), tol)) next
      r <- .id_pslq(c(list(t), lapply(constants, `[[`, "value")), prec, tol, M)
      s <- NULL
      if (!is.null(r) && all(vapply(r, function(u) abs(u) <= M, NA)) && r[[1L]] != 0) {
        s <- .id_pslqstring(r, constants)
      } else {
        q <- .id_pslq(list(.id_one, t, .id_pow(t, 2L, prec)), prec, tol, M)
        if (!is.null(q) && length(q) == 3L && q[[3L]] != 0 && max(abs(q[[1L]]), abs(q[[2L]]), abs(q[[3L]])) <= M)
          s <- .id_quadraticstring(t, q[[1L]], q[[2L]], q[[3L]], prec)
      }
      if (!is.null(s)) {
        ftn <- tr[[2L]]
        s <- if (cc$name == "1" && grepl("/$c", ftn, fixed = TRUE)) sub("/$c", "", gsub("$y", s, ftn, fixed = TRUE), fixed = TRUE)
             else gsub("$c", cc$name, gsub("$y", s, ftn, fixed = TRUE), fixed = TRUE)
        sols <- c(sols, s)
        if (!full) return(s)
      }
    }
  }
  if (!(x$m == 1 && x$e == 0)) {
    ilogs <- c(2, 3, 5, 7)
    lnI <- lapply(ilogs, function(i) .mp_log(.mp_from_q(gmp::as.bigq(i), prec), prec))
    logs <- list()
    for (cc in constants) {
      la <- .id_log(cc$value, prec)
      if (is.null(la)) next
      found <- any(vapply(lnI, function(li) !is.null(.id_findpoly(.mp_div(la, li, prec), 1L, prec)), NA))
      if (!found) logs[[length(logs) + 1L]] <- list(value = la, name = cc$name)
    }
    logs <- c(Map(function(v, i) list(value = v, name = as.character(i)), lnI, ilogs), logs)
    r <- .id_pslq(c(list(.mp_log(x, prec)), lapply(logs, `[[`, "value")), prec, tol, M)
    if (!is.null(r) && all(vapply(r, function(u) abs(u) <= M, NA)) && r[[1L]] != 0) {
      s <- .id_prodstring(r, logs)
      sols <- c(sols, s)
      if (!full) return(s)
    }
  }
  if (full && length(sols)) return(sols[order(nchar(sols))][[1L]])
  NULL
}

# ---------------------------------------------------------------------------
# nsimplify

#' Find simple exact forms of numbers
#'
#' `nsimplify(e)` replaces a number, or the floating point numbers in an
#' expression, by a simple exact form, as SymPy's `nsimplify()`: the value
#' is evaluated to 30 digits and identified (with mpmath's `identify()`) as
#' a rational, a quadratic irrational, or a simple formula in the given
#' `constants`; otherwise its Floats become rationals.
#'
#' @param e A number or expression.
#' @param constants A list of constants (such as `sym_pi`) to try.
#' @param tolerance The tolerance of the identification, or `NULL` for
#'   `10^-d` with `d` the smallest precision (at most 15 digits) of the
#'   Floats in `e`; with `rational = TRUE`, the tolerance of the rational
#'   approximation.
#' @param full Try harder, returning the shortest formula found.
#' @param rational `TRUE` converts the Floats to rationals only; `FALSE`
#'   returns the identified value or `e` itself; `NULL` falls back to
#'   rationals when no formula is found.
#' @param rational_conversion `"base10"` converts a Float through its
#'   decimal digits, `"exact"` through its exact binary value.
#' @return An expression.
#' @examples
#' nsimplify(0.333333333333)
#' nsimplify(1.414213562373095)
#' nsimplify(6.283185307179586, constants = list(sym_pi))
#' nsimplify(sym_pi, tolerance = 0.01)
#' @export
nsimplify <- function(e, constants = list(), tolerance = NULL, full = FALSE, rational = NULL,
                      rational_conversion = "base10") {
  expr <- as_sym(e)
  if (expr$type == "Integer") return(expr)
  floats <- Filter(function(z) z$type == "Float", c(list(expr), .all_subexprs(expr)))
  # a Float holding an infinite double and no multiprecision value; a large
  # finite Float (1e400) overflows the double but keeps its mpf
  infs <- Filter(function(z) is.null(z$extra$mp) && !is.finite(z$value), floats)
  if (length(infs)) expr <- xreplace(expr, lapply(infs, function(z) list(z, if (z$value > 0) S_Infinity else S_NegativeInfinity)))
  if (is_inf_const(expr) && !is_const(expr, "zoo")) return(expr)
  if (isTRUE(rational) || length(free_symbols(expr))) return(.real_to_rational(expr, tolerance, rational_conversion))
  floats <- Filter(function(z) z$type == "Float", c(list(expr), .all_subexprs(expr)))
  user_tol_null <- is.null(tolerance)
  if (is.null(tolerance))
    tolerance <- 10^(-min(c(15, vapply(floats, function(z) .prec_to_dps(.float_prec(z)), 0))))
  prec <- 30L
  bprec <- as.integer(prec * 3.33)
  cdict <- list()
  for (cst in constants) {
    cst <- as_sym(cst)
    v <- N(cst, prec)
    if (v$type != "Float") stop("ValueError: constants must be real-valued", call. = FALSE)
    cdict[[sstr(cst)]] <- list(value = .id_mpf(.float_mp(v), bprec), name = sstr(cst))
  }
  if (length(cdict)) cdict <- cdict[order(names(cdict))]
  exprval <- N(expr, prec, chop = TRUE)
  ri <- .re_im(exprval)
  re <- ri[[1L]]; im <- ri[[2L]]
  if (!(is_number_type(re) && is_number_type(im))) return(expr)
  wp <- .dps_to_prec(prec)
  nsimplify_real <- function(x) {
    xv <- if (x$type == "Float") .id_mpf(.float_mp(x), bprec) else .mp_from_q(x$value, bprec)
    s <- .id_identify(xv, unname(cdict), tolerance, 1000, full, wp)
    if (is.null(s)) stop("nsimplify: not identified", call. = FALSE)
    r <- sympify(s)
    if (!is_zero_num(x) && is_zero_num(r)) stop("nsimplify: identified as zero", call. = FALSE)
    if (isFALSE(ask_fact(r, "finite"))) stop("nsimplify: infinite", call. = FALSE)
    r
  }
  ok <- tryCatch({
    if (!is_zero_num(re)) re <- nsimplify_real(re)
    if (!is_zero_num(im)) im <- nsimplify_real(im)
    TRUE
  }, error = function(err) FALSE)
  if (!ok) {
    # a lone Float that failed with the default tolerance fails again inside
    # _real_to_rational, so that second identification is skipped
    skip <- if (expr$type == "Float" && user_tol_null) expr$key else NULL
    if (is.null(rational)) return(.real_to_rational(expr, rational_conversion = rational_conversion, skip = skip))
    return(expr)
  }
  rv <- .add(list(re, .mul(list(im, S_I))))
  if (!identical(rv$key, expr$key) || isFALSE(rational)) return(rv)
  .real_to_rational(expr, rational_conversion = rational_conversion)
}

# Fraction.limit_denominator: the closest rational with denominator at most
# max_den
.q_limit_denominator <- function(q, max_den) {
  q <- gmp::as.bigq(q); max_den <- gmp::as.bigz(max_den)
  if (gmp::denominator(q) <= max_den) return(q)
  p0 <- gmp::as.bigz(0); q0 <- gmp::as.bigz(1); p1 <- gmp::as.bigz(1); q1 <- gmp::as.bigz(0)
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  repeat {
    a <- n %/% d
    q2 <- q0 + a * q1
    if (q2 > max_den) break
    t <- p0 + a * p1; p0 <- p1; q0 <- q1; p1 <- t; q1 <- q2
    t <- n - a * d; n <- d; d <- t
  }
  k <- (max_den - q0) %/% q1
  b1 <- gmp::as.bigq(p0 + k * p1, q0 + k * q1)
  b2 <- gmp::as.bigq(p1, q1)
  if (abs(b2 - q) <= abs(b1 - q)) b2 else b1
}

# _real_to_rational: every Float of expr replaced by a rational
.real_to_rational <- function(expr, tolerance = NULL, rational_conversion = "base10", skip = NULL) {
  floats <- Filter(function(z) z$type == "Float", c(list(expr), .all_subexprs(expr)))
  floats <- floats[!duplicated(vapply(floats, `[[`, "", "key"))]
  reduce_num <- if (!is.null(tolerance) && tolerance < 1) ceiling(1 / tolerance) else NULL
  reps <- list()
  for (fl in floats) {
    q <- .mp_to_q(.float_mp(fl))
    exact <- function() make_rational(q)
    if (!is.null(reduce_num)) {
      r <- make_rational(.q_limit_denominator(q, reduce_num))
    } else if (!is.null(tolerance) && tolerance >= 1 && gmp::denominator(q) != 1) {
      r <- make_rational(.q_limit_denominator(gmp::as.bigq(tolerance * round(num_value_dbl(fl) / tolerance)), as.integer(tolerance)))
    } else {
      if (rational_conversion == "exact") { reps[[length(reps) + 1L]] <- list(fl, exact()); next }
      if (rational_conversion != "base10") stop("ValueError: rational_conversion must be 'base10' or 'exact'", call. = FALSE)
      r <- if (identical(fl$key, skip)) fl else nsimplify(fl, rational = FALSE)
      if (q != 0 && is_zero_num(r)) r <- exact()
      # the decimal digits of the Float, as Rational(str(fl/d))*d
      else if (!r$type %in% c("Integer", "Rational")) r <- make_rational(.decimal_to_bigq(.fmt_float_sym(fl, strip = TRUE)))
    }
    reps[[length(reps) + 1L]] <- list(fl, r)
  }
  if (!length(reps)) return(expr)
  xreplace(expr, reps)
}
