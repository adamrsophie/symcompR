# Real root isolation and refinement by continued fractions (Vincent,
# Akritas and Strzebonski), after sympy/polys/rootisolation.py. Integer
# polynomials are bigz vectors, highest degree first; Mobius transforms are
# quadruples (a, b, c, d) of bigz for the map x -> (a*x + b)/(c*x + d).

.zp <- function(v) gmp::as.bigz(v)
.zp_strip <- function(f) { nz <- which(f != 0); if (!length(nz)) .zp(integer(0)) else f[nz[1L]:length(f)] }
.zp_tc <- function(f) if (length(f)) f[length(f)] else .zp(0)
.zp_shift <- function(f, a) {
  n <- length(f)
  if (n < 2L || a == 0) return(f)
  for (i in seq_len(n - 1L)) for (j in 2:(n - i + 1L)) f[j] <- f[j] + a * f[j - 1L]
  f
}
.zp_reverse <- function(f) .zp_strip(rev(f))
.zp_mirror <- function(f) {
  n <- length(f) - 1L
  if (n < 0L) return(f)
  odd <- ((n - (seq_along(f) - 1L)) %% 2L) == 1L
  f[odd] <- -f[odd]
  f
}
.zp_scale <- function(f, a) {
  n <- length(f) - 1L
  if (n < 0L) return(f)
  f * a^(n - (seq_along(f) - 1L))
}
.zp_rshift <- function(f, k) if (k <= 0L) f else f[seq_len(max(0L, length(f) - k))]
.zp_sign_var <- function(f) {
  s <- ifelse(f[f != 0] > 0, 1, -1)
  if (length(s) < 2L) 0L else sum(s[-1L] != s[-length(s)])
}
.zp_mul <- function(f, g) {
  if (!length(f) || !length(g)) return(.zp(integer(0)))
  out <- .zp(rep(0, length(f) + length(g) - 1L))
  for (i in seq_along(f)) { idx <- i:(i + length(g) - 1L); out[idx] <- out[idx] + f[i] * g }
  out
}
.zp_add <- function(f, g) {
  n <- max(length(f), length(g))
  f <- c(.zp(rep(0, n - length(f))), f); g <- c(.zp(rep(0, n - length(g))), g)
  .zp_strip(f + g)
}
# dup_transform: q^n * f(p/q) for linear p and q
.zp_transform <- function(f, p, q) {
  if (!length(f)) return(f)
  n <- length(f) - 1L
  h <- f[1L]; Q <- list(.zp(1))
  for (i in seq_len(n)) Q[[i + 1L]] <- .zp_mul(Q[[i]], q)
  for (i in seq_len(n)) h <- .zp_add(.zp_mul(h, p), Q[[i + 1L]] * f[i + 1L])
  h
}

# ZZ.log(a, 2) as int(math.log(a, 2))
.zz_log2 <- function(a) {
  d <- as.double(a)
  if (is.finite(d)) floor(log(d) / log(2)) else .bitlen(a) - 1
}

.zp_upper_bound <- function(f) {
  n <- length(f); P <- c(); t <- rep(1, n)
  if (f[1L] < 0) f <- -f
  f <- rev(f)
  for (i in seq_len(n)) {
    if (f[i] >= 0) next
    a <- .zz_log2(-f[i]); QL <- list()
    if (i < n) for (j in (i + 1L):n) {
      if (f[j] <= 0) next
      q <- t[j] + a - .zz_log2(f[j])
      QL[[length(QL) + 1L]] <- c(q %/% (j - i), j)
    }
    if (!length(QL)) next
    # min over [q, j] pairs, lexicographic
    o <- order(vapply(QL, `[`, 0, 1L), vapply(QL, `[`, 0, 2L))[1L]
    q <- QL[[o]]
    t[q[2L]] <- t[q[2L]] + 1
    P <- c(P, q[1L])
  }
  if (!length(P)) return(NULL)
  k <- max(P) + 1
  if (k >= 0) gmp::as.bigq(gmp::as.bigz(2)^k) else gmp::as.bigq(1, gmp::as.bigz(2)^(-k))
}
.zp_lower_bound <- function(f) {
  b <- .zp_upper_bound(.zp_reverse(f))
  if (is.null(b)) NULL else 1 / b
}
.int_floor_q <- function(q) { q <- gmp::as.bigq(q); gmp::as.bigz(gmp::numerator(q) %/% gmp::denominator(q)) }

.mob_to_interval <- function(M) {
  s <- gmp::as.bigq(M[[1L]], M[[3L]]); t <- gmp::as.bigq(M[[2L]], M[[4L]])
  if (s <= t) list(s, t) else list(t, s)
}
.mob_from_interval <- function(s, t) {
  s <- gmp::as.bigq(s); t <- gmp::as.bigq(t)
  list(gmp::numerator(s), gmp::numerator(t), gmp::denominator(s), gmp::denominator(t))
}
.mob_width <- function(M) abs(gmp::as.bigq(M[[1L]], M[[3L]]) - gmp::as.bigq(M[[2L]], M[[4L]]))

# dup_step_refine_real_root
.zp_step_refine <- function(f, M, fast = FALSE) {
  a <- M[[1L]]; b <- M[[2L]]; c <- M[[3L]]; d <- M[[4L]]
  if (a == b && c == d) return(list(f = f, M = M))
  A <- .zp_lower_bound(f)
  A <- if (!is.null(A)) .int_floor_q(A) else .zp(0)
  if (fast && A > 16) { f <- .zp_scale(f, A); a <- A * a; c <- A * c; A <- .zp(1) }
  if (A >= 1) {
    f <- .zp_shift(f, A); b <- A * a + b; d <- A * c + d
    if (.zp_tc(f) == 0) return(list(f = f, M = list(b, b, d, d)))
  }
  g <- f; f <- .zp_shift(f, 1)
  a1 <- a; b1 <- a + b; c1 <- c; d1 <- c + d
  if (.zp_tc(f) == 0) return(list(f = f, M = list(b1, b1, d1, d1)))
  if (.zp_sign_var(f) == 1L) {
    M <- list(a1, b1, c1, d1)
  } else {
    f <- .zp_shift(.zp_reverse(g), 1)
    if (.zp_tc(f) == 0) f <- .zp_rshift(f, 1L)
    M <- list(b, a + b, d, c + d)
  }
  list(f = f, M = M)
}

# dup_inner_refine_real_root (mobius = TRUE form)
.zp_inner_refine <- function(f, M, eps = NULL, steps = NULL, disjoint = NULL, fast = FALSE) {
  while (M[[3L]] == 0) { r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M }
  if (!is.null(eps) && !is.null(steps)) {
    for (i in seq_len(steps)) {
      if (.mob_width(M) >= eps) { r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M } else break
    }
  } else {
    if (!is.null(eps)) while (.mob_width(M) >= eps) { r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M }
    if (!is.null(steps)) for (i in seq_len(steps)) { r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M }
  }
  if (!is.null(disjoint)) repeat {
    uv <- .mob_to_interval(M)
    if (uv[[2L]] <= disjoint || disjoint <= uv[[1L]]) break
    r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M
  }
  list(f = f, M = M)
}

# dup_refine_real_root over ZZ
.zp_refine <- function(f, s, t, eps = NULL, steps = NULL, disjoint = NULL, fast = FALSE) {
  s <- gmp::as.bigq(s); t <- gmp::as.bigq(t)
  if (s == t) return(list(s, t))
  if (s > t) { tmp <- s; s <- t; t <- tmp }
  negative <- FALSE
  if (s < 0) {
    if (t <= 0) { f <- .zp_mirror(f); tmp <- -t; t <- -s; s <- tmp; negative <- TRUE }
    else stop("Cannot refine a real root in (", as.character(s), ", ", as.character(t), ")", call. = FALSE)
  }
  if (negative && !is.null(disjoint)) disjoint <- if (disjoint < 0) -disjoint else NULL
  M <- .mob_from_interval(s, t)
  f <- .zp_transform(f, .zp_strip(c(M[[1L]], M[[2L]])), .zp_strip(c(M[[3L]], M[[4L]])))
  if (.zp_sign_var(f) != 1L)
    stop("there should be exactly one root in (", as.character(s), ", ", as.character(t), ") interval", call. = FALSE)
  r <- .zp_inner_refine(f, M, eps, steps, disjoint, fast)
  uv <- .mob_to_interval(r$M)
  if (negative) list(-uv[[2L]], -uv[[1L]]) else uv
}

# dup_inner_isolate_real_roots: list of list(f, M) for the positive roots
.zp_inner_isolate <- function(f, eps = NULL, fast = FALSE) {
  one <- .zp(1); zero <- .zp(0)
  k <- .zp_sign_var(f)
  if (k == 0L) return(list())
  if (k == 1L) return(list(.zp_inner_refine(f, list(one, zero, zero, one), eps = eps, fast = fast)))
  roots <- list(); stack <- list(list(one, zero, zero, one, f, k))
  while (length(stack)) {
    top <- stack[[length(stack)]]; stack <- stack[-length(stack)]
    a <- top[[1L]]; b <- top[[2L]]; c <- top[[3L]]; d <- top[[4L]]; f <- top[[5L]]; k <- top[[6L]]
    A <- .zp_lower_bound(f)
    A <- if (!is.null(A)) .int_floor_q(A) else zero
    if (fast && A > 16) { f <- .zp_scale(f, A); a <- A * a; c <- A * c; A <- one }
    if (A >= 1) {
      f <- .zp_shift(f, A); b <- A * a + b; d <- A * c + d
      if (.zp_tc(f) == 0) { roots[[length(roots) + 1L]] <- list(f = f, M = list(b, b, d, d)); f <- .zp_rshift(f, 1L) }
      k <- .zp_sign_var(f)
      if (k == 0L) next
      if (k == 1L) { roots[[length(roots) + 1L]] <- .zp_inner_refine(f, list(a, b, c, d), eps = eps, fast = fast); next }
    }
    f1 <- .zp_shift(f, 1)
    a1 <- a; b1 <- a + b; c1 <- c; d1 <- c + d; r <- 0L
    if (.zp_tc(f1) == 0) { roots[[length(roots) + 1L]] <- list(f = f1, M = list(b1, b1, d1, d1)); f1 <- .zp_rshift(f1, 1L); r <- 1L }
    k1 <- .zp_sign_var(f1)
    k2 <- k - k1 - r
    a2 <- b; b2 <- a + b; c2 <- d; d2 <- c + d
    f2 <- NULL
    if (k2 > 1L) {
      f2 <- .zp_shift(.zp_reverse(f), 1)
      if (.zp_tc(f2) == 0) f2 <- .zp_rshift(f2, 1L)
      k2 <- .zp_sign_var(f2)
    }
    if (k1 < k2) {
      tmp <- list(a1, b1, c1, d1, f1, k1)
      a1 <- a2; b1 <- b2; c1 <- c2; d1 <- d2; f1 <- f2; k1 <- k2
      a2 <- tmp[[1L]]; b2 <- tmp[[2L]]; c2 <- tmp[[3L]]; d2 <- tmp[[4L]]; f2 <- tmp[[5L]]; k2 <- tmp[[6L]]
    }
    if (!k1) next
    if (is.null(f1)) { f1 <- .zp_shift(.zp_reverse(f), 1); if (.zp_tc(f1) == 0) f1 <- .zp_rshift(f1, 1L) }
    if (k1 == 1L) roots[[length(roots) + 1L]] <- .zp_inner_refine(f1, list(a1, b1, c1, d1), eps = eps, fast = fast)
    else stack[[length(stack) + 1L]] <- list(a1, b1, c1, d1, f1, k1)
    if (!k2) next
    if (is.null(f2)) { f2 <- .zp_shift(.zp_reverse(f), 1); if (.zp_tc(f2) == 0) f2 <- .zp_rshift(f2, 1L) }
    if (k2 == 1L) roots[[length(roots) + 1L]] <- .zp_inner_refine(f2, list(a2, b2, c2, d2), eps = eps, fast = fast)
    else stack[[length(stack) + 1L]] <- list(a2, b2, c2, d2, f2, k2)
  }
  roots
}

.zp_discard_outside <- function(f, M, inf, sup, negative, fast) {
  repeat {
    uv <- .mob_to_interval(M); u <- uv[[1L]]; v <- uv[[2L]]
    if (negative) { tmp <- -v; v <- -u; u <- tmp }
    if ((is.null(inf) || u >= inf) && (is.null(sup) || v <= sup)) return(list(f = f, M = M))
    if ((!is.null(sup) && u > sup) || (!is.null(inf) && v < inf)) return(NULL)
    r <- .zp_step_refine(f, M, fast); f <- r$f; M <- r$M
  }
}

# positive (sgn = 1) or negative (sgn = -1) roots as list(f, M, neg)
.zp_isolate_side <- function(f, sgn, eps = NULL, inf = NULL, sup = NULL, fast = FALSE) {
  if (sgn > 0 && !is.null(sup) && sup < 0) return(list())
  if (sgn < 0 && !is.null(inf) && inf >= 0) return(list())
  roots <- .zp_inner_isolate(if (sgn > 0) f else .zp_mirror(f), eps, fast)
  if (!is.null(inf) || !is.null(sup))
    roots <- Filter(Negate(is.null), lapply(roots, function(r) .zp_discard_outside(r$f, r$M, inf, sup, sgn < 0, fast)))
  roots
}

.mob_disjoint <- function(M, N, strict = FALSE) {
  a1 <- M[[1L]]; b1 <- M[[2L]]; c1 <- M[[3L]]; d1 <- M[[4L]]
  a2 <- N[[1L]]; b2 <- N[[2L]]; c2 <- N[[3L]]; d2 <- N[[4L]]
  if (a1 * d1 == b1 * c1 && a2 * d2 == b2 * c2) return(TRUE)
  if (a1 * d1 > b1 * c1) { t1 <- a1; a1 <- b1; b1 <- t1; t1 <- c1; c1 <- d1; d1 <- t1 }
  if (a2 * d2 > b2 * c2) { t1 <- a2; a2 <- b2; b2 <- t1; t1 <- c2; c2 <- d2; d2 <- t1 }
  if (!strict) a2 * d1 >= c2 * b1 || b2 * c1 <= d2 * a1 else a2 * d1 > c2 * b1 || b2 * c1 < d2 * a1
}

# dup_isolate_real_roots: list of list(s, t, k) sorted, for f over ZZ
.zp_isolate_real_roots <- function(f, eps = NULL, inf = NULL, sup = NULL, fast = FALSE, sqf = FALSE) {
  f <- .zp_strip(f)
  if (length(f) <= 1L) return(list())
  out <- list()
  # _isolate_zero
  j <- 0L
  while (length(f) > 1L && .zp_tc(f) == 0) { f <- .zp_rshift(f, 1L); j <- j + 1L }
  if (j > 0L && (is.null(inf) || inf <= 0) && (is.null(sup) || 0 <= sup))
    out[[1L]] <- list(gmp::as.bigq(0), gmp::as.bigq(0), j)
  facs <- if (sqf) list(list(f, 1L)) else {
    x <- Dummy("x")
    fl <- sqf_list(.add(lapply(seq_along(f), function(i) .mul(list(.int(f[i]), .pow(x, .int(length(f) - i)))))), x)$factors
    lapply(fl, function(fk) list(.zp(vapply(.coeff_list(expand(fk[[1L]]), x), function(c) as.character(c$value), "")), as.integer(fk[[2L]])))
  }
  if (length(facs) == 1L) {
    k <- facs[[1L]][[2L]]; g <- facs[[1L]][[1L]]
    for (sg in c(-1, 1)) for (r in .zp_isolate_side(g, sg, eps, inf, sup, fast)) {
      uv <- .mob_to_interval(r$M)
      out[[length(out) + 1L]] <- if (sg < 0) list(-uv[[2L]], -uv[[1L]], k) else list(uv[[1L]], uv[[2L]], k)
    }
  } else {
    # _real_isolate_and_disjoin
    for (sg in c(-1, 1)) {
      I <- list()
      for (fk in facs) for (r in .zp_isolate_side(fk[[1L]], sg, eps, inf, sup, fast)) I[[length(I) + 1L]] <- list(f = r$f, M = r$M, k = fk[[2L]])
      for (i in seq_along(I)) {
        if (i < length(I)) for (jj in (i + 1L):length(I)) {
          while (!.mob_disjoint(I[[i]]$M, I[[jj]]$M)) {
            r1 <- .zp_inner_refine(I[[i]]$f, I[[i]]$M, steps = 1L, fast = fast); I[[i]]$f <- r1$f; I[[i]]$M <- r1$M
            r2 <- .zp_inner_refine(I[[jj]]$f, I[[jj]]$M, steps = 1L, fast = fast); I[[jj]]$f <- r2$f; I[[jj]]$M <- r2$M
          }
        }
      }
      for (r in I) {
        uv <- .mob_to_interval(r$M)
        out[[length(out) + 1L]] <- if (sg < 0) list(-uv[[2L]], -uv[[1L]], r$k) else list(uv[[1L]], uv[[2L]], r$k)
      }
    }
  }
  # sorted as tuples ((s, t), k)
  if (length(out) > 1L) {
    o <- order(vapply(out, function(z) as.double(z[[1L]]), 0), vapply(out, function(z) as.double(z[[2L]]), 0))
    out <- out[o]
  }
  out
}

# An expression as a primitive integer polynomial in x (dup_clear_denoms).
.zp_from_expr <- function(f, x) {
  cl <- .coeff_list(expand(as_sym(f)), x)
  q <- gmp::as.bigq(vapply(cl, function(c) {
    if (!(c$type %in% c("Integer", "Rational"))) stop("real root isolation needs rational coefficients", call. = FALSE)
    as.character(c$value)
  }, ""))
  den <- Reduce(gmp::lcm.bigz, lapply(seq_along(q), function(i) gmp::denominator(q[i])), gmp::as.bigz(1))
  .zp(gmp::numerator(q * den))
}
