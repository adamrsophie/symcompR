# Integration of rational functions, after sympy/integrals/rationaltools.py
# (ratint, ratint_ratpart, ratint_logpart, log_to_atan, log_to_real).

# Poly(p, x).cancel(Poly(q, x), include = False) over the field domain K:
# list(coeff, p, q) as dpolys over K
.dp_cancel <- function(F, G) {
  K0 <- F$K
  if (.dk_field(K0) && !(.dk(K0) %in% c("EX", "RR", "CC", "GF", "AlgField"))) {
    Kr <- .dom_get_ring(K0)
    cd <- function(P) {
      if (.dk(K0) == "QQ_I") {
        dens <- lapply(P$c, function(z) { ri <- .re_im(z); .dom_lcm(ZZ, .zs(gmp::denominator(gmp::as.bigq(ri[[1L]]$value))), .zs(gmp::denominator(gmp::as.bigq(ri[[2L]]$value)))) })
        common <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), dens, S_One)
      } else if (.dk(K0) == "QQ") {
        dens <- lapply(P$c, function(z) .zs(gmp::denominator(gmp::as.bigq(z$value))))
        common <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), dens, S_One)
      } else {
        dens <- lapply(P$c, function(z) expand(.as_numer_denom(z)[[2L]]))
        common <- Reduce(function(a, b) .mpoly_lcm(a, b, .dk_gens(K0)), dens, S_One)
      }
      list(common, .dp_convert(.dp_mul_ground(P, common), Kr))
    }
    a <- cd(F); b <- cd(G)
    cq <- a[[1L]]; f <- a[[2L]]; cp <- b[[1L]]; g <- b[[2L]]
    ig <- .dp_inner_gcd(f, g)
    p <- ig[[2L]]; q <- ig[[3L]]
    h <- .dom_gcd(Kr, cp, cq)
    cp <- .dom_quo(Kr, cp, h); cq <- .dom_quo(Kr, cq, h)
    p <- .dp_convert(p, K0); q <- .dp_convert(q, K0)
  } else {
    cp <- S_One; cq <- S_One
    ig <- .dp_inner_gcd(F, G)
    p <- ig[[2L]]; q <- ig[[3L]]
  }
  pn <- .dom_is_negative(K0, .dp_LC(p)); qn <- .dom_is_negative(K0, .dp_LC(q))
  if (pn && qn) { p <- .dp_neg(p); q <- .dp_neg(q) }
  else if (pn) { cp <- .neg(cp); p <- .dp_neg(p) }
  else if (qn) { cp <- .neg(cp); q <- .dp_neg(q) }
  list(.div(cp, cq), p, q)
}

# Poly.integrate
.dp_integrate <- function(f) {
  n <- length(f$c)
  if (!n) return(f)
  .dp_set(f, c(lapply(seq_len(n), function(i) .dn(f$K, .div(f$c[[i]], .int(n - i + 1L)))), list(S_Zero)))
}

# the univariate domain of several expressions in x (Poly._unify of each)
.dp_joint_domain <- function(es, x, field = FALSE) {
  cl <- unlist(lapply(es, function(e) .coeff_list(expand(e), x)), recursive = FALSE)
  .dp_domain(if (length(cl)) cl else list(S_Zero), field)
}

#' Integration of rational functions
#'
#' A port of SymPy's `ratint()` (sympy/integrals/rationaltools.py): the
#' Hermite reduction (`ratint_ratpart()`), the Lazard, Rioboo and Trager
#' logarithmic part (`ratint_logpart()`) and its conversion to real
#' functions (`log_to_real()`, `log_to_atan()`).
#'
#' @param f A rational function, or `list(p, q)`.
#' @param x The variable.
#' @param symbol The name of the variable of the roots in `RootSum` terms.
#' @param real `TRUE` to convert to real functions, `FALSE` for `RootSum`
#'   terms, `NULL` to decide from the atoms of `f`.
#' @return The antiderivative.
#' @examples
#' x <- Symbol("x")
#' ratint(1/(x^2 + 1), x)
#' ratint(36/(x^5 - 2*x^4 - 2*x^3 + 4*x^2 + x - 2), x)
#' @export
ratint <- function(f, x, symbol = "t", real = NULL) {
  x <- as_sym(x)
  if (is.list(f) && !is_sym(f)) { p <- as_sym(f[[1L]]); q <- as_sym(f[[2L]]) }
  else { f <- as_sym(f); nd <- .as_numer_denom(f); p <- nd[[1L]]; q <- nd[[2L]] }
  K <- .dp_joint_domain(list(p, q), x, field = TRUE)
  P <- .dp_from(p, x, K); Q <- .dp_from(q, x, K)
  cc <- .dp_cancel(P, Q)
  coeff <- cc[[1L]]; P <- cc[[2L]]; Q <- cc[[3L]]
  qr <- .dp_div(P, Q)
  result <- .dp_expr(.dp_integrate(qr[[1L]]))
  P <- qr[[2L]]
  if (.dp_zero(P)) return(.mul(list(coeff, result)))
  gh <- ratint_ratpart(.dp_expr(P), .dp_expr(Q), x)
  g <- gh[[1L]]; h <- gh[[2L]]
  nd <- .as_numer_denom(h)
  Pn <- .dp_from(nd[[1L]], x); Qn <- .dp_from(nd[[2L]], x)
  Kj <- .dp_joint_domain(list(nd[[1L]], nd[[2L]]), x)
  Pn <- .dp_convert(Pn, Kj); Qn <- .dp_convert(Qn, Kj)
  qr <- .dp_div(.dp_to_field(Pn), .dp_to_field(Qn))
  r <- .dp_to_ring(qr[[2L]])
  result <- .add(list(result, g, .dp_expr(.dp_integrate(qr[[1L]]))))
  if (!.dp_zero(r)) {
    t <- Dummy(if (is.character(symbol)) symbol else symbol$name)
    L <- .ratint_logpart_dp(r, Qn, x, t)
    if (is.null(real)) {
      ats <- if (is.list(f) && !is_sym(f)) c(.ri_atoms(p), .ri_atoms(q)) else .ri_atoms(f)
      ats <- Filter(function(a) !identical(a$key, x$key), ats)
      real <- all(vapply(ats, function(a) isTRUE(ask_fact(a, "extended_real")), NA))
    }
    eps <- list()
    for (hq in L) {
      hh <- .dp_primitive(hq[[1L]])[[2L]]; qq <- hq[[2L]]
      R <- if (real) log_to_real(hh, qq, x, t) else NULL
      if (is.null(R)) R <- RootSum(.dp_expr(qq), Lambda(t, .mul(list(t, fcall("log", list(.dp_expr(hh)))))), t, quadratic = TRUE)
      eps[[length(eps) + 1L]] <- R
    }
    result <- .add(c(list(result), eps))
  }
  .mul(list(coeff, result))
}

# f.atoms(): symbols and numbers
.ri_atoms <- function(e) {
  out <- list()
  rec <- function(z) {
    if (!length(z$args)) { out[[z$key]] <<- z; return(invisible()) }
    for (a in z$args) rec(a)
  }
  rec(e)
  unname(out)
}

#' @rdname ratint
#' @param g The denominator.
#' @export
ratint_ratpart <- function(f, g, x) {
  x <- as_sym(x)
  K <- .dp_joint_domain(list(f, g), x, field = TRUE)
  F <- .dp_from(f, x, K); G <- .dp_from(g, x, K)
  ig <- .dp_inner_gcd(G, .dp_diff(G))
  u <- ig[[1L]]; v <- ig[[2L]]
  n <- .dp_deg(u); m <- .dp_deg(v)
  A_c <- lapply(seq_len(max(n, 0)) - 1L, function(i) Dummy(paste0("a", n - i)))
  B_c <- lapply(seq_len(max(m, 0)) - 1L, function(i) Dummy(paste0("b", m - i)))
  Aexpr <- .add(lapply(seq_along(A_c), function(i) .mul(list(A_c[[i]], .pow(x, .int(n - i))))))
  Bexpr <- .add(lapply(seq_along(B_c), function(i) .mul(list(B_c[[i]], .pow(x, .int(m - i))))))
  w <- .dp_quo(.dp_mul(.dp_diff(u), v), u)
  H <- expand(.add(list(.dp_expr(F), .neg(.mul(list(.diff1(Aexpr, x), .dp_expr(v)))), .mul(list(Aexpr, .dp_expr(w))),
                        .neg(.mul(list(Bexpr, .dp_expr(u)))))))
  eqs <- Filter(function(e) !is_zero_num(e), .coeff_list(H, x))
  unk <- c(A_c, B_c)
  sol <- if (length(unk)) .ri_solve_linear(eqs, unk) else list()
  A <- subs(Aexpr, lapply(seq_along(unk), function(i) list(unk[[i]], sol[[i]])))
  B <- subs(Bexpr, lapply(seq_along(unk), function(i) list(unk[[i]], sol[[i]])))
  list(.cancel_in(.div(A, .dp_expr(u)), x), .cancel_in(.div(B, .dp_expr(v)), x))
}

# the unique solution of a linear system (Gauss-Jordan elimination)
.ri_solve_linear <- function(eqs, unk) {
  n <- length(unk)
  rows <- lapply(eqs, function(e) {
    cf <- lapply(unk, function(u) cancel(.diff1(e, u)))
    rest <- cancel(.neg(subs(e, lapply(unk, function(u) list(u, S_Zero)))))
    c(cf, list(rest))
  })
  sol <- vector("list", n); used <- rep(FALSE, length(rows))
  piv <- integer(0)
  for (j in seq_len(n)) {
    k <- which(!used & vapply(rows, function(r) !is_zero_num(r[[j]]), NA))[1L]
    if (is.na(k)) next
    used[k] <- TRUE; piv[j] <- k
    pr <- rows[[k]]; pv <- pr[[j]]
    pr <- lapply(pr, function(z) cancel(.div(z, pv)))
    rows[[k]] <- pr
    for (i in seq_along(rows)) if (i != k && !is_zero_num(rows[[i]][[j]])) {
      fct <- rows[[i]][[j]]
      rows[[i]] <- Map(function(a, b) cancel(.sub(a, .mul(list(fct, b)))), rows[[i]], pr)
    }
  }
  for (j in seq_len(n)) sol[[j]] <- if (length(piv) >= j && !is.na(piv[j])) rows[[piv[j]]][[n + 1L]] else S_Zero
  sol
}

# cancel(f, x): SymPy's cancel with x as the only generator
.cancel_in <- function(f, x) {
  f <- signsimp(as_sym(f))
  if (is_number_type(f) || !.is_expr_node(f)) return(f)
  f <- factor_terms(f, radical = TRUE)
  nd <- .as_numer_denom(f)
  p <- expand(nd[[1L]]); q <- expand(nd[[2L]])
  if (!.has1(p, x) && !.has1(q, x) && FALSE) return(expand(f))
  K <- .dp_joint_domain(list(p, q), x)
  F <- .dp_from(p, x, K); G <- .dp_from(q, x, K)
  if (.dk_field(K) && !(.dk(K) %in% c("EX", "RR", "CC", "GF", "AlgField"))) {
    cc <- .dp_cancel(F, G)
    P <- .dp_mul_ground(cc[[2L]], .as_numer_denom(cc[[1L]])[[1L]])
    Q <- .dp_mul_ground(cc[[3L]], .as_numer_denom(cc[[1L]])[[2L]])
    # PolyElement.cancel multiplies the cofactors of the denominators back
  } else {
    ig <- .dp_inner_gcd(F, G)
    P <- ig[[2L]]; Q <- ig[[3L]]
    if (.dom_is_negative(K, .dp_LC(Q))) { P <- .dp_neg(P); Q <- .dp_neg(Q) }
  }
  .div(.dp_expr(P), .dp_expr(Q))
}

#' @rdname ratint
#' @param t The variable of the roots.
#' @export
ratint_logpart <- function(f, g, x, t = NULL) {
  x <- as_sym(x)
  if (is.null(t)) t <- Dummy("t")
  K <- .dp_joint_domain(list(f, g), x)
  L <- .ratint_logpart_dp(.dp_from(f, x, K), .dp_from(g, x, K), x, t)
  lapply(L, function(hq) list(.dp_expr(hq[[1L]]), .dp_expr(hq[[2L]])))
}

# ratint_logpart on dpolys: list(list(h (dpoly in x), q (dpoly in t)))
.ratint_logpart_dp <- function(f, g, x, t) {
  fe <- .dp_expr(f); ge <- .dp_expr(g)
  be <- expand(.sub(fe, .mul(list(.diff1(ge, x), t))))
  Kb <- .dp_joint_domain(list(fe, ge, be), x)
  a <- .dp_convert(.dp_from(ge, x, Kb), Kb); b <- .dp_from(be, x, Kb)
  pr <- .dp_prs_resultant(a, b)
  res <- expand(pr[[1L]]); R <- pr[[2L]]
  Rmap <- list()
  for (r in R) Rmap[[as.character(.dp_deg(r))]] <- r
  include_sign <- function(c, K, sqf) {
    if (isTRUE(ask_fact(c, "extended_real")) && .is_true(Lt(c, 0))) {
      h <- sqf[[1L]][[1L]]
      sqf[[1L]][[1L]] <- .dp_mul_ground(h, c)
    }
    sqf
  }
  resp <- .dp_from(res, t)
  sl <- .dp_sqf_list(resp)
  C <- sl[[1L]]; res_sqf <- include_sign(sl[[1L]], resp$K, sl[[2L]])
  H <- list()
  for (qi in res_sqf) {
    q <- .dp_primitive(qi[[1L]])[[2L]]; i <- qi[[2L]]
    if (.dp_deg(g) == i) { H[[length(H) + 1L]] <- list(g, q); next }
    h <- Rmap[[as.character(i)]]
    hlc <- .dp_from(.dp_LC(h), t, field = TRUE)
    hs <- .dp_sqf_list(hlc, all = TRUE)
    hs2 <- include_sign(hs[[1L]], hlc$K, hs[[2L]])
    qf <- .dp_to_field(.dp_convert(q, .dp_joint_domain(list(.dp_expr(q), .dp_expr(hlc)), t, field = TRUE)))
    for (aj in hs2) {
      aF <- .dp_convert(aj[[1L]], qf$K)
      gq <- .dp_gcd(aF, .dp_convert(qf, qf$K))
      dv <- .dp_expr(.dp_pow(gq, aj[[2L]]))
      h <- .dp_set(h, lapply(h$c, function(z) expand(cancel(.div(z, dv)))))
    }
    hlcF <- .dp_convert(hlc, qf$K)
    inv <- .dp_invert(hlcF, qf)
    hc <- Filter(function(z) !is_zero_num(z), h$c)
    hm <- which(vapply(h$c, function(z) !is_zero_num(z), NA))
    coeffs <- list(S_One)
    for (cf in hc[-1L]) {
      cp <- .dp_from(cf, t, qf$K)
      T <- .dp_rem(.dp_mul(inv, cp), qf)
      coeffs[[length(coeffs) + 1L]] <- .dp_expr(T)
    }
    deg <- length(h$c) - 1L
    hexpr <- .add(lapply(seq_along(hm), function(k) .mul(list(coeffs[[k]], .pow(x, .int(deg - hm[k] + 1L))))))
    H[[length(H) + 1L]] <- list(.dp_from(hexpr, x), q)
  }
  H
}

#' @rdname ratint
#' @export
log_to_atan <- function(f, g) {
  if (!inherits(f, "dpoly")) stop("log_to_atan() takes dpolys")
  if (.dp_deg(f) < .dp_deg(g)) { tmp <- .dp_neg(g); g <- f; f <- tmp }
  f <- .dp_to_field(f); g <- .dp_to_field(g)
  K <- .dp_joint_domain(list(.dp_expr(f), .dp_expr(g)), f$x, field = TRUE)
  f <- .dp_convert(f, K); g <- .dp_convert(g, K)
  pq <- .dp_div(f, g)
  if (.dp_zero(pq[[2L]])) return(.mul(list(S_Two, fcall("atan", list(.dp_expr(pq[[1L]]))))))
  st <- .dp_gcdex(g, .dp_neg(f))
  s <- st[[1L]]; t <- st[[2L]]; h <- st[[3L]]
  u <- .dp_quo(.dp_add(.dp_mul(f, s), .dp_mul(g, t)), h)
  .add(list(.mul(list(S_Two, fcall("atan", list(.dp_expr(u))))), log_to_atan(s, t)))
}

# _roots_real_complex
.roots_real_complex <- function(q, t) {
  rs <- roots(.dp_expr(q), t)
  mult <- attr(rs, "multiplicity")
  rl <- unclass(rs)[seq_along(rs)]
  if (sum(mult) != .dp_deg(q)) return(NULL)
  keys <- vapply(rl, `[[`, "", "key")
  reals <- list(); complexes <- list()
  remaining <- seq_along(rl)
  while (length(remaining)) {
    i <- remaining[length(remaining)]; remaining <- remaining[-length(remaining)]
    r <- rl[[i]]
    rc <- fcall("conjugate", list(r))
    j <- match(rc$key, keys)
    if (!identical(r$key, rc$key) && !is.na(j)) {
      remaining <- setdiff(remaining, j)
      ri <- .as_real_imag(r)
      o <- .bool_print_order(list(ri[[2L]], .neg(ri[[2L]])))
      complexes[[length(complexes) + 1L]] <- list(ri[[1L]], o[[2L]], mult[i])
    } else reals[[length(reals) + 1L]] <- list(r, mult[i])
  }
  list(reals, complexes)
}

#' @rdname ratint
#' @param h,q Polynomials (dpolys) in `x` and `t`.
#' @export
log_to_real <- function(h, q, x, t) {
  rs <- .roots_real_complex(q, t)
  if (is.null(rs)) return(NULL)
  reals <- rs[[1L]]; complexes <- rs[[2L]]
  u <- Dummy("u"); v <- Dummy("v")
  H <- expand(subs(.dp_expr(h), t, .add(list(u, .mul(list(S_I, v))))))
  cm <- collect(H, S_I, evaluate = FALSE)
  a <- .ri_collect_get(cm, S_One); b <- .ri_collect_get(cm, S_I)
  result <- list()
  for (cp in complexes) {
    ru <- cp[[1L]]; rv <- cp[[2L]]
    A <- .dp_from(subs(a, list(list(u, ru), list(v, rv))), x)
    B <- .dp_from(subs(b, list(list(u, ru), list(v, rv))), x)
    Kab <- .dp_joint_domain(list(.dp_expr(A), .dp_expr(B)), x)
    A2 <- .dp_convert(A, Kab); B2 <- .dp_convert(B, Kab)
    AB <- .dp_expr(.dp_add(.dp_mul(A2, A2), .dp_mul(B2, B2)))
    result[[length(result) + 1L]] <- .add(list(.mul(list(ru, fcall("log", list(AB)))), .mul(list(rv, log_to_atan(A, B)))))
  }
  for (rr in reals) {
    r <- rr[[1L]]
    result[[length(result) + 1L]] <- .mul(list(r, fcall("log", list(subs(.dp_expr(h), t, r)))))
  }
  .add(result)
}

.ri_collect_get <- function(cm, k) { v <- cm[[.str(k)]]; if (is.null(v)) S_Zero else v }
