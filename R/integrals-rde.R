# The Risch differential equation, after sympy/integrals/rde.py.

order_at <- function(a, p, t) {
  if (P_is_zero(a)) return(Inf)
  if (identical(P_expr(p)$key, as_sym(t)$key) && P_deg(p) == 1L) {
    # a.as_poly(t).ET()[0][0]: the exponent of the last nonzero term
    ms <- P_monoms(.P(P_expr(a), t))
    return(ms[length(ms)])
  }
  power_list <- list()
  p1 <- p
  r <- P_rem(a, p1)
  tracks <- 1
  while (P_is_zero(r)) {
    power_list[[length(power_list) + 1L]] <- list(p1, tracks)
    p1 <- P_mul(p1, p1)
    tracks <- tracks * 2
    r <- P_rem(a, p1)
  }
  n <- 0
  product <- .P(S_One, t)
  while (length(power_list)) {
    final <- power_list[[length(power_list)]]; power_list[[length(power_list)]] <- NULL
    productf <- P_mul(product, final[[1L]])
    r <- P_rem(a, productf)
    if (P_is_zero(r)) { n <- n + final[[2L]]; product <- productf }
  }
  n
}

order_at_oo <- function(a, d, t) {
  if (P_is_zero(a)) return(Inf)
  P_deg(d) - P_deg(a)
}

weak_normalizer <- function(a, d, DE, z = NULL) {
  if (is.null(z)) z <- Dummy("z")
  sp <- splitfactor(d, DE)
  dn <- sp[[1L]]
  g <- P_gcd(dn, P_diff(dn))
  d_sqf_part <- P_quo(dn, g)
  d1 <- P_quo(d_sqf_part, P_gcd(d_sqf_part, g))
  ab <- gcdex_diophantine(.P(P_expr(P_quo(d, d1)), DE$t), .P(P_expr(d1), DE$t), .P(P_expr(a), DE$t))
  r <- P_resultant(.P(P_expr(P_sub(a, P_mul(.P(z, DE$t), derivation(d1, DE)))), DE$t), .P(P_expr(d1), DE$t))
  r <- .P(r, z)
  if (!.has1(P_expr(r), z)) return(list(.P(S_One, DE$t), list(a, d)))
  rr <- real_roots(P_expr(r), z)
  N <- Filter(function(i) i$type == "Integer" && i$value > 0, rr)
  q <- .P(S_One, DE$t)
  for (n in N) q <- P_mul(q, P_gcd(P_sub(a, P_mul(.P(n, DE$t), derivation(d1, DE))), d1))
  dq <- derivation(q, DE)
  sn <- P_sub(P_mul(q, a), P_mul(d, dq))
  sd <- P_mul(q, d)
  cc <- P_cancel(sn, sd, include = TRUE)
  list(q, list(cc[[1L]], cc[[2L]]))
}

normal_denom <- function(fa, fd, ga, gd, DE) {
  dn <- splitfactor(fd, DE)[[1L]]
  en <- splitfactor(gd, DE)[[1L]]
  p <- P_gcd(dn, en)
  h <- P_quo(P_gcd(en, P_diff(en)), P_gcd(p, P_diff(p)))
  a <- P_mul(dn, h)
  c <- P_mul(a, h)
  if (!P_is_zero(P_div(c, en)[[2L]])) .nonelem("en does not divide dn*h**2 in normal_denom().")
  ca <- P_mul(c, ga)
  cc <- P_cancel(ca, gd, include = TRUE)
  ba <- P_sub(P_mul(a, fa), P_mul(P_mul(dn, derivation(h, DE)), fd))
  bb <- P_cancel(ba, fd, include = TRUE)
  list(a, list(bb[[1L]], bb[[2L]]), list(cc[[1L]], cc[[2L]]), h)
}

.special_denom_cancel_bound <- function(a, ba, bd, n, DE, case) {
  if (case == "exp") {
    dcoeff <- P_quo(DE$d, .P(DE$t, DE$t))
    .with_dec(DE, {
      al <- frac_in(.div(.div(.neg(P_eval(ba, S_Zero)), P_eval(bd, S_Zero)), P_eval(a, S_Zero)), DE$t)
      et <- frac_in(dcoeff, DE$t)
      A <- parametric_log_deriv(al[[1L]], al[[2L]], et[[1L]], et[[2L]], DE)
      if (!is.null(A) && is_one_num(as_sym(A[[1L]]))) n <- min(n, as.double(gmp::as.bigq(as_sym(A[[2L]])$value)))
    })
  } else if (case == "tan") {
    .ri_notimpl("the hypertangent case of special_denom")
  }
  n
}

special_denom <- function(a, ba, bd, ca, cd, DE, case = "auto") {
  if (case == "auto") case <- DE$case
  if (case == "exp") p <- .P(DE$t, DE$t)
  else if (case == "tan") p <- .P(.add(list(.sq(DE$t), S_One)), DE$t)
  else if (case %in% c("primitive", "base")) {
    B <- P_quo(P_to_field(ba), bd)
    C <- P_quo(P_to_field(ca), cd)
    return(list(a, B, C, .P(S_One, DE$t)))
  } else if (case %in% c("other_linear", "other_nonlinear")) .ri_notimpl(paste0("The ", case, " case is not implemented in special_denom()."))
  else stop("bad case")
  nb <- order_at(ba, p, DE$t) - order_at(bd, p, DE$t)
  nc <- order_at(ca, p, DE$t) - order_at(cd, p, DE$t)
  n <- min(0, nc - min(0, nb))
  if (!nb) n <- .special_denom_cancel_bound(a, ba, bd, n, DE, case)
  N <- max(0, -nb, n - nc)
  pN <- P_pow(p, N)
  pn <- P_pow(p, -n)
  A <- P_mul(a, pN)
  B <- P_add(P_mul(ba, P_quo(pN, bd)), P_mul(P_mul(P_mul(.P(.int(n), DE$t), a), P_quo(derivation(p, DE), p)), pN))
  C <- P_quo(P_mul(P_mul(ca, pN), pn), cd)
  list(A, B, C, pn)
}

bound_degree <- function(a, b, cQ, DE, case = "auto", parametric = FALSE) {
  if (case == "auto") case <- DE$case
  da <- P_deg(a); db <- P_deg(b)
  dc <- if (parametric) max(vapply(cQ, P_deg, 0)) else P_deg(cQ)
  alpha <- cancel(.div(.neg(P_LC(.P(P_expr(b), DE$t))), P_LC(.P(P_expr(a), DE$t))))
  ival <- function(e) as.double(gmp::as.bigq(e$value))
  if (case == "base") {
    n <- max(0, dc - max(db, da - 1))
    if (db == da - 1 && alpha$type == "Integer") n <- max(0, ival(alpha), dc - db)
  } else if (case == "primitive") {
    n <- if (db > da) max(0, dc - db) else max(0, dc - da + 1)
    et <- frac_in(DE$d, DE$T[[.de_idx(DE) - 1L]])
    t1 <- DE$t
    .with_dec(DE, {
      al <- frac_in(alpha, DE$t)
      if (db == da - 1) {
        A <- tryCatch(limited_integrate(al[[1L]], al[[2L]], list(list(et[[1L]], et[[2L]])), DE), NonElementaryIntegralException = function(e) NULL)
        if (!is.null(A)) {
          m <- A[[2L]]
          if (length(m) != 1L) stop("Length of m should be 1")
          m <- as_sym(m[[1L]])
          if (m$type == "Integer") n <- max(n, ival(m))
        }
      } else if (db == da) {
        A <- is_log_deriv_k_t_radical_in_field(al[[1L]], al[[2L]], DE)
        if (!is.null(A)) {
          aa <- A[[1L]]; z <- A[[2L]]
          if (is_one_num(as_sym(aa))) {
            beta <- .div(.neg(P_LC(P_add(P_mul(a, .P(P_expr(derivation(z, DE)), t1)), P_mul(b, .P(z, t1))))), .mul(list(z, P_LC(a))))
            bt <- frac_in(beta, DE$t)
            A <- tryCatch(limited_integrate(bt[[1L]], bt[[2L]], list(list(et[[1L]], et[[2L]])), DE), NonElementaryIntegralException = function(e) NULL)
            if (!is.null(A)) {
              m <- as_sym(A[[2L]][[1L]])
              if (m$type == "Integer") n <- max(n, ival(m))
            }
          }
        }
      }
    })
  } else if (case == "exp") {
    n <- max(0, dc - max(db, da))
    if (da == db) {
      et <- frac_in(P_quo(DE$d, .P(DE$t, DE$t)), DE$T[[.de_idx(DE) - 1L]])
      .with_dec(DE, {
        al <- frac_in(alpha, DE$t)
        A <- parametric_log_deriv(al[[1L]], al[[2L]], et[[1L]], et[[2L]], DE)
        if (!is.null(A) && is_one_num(as_sym(A[[1L]]))) n <- max(n, ival(as_sym(A[[2L]])))
      })
    }
  } else if (case %in% c("tan", "other_nonlinear")) {
    delta <- P_deg(DE$d); lam <- P_LC(DE$d)
    alpha <- cancel(.div(alpha, lam))
    n <- max(0, dc - max(da + delta - 1, db))
    if (db == da + delta - 1 && alpha$type == "Integer") n <- max(0, ival(alpha), dc - db)
  } else if (case == "other_linear") .ri_notimpl("The other_linear case is not implemented in bound_degree().")
  else stop("bad case")
  n
}

spde <- function(a, b, c, n, DE) {
  zero <- .P(S_Zero, DE$t)
  alpha <- .P(S_One, DE$t); beta <- .P(S_Zero, DE$t)
  repeat {
    if (P_is_zero(c)) return(list(zero, zero, 0, zero, beta))
    if (n < 0) .nonelem("n became negative in spde().")
    g <- P_gcd(a, b)
    if (!P_is_zero(P_rem(c, g))) .nonelem("gcd(a, b) does not divide c in spde().")
    a <- P_quo(a, g); b <- P_quo(b, g); c <- P_quo(c, g)
    if (P_deg(a) == 0) {
      b <- P_quo(P_to_field(b), a)
      c <- P_quo(P_to_field(c), a)
      return(list(b, c, n, alpha, beta))
    }
    rz <- gcdex_diophantine(b, a, c)
    r <- rz[[1L]]; z <- rz[[2L]]
    b <- P_add(b, derivation(a, DE))
    c <- P_sub(z, derivation(r, DE))
    n <- n - P_deg(a)
    beta <- P_add(beta, P_mul(alpha, r))
    alpha <- P_mul(alpha, a)
  }
}

no_cancel_b_large <- function(b, c, n, DE) {
  q <- .P(S_Zero, DE$t)
  while (!P_is_zero(c)) {
    m <- P_deg(c) - P_deg(b)
    if (!(0 <= m && m <= n)) .nonelem("deg(c) - deg(b) is not between 0 and n in no_cancel_b_large().")
    p <- .P(.mul(list(.div(P_LC(.P(P_expr(c), DE$t)), P_LC(.P(P_expr(b), DE$t))), .pow(DE$t, .int(m)))), DE$t)
    q <- P_add(q, p)
    n <- m - 1
    c <- P_sub(P_sub(c, derivation(p, DE)), P_mul(b, p))
  }
  q
}

no_cancel_b_small <- function(b, c, n, DE) {
  q <- .P(S_Zero, DE$t)
  while (!P_is_zero(c)) {
    m <- if (n == 0) 0L else P_deg(c) - P_deg(DE$d) + 1L
    if (!(0 <= m && m <= n)) .nonelem("m is not between 0 and n in no_cancel_b_small().")
    if (m > 0) p <- .P(.mul(list(.div(P_LC(.P(P_expr(c), DE$t)), .mul(list(.int(m), P_LC(.P(P_expr(DE$d), DE$t))))), .pow(DE$t, .int(m)))), DE$t)
    else {
      if (P_deg(b) != P_deg(c)) .nonelem("deg(b) != deg(c) in the m == 0 case of no_cancel_b_small().")
      if (P_deg(b) == 0) return(list(q, .P(P_expr(b), DE$T[[.de_idx(DE) - 1L]]), .P(P_expr(c), DE$T[[.de_idx(DE) - 1L]])))
      p <- .P(.div(P_LC(.P(P_expr(c), DE$t)), P_LC(.P(P_expr(b), DE$t))), DE$t)
    }
    q <- P_add(q, p)
    n <- m - 1
    c <- P_sub(P_sub(c, derivation(p, DE)), P_mul(b, p))
  }
  q
}

no_cancel_equal <- function(b, c, n, DE) {
  q <- .P(S_Zero, DE$t)
  lc <- cancel(.div(.neg(P_LC(.P(P_expr(b), DE$t))), P_LC(.P(P_expr(DE$d), DE$t))))
  M <- if (lc$type == "Integer" && lc$value > 0) as.double(gmp::as.bigq(lc$value)) else -1
  while (!P_is_zero(c)) {
    m <- max(M, P_deg(c) - P_deg(DE$d) + 1)
    if (!(0 <= m && m <= n)) .nonelem("m is not between 0 and n in no_cancel_equal().")
    u <- cancel(.add(list(.mul(list(.int(m), P_LC(.P(P_expr(DE$d), DE$t)))), P_LC(.P(P_expr(b), DE$t)))))
    if (isTRUE(ask_fact(u, "zero"))) return(list(q, m, c))
    if (m > 0) p <- .P(.mul(list(.div(P_LC(.P(P_expr(c), DE$t)), u), .pow(DE$t, .int(m)))), DE$t)
    else if (P_deg(c) != P_deg(DE$d) - 1) .nonelem("deg(c) != deg(D) - 1 in the m == 0 case of no_cancel_equal().")
    else p <- .P(.div(P_LC(.P(P_expr(c), DE$t)), P_LC(.P(P_expr(b), DE$t))), DE$t)
    q <- P_add(q, p)
    n <- m - 1
    c <- P_sub(P_sub(c, derivation(p, DE)), P_mul(b, p))
  }
  q
}

cancel_primitive <- function(b, c, n, DE) {
  bb <- .with_dec(DE, {
    bf <- frac_in(b, DE$t)
    A <- is_log_deriv_k_t_radical_in_field(bf[[1L]], bf[[2L]], DE)
    if (!is.null(A) && is_one_num(as_sym(A[[1L]]))) .ri_notimpl("is_deriv_in_field() is required to  solve this problem.")
    bf
  })
  if (P_is_zero(c)) return(c)
  if (n < P_deg(c)) .nonelem("n < deg(c) in cancel_primitive().")
  q <- .P(S_Zero, DE$t)
  while (!P_is_zero(c)) {
    m <- P_deg(c)
    if (n < m) .nonelem("n < deg(c) in cancel_primitive().")
    s <- .with_dec(DE, {
      a2 <- frac_in(P_LC(c), DE$t)
      rischDE(bb[[1L]], bb[[2L]], a2[[1L]], a2[[2L]], DE)
    })
    stm <- .P(.mul(list(.div(P_expr(s[[1L]]), P_expr(s[[2L]])), .pow(DE$t, .int(m)))), DE$t)
    q <- P_add(q, stm)
    n <- m - 1
    c <- P_sub(c, P_add(P_mul(b, stm), derivation(stm, DE)))
  }
  q
}

cancel_exp <- function(b, c, n, DE) {
  eta <- P_expr(P_quo(DE$d, .P(DE$t, DE$t)))
  ef <- .with_dec(DE, {
    et <- frac_in(eta, DE$t)
    bf <- frac_in(b, DE$t)
    A <- parametric_log_deriv(bf[[1L]], bf[[2L]], et[[1L]], et[[2L]], DE)
    if (!is.null(A) && is_one_num(as_sym(A[[1L]]))) .ri_notimpl("is_deriv_in_field() is required to solve this problem.")
    et
  })
  etaa <- ef[[1L]]; etad <- ef[[2L]]
  if (P_is_zero(c)) return(c)
  if (n < P_deg(c)) .nonelem("n < deg(c) in cancel_exp().")
  q <- .P(S_Zero, DE$t)
  while (!P_is_zero(c)) {
    m <- P_deg(c)
    if (n < m) .nonelem("n < deg(c) in cancel_exp().")
    a1 <- P_expr(b)
    s <- .with_dec(DE, {
      a1f <- frac_in(a1, DE$t)
      a1a <- P_add(P_mul(a1f[[1L]], etad), P_mul(P_mul(etaa, a1f[[2L]]), .P(.int(m), DE$t)))
      a1d <- P_mul(a1f[[2L]], etad)
      a2 <- frac_in(P_LC(c), DE$t)
      rischDE(a1a, a1d, a2[[1L]], a2[[2L]], DE)
    })
    stm <- .P(.mul(list(.div(P_expr(s[[1L]]), P_expr(s[[2L]])), .pow(DE$t, .int(m)))), DE$t)
    q <- P_add(q, stm)
    n <- m - 1
    c <- P_sub(c, P_add(P_mul(b, stm), derivation(stm, DE)))
  }
  q
}

solve_poly_rde <- function(b, cQ, n, DE, parametric = FALSE) {
  dd <- P_deg(DE$d)
  if (!P_is_zero(b) && (DE$case == "base" || P_deg(b) > max(0, dd - 1))) {
    if (parametric) return(prde_no_cancel_b_large(b, cQ, n, DE))
    return(no_cancel_b_large(b, cQ, n, DE))
  }
  if ((P_is_zero(b) || P_deg(b) < dd - 1) && (DE$case == "base" || dd >= 2)) {
    if (parametric) return(prde_no_cancel_b_small(b, cQ, n, DE))
    R <- no_cancel_b_small(b, cQ, n, DE)
    if (inherits(R, "dpoly")) return(R)
    h <- R[[1L]]
    y <- .with_dec(DE, {
      b0 <- .P(P_expr(R[[2L]]), DE$t); c0 <- .P(P_expr(R[[3L]]), DE$t)
      .P(P_expr(solve_poly_rde(b0, c0, n, DE)), DE$t)
    })
    return(P_add(h, y))
  }
  if (dd >= 2 && P_deg(b) == dd - 1 && .is_true(Gt(.int(n), .div(.neg(P_LC(b)), P_LC(DE$d))))) {
    if (parametric) .ri_notimpl("prde_no_cancel_b_equal() is not yet implemented.")
    R <- no_cancel_equal(b, cQ, n, DE)
    if (inherits(R, "dpoly")) return(R)
    y <- solve_poly_rde(b, R[[3L]], R[[2L]], DE)
    return(P_add(R[[1L]], y))
  }
  if (P_is_zero(b)) .ri_notimpl("Remaining cases for Poly (P)RDE are not yet implemented (is_deriv_in_field() required).")
  if (DE$case == "exp") {
    if (parametric) .ri_notimpl("Parametric RDE cancellation hyperexponential case is not yet implemented.")
    return(cancel_exp(b, cQ, n, DE))
  }
  if (DE$case == "primitive") {
    if (parametric) .ri_notimpl("Parametric RDE cancellation primitive case is not yet implemented.")
    return(cancel_primitive(b, cQ, n, DE))
  }
  .ri_notimpl(paste0("Other Poly (P)RDE cancellation cases are not yet implemented (", DE$case, ")."))
}

rischDE <- function(fa, fd, ga, gd, DE) {
  wn <- weak_normalizer(fa, fd, DE)
  q <- wn[[1L]]; fa <- wn[[2L]][[1L]]; fd <- wn[[2L]][[2L]]
  gg <- P_cancel(P_mul(ga, q), gd, include = TRUE); ga <- gg[[1L]]; gd <- gg[[2L]]
  nd <- normal_denom(fa, fd, ga, gd, DE)
  a <- nd[[1L]]; ba <- nd[[2L]][[1L]]; bd <- nd[[2L]][[2L]]; ca <- nd[[3L]][[1L]]; cd <- nd[[3L]][[2L]]; hn <- nd[[4L]]
  sd <- special_denom(a, ba, bd, ca, cd, DE)
  A <- sd[[1L]]; B <- sd[[2L]]; C <- sd[[3L]]; hs <- sd[[4L]]
  n <- tryCatch(bound_degree(A, B, C, DE), NotImplementedError = function(e) Inf)
  sp <- spde(A, B, C, n, DE)
  B <- sp[[1L]]; C <- sp[[2L]]; m <- sp[[3L]]; alpha <- sp[[4L]]; beta <- sp[[5L]]
  y <- if (P_is_zero(C)) C else solve_poly_rde(B, C, m, DE)
  list(P_add(P_mul(alpha, y), beta), P_mul(P_mul(q, hn), hs))
}
