# Parametric Risch differential equations, after sympy/integrals/prde.py,
# with the PolyMatrix operations it uses (entries are expressions kept in
# lowest terms).

# ---------------------------------------------------------------------------
# matrices of expressions: list(nr, nc, e) with e row-major

.pm <- function(nr, nc, e = NULL) {
  if (is.null(e)) e <- rep(list(S_Zero), nr * nc)
  list(nr = as.integer(nr), nc = as.integer(nc), e = lapply(e, function(z) if (inherits(z, "dpoly")) P_expr(z) else as_sym(z)))
}
.pm_fun <- function(nr, nc, f) .pm(nr, nc, unlist(lapply(seq_len(nr) - 1L, function(i) lapply(seq_len(nc) - 1L, function(j) f(i, j))), recursive = FALSE))
.pm_get <- function(M, i, j) M$e[[(i - 1L) * M$nc + j]]
.pm_set <- function(M, i, j, v) { M$e[[(i - 1L) * M$nc + j]] <- v; M }
.pm_zeros <- function(nr, nc) .pm(nr, nc)
.pm_eye <- function(n) .pm_fun(n, n, function(i, j) if (i == j) S_One else S_Zero)
.pm_row <- function(M, i) .pm(1L, M$nc, M$e[(i - 1L) * M$nc + seq_len(M$nc)])
.pm_col <- function(M, j) .pm(M$nr, 1L, lapply(seq_len(M$nr), function(i) .pm_get(M, i, j)))
.pm_row_join <- function(A, B) {
  if (A$nr == 0L && A$nc == 0L) return(B)
  .pm_fun(A$nr, A$nc + B$nc, function(i, j) if (j < A$nc) .pm_get(A, i + 1L, j + 1L) else .pm_get(B, i + 1L, j - A$nc + 1L))
}
.pm_col_join <- function(A, B) {
  if (A$nr == 0L && A$nc == 0L) return(B)
  .pm(A$nr + B$nr, max(A$nc, B$nc), c(A$e, B$e))
}
.pm_sub <- function(A, B) .pm(A$nr, A$nc, Map(function(a, b) cancel(.sub(a, b)), A$e, B$e))
.pm_mul <- function(A, B) .pm_fun(A$nr, B$nc, function(i, j) cancel(.add(lapply(seq_len(A$nc), function(k) .mul(list(.pm_get(A, i + 1L, k), .pm_get(B, k, j + 1L)))))))
.pm_cols_range <- function(M, js) .pm_fun(M$nr, length(js), function(i, j) .pm_get(M, i + 1L, js[j + 1L]))
.pm_is_empty <- function(M) M$nr == 0L || M$nc == 0L

# cancel(), skipped for an expanded polynomial with rational coefficients
# (sums of rational multiples of products of symbol powers), which cancel()
# returns unchanged
.cancel_poly <- function(z) {
  atom <- function(f) f$type %in% c("Integer", "Rational", "Symbol", "Dummy") ||
    (f$type == "Pow" && f$args[[1L]]$type %in% c("Symbol", "Dummy") &&
       f$args[[2L]]$type == "Integer" && f$args[[2L]]$value > 0)
  mono <- function(t) if (t$type == "Mul") all(vapply(t$args, atom, NA)) else atom(t)
  ok <- if (z$type == "Add") all(vapply(z$args, mono, NA)) else mono(z)
  if (ok) z else cancel(z)
}

# A sparse row is list(j = increasing column indices, v = nonzero entries).
.sp_row <- function(entries) {
  v <- lapply(entries, function(z) if (is_zero_num(z)) z else .cancel_poly(z))
  nz <- which(!vapply(v, is_zero_num, NA))
  list(j = nz, v = v[nz])
}

# a - f*b for sparse rows
.sp_axpy <- function(a, f, b) {
  cols <- sort(unique(c(a$j, b$j)))
  ia <- match(cols, a$j); ib <- match(cols, b$j)
  v <- vector("list", length(cols))
  for (t in seq_along(cols)) {
    v[[t]] <- if (is.na(ib[t])) a$v[[ia[t]]]
      else if (is.na(ia[t])) .cancel_poly(.neg(.mul(list(f, b$v[[ib[t]]]))))
      else .cancel_poly(.sub(a$v[[ia[t]]], .mul(list(f, b$v[[ib[t]]]))))
  }
  keep <- !vapply(v, is_zero_num, NA)
  list(j = cols[keep], v = v[keep])
}

# Reduced row echelon form of sparse rows, one row at a time as in
# sdm_irref: list(rows, piv) with the nonzero rows ordered by pivot column;
# NULL as soon as column `stop` becomes a pivot. The pivot rows stay reduced
# against each other, so clearing a row's pivot columns adds no others.
.sp_rref <- function(rows, stop = NULL) {
  P <- list(); pc <- integer()
  for (r in rows) {
    for (h in which(pc %in% r$j)) {
      k <- match(pc[h], r$j)
      r <- .sp_axpy(r, r$v[[k]], P[[h]])
    }
    if (!length(r$j)) next
    c0 <- r$j[1L]
    if (!is.null(stop) && c0 == stop) return(NULL)
    pv <- r$v[[1L]]
    r$v <- c(list(S_One), lapply(r$v[-1L], function(z) .cancel_poly(.div(z, pv))))
    for (h in seq_along(P)) {
      k <- match(c0, P[[h]]$j)
      if (!is.na(k)) P[[h]] <- .sp_axpy(P[[h]], P[[h]]$v[[k]], r)
    }
    P[[length(P) + 1L]] <- r; pc <- c(pc, c0)
  }
  o <- order(pc)
  list(rows = P[o], piv = pc[o])
}

# reduced row echelon form: list(R, pivots)
.pm_rref <- function(M) {
  nr <- M$nr; nc <- M$nc
  s <- .sp_rref(lapply(seq_len(nr), function(i) .sp_row(M$e[(i - 1L) * nc + seq_len(nc)])))
  e <- rep(list(S_Zero), nr * nc)
  for (k in seq_along(s$rows)) e[(k - 1L) * nc + s$rows[[k]]$j] <- s$rows[[k]]$v
  M$e <- e
  list(M, s$piv)
}
.pm_nullspace <- function(M) {
  rr <- .pm_rref(M); R <- rr[[1L]]; piv <- rr[[2L]]
  free <- setdiff(seq_len(M$nc), piv)
  lapply(free, function(fj) {
    v <- rep(list(S_Zero), M$nc)
    v[[fj]] <- S_One
    for (k in seq_along(piv)) v[[piv[k]]] <- cancel(.neg(.pm_get(R, k, fj)))
    .pm(M$nc, 1L, v)
  })
}

# ---------------------------------------------------------------------------

prde_normal_denom <- function(fa, fd, G, DE) {
  dn <- splitfactor(fd, DE)[[1L]]
  gd <- Reduce(function(i, j) P_lcm(i, j), lapply(G, `[[`, 2L), .P(S_One, DE$t))
  en <- splitfactor(gd, DE)[[1L]]
  p <- P_gcd(dn, en)
  h <- P_quo(P_gcd(en, P_diff(en)), P_gcd(p, P_diff(p)))
  a <- P_mul(dn, h)
  c <- P_mul(a, h)
  ba <- P_sub(P_mul(a, fa), P_mul(P_mul(dn, derivation(h, DE)), fd))
  bb <- P_cancel(ba, fd, include = TRUE)
  G <- lapply(G, function(AD) P_cancel(P_mul(c, AD[[1L]]), AD[[2L]], include = TRUE))
  list(a, list(bb[[1L]], bb[[2L]]), G, h)
}

prde_special_denom <- function(a, ba, bd, G, DE, case = "auto") {
  if (case == "auto") case <- DE$case
  if (case == "exp") p <- .P(DE$t, DE$t)
  else if (case == "tan") p <- .P(.add(list(.sq(DE$t), S_One)), DE$t)
  else if (case %in% c("primitive", "base")) return(list(a, P_quo(ba, bd), G, .P(S_One, DE$t)))
  else .ri_notimpl(paste0("The ", case, " case is not implemented in prde_special_denom()."))
  nb <- order_at(ba, p, DE$t) - order_at(bd, p, DE$t)
  nc <- min(vapply(G, function(g) order_at(g[[1L]], p, DE$t) - order_at(g[[2L]], p, DE$t), 0))
  n <- min(0, nc - min(0, nb))
  if (!nb) n <- .special_denom_cancel_bound(a, ba, bd, n, DE, case)
  N <- max(0, -nb)
  pN <- P_pow(p, N); pn <- P_pow(p, -n)
  A <- P_mul(a, pN)
  B <- P_add(P_mul(ba, P_quo(pN, bd)), P_mul(P_mul(P_mul(.P(.int(n), DE$t), a), P_quo(derivation(p, DE), p)), pN))
  G <- lapply(G, function(g) P_cancel(P_mul(P_mul(g[[1L]], pN), pn), g[[2L]], include = TRUE))
  list(A, B, G, pn)
}

prde_linear_constraints <- function(a, b, G, DE) {
  m <- length(G)
  d <- Reduce(function(i, j) P_lcm(i, j), lapply(G, `[[`, 2L))
  d <- .P(P_expr(d), d$x, field = TRUE)
  Q <- lapply(G, function(g) P_div(P_mul(g[[1L]], P_quo(d, g[[2L]])), d))
  if (!all(vapply(Q, function(q) P_is_zero(q[[2L]]), NA))) {
    N <- max(vapply(Q, function(q) P_deg(q[[2L]]), 0))
    M <- .pm_fun(N + 1L, m, function(i, j) P_nth(Q[[j + 1L]][[2L]], i))
  } else M <- .pm(0L, m, list())
  list(lapply(Q, `[[`, 1L), M)
}

poly_linear_constraints <- function(p, d) {
  m <- length(p)
  qr <- lapply(p, function(pi) P_div(pi, d))
  if (!all(vapply(qr, function(z) P_is_zero(z[[2L]]), NA))) {
    n <- max(vapply(qr, function(z) P_deg(z[[2L]]), 0))
    M <- .pm_fun(n + 1L, m, function(i, j) P_nth(qr[[j + 1L]][[2L]], i))
  } else M <- .pm(0L, m, list())
  list(lapply(qr, `[[`, 1L), M)
}

constant_system <- function(A, u, DE) {
  if (.pm_is_empty(A)) return(list(A, u))
  Au <- .pm_rref(.pm_row_join(A, u))[[1L]]
  A <- .pm_cols_range(Au, seq_len(Au$nc - 1L)); u <- .pm_col(Au, Au$nc)
  D <- function(x) derivation(x, DE, basic = TRUE)
  nr0 <- A$nr; nc0 <- A$nc
  for (j in seq_len(nc0)) for (i in seq_len(nr0)) {
    aij <- .pm_get(A, i, j)
    if (.has_any_of(aij, DE$T)) {
      Ri <- .pm_row(A, i)
      DAij <- D(aij)
      Rm1 <- .pm(1L, Ri$nc, lapply(Ri$e, function(z) cancel(.div(D(z), DAij))))
      um1 <- cancel(.div(D(.pm_get(u, i, 1L)), DAij))
      Aj <- .pm_col(A, j)
      A <- .pm_sub(A, .pm_mul(Aj, Rm1))
      u <- .pm_sub(u, .pm_mul(Aj, .pm(1L, 1L, list(um1))))
      A <- .pm_col_join(A, Rm1)
      u <- .pm_col_join(u, .pm(1L, 1L, list(um1)))
    }
  }
  list(A, u)
}

prde_spde <- function(a, b, Q, n, DE) {
  RZ <- lapply(Q, function(qi) gcdex_diophantine(b, a, qi))
  R <- lapply(RZ, `[[`, 1L); Z <- lapply(RZ, `[[`, 2L)
  A <- a
  B <- P_add(b, derivation(a, DE))
  Qq <- Map(function(ri, zi) P_sub(zi, derivation(ri, DE)), R, Z)
  list(A, B, Qq, R, n - P_deg(a))
}

prde_no_cancel_b_large <- function(b, Q, n, DE) {
  db <- P_deg(b); m <- length(Q)
  H <- rep(list(.P(S_Zero, DE$t)), m)
  for (N in n:0) for (i in seq_len(m)) {
    si <- .div(P_nth(Q[[i]], N + db), P_LC(b))
    sitn <- .P(.mul(list(si, .pow(DE$t, .int(N)))), DE$t)
    H[[i]] <- P_add(H[[i]], sitn)
    Q[[i]] <- P_sub(P_sub(Q[[i]], derivation(sitn, DE)), P_mul(b, sitn))
  }
  dc <- if (all(vapply(Q, P_is_zero, NA))) -1L else max(vapply(Q, P_deg, 0))
  M <- .pm_fun(dc + 1L, m, function(i, j) P_nth(Q[[j + 1L]], i))
  A <- constant_system(M, .pm_zeros(dc + 1L, 1L), DE)[[1L]]
  c <- .pm_eye(m)
  A <- .pm_col_join(.pm_row_join(A, .pm_zeros(A$nr, m)), .pm_row_join(c, .pm(m, m, lapply(c$e, .neg))))
  list(H, A)
}

prde_no_cancel_b_small <- function(b, Q, n, DE) {
  m <- length(Q)
  H <- rep(list(.P(S_Zero, DE$t)), m)
  if (n >= 1) for (N in n:1) for (i in seq_len(m)) {
    si <- .div(P_nth(Q[[i]], N + P_deg(DE$d) - 1L), .mul(list(.int(N), P_LC(DE$d))))
    sitn <- .P(.mul(list(si, .pow(DE$t, .int(N)))), DE$t)
    H[[i]] <- P_add(H[[i]], sitn)
    Q[[i]] <- P_sub(P_sub(Q[[i]], derivation(sitn, DE)), P_mul(b, sitn))
  }
  if (P_deg(b) > 0) {
    for (i in seq_len(m)) {
      si <- .P(.div(P_nth(Q[[i]], P_deg(b)), P_LC(b)), DE$t)
      H[[i]] <- P_add(H[[i]], si)
      Q[[i]] <- P_sub(P_sub(Q[[i]], derivation(si, DE)), P_mul(b, si))
    }
    dc <- if (all(vapply(Q, P_is_zero, NA))) -1L else max(vapply(Q, P_deg, 0))
    M <- .pm_fun(dc + 1L, m, function(i, j) P_nth(Q[[j + 1L]], i))
    A <- constant_system(M, .pm_zeros(dc + 1L, 1L), DE)[[1L]]
    c <- .pm_eye(m)
    A <- .pm_col_join(.pm_row_join(A, .pm_zeros(A$nr, m)), .pm_row_join(c, .pm(m, m, lapply(c$e, .neg))))
    return(list(H, A))
  }
  t <- DE$t
  if (DE$case != "base") {
    fb <- .with_dec(DE, {
      t0 <- DE$t
      bf <- frac_in(b, t0, field = TRUE)
      Q0 <- lapply(Q, function(qi) frac_in(P_TC(qi), t0, field = TRUE))
      param_rischDE(bf[[1L]], bf[[2L]], Q0, DE)
    })
    f <- lapply(fb[[1L]], function(fad) .P(.div(P_expr(fad[[1L]]), P_expr(fad[[2L]])), t, field = TRUE))
    B <- fb[[2L]]
  } else {
    f <- list(.P(S_One, t, field = TRUE))
    B <- .pm(1L, m + 1L, c(lapply(Q, P_TC), list(S_Zero)))
  }
  d <- max(vapply(Q, P_deg, 0))
  if (d > 0) {
    M <- .pm_fun(d, m, function(i, j) P_nth(Q[[j + 1L]], i + 1L))
    A <- constant_system(M, .pm_zeros(d, 1L), DE)[[1L]]
  } else A <- .pm(0L, m, list())
  r <- length(f)
  I <- .pm_eye(m)
  A <- .pm_row_join(A, .pm_zeros(A$nr, r + m))
  B <- .pm_row_join(B, .pm_zeros(B$nr, m))
  C <- .pm_row_join(.pm_row_join(I, .pm_zeros(m, r)), .pm(m, m, lapply(I$e, .neg)))
  list(c(f, H), .pm_col_join(.pm_col_join(A, B), C))
}

prde_cancel_liouvillian <- function(b, Q, n, DE) {
  H <- list()
  if (DE$case == "primitive") bb <- .with_dec(DE, frac_in(b, DE$t, field = TRUE))
  M <- NULL
  for (i in n:0) {
    if (DE$case == "exp") {
      bb <- .with_dec(DE, frac_in(P_add(b, .P(.mul(list(.int(i), .div(P_expr(derivation(DE$t, DE)), DE$t))), b$x)), DE$t, field = TRUE))
    }
    fA <- .with_dec(DE, {
      Qy <- lapply(Q, function(q) frac_in(P_nth(q, i), DE$t, field = TRUE))
      param_rischDE(bb[[1L]], bb[[2L]], Qy, DE)
    })
    fi <- lapply(fA[[1L]], function(fad) .P(.div(P_expr(fad[[1L]]), P_expr(fad[[2L]])), DE$t, field = TRUE))
    Ai <- fA[[2L]]
    ri <- length(fi)
    M <- if (i == n) Ai else .pm_col_join(Ai, .pm_row_join(M, .pm_zeros(M$nr, ri)))
    Fi <- list(); hi <- list()
    for (j in seq_len(ri)) {
      hji <- P_mul(fi[[j]], .P(.pow(DE$t, .int(i)), fi[[j]]$x))
      hi[[j]] <- hji
      Fi[[j]] <- P_neg(P_sub(derivation(hji, DE), P_mul(b, hji)))
    }
    H <- c(H, hi)
    Q <- c(Q, Fi)
  }
  list(H, M)
}

param_poly_rischDE <- function(a, b, q, n, DE) {
  m <- length(q)
  if (n < 0) {
    if (all(vapply(q, P_is_zero, NA))) return(list(list(), .pm_zeros(1L, m)))
    N <- max(vapply(q, P_deg, 0))
    M <- .pm_fun(N + 1L, m, function(i, j) P_nth(q[[j + 1L]], i))
    A <- constant_system(M, .pm_zeros(M$nr, 1L), DE)[[1L]]
    return(list(list(), A))
  }
  if (P_is_ground(a)) {
    a <- P_LC(a)
    b <- P_exquo_ground(P_to_field(b), a)
    q <- lapply(q, function(qi) P_exquo_ground(P_to_field(qi), a))
    dd <- P_deg(DE$d)
    if (!P_is_zero(b) && (DE$case == "base" || P_deg(b) > max(0, dd - 1))) return(prde_no_cancel_b_large(b, q, n, DE))
    if ((P_is_zero(b) || P_deg(b) < dd - 1) && (DE$case == "base" || dd >= 2)) return(prde_no_cancel_b_small(b, q, n, DE))
    if (dd >= 2 && P_deg(b) == dd - 1 && .is_true(Gt(.int(n), .div(.neg(P_LC(b)), P_LC(DE$d))))) .ri_notimpl("prde_no_cancel_b_equal() is not yet implemented.")
    if (DE$case %in% c("primitive", "exp")) return(prde_cancel_liouvillian(b, q, n, DE))
    .ri_notimpl("non-linear and hypertangent cases have not yet been implemented")
  }
  alpha <- .P(S_One, a$x, domain = a$K); beta <- rep(list(.P(S_Zero, a$x, domain = a$K)), m)
  d <- NULL
  while (n >= 0) {
    sp <- prde_spde(a, b, q, n, DE)
    a <- sp[[1L]]; b <- sp[[2L]]; q <- sp[[3L]]; r <- sp[[4L]]; n <- sp[[5L]]
    beta <- Map(function(bi, ri) P_add(bi, P_mul(alpha, ri)), beta, r)
    alpha <- P_mul(alpha, a)
    d <- P_gcd(a, b)
    if (!P_is_ground(d)) break
  }
  qm <- poly_linear_constraints(q, d)
  qq <- qm[[1L]]; M <- qm[[2L]]
  A <- constant_system(M, .pm_zeros(M$nr, 1L), DE)[[1L]]
  V <- .pm_nullspace(A)
  if (!length(V)) return(list(list(), .pm_eye(m)))
  r <- lapply(V, function(vj) .P(.add(Map(function(qk, vk) .mul(list(P_expr(qk), vk)), qq, vj$e)), DE$t))
  f <- lapply(V, function(vj) .P(.add(Map(function(bk, vk) .mul(list(P_expr(bk), vk)), beta, vj$e)), DE$t))
  gB <- param_poly_rischDE(P_quo(a, d), P_quo(b, d), r, n, DE)
  g <- gB[[1L]]; B <- gB[[2L]]
  h <- c(f, lapply(g, function(gk) P_mul(alpha, gk)))
  A <- .pm(m, m, lapply(.pm_eye(m)$e, .neg))
  for (vj in V) A <- .pm_row_join(A, vj)
  A <- .pm_row_join(A, .pm_zeros(m, length(g)))
  A <- .pm_col_join(A, .pm_row_join(.pm_zeros(B$nr, m), B))
  list(h, A)
}

param_rischDE <- function(fa, fd, G, DE) {
  m <- length(G)
  wn <- weak_normalizer(fa, fd, DE)
  q <- wn[[1L]]; fa <- wn[[2L]][[1L]]; fd <- wn[[2L]][[2L]]
  gamma <- q
  G <- lapply(G, function(g) P_cancel(P_mul(q, g[[1L]]), g[[2L]], include = TRUE))
  nd <- prde_normal_denom(fa, fd, G, DE)
  a <- nd[[1L]]; ba <- nd[[2L]][[1L]]; bd <- nd[[2L]][[2L]]; G <- nd[[3L]]; hn <- nd[[4L]]
  gamma <- P_mul(gamma, hn)
  sd <- prde_special_denom(a, ba, bd, G, DE)
  A <- sd[[1L]]; B <- sd[[2L]]; G <- sd[[3L]]; hs <- sd[[4L]]
  gamma <- P_mul(gamma, hs)
  g <- P_gcd(A, B)
  a <- P_quo(A, g); b <- P_quo(B, g)
  g <- lapply(G, function(gg) P_cancel(gg[[1L]], P_mul(gg[[2L]], g), include = TRUE))
  lc <- prde_linear_constraints(a, b, g, DE)
  qv <- lc[[1L]]; M <- lc[[2L]]
  M <- constant_system(M, .pm_zeros(M$nr, 1L), DE)[[1L]]
  V <- .pm_nullspace(M)
  if (!length(V)) return(list(list(), .pm_eye(m)))
  r <- lapply(V, function(vj) .P(.add(Map(function(qk, vk) .mul(list(P_expr(qk), vk)), qv, vj$e)), DE$t))
  n <- tryCatch(bound_degree(a, b, r, DE, parametric = TRUE), NotImplementedError = function(e) 5)
  hB <- param_poly_rischDE(a, b, r, n, DE)
  h <- hB[[1L]]; B <- hB[[2L]]
  A <- .pm(m, m, lapply(.pm_eye(m)$e, .neg))
  for (vj in V) A <- .pm_row_join(A, vj)
  A <- .pm_row_join(A, .pm_zeros(m, length(h)))
  A <- .pm_col_join(A, .pm_row_join(.pm_zeros(B$nr, m), B))
  W <- .pm_nullspace(A)
  v <- length(h)
  elements <- lapply(W, function(wl) c(wl$e[seq_len(m)], if (v) wl$e[(length(wl$e) - v + 1L):length(wl$e)] else list()))
  Mm <- .pm(length(W), m + v, unlist(elements, recursive = FALSE))
  N <- .pm_nullspace(Mm)
  C <- .pm(length(N), m + v, unlist(lapply(N, `[[`, "e"), recursive = FALSE))
  list(lapply(h, function(hk) P_cancel(hk, gamma, include = TRUE)), C)
}

limited_integrate <- function(fa, fd, G, DE) {
  fa <- P_mul(fa, .P(.inv(P_LC(fd)), DE$t)); fd <- P_monic(fd)
  Fa <- .P(S_Zero, DE$t); Fd <- .P(S_One, DE$t)
  G <- c(list(list(fa, fd)), G)
  hA <- param_rischDE(Fa, Fd, G, DE)
  h <- hA[[1L]]; A <- hA[[2L]]
  V <- Filter(function(v) !is_zero_num(v$e[[1L]]), .pm_nullspace(A))
  if (!length(V)) return(NULL)
  c0 <- V[[1L]]$e[[1L]]
  v <- lapply(V[[1L]]$e, function(z) cancel(.div(z, .neg(c0))))
  r <- length(h)
  m <- length(v) - r - 1L
  C <- v[seq_len(m) + 1L]
  y <- .neg(.add(lapply(seq_len(r), function(i) .div(.mul(list(v[[m + 1L + i]], P_expr(h[[i]][[1L]]))), P_expr(h[[i]][[2L]])))))
  nd <- .as_numer_denom(y)
  Ya <- .P(nd[[1L]], DE$t); Yd <- .P(nd[[2L]], DE$t)
  list(list(P_mul(Ya, .P(.inv(P_LC(Yd)), DE$t)), P_monic(Yd)), C)
}

parametric_log_deriv_heu <- function(fa, fd, wa, wd, DE, c1 = NULL) {
  f <- .div(P_expr(fa), P_expr(fd)); w <- .div(P_expr(wa), P_expr(wd))
  if (f$type %in% c("Integer", "Rational") && w$type %in% c("Integer", "Rational")) {
    if (is_zero_num(f) && is_zero_num(w)) return(list(S_One, S_Zero, .P(S_One, DE$t)))
    fq <- gmp::as.bigq(f$value); wq <- gmp::as.bigq(w$value)
    xx <- gmp::numerator(fq) * gmp::denominator(wq); yy <- gmp::numerator(wq) * gmp::denominator(fq)
    g <- gmp::gcd.bigz(xx, yy); xx <- xx %/% g; yy <- yy %/% g
    n <- yy; m <- xx
    if (n == 0) return(NULL)
    if (n < 0) { n <- -n; m <- -m }
    return(list(.zs(n), .zs(m), .P(S_One, DE$t)))
  }
  if (is.null(c1)) c1 <- Dummy("c1")
  pa <- P_div(fa, fd); qb <- P_div(wa, wd)
  p <- pa[[1L]]; a <- pa[[2L]]; q <- qb[[1L]]; b <- qb[[2L]]
  B <- max(0, P_deg(derivation(DE$t, DE)) - 1)
  C <- max(P_deg(p), P_deg(q))
  solve_c1 <- function(eqs) {
    eqs <- Filter(function(e) !is_zero_num(cancel(e)), eqs)
    if (!length(eqs)) return(NULL)
    vals <- lapply(eqs, function(e) { sl <- solve(e, c1); if (length(sl) == 1L) sl[[1L]] else NULL })
    if (any(vapply(vals, is.null, NA))) return(NULL)
    if (length(unique(vapply(vals, function(v) cancel(v)$key, ""))) != 1L) return(NULL)
    cancel(vals[[1L]])
  }
  if (P_deg(q) > B) {
    s <- solve_c1(lapply((B + 1):C, function(i) .sub(P_nth(p, i), .mul(list(c1, P_nth(q, i))))))
    if (is.null(s) || !(s$type %in% c("Integer", "Rational"))) return(NULL)
    MN <- .as_numer_denom(s)
    nfmwa <- P_sub(P_mul(P_mul(.P(MN[[2L]], q$x), fa), wd), P_mul(P_mul(.P(MN[[1L]], q$x), wa), fd))
    nfmwd <- P_mul(fd, wd)
    Qv <- is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE, "auto")
    if (is.null(Qv)) return(NULL)
    Q <- as_sym(Qv[[1L]]); v <- Qv[[2L]]
    if (is_zero_num(Q) || is_zero_num(as_sym(v))) return(NULL)
    return(list(.mul(list(Q, MN[[2L]])), .mul(list(Q, MN[[1L]])), v))
  }
  if (P_deg(p) > B) return(NULL)
  c <- lcm(P_LC(.P(P_expr(fd), DE$t)), P_LC(.P(P_expr(wd), DE$t)))
  l <- P_mul(P_lcm(P_monic(fd), P_monic(wd)), .P(c, DE$t))
  sp <- splitfactor(l, DE)
  ln <- sp[[1L]]; ls <- sp[[2L]]
  z <- P_mul(ls, P_gcd(ln, P_diff(ln)))
  if (P_deg(z) < 1) {
    if (DE$case == "base") {
      if (P_is_zero(q)) {
        if (!P_is_zero(p)) return(NULL)
        .ri_notimpl("parametric_log_deriv_heu() heuristic failed: the full parametric logarithmic derivative problem (using the structure theorems) is not yet implemented.")
      }
      if (P_is_zero(p)) cc <- S_Zero
      else {
        if (!identical(P_expr(P_mul(p, .P(P_LC(q), DE$t)))$key, P_expr(P_mul(q, .P(P_LC(p), DE$t)))$key)) return(NULL)
        cc <- cancel(.div(P_LC(p), P_LC(q)))
        if (!(cc$type %in% c("Integer", "Rational"))) return(NULL)
      }
      MN <- .as_numer_denom(cc)
      nfmwa <- P_sub(P_mul(P_mul(.P(MN[[2L]], DE$t), fa), wd), P_mul(P_mul(.P(MN[[1L]], DE$t), wa), fd))
      nfmwd <- P_mul(fd, wd)
      Qv <- is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE)
      if (is.null(Qv)) return(NULL)
      Q <- as_sym(Qv[[1L]]); v <- Qv[[2L]]
      if (is_zero_num(Q) || is_zero_num(as_sym(v))) return(NULL)
      return(list(.mul(list(Q, MN[[2L]])), .mul(list(Q, MN[[1L]])), v))
    }
    .ri_notimpl("parametric_log_deriv_heu() heuristic failed: the full parametric logarithmic derivative problem (using the structure theorems) is not yet implemented.")
  }
  ur1 <- P_div(P_mul(fa, P_quo(l, fd)), z); ur2 <- P_div(P_mul(wa, P_quo(l, wd)), z)
  r1 <- ur1[[2L]]; r2 <- ur2[[2L]]
  s <- solve_c1(lapply(seq_len(P_deg(z)) - 1L, function(i) .sub(P_nth(r1, i), .mul(list(c1, P_nth(r2, i))))))
  if (is.null(s) || !(s$type %in% c("Integer", "Rational"))) return(NULL)
  MN <- .as_numer_denom(s)
  nfmwa <- P_sub(P_mul(P_mul(.P(MN[[2L]], DE$t), fa), wd), P_mul(P_mul(.P(MN[[1L]], DE$t), wa), fd))
  nfmwd <- P_mul(fd, wd)
  Qv <- is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE)
  if (is.null(Qv)) return(NULL)
  Q <- as_sym(Qv[[1L]]); v <- Qv[[2L]]
  if (is_zero_num(Q) || is_zero_num(as_sym(v))) return(NULL)
  list(.mul(list(Q, MN[[2L]])), .mul(list(Q, MN[[1L]])), v)
}

parametric_log_deriv <- function(fa, fd, wa, wd, DE) parametric_log_deriv_heu(fa, fd, wa, wd, DE)

.structure_check <- function(DE) {
  if (length(DE$exts) != length(DE$D)) {
    prim <- which(unlist(DE$cases) == "primitive")
    if (any(unlist(DE$cases) == "tan") || length(setdiff(prim, .de_indices(DE, "log"))))
      .ri_notimpl("Real version of the structure theorems with hypertangent support is not yet implemented.")
    .ri_notimpl("Nonelementary extensions not supported in the structure theorems.")
  }
}
.structure_parts <- function(DE) {
  E_part <- lapply(.de_indices(DE, "exp"), function(i) P_expr(P_quo(DE$D[[i]], .P(DE$T[[i]], DE$T[[i]]))))
  L_part <- lapply(.de_indices(DE, "log"), function(i) P_expr(DE$D[[i]]))
  c(E_part, L_part)
}

is_deriv_k <- function(fa, fd, DE) {
  dfa <- P_sub(P_mul(fd, derivation(fa, DE)), P_mul(fa, derivation(fd, DE))); dfd <- P_mul(fd, fa)
  cc <- P_cancel(dfa, dfd, include = TRUE); dfa <- cc[[1L]]; dfd <- cc[[2L]]
  .structure_check(DE)
  parts <- .structure_parts(DE)
  lhs <- .pm(1L, length(parts), parts)
  rhs <- .pm(1L, 1L, list(.div(P_expr(dfa), P_expr(dfd))))
  Au <- constant_system(lhs, rhs, DE)
  A <- Au[[1L]]; u <- Au[[2L]]$e
  if (.pm_is_empty(A) || !all(vapply(u, function(i) is_zero_num(derivation(i, DE, basic = TRUE)), NA))) return(NULL)
  if (!all(vapply(u, function(i) i$type %in% c("Integer", "Rational"), NA))) .ri_notimpl("Cannot work with non-rational coefficients in this case.")
  terms <- c(lapply(.de_indices(DE, "exp"), function(i) DE$extargs[[i]]), lapply(.de_indices(DE, "log"), function(i) DE$T[[i]]))
  ans <- Map(function(a, b) list(a, b), terms, u)
  result <- .add(lapply(ans, function(p) .mul(list(p[[1L]], p[[2L]]))))
  argterms <- c(lapply(.de_indices(DE, "exp"), function(i) DE$T[[i]]), lapply(.de_indices(DE, "log"), function(i) DE$extargs[[i]]))
  l <- list(); ld <- list()
  for (k in seq_along(argterms)) {
    nd <- .as_numer_denom(argterms[[k]]); j <- u[[k]]
    sq <- sqf_list(nd[[1L]])
    l[[k]] <- .mul(c(list(.pow(sq$coeff, j)), lapply(sq$factors, function(be) .pow(be[[1L]], .mul(list(.int(be[[2L]]), j))))))
    sq <- sqf_list(nd[[2L]])
    ld[[k]] <- .mul(c(list(.pow(sq$coeff, j)), lapply(sq$factors, function(be) .pow(be[[1L]], .mul(list(.int(be[[2L]]), j))))))
  }
  const <- cancel(.mul(list(.div(.div(P_expr(fa), P_expr(fd)), .mul(l)), .mul(ld))))
  list(ans, result, const)
}

is_log_deriv_k_t_radical <- function(fa, fd, DE, Df = TRUE) {
  if (Df) {
    cc <- P_cancel(P_sub(P_mul(fd, derivation(fa, DE)), P_mul(fa, derivation(fd, DE))), P_pow(fd, 2L), include = TRUE)
    dfa <- cc[[1L]]; dfd <- cc[[2L]]
  } else { dfa <- fa; dfd <- fd }
  .structure_check(DE)
  parts <- .structure_parts(DE)
  lhs <- .pm(1L, length(parts), parts)
  rhs <- .pm(1L, 1L, list(.div(P_expr(dfa), P_expr(dfd))))
  Au <- constant_system(lhs, rhs, DE)
  A <- Au[[1L]]; u <- Au[[2L]]$e
  if (.pm_is_empty(A) || !all(vapply(u, function(i) is_zero_num(derivation(i, DE, basic = TRUE)), NA))) return(NULL)
  if (!all(vapply(u, function(i) i$type %in% c("Integer", "Rational"), NA))) .ri_notimpl("Cannot work with non-rational coefficients in this case.")
  n <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), lapply(u, function(i) .as_numer_denom(i)[[2L]]), S_One)
  u <- lapply(u, function(i) .mul(list(i, n)))
  terms <- c(lapply(.de_indices(DE, "exp"), function(i) DE$T[[i]]), lapply(.de_indices(DE, "log"), function(i) DE$extargs[[i]]))
  ans <- Map(function(a, b) list(a, b), terms, u)
  result <- .mul(lapply(ans, function(p) .pow(p[[1L]], p[[2L]])))
  argterms <- c(lapply(.de_indices(DE, "exp"), function(i) DE$extargs[[i]]), lapply(.de_indices(DE, "log"), function(i) DE$T[[i]]))
  const <- cancel(.sub(.div(P_expr(fa), P_expr(fd)), .add(Map(function(i, j) .mul(list(i, .div(j, n))), argterms, u))))
  list(ans, result, n, const)
}

is_log_deriv_k_t_radical_in_field <- function(fa, fd, DE, case = "auto", z = NULL) {
  cc <- P_cancel(fa, fd, include = TRUE); fa <- cc[[1L]]; fd <- cc[[2L]]
  if (is.null(z)) z <- Dummy("z")
  rr <- residue_reduce(fa, fd, DE, z = z)
  H <- rr[[1L]]; b <- rr[[2L]]
  if (!b) return(NULL)
  rts <- lapply(H, function(hi) real_roots(P_expr(hi[[1L]]), z))
  if (!all(vapply(seq_along(H), function(k) length(rts[[k]]) == P_deg(H[[k]][[1L]]) &&
                    all(vapply(rts[[k]], function(r) r$type %in% c("Integer", "Rational"), NA)), NA))) return(NULL)
  residueterms <- unlist(lapply(seq_along(H), function(j) lapply(rts[[j]], function(i) list(subs(P_expr(H[[j]][[2L]]), z, i), i))), recursive = FALSE)
  p <- cancel(.sub(.div(P_expr(fa), P_expr(fd)), residue_reduce_derivation(H, DE, z)))
  p <- P_or_null(p, DE$t)
  if (is.null(p)) return(NULL)
  if (P_deg(p) >= max(1, P_deg(DE$d))) return(NULL)
  if (case == "auto") case <- DE$case
  if (case == "exp") {
    ww <- P_cancel(derivation(DE$t, DE), .P(DE$t, DE$t), include = TRUE)
    A <- .with_dec(DE, {
      pp <- frac_in(p, DE$t, cancel = TRUE)
      w2 <- frac_in(ww, DE$t)
      parametric_log_deriv(pp[[1L]], pp[[2L]], w2[[1L]], w2[[2L]], DE)
    })
    if (is.null(A)) return(NULL)
    n <- A[[1L]]; e <- A[[2L]]; u <- A[[3L]]
    u <- .mul(list(if (inherits(u, "dpoly")) P_expr(u) else u, .pow(DE$t, e)))
  } else if (case == "primitive") {
    A <- .with_dec(DE, {
      pp <- frac_in(p, DE$t)
      is_log_deriv_k_t_radical_in_field(pp[[1L]], pp[[2L]], DE, case = "auto")
    })
    if (is.null(A)) return(NULL)
    n <- A[[1L]]; u <- A[[2L]]
  } else if (case == "base") {
    if (!.dp_is_sqf(fd) || P_deg(fa) >= P_deg(fd)) return(NULL)
    n <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), lapply(residueterms, function(r) .as_numer_denom(r[[2L]])[[2L]]), S_One)
    u <- .mul(lapply(residueterms, function(r) .pow(r[[1L]], .mul(list(r[[2L]], n)))))
    return(list(n, u))
  } else if (case == "tan") .ri_notimpl("The hypertangent case is not yet implemented for is_log_deriv_k_t_radical_in_field()")
  else stop("bad case")
  n <- as_sym(n)
  common <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), c(lapply(residueterms, function(r) .as_numer_denom(r[[2L]])[[2L]]), list(n)), S_One)
  residueterms <- lapply(residueterms, function(r) list(r[[1L]], .mul(list(r[[2L]], common))))
  m <- .div(common, n)
  u <- cancel(.mul(c(list(.pow(if (inherits(u, "dpoly")) P_expr(u) else u, m)), lapply(residueterms, function(r) .pow(r[[1L]], r[[2L]])))))
  list(common, u)
}

.dp_is_sqf <- function(f) P_deg(P_gcd(f, P_diff(f))) <= 0

real_imag <- function(ba, bd, gen) {
  .ri_notimpl("real_imag (hypertangent case)")
}
