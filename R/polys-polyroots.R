# Roots of univariate polynomials, after sympy/polys/polyroots.py: roots(),
# roots_linear, roots_quadratic, roots_cubic, roots_quartic, roots_binomial,
# roots_cyclotomic, preprocess_roots, _integer_basis, and the helpers they
# use (dup_decompose, dup_cyclotomic_p, to_rational_coeffs, expand_2arg).
# Results keep SymPy's dict order.

.dk_numerical <- function(K) .dk(K) %in% c("ZZ", "QQ", "RR", "CC", "ZZ_I", "QQ_I", "AlgField")
.dk_composite <- function(K) .dk(K) %in% c("PolyRing", "FracField")

# expand_2arg: distribute a Number coefficient over an Add, unevaluated
expand_2arg <- function(e) .map_bottom_up(as_sym(e), function(z) {
  if (z$type == "Mul") {
    cm <- as_coeff_mul_split(z)
    c <- cm[[1L]]; r <- cm[[2L]]
    if (is_number_type(c) && r$type == "Add") return(.add(lapply(r$args, function(ri) .mul(list(c, ri)))))
  }
  z
})

.pr_all_coeffs <- function(f) lapply(f$c, .gd_out)
.pr_nth <- function(f, k) { n <- .dp_deg(f); if (k > n) S_Zero else .gd_out(f$c[[n - k + 1L]]) }
.pr_length <- function(f) sum(!vapply(f$c, is_zero_num, NA))

.pr_dom_simplify <- function(f, r) {
  dom <- f$K
  if (.dk_numerical(dom)) return(r)
  if (.dk_composite(dom)) sym_factor(r) else simplify(r)
}

roots_linear <- function(f) {
  r <- .neg(.div(.pr_nth(f, 0L), .pr_nth(f, 1L)))
  list(.pr_dom_simplify(f, r))
}

roots_quadratic <- function(f) {
  cf <- .pr_all_coeffs(f)
  a <- cf[[1L]]; b <- cf[[2L]]; c <- cf[[3L]]
  dom <- f$K
  .sqrt_sq <- function(d) {
    co <- list(); other <- list()
    for (di in (if (d$type == "Mul") .sympy_args(d) else list(d))) {
      if (di$type == "Pow" && di$args[[2L]]$type == "Integer" && .int_val(di$args[[2L]]) %% 2 == 0)
        co[[length(co) + 1L]] <- .pow(di$args[[1L]], .int(.int_val(di$args[[2L]]) %/% 2))
      else other[[length(other) + 1L]] <- di
    }
    if (length(co)) return(.mul(list(.mul(co), .sqrt(.mul(other)))))
    .sqrt(d)
  }
  simp <- function(e) if (.dk_composite(dom)) sym_factor(e) else simplify(e)
  if (is_zero_num(c)) {
    r0 <- S_Zero; r1 <- .neg(.div(b, a))
    if (!.dk_numerical(dom)) r1 <- simp(r1)
    else if (isTRUE(ask_fact(r1, "negative"))) { t <- r0; r0 <- r1; r1 <- t }
  } else if (is_zero_num(b)) {
    r <- .neg(.div(c, a))
    if (!.dk_numerical(dom)) r <- simp(r)
    R <- .sqrt_sq(r)
    r0 <- .neg(R); r1 <- R
  } else {
    d <- .sub(.sq(b), .mul(list(.int(4), a, c)))
    A <- .mul(list(S_Two, a))
    B <- .div(.neg(b), A)
    if (!.dk_numerical(dom)) { d <- simp(d); B <- simp(B) }
    D <- factor_terms(.div(.sqrt_sq(d), A))
    r0 <- .sub(B, D); r1 <- .add(list(B, D))
    if (isTRUE(ask_fact(a, "negative"))) { t <- r0; r0 <- r1; r1 <- t }
    else if (!.dk_numerical(dom)) { r0 <- expand_2arg(r0); r1 <- expand_2arg(r1) }
  }
  list(r0, r1)
}

.pr_roots_list <- function(cf, ...) roots(lapply(cf, as_sym), multiple = TRUE, ...)

roots_cubic <- function(f, trig = FALSE) {
  if (trig) {
    cf <- .pr_all_coeffs(f)
    a <- cf[[1L]]; b <- cf[[2L]]; c <- cf[[3L]]; d <- cf[[4L]]
    p <- .div(.sub(.mul(list(.int(3), a, c)), .sq(b)), .mul(list(.int(3), .sq(a))))
    q <- .div(.add(list(.mul(list(S_Two, .pow(b, .int(3)))), .mul(list(.int(-9), a, b, c)), .mul(list(.int(27), .sq(a), d)))), .mul(list(.int(27), .pow(a, .int(3)))))
    D <- .add(list(.mul(list(.int(18), a, b, c, d)), .mul(list(.int(-4), .pow(b, .int(3)), d)), .mul(list(.sq(b), .sq(c))),
                   .mul(list(.int(-4), a, .pow(c, .int(3)))), .mul(list(.int(-27), .sq(a), .sq(d)))))
    if (.is_true(Gt(D, 0))) {
      rv <- lapply(0:2, function(k) .mul(list(S_Two, .sqrt(.div(.neg(p), .int(3))),
        fcall("cos", list(.sub(.div(fcall("acos", list(.mul(list(.div(q, p), .sqrt(.div(.int(-3), p)), Rational(3, 2))))), .int(3)),
                               .mul(list(.int(k), S_Pi, Rational(2, 3)))))))))
      return(lapply(rv, function(i) .sub(i, .div(.div(b, .int(3)), a))))
    }
  }
  cf <- .pr_all_coeffs(.dp_monic(f))
  a <- cf[[2L]]; b <- cf[[3L]]; c <- cf[[4L]]
  if (is_zero_num(c)) {
    xx <- .pr_roots_list(list(S_One, a, b))
    return(list(xx[[1L]], S_Zero, xx[[2L]]))
  }
  p <- .sub(b, .div(.sq(a), .int(3)))
  q <- .add(list(c, .neg(.div(.mul(list(a, b)), .int(3))), .div(.mul(list(S_Two, .pow(a, .int(3)))), .int(27))))
  pon3 <- .div(p, .int(3)); aon3 <- .div(a, .int(3))
  u1 <- NULL
  if (is_zero_num(p)) {
    if (is_zero_num(q)) return(rep(list(.neg(aon3)), 3L))
    u1 <- if (isTRUE(ask_fact(q, "positive"))) .neg(root(q, 3)) else root(.neg(q), 3)
  } else if (is_zero_num(q)) {
    yy <- .pr_roots_list(list(S_One, S_Zero, p))
    return(lapply(list(yy[[1L]], S_Zero, yy[[2L]]), function(tmp) .sub(tmp, aon3)))
  } else if (isTRUE(ask_fact(q, "real")) && isTRUE(ask_fact(q, "negative"))) {
    u1 <- .neg(root(.add(list(.div(.neg(q), S_Two), .sqrt(.add(list(.div(.sq(q), .int(4)), .pow(pon3, .int(3))))))), 3))
  }
  coeff <- .div(.mul(list(S_I, .sqrt(.int(3)))), S_Two)
  if (is.null(u1)) {
    us <- list(S_One, .add(list(Rational(-1, 2), coeff)), .sub(Rational(-1, 2), coeff))
    bb <- a; cc <- b; dd <- c
    D0 <- .sub(.sq(bb), .mul(list(.int(3), cc)))
    D1 <- .add(list(.mul(list(S_Two, .pow(bb, .int(3)))), .mul(list(.int(-9), bb, cc)), .mul(list(.int(27), dd))))
    C <- root(.div(.add(list(D1, .sqrt(.sub(.sq(D1), .mul(list(.int(4), .pow(D0, .int(3)))))))), S_Two), 3)
    return(lapply(us, function(uk) .neg(.div(.add(list(bb, .mul(list(uk, C)), .div(.div(D0, C), uk))), .int(3)))))
  }
  u2 <- .mul(list(u1, .add(list(Rational(-1, 2), coeff))))
  u3 <- .mul(list(u1, .sub(Rational(-1, 2), coeff)))
  if (is_zero_num(p)) return(lapply(list(u1, u2, u3), function(u) .sub(u, aon3)))
  lapply(list(u1, u2, u3), function(u) .add(list(.neg(u), .div(pon3, u), .neg(aon3))))
}

.roots_quartic_euler <- function(p, q, r, a) {
  x <- Dummy("x")
  eq <- .add(list(.mul(list(.int(64), .pow(x, .int(3)))), .mul(list(.int(32), p, .sq(x))),
                  .mul(list(.sub(.mul(list(.int(4), .sq(p))), .mul(list(.int(16), r))), x)), .neg(.sq(q))))
  xs <- roots(.dp_from(eq, x), x, cubics = FALSE)
  xs <- Filter(function(sol) isTRUE(ask_fact(sol, "rational")) && isTRUE(ask_fact(sol, "nonzero")), unclass(xs)[seq_along(xs)])
  if (!length(xs)) return(NULL)
  R <- xs[[which.max(vapply(xs, function(v) as.double(gmp::as.bigq(v$value)), 0))]]
  c1 <- .sqrt(R)
  B <- .div(.mul(list(.neg(q), c1)), .mul(list(.int(4), R)))
  A <- .sub(.neg(R), .div(p, S_Two))
  c2 <- .sqrt(.add(list(A, B))); c3 <- .sqrt(.sub(A, B))
  list(.sub(.sub(c1, c2), a), .sub(.sub(.neg(c1), c3), a), .sub(.add(list(.neg(c1), c3)), a), .sub(.add(list(c1, c2)), a))
}

roots_quartic <- function(f) {
  cf <- .pr_all_coeffs(.dp_monic(f))
  a <- cf[[2L]]; b <- cf[[3L]]; c <- cf[[4L]]; d <- cf[[5L]]
  if (is_zero_num(d)) return(c(list(S_Zero), .pr_roots_list(list(S_One, a, b, c))))
  if (!is_zero_num(a) && identical(.sq(.div(c, a))$key, d$key)) {
    x <- f$x; m <- .div(c, a)
    g <- .dp_from(.add(list(.sq(x), .mul(list(a, x)), b, .mul(list(.int(-2), m)))), x)
    zz <- roots_quadratic(g)
    h1 <- .dp_from(.add(list(.sq(x), .neg(.mul(list(zz[[1L]], x))), m)), x)
    h2 <- .dp_from(.add(list(.sq(x), .neg(.mul(list(zz[[2L]], x))), m)), x)
    return(c(roots_quadratic(h1), roots_quadratic(h2)))
  }
  a2 <- .sq(a)
  e <- .sub(b, .div(.mul(list(.int(3), a2)), .int(8)))
  ff <- .mexpand(.add(list(c, .mul(list(a, .sub(.div(a2, .int(8)), .div(b, S_Two)))))))
  aon4 <- .div(a, .int(4))
  g <- .mexpand(.sub(d, .mul(list(aon4, .add(list(.mul(list(a, .sub(.div(.mul(list(.int(3), a2)), .int(64)), .div(b, .int(4))))), c))))))
  if (isTRUE(ask_fact(ff, "zero"))) {
    yy <- lapply(.pr_roots_list(list(S_One, e, g)), .sqrt)
    return(lapply(list(.neg(yy[[1L]]), .neg(yy[[2L]]), yy[[1L]], yy[[2L]]), function(tmp) .sub(tmp, aon4)))
  }
  if (isTRUE(ask_fact(g, "zero"))) {
    y <- c(list(S_Zero), .pr_roots_list(list(S_One, S_Zero, e, ff)))
    return(lapply(y, function(tmp) .sub(tmp, aon4)))
  }
  sols <- .roots_quartic_euler(e, ff, g, aon4)
  if (!is.null(sols)) return(sols)
  p <- .sub(.neg(.div(.sq(e), .int(12))), g)
  q <- .add(list(.neg(.div(.pow(e, .int(3)), .int(108))), .div(.mul(list(e, g)), .int(3)), .neg(.div(.sq(ff), .int(8)))))
  TH <- Rational(1, 3)
  ans <- function(y) {
    w <- .sqrt(.add(list(e, .mul(list(S_Two, y)))))
    arg1 <- .add(list(.mul(list(.int(3), e)), .mul(list(S_Two, y))))
    arg2 <- .div(.mul(list(S_Two, ff)), w)
    out <- list()
    for (sg in c(-1L, 1L)) {
      rt <- .sqrt(.neg(.add(list(arg1, .mul(list(.int(sg), arg2))))))
      for (tg in c(-1L, 1L)) out[[length(out) + 1L]] <- .sub(.div(.sub(.mul(list(.int(sg), w)), .mul(list(.int(tg), rt))), S_Two), aon4)
    }
    out
  }
  p <- .mexpand(p)
  y1 <- .sub(.mul(list(e, Rational(-5, 6))), .pow(q, TH))
  if (isTRUE(ask_fact(p, "zero"))) return(ans(y1))
  rt <- .sqrt(.add(list(.div(.sq(q), .int(4)), .div(.pow(p, .int(3)), .int(27)))))
  r <- .add(list(.div(.neg(q), S_Two), rt))
  u <- .pow(r, TH)
  y2 <- .add(list(.mul(list(e, Rational(-5, 6))), u, .neg(.div(.div(p, u), .int(3)))))
  if (identical(ask_fact(p, "zero"), FALSE)) return(ans(y2))
  Map(function(a1, a2) Piecewise(list(a1, Eq(p, 0)), list(a2, S_true)), ans(y1), ans(y2))
}

roots_binomial <- function(f) {
  n <- .dp_deg(f)
  a <- .pr_nth(f, n); b <- .pr_nth(f, 0L)
  base <- .neg(cancel(.div(b, a)))
  cm <- as_coeff_mul_split(sym_factor(base))
  coeff <- cm[[1L]]; rest <- cm[[2L]]
  minus <- num_sign(coeff) < 0
  terms <- list(if (minus) root(.neg(coeff), n) else root(coeff, n))
  pd <- .as_powers_list(rest)
  for (pr in pd) {
    be <- .as_base_exp(pr[[1L]])
    terms[[length(terms) + 1L]] <- .pow(root(be[[1L]], n), .mul(list(be[[2L]], pr[[2L]])))
  }
  alpha <- .mul(terms)
  if (!length(free_symbols(alpha))) alpha <- expand_complex(alpha)
  even <- n %% 2L == 0L
  zetas <- list()
  ez <- function(k, odd) fcall("exp", list(.div(.mul(list(.int(if (odd) 2L * k + 1L else 2L * k), S_I, S_Pi)), .int(n))))
  if (even && !minus) {
    zetas <- list(S_NegativeOne, S_One)
    for (k in rev(seq_len(n %/% 2L - 1L))) { z <- ez(k, FALSE); zetas <- c(zetas, list(fcall("conjugate", list(z)), z)) }
  } else if (even && minus) {
    for (k in rev(seq_len(n %/% 2L) - 1L)) { z <- ez(k, TRUE); zetas <- c(zetas, list(fcall("conjugate", list(z)), z)) }
  } else if (!even && !minus) {
    zetas <- list(S_One)
    for (k in rev(seq_len(n %/% 2L))) { z <- ez(k, FALSE); zetas <- c(zetas, list(fcall("conjugate", list(z)), z)) }
  } else {
    zetas <- list(S_NegativeOne)
    for (k in rev(seq_len(n %/% 2L) - 1L)) { z <- ez(k, TRUE); zetas <- c(zetas, list(fcall("conjugate", list(z)), z)) }
  }
  zetas <- lapply(zetas, expand_complex)
  lapply(zetas, function(zeta) expand(.mul(list(alpha, zeta)), power_base = FALSE))
}

# Expr.as_powers_dict, in the order of the arguments
.as_powers_list <- function(e) {
  out <- list(); keys <- character()
  add <- function(b, ex) {
    k <- b$key
    i <- match(k, keys)
    if (is.na(i)) { keys <<- c(keys, k); out[[length(out) + 1L]] <<- list(b, ex) }
    else out[[i]][[2L]] <<- .add(list(out[[i]][[2L]], ex))
  }
  for (t in (if (e$type == "Mul") .sympy_args(e) else list(e))) {
    be <- .as_base_exp(t)
    add(be[[1L]], be[[2L]])
  }
  out
}
.as_base_exp <- function(e) {
  if (e$type == "Pow") return(list(e$args[[1L]], e$args[[2L]]))
  if (e$type == "Function" && e$name == "exp") return(list(fcall("exp", list(S_One)), e$args[[1L]]))
  list(e, S_One)
}

.inv_totient_estimate <- function(m) {
  ds <- as.integer(as.character(divisors(m)))
  primes <- Filter(function(p) isTRUE(isprime(p)), ds + 1L)
  a <- 1; b <- 1
  for (p in primes) { a <- a * p; b <- b * (p - 1) }
  L <- m
  U <- ceiling(m * (a / b))
  P <- 2; p <- 2; primes <- integer()
  while (P <= U) { p <- as.integer(as.character(nextprime(p))); primes <- c(primes, p); P <- P * p }
  P <- P %/% p
  b <- 1
  for (p2 in primes[-length(primes)]) b <- b * (p2 - 1)
  U <- ceiling(m * (P / b))
  c(L, U)
}

roots_cyclotomic <- function(f, factor = FALSE) {
  LU <- .inv_totient_estimate(.dp_deg(f))
  fe <- .dp_expr(f)
  n <- NA
  for (k in LU[1L]:LU[2L]) {
    g <- cyclotomic_poly(k, f$x)
    if (identical(expand(g)$key, fe$key)) { n <- k; break }
  }
  if (is.na(n)) stop("failed to find index of a cyclotomic polynomial")
  h <- n %/% 2L
  ks <- Filter(function(i) .igcd2(i, n) == 1L, seq_len(n))
  key1 <- vapply(ks, function(v) if (v <= h) v else abs(v - n), 0)
  key2 <- vapply(ks, function(v) if (v <= h) -1 else 1, 0)
  ks <- ks[order(key1, key2)]
  d <- .div(.mul(list(S_Two, S_I, S_Pi)), .int(n))
  lapply(rev(ks), function(k) expand_complex(fcall("exp", list(.mul(list(.int(k), d))))))
}
.igcd2 <- function(a, b) { while (b) { t <- a %% b; a <- b; b <- t }; a }

# dup_cyclotomic_p
.dp_is_cyclotomic <- function(f, irreducible = FALSE) {
  K <- f$K
  if (.dk(K) == "QQ") {
    if (!all(vapply(f$c, function(z) z$type == "Integer", NA))) return(FALSE)
  } else if (.dk(K) != "ZZ") return(FALSE)
  cf <- vapply(f$c, function(z) as.numeric(gmp::as.bigq(z$value)), 0)
  .cyclotomic_p_int(cf, irreducible)
}
.cyclotomic_p_int <- function(cf, irreducible = FALSE) {
  strip <- function(v) { while (length(v) && v[1L] == 0) v <- v[-1L]; v }
  sqr <- function(v) { if (!length(v)) return(v); n <- length(v); out <- numeric(2L * n - 1L); for (i in seq_len(n)) for (j in seq_len(n)) out[i + j - 1L] <- out[i + j - 1L] + v[i] * v[j]; out }
  sub <- function(u, v) { lu <- length(u); lv <- length(v); if (lu < lv) u <- c(rep(0, lv - lu), u) else if (lv < lu) v <- c(rep(0, lu - lv), v); strip(u - v) }
  lc <- cf[1L]; tc <- cf[length(cf)]
  if (lc != 1 || (tc != -1 && tc != 1)) return(FALSE)
  if (!irreducible) {
    fl <- factor_list(.coeffs_to_expr(lapply(cf, function(v) .int(v)), Symbol("x")))
    if (!is_one_num(fl$coeff) || length(fl$factors) != 1L || fl$factors[[1L]][[2L]] != 1) return(FALSE)
  }
  n <- length(cf) - 1L
  rev_cf <- rev(cf)  # rev_cf[i + 1] is the coefficient of x^i
  g <- numeric(); h <- numeric()
  for (i in seq(n, 0, by = -2)) g <- c(rev_cf[i + 1L], g)
  if (n >= 1) for (i in seq(n - 1L, 0, by = -2)) h <- c(rev_cf[i + 1L], h)
  g <- sqr(strip(g)); h <- sqr(strip(h))
  F <- sub(g, if (length(h)) c(h, 0) else h)
  if (length(F) && F[1L] < 0) F <- -F
  if (identical(F, cf)) return(TRUE)
  mir <- cf * (-1)^(n - (seq_along(cf) - 1L))
  if (length(mir) && mir[1L] < 0) mir <- -mir
  if (identical(F, mir) && .cyclotomic_p_int(mir)) return(TRUE)
  Fe <- .coeffs_to_expr(lapply(F, function(v) .int(v)), Symbol("x"))
  G <- .coeff_list(expand(sqf_part(Fe)), Symbol("x"))
  Gv <- vapply(G, function(z) as.numeric(gmp::as.bigq(z$value)), 0)
  if (identical(sqr(Gv), F) && .cyclotomic_p_int(Gv)) return(TRUE)
  FALSE
}

# dup_decompose over a field
.dp_decompose <- function(f) {
  right <- function(f, s) {
    n <- .dp_deg(f); lc <- .dp_LC(f); K <- f$K
    fd <- list(); for (i in seq_along(f$c)) if (!is_zero_num(f$c[[i]])) fd[[as.character(n - i + 1L)]] <- f$c[[i]]
    g <- list(); g[[as.character(s)]] <- S_One
    r <- n %/% s
    for (i in seq_len(s - 1L)) {
      coeff <- S_Zero
      for (j in 0:(i - 1L)) {
        a <- fd[[as.character(n + j - i)]]; b <- g[[as.character(s - j)]]
        if (is.null(a) || is.null(b)) next
        coeff <- .dn(K, .add(list(coeff, .mul(list(.int(i - r * j), a, b)))))
      }
      g[[as.character(s - i)]] <- .dom_quo(K, coeff, .dn(K, .mul(list(.int(i * r), lc))))
    }
    cl <- lapply(s:0, function(k) { v <- g[[as.character(k)]]; if (is.null(v)) S_Zero else v })
    .dp_set(f, cl)
  }
  left <- function(f, h) {
    g <- list(); i <- 0L
    while (!.dp_zero(f)) {
      qr <- .dp_div(f, h)
      if (.dp_deg(qr[[2L]]) > 0) return(NULL)
      g[[as.character(i)]] <- .dp_LC(qr[[2L]])
      f <- qr[[1L]]; i <- i + 1L
    }
    if (!length(g)) return(.dp_set(f, list()))
    m <- max(as.integer(names(g)))
    .dp_set(f, lapply(m:0, function(k) { v <- g[[as.character(k)]]; if (is.null(v)) S_Zero else v }))
  }
  dec1 <- function(f) {
    df <- .dp_deg(f)
    if (df < 3L) return(NULL)
    for (s in 2:(df - 1L)) {
      if (df %% s != 0L) next
      h <- right(f, s)
      g <- left(f, h)
      if (!is.null(g)) return(list(g, h))
    }
    NULL
  }
  F <- list()
  repeat {
    r <- dec1(f)
    if (is.null(r)) break
    f <- r[[1L]]; F <- c(list(r[[2L]]), F)
  }
  c(list(f), F)
}

# to_rational_coeffs
.to_rational_coeffs <- function(f) {
  has_square_roots <- function(p) {
    has_sq <- FALSE
    for (y in Filter(function(z) !is_zero_num(z), p$c)) for (xx in (if (y$type == "Add") y$args else list(y))) {
      fs <- .as_powers_list(xx)
      r <- unlist(lapply(fs, function(pr) {
        b <- pr[[1L]]; wx <- pr[[2L]]
        if (!length(free_symbols(b)) && is_number_type(b) && wx$type %in% c("Integer", "Rational")) {
          q <- as.numeric(gmp::denominator(gmp::as.bigq(wx$value)))
          if (q >= 2) return(q)
        }
        NULL
      }))
      if (!length(r)) next
      if (min(r) == 2) has_sq <- TRUE
      if (max(r) > 2) return(FALSE)
    }
    has_sq
  }
  if (.dk(f$K) != "EX" || !has_square_roots(f)) return(NULL)
  f1 <- .dp_monic(f)
  n <- .dp_deg(f); lc <- .dp_LC(f)
  coeffs <- lapply(f1$c[-1L], simplify)
  if (length(coeffs) > 1L && !is_zero_num(coeffs[[length(coeffs) - 1L]])) {
    rescale1 <- simplify(.div(coeffs[[length(coeffs) - 1L]], coeffs[[length(coeffs)]]))
    coeffs1 <- list(); ok <- TRUE
    for (i in seq_along(coeffs)) {
      cx <- simplify(.mul(list(coeffs[[i]], .pow(rescale1, .int(i)))))
      if (!isTRUE(ask_fact(cx, "rational"))) { ok <- FALSE; break }
      coeffs1[[i]] <- cx
    }
    if (ok) {
      rescale_x <- simplify(.inv(rescale1))
      x <- f$x
      v <- c(list(.pow(x, .int(n))), lapply(seq_len(n), function(i) .mul(list(coeffs1[[i]], .pow(x, .int(n - i))))))
      return(list(lc, rescale_x, NULL, .dp_from(.add(v), x)))
    }
  }
  c <- simplify(coeffs[[1L]])
  if (c$type == "Add" && !isTRUE(ask_fact(c, "rational"))) {
    nonrat <- Filter(function(z) !isTRUE(ask_fact(z, "rational")), c$args)
    alpha <- .div(.neg(.add(nonrat)), .int(n))
    f2 <- .dp_from(subs(.dp_expr(f1), f1$x, .add(list(f1$x, alpha))), f1$x)
    return(list(NULL, NULL, alpha, f2))
  }
  NULL
}

# _integer_basis
.pr_integer_basis <- function(poly) {
  n <- .dp_deg(poly)
  idx <- which(!vapply(poly$c, is_zero_num, NA))
  monoms <- n - (idx - 1L)
  coeffs <- lapply(poly$c[idx], function(z) abs(gmp::numerator(gmp::as.bigq(z$value))))
  if (coeffs[[1L]] < coeffs[[length(coeffs)]]) {
    coeffs <- rev(coeffs)
    m0 <- monoms[1L]
    monoms <- m0 - rev(monoms)
  } else return(NULL)
  monoms <- monoms[-length(monoms)]
  coeffs <- coeffs[-length(coeffs)]
  if (length(monoms) == 1L) {
    r <- .pow(.zs(coeffs[[1L]]), .div(S_One, .int(monoms[1L])))
    if (r$type == "Integer") return(gmp::numerator(gmp::as.bigq(r$value)))
    return(NULL)
  }
  g <- Reduce(gmp::gcd.bigz, coeffs)
  divs <- rev(as.list(divisors(.zs(g))))
  divs <- divs[-length(divs)]
  for (dv in divs) {
    dz <- gmp::as.bigz(as.character(dv))
    if (all(vapply(seq_along(monoms), function(i) gmp::mod.bigz(coeffs[[i]], dz^monoms[i]) == 0, NA))) return(dz)
  }
  NULL
}

# preprocess_roots: list(coeff, poly)
.pr_preprocess_roots <- function(poly) {
  coeff <- S_One
  K <- poly$K
  if (.dk_field(K)) {
    if (.dk(K) %in% c("QQ", "QQ_I")) {
      den1 <- function(z) { ri <- .re_im(z); .dom_lcm(ZZ, .zs(gmp::denominator(gmp::as.bigq(ri[[1L]]$value))), .zs(gmp::denominator(gmp::as.bigq(ri[[2L]]$value)))) }
      common <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), lapply(poly$c, den1), S_One)
      poly <- .dp_convert(.dp_mul_ground(poly, common), .dom_get_ring(K))
    } else if (.dk(K) == "FracField") {
      common <- Reduce(function(a, b) .mpoly_lcm(a, b, .dk_dgens(K)), lapply(poly$c, function(z) expand(.as_numer_denom(z)[[2L]])), S_One)
      poly <- .dp_convert(.dp_mul_ground(poly, common), .dom_get_ring(K))
    } else if (.dk(K) == "EX") {
      common <- Reduce(function(a, b) lcm(a, b), lapply(poly$c, function(z) .as_numer_denom(z)[[2L]]), S_One)
      poly <- .dp_mul_ground(poly, common)
    }
  }
  poly <- .dp_primitive(poly)[[2L]]
  # retract
  poly <- .dp_from(.dp_expr(poly), poly$x, .dp_domain(poly$c))
  if (.dk(poly$K) == "PolyRing" && all(vapply(poly$c, function(c) is_zero_num(c) || .ri_is_term(c), NA))) {
    gens <- .dk_gens(poly$K); dgens <- .dk_dgens(poly$K)
    n <- .dp_deg(poly)
    nz <- which(!vapply(poly$c, is_zero_num, NA))
    base <- n - (nz - 1L)
    for (gi in seq_along(gens)) {
      gen <- gens[[gi]]; dgen <- dgens[[gi]]
      strip <- vapply(poly$c[nz], function(c) degree(c, dgen), 0)
      bb <- base
      rev_ <- FALSE
      if (strip[1L] < strip[length(strip)]) { strip <- rev(strip); rev_ <- TRUE }
      ratio <- NULL; ok <- TRUE
      for (k in seq_along(bb)) {
        a <- bb[k]; b <- strip[k]
        if (!a && !b) next
        if (!a || !b || b %% a != 0) { ok <- FALSE; break }
        rr <- b %/% a
        if (is.null(ratio)) ratio <- rr else if (ratio != rr) { ok <- FALSE; break }
      }
      if (ok && !is.null(ratio)) {
        if (rev_) ratio <- -ratio
        poly <- .dp_from(.gd_out(subs(.dp_expr_raw(poly), dgen, S_One)), poly$x)
        coeff <- .mul(list(coeff, .pow(gen, .int(-ratio))))
      }
    }
    poly <- .dp_from(.dp_expr(poly), poly$x)
  }
  if (.dk(poly$K) == "ZZ") {
    basis <- .pr_integer_basis(poly)
    if (!is.null(basis)) {
      n <- .dp_deg(poly)
      poly <- .dp_set(poly, lapply(seq_along(poly$c), function(i) {
        k <- n - i + 1L
        .zs(gmp::numerator(gmp::as.bigq(poly$c[[i]]$value)) %/% (basis^(n - k)))
      }))
      coeff <- .mul(list(coeff, .zs(basis)))
    }
  }
  list(coeff, poly)
}
.ri_is_term <- function(c) !(expand(c)$type == "Add")

.roots_sympy <- function(f, x = NULL, auto = TRUE, cubics = TRUE, trig = FALSE, quartics = TRUE, quintics = FALSE,
                         multiple = FALSE, filter = NULL, predicate = NULL, strict = FALSE) {
  if (is.list(f) && !is_sym(f) && !inherits(f, "dpoly")) {
    x <- Dummy("x")
    cl <- lapply(f, as_sym)
    f <- dpoly(lapply(cl, function(z) .dn(.dp_domain(cl, field = TRUE), z)), x, .dp_domain(cl, field = TRUE))
  } else if (inherits(f, "dpoly")) {
    x <- f$x
  } else {
    if (inherits(f, "Poly")) { if (is.null(x)) x <- f$gens[[1L]]; f <- as_expr(f) }
    f <- as_sym(f)
    if (is.null(x)) {
      fs <- free_symbols(f)
      if (!length(fs)) return(if (multiple) list() else structure(list(), multiplicity = integer()))
      if (length(fs) != 1L) stop("multivariate polynomials are not supported", call. = FALSE)
      x <- fs[[1L]]
    }
    x <- as_sym(x)
    fe <- f
    f <- .dp_from(fe, x)
    n <- .dp_deg(f)
    if (.pr_length(f) == 2L && n > 2L) {
      ex <- expand(fe)
      terms <- if (ex$type == "Add") ex$args else list(ex)
      con <- .add(Filter(function(z) !.has1(z, x), terms)); dep <- .add(Filter(function(z) .has1(z, x), terms))
      fcon <- .neg(sym_factor(.neg(con)))
      if (!identical(fcon$key, con$key)) {
        con <- fcon
        bases <- list()
        for (i in (if (con$type == "Mul") .sympy_args(con) else list(con))) {
          if (i$type == "Pow") {
            be <- .as_base_exp(i)
            if (be[[2L]]$type == "Integer" && be[[1L]]$type == "Add") bases[[length(bases) + 1L]] <- list(be[[1L]], Dummy("Dummy", positive = TRUE))
          }
        }
        if (length(bases)) {
          rv <- .roots_sympy(subs(.add(list(dep, con)), bases), x, auto, cubics, trig, quartics, quintics, FALSE, filter, predicate)
          back <- lapply(bases, function(b) list(b[[2L]], b[[1L]]))
          out <- lapply(unclass(rv)[seq_along(rv)], function(k) factor_terms(subs(k, back)))
          mul <- attr(rv, "multiplicity")
          if (multiple) return(unlist(Map(function(r, m) rep(list(r), m), out, mul), recursive = FALSE))
          return(structure(out, multiplicity = mul))
        }
      }
    }
  }
  res_keys <- character(); res_vals <- list(); res_mult <- integer()
  upd <- function(r, k) {
    i <- match(r$key, res_keys)
    if (is.na(i)) { res_keys <<- c(res_keys, r$key); res_vals[[length(res_vals) + 1L]] <<- r; res_mult <<- c(res_mult, as.integer(k)) }
    else res_mult[i] <<- res_mult[i] + as.integer(k)
  }
  d0 <- Dummy("x0")
  orig_x <- x
  f <- .dp_set(dpoly(f$c, d0, f$K), f$c)
  # terms_gcd
  k0 <- 0L
  while (length(f$c) > 1L && is_zero_num(f$c[[length(f$c)]])) { f$c <- f$c[-length(f$c)]; k0 <- k0 + 1L }
  zeros <- if (k0) list(list(S_Zero, k0)) else list()
  if (.dp_zero(f)) return(if (multiple) list() else structure(list(), multiplicity = integer()))
  pp <- .pr_preprocess_roots(f)
  coeff <- pp[[1L]]; f <- pp[[2L]]
  if (auto && !.dk_field(f$K)) f <- .dp_to_field(f)
  if (.dk(f$K) == "QQ_I") f <- .dp_convert(f, EX)
  rescale_x <- NULL; translate_x <- NULL
  try_heuristics <- function(f) {
    if (.dp_deg(f) <= 0) return(list())
    if (.pr_length(f) == 1L) return(rep(list(S_Zero), .dp_deg(f)))
    if (.pr_length(f) == 2L) {
      if (.dp_deg(f) == 1L) return(lapply(roots_linear(f), cancel))
      return(roots_binomial(f))
    }
    result <- list()
    for (i in c(-1L, 1L)) {
      if (is_zero_num(.dp_eval(f, .int(i)))) {
        f <- .dp_quo(f, .dp_from(.sub(f$x, .int(i)), f$x, f$K))
        result <- c(result, list(.int(i)))
        break
      }
    }
    n <- .dp_deg(f)
    if (n == 1L) return(c(result, lapply(roots_linear(f), cancel)))
    if (n == 2L) return(c(result, lapply(roots_quadratic(f), cancel)))
    if (.dp_is_cyclotomic(f)) return(c(result, roots_cyclotomic(f)))
    if (n == 3L && isTRUE(cubics)) return(c(result, roots_cubic(f, trig = trig)))
    if (n == 4L && isTRUE(quartics)) return(c(result, roots_quartic(f)))
    if (n == 5L && isTRUE(quintics)) return(c(result, roots_quintic(f)))
    result
  }
  try_decompose <- function(f) {
    facs <- .dp_decompose(f)
    rts <- try_heuristics(facs[[1L]])
    for (cf in facs[-1L]) {
      previous <- rts; rts <- list()
      for (r in previous) {
        ge <- .sub(.dp_expr(cf), r)
        g <- .dp_from(ge, cf$x, .dp_domain(.coeff_list(expand(ge), cf$x), field = TRUE))
        rts <- c(rts, try_heuristics(g))
      }
    }
    rts
  }
  if (.dp_deg(f) > 0) {
    dom <- f$K
    if (.dk(dom) %in% c("RR", "CC")) {
      for (r in nroots(.dp_expr(f), f$x)) upd(r, 1L)
    } else if (.dp_deg(f) == 1L) {
      upd(roots_linear(f)[[1L]], 1L)
    } else if (.pr_length(f) == 2L) {
      for (r in (if (.dp_deg(f) == 2L) roots_quadratic(f) else roots_binomial(f))) upd(r, 1L)
    } else {
      fl <- factor_list(.dp_expr(f))
      factors <- fl$factors
      if (length(factors) == 1L && .dp_deg(f) == 2L) {
        for (r in roots_quadratic(f)) upd(r, 1L)
      } else if (length(factors) == 1L && factors[[1L]][[2L]] == 1) {
        if (.dk(f$K) == "EX") {
          rr <- .to_rational_coeffs(f)
          if (!is.null(rr)) {
            if (is.null(rr[[1L]])) { translate_x <- rr[[3L]]; f <- rr[[4L]] }
            else { rescale_x <- rr[[2L]]; f <- rr[[4L]] }
            sub <- .roots_sympy(f)
            if (length(sub)) { for (j in seq_along(sub)) upd(sub[[j]], attr(sub, "multiplicity")[j]) }
            else for (r in try_decompose(f)) upd(r, 1L)
          } else for (r in try_heuristics(f)) upd(r, 1L)
        } else for (r in try_decompose(f)) upd(r, 1L)
      } else {
        for (fk in factors) {
          g <- .dp_from(fk[[1L]], f$x, field = TRUE)
          for (r in try_heuristics(g)) upd(r, fk[[2L]])
        }
      }
    }
  }
  if (!is_one_num(coeff)) {
    old <- list(res_vals, res_mult); res_keys <- character(); res_vals <- list(); res_mult <- integer()
    for (j in seq_along(old[[1L]])) upd(.mul(list(coeff, old[[1L]][[j]])), old[[2L]][j])
  }
  keepf <- function(test) {
    keep <- vapply(res_vals, test, NA)
    res_vals <<- res_vals[keep]; res_mult <<- res_mult[keep]; res_keys <<- res_keys[keep]
  }
  if (!is.null(filter) && !identical(filter, "C")) {
    test <- switch(filter,
      Z = function(r) r$type == "Integer",
      Q = function(r) r$type %in% c("Integer", "Rational"),
      R = function(r) all(vapply(.as_numer_denom(r), function(a) isTRUE(ask_fact(a, "real")), NA)),
      I = function(r) isTRUE(ask_fact(r, "imaginary")),
      real = function(r) isTRUE(ask_fact(r, "real")),
      rational = function(r) isTRUE(ask_fact(r, "rational")),
      complex = function(r) TRUE,
      stop("Invalid filter: ", filter, call. = FALSE))
    keepf(test)
  }
  if (!is.null(predicate)) keepf(predicate)
  if (!is.null(rescale_x)) { res_vals <- lapply(res_vals, function(k) .mul(list(k, rescale_x))); res_keys <- vapply(res_vals, `[[`, "", "key") }
  if (!is.null(translate_x)) { res_vals <- lapply(res_vals, function(k) .add(list(k, translate_x))); res_keys <- vapply(res_vals, `[[`, "", "key") }
  for (z in zeros) {
    i <- match(z[[1L]]$key, res_keys)
    if (is.na(i)) { res_keys <- c(res_keys, z[[1L]]$key); res_vals[[length(res_vals) + 1L]] <- z[[1L]]; res_mult <- c(res_mult, z[[2L]]) }
    else res_mult[i] <- z[[2L]]
  }
  # the dummy generator back to x
  res_vals <- lapply(res_vals, function(r) subs(r, d0, orig_x))
  if (strict && sum(res_mult) < .dp_deg(f)) stop("UnsolvableFactorError", call. = FALSE)
  if (!multiple) return(structure(res_vals, multiplicity = res_mult))
  o <- .bool_print_order(res_vals)
  idx <- match(vapply(o, `[[`, "", "key"), vapply(res_vals, `[[`, "", "key"))
  unlist(lapply(idx, function(i) rep(list(res_vals[[i]]), res_mult[i])), recursive = FALSE)
}
