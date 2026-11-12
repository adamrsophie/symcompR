# Polynomial options of SymPy's polytools functions: modulus = p (finite
# fields, after sympy/polys/galoistools.py), and extension = a, gaussian =
# TRUE or domain = QQ<a> (algebraic number fields, after dup_ext_factor in
# sympy/polys/factortools.py and the ANP arithmetic of polyclasses.py).

.POLY_OPTS <- c("modulus", "extension", "gaussian", "domain", "symmetric", "greedy", "expand",
                "frac", "polys", "strict", "auto", "order", "wrt", "sort", "field", "composite",
                "include", "all", "split", "fraction")

# Split ... into generators and named options.
.split_poly_opts <- function(dots) {
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  is_opt <- nm %in% .POLY_OPTS
  list(gens = unname(dots[!is_opt]), opts = dots[is_opt])
}

# R/polys-options.R
# The ground field asked for: list(kind = "QQ"), list(kind = "GF", p, symmetric)
# or list(kind = "NF", exts); the Gaussian domains add gauss = "ZZ_I" or
# "QQ_I" (gaussian = TRUE is QQ_I, extension = I the algebraic field QQ<I>).
.poly_field <- function(opts, exprs = list()) {
  if (!is.null(opts$modulus)) {
    p <- gmp::as.bigz(as.character(as_sym(opts$modulus)$value))
    return(.gf_field(p, !isFALSE(opts$symmetric)))
  }
  d <- opts$domain
  if (!is.null(d)) {
    if (is.character(d)) {
      if (d %in% c("ZZ_I", "QQ_I")) return(list(kind = "NF", exts = list(S_I), gauss = d))
      m <- regmatches(d, regexec("^(GF|FF)\\((\\d+)\\)$", d))[[1L]]
      if (length(m)) return(.gf_field(gmp::as.bigz(m[3L]), !isFALSE(opts$symmetric)))
      if (d %in% c("ZZ", "QQ")) return(list(kind = "QQ"))
      m <- regmatches(d, regexec("^QQ<(.*)>$", d))[[1L]]
      if (length(m)) return(list(kind = "NF", exts = lapply(strsplit(m[2L], ",", fixed = TRUE)[[1L]], sympify)))
      stop("unsupported domain ", d, call. = FALSE)
    }
    if (inherits(d, "sym_domain")) {
      dd <- attr(d, "dom")
      if (dd$kind %in% c("ZZ", "QQ")) return(list(kind = "QQ"))
      if (dd$kind == "GF") return(.gf_field(gmp::as.bigz(dd$mod), !isFALSE(opts$symmetric)))
      if (dd$kind %in% c("ZZ_I", "QQ_I")) return(list(kind = "NF", exts = list(S_I), gauss = dd$kind))
      if (dd$kind == "Algebraic") return(list(kind = "NF", exts = dd$ext))
      stop("unsupported domain ", format(d), call. = FALSE)
    }
  }
  # Gaussian.postprocess: gaussian = True is the domain QQ_I
  if (isTRUE(opts$gaussian)) return(list(kind = "NF", exts = list(S_I), gauss = "QQ_I"))
  ext <- opts$extension
  if (isTRUE(ext)) {
    # extension = TRUE: the algebraic numbers among the coefficients
    al <- .algebraic_atoms(exprs)
    if (!length(al)) return(list(kind = "QQ"))
    return(list(kind = "NF", exts = al))
  }
  if (!is.null(ext) && !isFALSE(ext)) {
    exts <- if (is_sym(ext)) list(ext) else lapply(as.list(ext), as_sym)
    exts <- Filter(function(a) !is_exact_number(a), exts)
    if (!length(exts)) return(list(kind = "QQ"))
    return(list(kind = "NF", exts = exts))
  }
  list(kind = "QQ")
}

.gf_field <- function(p, symmetric = TRUE) {
  if (p < 2 || !gmp::isprime(p)) stop("modulus must be a positive prime number, got ", as.character(p), call. = FALSE)
  if (p >= 2^25) stop("moduli at or above 2^25 are not supported", call. = FALSE)
  list(kind = "GF", p = as.numeric(p), symmetric = symmetric)
}

# Algebraic irrationals (I and rational powers of rationals) in expressions.
.algebraic_atoms <- function(exprs) {
  out <- list()
  walk <- function(z) {
    if (is_const(z, "I")) { out[[z$key]] <<- z; return(invisible()) }
    if (z$type == "Pow" && is_exact_number(z$args[[1L]]) && z$args[[2L]]$type == "Rational") {
      out[[z$key]] <<- z; return(invisible())
    }
    for (a in z$args) walk(a)
  }
  for (e in exprs) walk(expand(as_sym(e)))
  unname(out)
}

# ---------------------------------------------------------------------------
# Finite fields

# gf_sqf_list: square-free decomposition of f over GF(p), f a double vector
# (highest degree first); returns list(lc, factors = list(list(g, k))).
.gf_sqf_list <- function(f, p) {
  f <- .gf_trim(f %% p)
  if (length(f) <= 1L) return(list(lc = if (length(f)) f[1L] else 0, factors = list()))
  lc <- f[1L]; f <- .gf_monic(f, p)
  n <- 1; sqf <- FALSE; factors <- list()
  repeat {
    Fd <- .gf_diff(f, p)
    if (length(Fd)) {
      g <- .gf_gcd(f, Fd, p)
      h <- .gf_divmod(f, g, p)$q
      i <- 1
      while (!(length(h) == 1L && h[1L] == 1)) {
        G <- .gf_gcd(g, h, p)
        H <- .gf_divmod(h, G, p)$q
        if (length(H) > 1L) factors[[length(factors) + 1L]] <- list(.gf_monic(H, p), i * n)
        g <- .gf_divmod(g, G, p)$q; h <- G; i <- i + 1
      }
      if (length(g) == 1L && g[1L] == 1) sqf <- TRUE else f <- g
    }
    if (sqf) break
    # f = g(x^p): take the p-th root
    d <- (length(f) - 1L) %/% p
    f <- f[(0:d) * p + 1L]
    n <- n * p
  }
  list(lc = lc, factors = factors)
}

# gf_edf_zassenhaus for p = 2: the trace map replaces r^((p^n - 1)/2).
.gf_edf2 <- function(f, n) {
  # gf_edf_zassenhaus for p = 2: the trace map of t = x, x^3, x^5, ...
  factors <- list(f)
  d <- length(f) - 1L
  if (d <= n) return(factors)
  N <- d %/% n
  t <- c(1, 0)
  while (length(factors) < N) {
    h <- t; r <- t
    for (i in seq_len(n - 1L)) { r <- .gf_powmod(r, 2, f, 2); h <- .gf_add(h, r, 2) }
    g <- .gf_gcd(f, h, 2)
    t <- c(t, 0, 0)
    if (length(g) > 1L && !(length(g) == length(f) && all(g == f)))
      factors <- c(.gf_edf2(g, n), .gf_edf2(.gf_divmod(f, g, 2)$q, n))
  }
  factors
}

# gf_factor: lc and the monic irreducible factors with multiplicities, in
# _sort_factors order (degree, multiplicity, coefficients in [0, p)).
.gf_factor <- function(f, p) {
  s <- .gf_sqf_list(f, p)
  out <- list()
  for (gk in s$factors) {
    g <- gk[[1L]]
    irr <- if (p == 2) unlist(lapply(.gf_ddf(g, 2), function(dd) .gf_edf2(dd[[1L]], dd[[2L]])), recursive = FALSE)
           else .gf_factor_sqf(g, p)
    for (h in irr) out[[length(out) + 1L]] <- list(.gf_monic(h, p), gk[[2L]])
  }
  list(lc = s$lc, factors = .gf_sort(out))
}

.gf_sort <- function(fs) {
  if (length(fs) < 2L) return(fs)
  key <- function(a, b) {
    fa <- a[[1L]]; fb <- b[[1L]]
    if (length(fa) != length(fb)) return(length(fa) - length(fb))
    if (a[[2L]] != b[[2L]]) return(a[[2L]] - b[[2L]])
    d <- which(fa != fb)
    if (!length(d)) 0 else fa[d[1L]] - fb[d[1L]]
  }
  # insertion sort with the comparator (lists are short)
  out <- list()
  for (f in fs) {
    i <- length(out) + 1L
    while (i > 1L && key(out[[i - 1L]], f) > 0) i <- i - 1L
    out <- append(out, list(f), i - 1L)
  }
  out
}

# An integer of GF(p) as SymPy prints it: symmetric in (-p/2, p/2].
.gf_int <- function(a, p, symmetric = TRUE) if (symmetric && a > p %/% 2) a - p else a

.gf_to_expr <- function(f, x, p, symmetric = TRUE) {
  n <- length(f) - 1L
  if (n < 0L) return(S_Zero)
  .add(lapply(seq_along(f), function(i) .mul(list(.int(.gf_int(f[i], p, symmetric)), .pow(x, .int(n - i + 1L))))))
}

# An expression as a polynomial in x over GF(p); rationals map to p-adic
# residues (1/2 is the inverse of 2).
.gf_from_expr <- function(e, x, p) {
  cl <- .coeff_list(expand(e), x)
  vapply(cl, function(c) {
    if (!(c$type %in% c("Integer", "Rational"))) stop("expression is not a polynomial over GF(", p, "): ", .str(e), call. = FALSE)
    q <- gmp::as.bigq(c$value)
    nu <- as.numeric(gmp::mod.bigz(gmp::numerator(q), p)); de <- as.numeric(gmp::mod.bigz(gmp::denominator(q), p))
    (nu * .gf_inv(de, p)) %% p
  }, 0)
}

.gf_univariate_gen <- function(exprs, gens) {
  g <- .gens_of(exprs, gens)
  if (length(g) > 1L) stop("multivariate polynomials over finite fields are not supported", call. = FALSE)
  if (!length(g)) return(NULL)
  g[[1L]]
}

# ---------------------------------------------------------------------------
# Algebraic number fields Q(theta). Elements are bigq vectors (highest power
# first, as ANP) reduced modulo the monic minimal polynomial of theta.

.nf_field <- function(exts) {
  exts <- lapply(exts, as_sym)
  t <- Dummy("t")
  if (length(exts) == 1L) {
    theta <- exts[[1L]]
  } else {
    pe <- primitive_element(exts, t)
    theta <- expand(.add(Map(function(c, a) .mul(list(as_sym(c), a)), pe[[2L]], exts)))
  }
  mp <- minimal_polynomial(theta, t)
  m <- up_monic(up_from_qpoly(expr_to_qpoly(mp, list(t), strict = TRUE), 1L))
  list(theta = theta, mod = m, d = length(m) - 1L, exts = exts,
       pow = lapply(seq_len(length(m) - 1L) - 1L, function(i) expand(.pow(theta, .int(i)))))
}

.nfe_red <- function(a, K) up_divmod(up_trim(a), K$mod)$r
.nfe_mul <- function(a, b, K) if (!length(a) || !length(b)) .bq0() else .nfe_red(up_mul(a, b), K)
.nfe_inv <- function(a, K) {
  # extended Euclid in Q[t] modulo the minimal polynomial
  r0 <- K$mod; r1 <- up_trim(a); s0 <- .bq0(); s1 <- gmp::as.bigq(1)
  while (length(r1)) {
    d <- up_divmod(r0, r1)
    r0 <- r1; r1 <- d$r
    tmp <- up_sub(s0, up_mul(d$q, s1)); s0 <- s1; s1 <- tmp
  }
  if (length(r0) != 1L) stop("element is not invertible", call. = FALSE)
  .nfe_red(s0 / r0[1L], K)
}
.nfe_const <- function(q) up_trim(gmp::as.bigq(q))

.nfe_from_expr <- function(c, K) {
  c <- expand(as_sym(c))
  if (is_exact_number(c)) return(.nfe_const(c$value))
  y <- Dummy("y")
  p <- .to_field_poly(c, K$theta, y, K$d)
  if (is.null(p)) stop(.str(c), " is not in Q(", .str(K$theta), ")", call. = FALSE)
  .nfe_red(up_from_qpoly(expr_to_qpoly(expand(p), list(y), strict = TRUE), 1L), K)
}

.nfe_to_expr <- function(a, K) {
  a <- up_trim(a)
  if (!length(a)) return(S_Zero)
  n <- length(a)
  expand(.add(lapply(seq_len(n), function(i) .mul(list(make_rational(a[i]), K$pow[[n - i + 1L]])))))
}

# Polynomials over K: lists of elements, highest degree first.
.nfp_trim <- function(f) { while (length(f) && !length(up_trim(f[[1L]]))) f <- f[-1L]; f }
.nfp_monic <- function(f, K) { f <- .nfp_trim(f); if (!length(f)) return(f); li <- .nfe_inv(f[[1L]], K); lapply(f, .nfe_mul, b = li, K = K) }
.nfp_add <- function(f, g) {
  n <- max(length(f), length(g))
  f <- c(rep(list(.bq0()), n - length(f)), f); g <- c(rep(list(.bq0()), n - length(g)), g)
  .nfp_trim(Map(function(a, b) up_add(a, b), f, g))
}
.nfp_neg <- function(f) lapply(f, function(a) -a)
.nfp_scale <- function(f, c, K) .nfp_trim(lapply(f, .nfe_mul, b = c, K = K))
.nfp_mul <- function(f, g, K) {
  if (!length(f) || !length(g)) return(list())
  out <- rep(list(.bq0()), length(f) + length(g) - 1L)
  for (i in seq_along(f)) for (j in seq_along(g)) out[[i + j - 1L]] <- up_add(out[[i + j - 1L]], .nfe_mul(f[[i]], g[[j]], K))
  .nfp_trim(out)
}
.nfp_divmod <- function(f, g, K) {
  g <- .nfp_trim(g); f <- .nfp_trim(f)
  if (!length(g)) stop("polynomial division by zero", call. = FALSE)
  dg <- length(g) - 1L
  if (length(f) - 1L < dg) return(list(q = list(), r = f))
  li <- .nfe_inv(g[[1L]], K)
  q <- rep(list(.bq0()), length(f) - dg)
  r <- f
  for (i in seq_along(q)) {
    cf <- .nfe_mul(r[[i]], li, K)
    q[[i]] <- cf
    if (length(cf)) for (j in 0:dg) r[[i + j]] <- up_sub(r[[i + j]], .nfe_mul(g[[j + 1L]], cf, K))
  }
  rr <- if (dg > 0) r[(length(r) - dg + 1L):length(r)] else list()
  list(q = .nfp_trim(q), r = .nfp_trim(rr))
}
.nfp_gcd <- function(f, g, K) {
  f <- .nfp_trim(f); g <- .nfp_trim(g)
  while (length(g)) { r <- .nfp_divmod(f, g, K)$r; f <- g; g <- r }
  .nfp_monic(f, K)
}
.nfp_diff <- function(f) {
  n <- length(f) - 1L
  if (n <= 0L) return(list())
  .nfp_trim(lapply(seq_len(n), function(i) f[[i]] * (n - i + 1L)))
}
# f(x + a): Taylor shift by Horner steps
.nfp_shift <- function(f, a, K) {
  n <- length(f)
  if (n < 2L) return(f)
  for (i in seq_len(n - 1L)) for (j in 2:(n - i + 1L)) f[[j]] <- up_add(f[[j]], .nfe_mul(f[[j - 1L]], a, K))
  f
}
.nfp_from_expr <- function(e, x, K) lapply(.coeff_list(expand(e), x), .nfe_from_expr, K = K)
.nfp_to_expr <- function(f, x, K) {
  n <- length(f) - 1L
  if (n < 0L) return(S_Zero)
  .add(lapply(seq_along(f), function(i) .mul(list(.nfe_to_expr(f[[i]], K), .pow(x, .int(n - i + 1L))))))
}

# dup_sqf_norm: the smallest s >= 0 with Res_t(m(t), f(x - s*t, t))
# square-free; returns s, the shifted f and the norm (bigq vector).
.nfp_sqf_norm <- function(f, K) {
  t <- Dummy("t"); x <- Dummy("x")
  mt <- up_to_expr(K$mod, t)
  s <- 0L
  repeat {
    Fx <- .add(lapply(seq_along(f), function(i) .mul(list(up_to_expr(f[[i]], t), .pow(x, .int(length(f) - i))))))
    r <- up_from_qpoly(expr_to_qpoly(expand(resultant(mt, Fx, t)), list(x), strict = TRUE), 1L)
    if (up_deg(up_gcd(r, up_diff(r))) <= 0) return(list(s = s, f = f, r = r))
    f <- .nfp_shift(f, gmp::as.bigq(c(-1, 0)), K)
    s <- s + 1L
  }
}

# dup_ext_factor: lc and the monic irreducible factors over K, each with its
# multiplicity, for f a polynomial over K.
.nf_factor <- function(f, K) {
  f <- .nfp_trim(f)
  n <- length(f) - 1L
  if (n <= 0L) return(list(lc = if (length(f)) f[[1L]] else .bq0(), factors = list()))
  lc <- f[[1L]]
  f <- .nfp_monic(f, K)
  if (n == 1L) return(list(lc = lc, factors = list(list(f, 1L))))
  Ff <- f
  g0 <- .nfp_gcd(f, .nfp_diff(f), K)
  if (length(g0) > 1L) f <- .nfp_monic(.nfp_divmod(f, g0, K)$q, K)
  sn <- .nfp_sqf_norm(f, K)
  rx <- Dummy("x")
  fl <- factor_list(up_to_expr(sn$r, rx), rx)$factors
  if (length(fl) == 1L) {
    facs <- list(f)
  } else {
    H <- .nfe_red(gmp::as.bigq(c(sn$s, 0)), K)
    g <- sn$f
    facs <- lapply(fl, function(fk) {
      h <- lapply(.coeff_list(expand(fk[[1L]]), rx), function(c) .nfe_const(c$value))
      h <- .nfp_gcd(h, g, K)
      g <<- .nfp_monic(.nfp_divmod(g, h, K)$q, K)
      .nfp_shift(h, H, K)
    })
  }
  # dup_trial_division: multiplicities in the original polynomial
  out <- lapply(facs, function(h) {
    k <- 0L; q <- Ff
    repeat { d <- .nfp_divmod(q, h, K); if (length(d$r)) break; q <- d$q; k <- k + 1L }
    list(h, k)
  })
  list(lc = lc, factors = out)
}

# _sort_factors order over K: degree, multiplicity, then the coefficients'
# ANP representations compared as lists.
.nf_sort <- function(fs) {
  cmp_q <- function(a, b) {
    for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(if (a[i] < b[i]) -1 else 1)
    sign(length(a) - length(b))
  }
  key <- function(a, b) {
    fa <- a[[1L]]; fb <- b[[1L]]
    if (length(fa) != length(fb)) return(length(fa) - length(fb))
    if (a[[2L]] != b[[2L]]) return(a[[2L]] - b[[2L]])
    for (i in seq_along(fa)) { r <- cmp_q(up_trim(fa[[i]]), up_trim(fb[[i]])); if (r != 0) return(r) }
    0
  }
  out <- list()
  for (f in fs) {
    i <- length(out) + 1L
    while (i > 1L && key(out[[i - 1L]], f) > 0) i <- i - 1L
    out <- append(out, list(f), i - 1L)
  }
  out
}

# dup_sqf_list over K (Yun's algorithm with monic gcds).
.nf_sqf_list <- function(f, K) {
  f <- .nfp_trim(f)
  if (length(f) <= 1L) return(list(lc = if (length(f)) f[[1L]] else .bq0(), factors = list()))
  lc <- f[[1L]]; f <- .nfp_monic(f, K)
  out <- list(); i <- 1L
  h <- .nfp_diff(f)
  g <- .nfp_gcd(f, h, K)
  p <- .nfp_divmod(f, g, K)$q
  q <- .nfp_divmod(h, g, K)$q
  repeat {
    d <- .nfp_add(q, .nfp_neg(.nfp_diff(p)))
    if (!length(d)) { out[[length(out) + 1L]] <- list(.nfp_monic(p, K), i); break }
    h <- .nfp_gcd(p, d, K)
    p <- .nfp_divmod(p, h, K)$q
    q <- .nfp_divmod(d, h, K)$q
    if (length(h) > 1L) out[[length(out) + 1L]] <- list(h, i)
    i <- i + 1L
  }
  list(lc = lc, factors = out)
}

.sort_by_method <- function(fs, method, sorter) {
  if (method != "sqf" || length(fs) < 2L) return(sorter(fs))
  ks <- vapply(fs, function(f) as.numeric(f[[2L]]), 0)
  unlist(lapply(sort(unique(ks)), function(k) sorter(fs[ks == k])), recursive = FALSE)
}

# factor_list or sqf_list of a polynomial expression over a finite or
# algebraic field: list(coeff, factors = list(list(expr, k))).
.field_factor_list <- function(e, Fd, gens, method = "factor") {
  e <- expand(as_sym(e))
  if (Fd$kind == "GF") {
    x <- .gf_univariate_gen(list(e), gens)
    if (is.null(x)) return(list(coeff = e, factors = list()))
    p <- Fd$p
    v <- .gf_from_expr(e, x, p)
    r <- if (method == "sqf") .gf_sqf_list(v, p) else .gf_factor(v, p)
    fs <- .sort_by_method(r$factors, method, .gf_sort)
    return(list(coeff = .int(.gf_int(r$lc, p, Fd$symmetric)),
                factors = lapply(fs, function(fk) list(.gf_to_expr(fk[[1L]], x, p, Fd$symmetric), .int(fk[[2L]])))))
  }
  K <- .kd_field(Fd$exts)
  gs <- .kd_gens(e, gens)
  if (!length(gs)) return(list(coeff = e, factors = list()))
  if (method == "sqf") {
    if (length(gs) > 1L) stop("multivariate square-free decomposition over algebraic fields is not supported", call. = FALSE)
    x <- gs[[1L]]
    f <- .nfp_from_expr(e, x, K)
    # dup_terms_gcd: the power of x first
    j <- 0L
    while (length(f) > 1L && !length(up_trim(f[[length(f)]]))) { f <- f[-length(f)]; j <- j + 1L }
    r <- .nf_sqf_list(f, K)
    fs <- r$factors
    if (j > 0L) {
      xf <- list(.nfe_const(1), .bq0())
      hit <- which(vapply(fs, function(fk) fk[[2L]] == j, TRUE))
      if (length(hit)) fs[[hit[1L]]][[1L]] <- .nfp_mul(fs[[hit[1L]]][[1L]], xf, K)
      else fs <- c(list(list(xf, j)), fs)
    }
    fs <- .sort_by_method(fs, method, .nf_sort)
    return(list(coeff = .nfe_to_expr(r$lc, K),
                factors = lapply(fs, function(fk) list(.nfp_to_expr(fk[[1L]], x, K), .int(fk[[2L]])))))
  }
  u <- length(gs) - 1L
  r <- .kd_factor_list(.kd_from_expr(e, gs, K), u, K, if (is.null(Fd$gauss)) "alg" else Fd$gauss)
  list(coeff = .nfe_to_expr(r$coeff, K),
       factors = lapply(r$factors, function(fk) list(.kd_to_expr(fk[[1L]], gs, u, K), .int(fk[[2L]]))))
}

# factor() with domain options: numerator and denominator factored over the field.
.field_factor <- function(f, Fd, gens, method = "factor") {
  # _symbolic_factor_list: every factor of together(f) on its own
  if (method == "factor" && Fd$kind == "NF") return(.symbolic_factor(as_sym(f), gens, Fd))
  nd <- .as_numer_denom(together(as_sym(f)))
  n <- .field_factor_list(nd[[1L]], Fd, gens, method)
  d <- .field_factor_list(nd[[2L]], Fd, gens, method)
  body <- .mul(c(lapply(n$factors, function(p) .pow(p[[1L]], p[[2L]])),
                 lapply(d$factors, function(p) .pow(p[[1L]], .neg(p[[2L]])))))
  .keep_coeff(.div(n$coeff, d$coeff), body)
}

# gcd, lcm and division over a finite or algebraic field (univariate).
.field_binary <- function(f, g, Fd, gens, op) {
  f <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  g <- as_sym(if (inherits(g, "Poly")) as_expr(g) else g)
  if (Fd$kind == "GF") return(.gf_binary(f, g, Fd, gens, op))
  K <- .nf_field(Fd$exts)
  gs <- Filter(function(z) !(z$type == "Pow" && is_exact_number(z$args[[1L]])) && !is_const(z, "I"), .gens_of(list(f, g), gens))
  if (length(gs) != 1L) stop("only univariate polynomials over algebraic fields are supported here", call. = FALSE)
  x <- gs[[1L]]
  a <- .nfp_from_expr(f, x, K); b <- .nfp_from_expr(g, x, K)
  ex <- function(v) .nfp_to_expr(v, x, K)
  switch(op,
    gcd = ex(.nfp_gcd(a, b, K)),
    lcm = if (!length(.nfp_trim(a)) || !length(.nfp_trim(b))) S_Zero
          else ex(.nfp_monic(.nfp_divmod(.nfp_mul(a, b, K), .nfp_gcd(a, b, K), K)$q, K)),
    div = { d <- .nfp_divmod(a, b, K); list(q = ex(d$q), r = ex(d$r)) })
}

.gf_binary <- function(f, g, Fd, gens, op) {
  x <- .gf_univariate_gen(list(f, g), gens)
  p <- Fd$p
  a <- .gf_from_expr(f, x, p); b <- .gf_from_expr(g, x, p)
  ex <- function(v) .gf_to_expr(v, x, p, Fd$symmetric)
  switch(op,
    gcd = ex(.gf_gcd(a, b, p)),
    lcm = if (!length(a) || !length(b)) S_Zero else ex(.gf_monic(.gf_divmod(.gf_mul(a, b, p), .gf_gcd(a, b, p), p)$q, p)),
    div = { d <- .gf_divmod(a, b, p); list(q = ex(d$q), r = ex(d$r)) })
}

# The factor of fs (polynomials in x) that vanishes at v, after
# _choose_factor: numerical evaluation at increasing precision.
.choose_factor <- function(fs, x, v) {
  if (length(fs) == 1L) return(fs[[1L]])
  for (dps in c(15, 30, 60, 120)) {
    vals <- vapply(fs, function(f) {
      z <- tryCatch(N(subs(f, x, v), dps), error = function(err) NULL)
      if (is.null(z)) Inf else Mod(as.complex(z))
    }, 0)
    tol <- 10^(-dps / 2)
    hit <- which(vals < tol)
    if (length(hit) == 1L) return(fs[[hit]])
  }
  stop("cannot decide which factor has ", .str(v), " as a root", call. = FALSE)
}

#' Primitive elements and number fields
#'
#' `primitive_element(extension, x)` finds a primitive element of the field
#' generated by the algebraic numbers in `extension`, as SymPy's
#' `primitive_element()`: it returns its minimal polynomial and the integer
#' coefficients `c` with `theta = sum(c * extension)`; with `ex = TRUE` also
#' the coordinates of each extension in the power basis of `theta` (highest
#' power first). `to_number_field(extension, theta)` expresses the primitive
#' element of `extension` in the field of `theta` as an [AlgebraicNumber()],
#' and `field_isomorphism(a, b)` returns the coordinates of `a` in the
#' power basis of `b` (highest power first), or `NULL` when `a` is not in
#' `Q(b)`.
#'
#' @param extension An algebraic number or a list of them.
#' @param x The variable of the minimal polynomial (a dummy by default).
#' @param ex Also return the coordinates of the extensions.
#' @param theta A generator of the target field.
#' @param gen The variable of the minimal polynomial of the result.
#' @param alias A symbol (or name) printed for the generator.
#' @param a,b Algebraic numbers.
#' @return A list, an `AlgebraicNumber`, or `NULL`.
#' @examples
#' x <- Symbol("x")
#' primitive_element(list(sqrt(Integer(2)), sqrt(Integer(3))), x)
#' field_isomorphism(sqrt(Integer(2)), sqrt(Integer(2)) + sqrt(Integer(3)))
#' to_number_field(sqrt(Integer(2)), sqrt(Integer(2)) + sqrt(Integer(3)))
#' @name primitive_element
NULL

#' @rdname primitive_element
#' @export
primitive_element <- function(extension, x = NULL, ex = FALSE) {
  exts <- if (is_sym(extension)) list(extension) else lapply(as.list(extension), as_sym)
  if (!length(exts)) stop("Cannot compute primitive element for empty extension", call. = FALSE)
  x <- if (is.null(x)) Dummy("x") else as_sym(x)
  gen <- exts[[1L]]; coeffs <- list(S_One)
  g <- minimal_polynomial(gen, x)
  for (e in exts[-1L]) {
    if (is_exact_number(e)) { coeffs[[length(coeffs) + 1L]] <- S_Zero; next }
    K <- .nf_field(list(e))
    fl <- .field_factor_list(g, list(kind = "NF", exts = list(e)), list(x), "factor")
    gk <- .choose_factor(lapply(fl$factors, `[[`, 1L), x, gen)
    sn <- .nfp_sqf_norm(.nfp_from_expr(gk, x, K), K)
    g <- up_to_expr(up_primitive(sn$r)$poly, x)
    gen <- .add(list(gen, .mul(list(.int(sn$s), e))))
    coeffs[[length(coeffs) + 1L]] <- .int(sn$s)
  }
  # _canonicalize: primitive with a positive leading coefficient
  g <- expand(g)
  pr <- primitive(g, x)$primitive
  if (could_extract_minus_sign(LC(pr, x))) pr <- expand(.neg(pr))
  if (!ex) return(list(pr, coeffs))
  theta <- expand(.add(Map(function(c, a) .mul(list(c, a)), coeffs, exts)))
  d <- as.integer(degree(pr, x))
  y <- Dummy("y")
  H <- lapply(exts, function(a) {
    p <- .to_field_poly(a, theta, y, d)
    lapply(.coeff_list(expand(p), y), function(c) c)
  })
  list(pr, coeffs, H)
}
