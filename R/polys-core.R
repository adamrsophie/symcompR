# Sparse multivariate polynomials over QQ, after sympy/polys/rings.py.
#
# A qpoly is list(e = integer exponent matrix (terms x variables),
# c = bigq coefficient vector), terms sorted by a monomial order with the
# leading term first, no zero coefficients and no repeated monomials.

.ORDER <- new.env()
.ORDER$current <- "lex"

.bq0 <- function() gmp::as.bigq(integer(0))

qp <- function(e, c) structure(list(e = e, c = c), class = "qpoly")

qp_zero <- function(n) qp(matrix(0L, 0L, n), .bq0())

qp_const <- function(v, n) {
  v <- gmp::as.bigq(v)
  if (v == 0) return(qp_zero(n))
  qp(matrix(0L, 1L, n), v)
}

qp_var <- function(i, n, k = 1L) {
  e <- matrix(0L, 1L, n)
  e[1L, i] <- as.integer(k)
  qp(e, gmp::as.bigq(1))
}

qp_nterms <- function(p) length(p$c)
qp_nvars <- function(p) ncol(p$e)
qp_is_zero <- function(p) length(p$c) == 0L

.mono_order_index <- function(e, order = .ORDER$current) {
  n <- ncol(e)
  if (nrow(e) < 2L) return(seq_len(nrow(e)))
  if (n == 0L) return(seq_len(nrow(e)))
  cols <- lapply(seq_len(n), function(j) e[, j])
  if (order == "lex") {
    return(do.call(base::order, c(lapply(cols, `-`), list(method = "radix"))))
  }
  deg <- rowSums(e)
  if (order == "grlex") {
    return(do.call(base::order, c(list(-deg), lapply(cols, `-`), list(method = "radix"))))
  }
  # grevlex
  do.call(base::order, c(list(-deg), rev(cols), list(method = "radix")))
}

.mono_keys <- function(e) {
  if (ncol(e) == 0L) return(rep("", nrow(e)))
  if (ncol(e) == 1L) return(as.character(e[, 1L]))
  do.call(paste, c(lapply(seq_len(ncol(e)), function(j) e[, j]), sep = ","))
}

# Combine duplicates, drop zeros, sort.
qp_norm <- function(e, c, order = .ORDER$current) {
  if (!length(c)) return(qp(matrix(0L, 0L, ncol(e)), .bq0()))
  storage.mode(e) <- "integer"
  keys <- .mono_keys(e)
  if (anyDuplicated(keys)) {
    o <- base::order(keys, method = "radix")
    keys <- keys[o]
    e <- e[o, , drop = FALSE]
    c <- c[o]
    r <- rle(keys)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    cs <- cumsum(c)
    sums <- cs[ends]
    multi <- which(starts > 1L)
    if (length(multi)) sums[multi] <- sums[multi] - cs[starts[multi] - 1L]
    e <- e[starts, , drop = FALSE]
    c <- sums
  }
  nz <- !(c == 0)
  if (!all(nz)) {
    e <- e[nz, , drop = FALSE]
    c <- c[nz]
  }
  if (length(c) > 1L) {
    o <- .mono_order_index(e, order)
    e <- e[o, , drop = FALSE]
    c <- c[o]
  }
  qp(e, c)
}

qp_add <- function(a, b) {
  if (qp_is_zero(a)) return(b)
  if (qp_is_zero(b)) return(a)
  qp_norm(rbind(a$e, b$e), c(a$c, b$c))
}

qp_neg <- function(a) qp(a$e, -a$c)

qp_sub <- function(a, b) {
  if (qp_is_zero(b)) return(a)
  qp_add(a, qp_neg(b))
}

qp_scale <- function(a, k) {
  k <- gmp::as.bigq(k)
  if (k == 0) return(qp_zero(qp_nvars(a)))
  qp(a$e, a$c * k)
}

qp_mul_term <- function(a, mono, coef) {
  if (qp_is_zero(a)) return(a)
  e <- a$e + matrix(mono, nrow(a$e), length(mono), byrow = TRUE)
  qp(e, a$c * coef)
}

qp_mul <- function(a, b) {
  na <- length(a$c)
  nb <- length(b$c)
  if (!na || !nb) return(qp_zero(qp_nvars(a)))
  if (na == 1L) return(qp_norm(b$e + matrix(a$e[1L, ], nb, ncol(b$e), byrow = TRUE), b$c * a$c))
  if (nb == 1L) return(qp_norm(a$e + matrix(b$e[1L, ], na, ncol(a$e), byrow = TRUE), a$c * b$c))
  ia <- rep(seq_len(na), each = nb)
  ib <- rep(seq_len(nb), times = na)
  e <- a$e[ia, , drop = FALSE] + b$e[ib, , drop = FALSE]
  qp_norm(e, a$c[ia] * b$c[ib])
}

qp_pow <- function(a, n) {
  n <- as.integer(n)
  if (n < 0L) stop("negative power of polynomial")
  res <- qp_const(1, qp_nvars(a))
  base <- a
  while (n > 0L) {
    if (n %% 2L == 1L) res <- qp_mul(res, base)
    n <- n %/% 2L
    if (n > 0L) base <- qp_mul(base, base)
  }
  res
}

qp_equal <- function(a, b) {
  if (length(a$c) != length(b$c)) return(FALSE)
  if (!length(a$c)) return(TRUE)
  all(a$e == b$e) && all(a$c == b$c)
}

qp_lc <- function(p) if (length(p$c)) p$c[1L] else gmp::as.bigq(0)
qp_lm <- function(p) p$e[1L, ]
qp_degree <- function(p, i = 1L) if (length(p$c)) max(p$e[, i]) else -Inf
qp_total_degree <- function(p) if (length(p$c)) max(rowSums(p$e)) else -Inf

qp_is_const <- function(p) !length(p$c) || (length(p$c) == 1L && all(p$e[1L, ] == 0L))

qp_monic <- function(p) if (qp_is_zero(p)) p else qp(p$e, p$c / p$c[1L])

# Re-sort under another monomial order.
qp_reorder <- function(p, order) {
  if (length(p$c) < 2L) return(p)
  o <- .mono_order_index(p$e, order)
  qp(p$e[o, , drop = FALSE], p$c[o])
}

# Polynomial with variables permuted / embedded: new column j takes old
# column map[j] (0 means a new variable).
qp_remap <- function(p, map, n) {
  e <- matrix(0L, nrow(p$e), n)
  for (j in seq_len(n)) if (map[j] > 0L) e[, j] <- p$e[, map[j]]
  qp_norm(e, p$c)
}

# Evaluate variable i at rational value v.
qp_eval_var <- function(p, i, v) {
  if (qp_is_zero(p)) return(p)
  v <- gmp::as.bigq(v)
  pw <- p$e[, i]
  e <- p$e
  e[, i] <- 0L
  if (v == 0) {
    keep <- pw == 0L
    return(qp_norm(e[keep, , drop = FALSE], p$c[keep]))
  }
  maxp <- max(pw)
  powers <- gmp::as.bigq(rep(1, maxp + 1L))
  if (maxp > 0L) for (k in seq_len(maxp)) powers[k + 1L] <- powers[k] * v
  qp_norm(e, p$c * powers[pw + 1L])
}

# Substitute polynomial q for variable i.
qp_subs_var <- function(p, i, q) {
  if (qp_is_zero(p)) return(p)
  pw <- p$e[, i]
  maxp <- max(pw)
  if (maxp == 0L) return(p)
  n <- qp_nvars(p)
  pows <- vector("list", maxp + 1L)
  pows[[1L]] <- qp_const(1, n)
  for (k in seq_len(maxp)) pows[[k + 1L]] <- qp_mul(pows[[k]], q)
  base <- p$e
  base[, i] <- 0L
  res <- qp_zero(n)
  for (k in 0:maxp) {
    w <- which(pw == k)
    if (!length(w)) next
    part <- qp_norm(base[w, , drop = FALSE], p$c[w])
    res <- qp_add(res, qp_mul(part, pows[[k + 1L]]))
  }
  res
}

qp_diff <- function(p, i) {
  if (qp_is_zero(p)) return(p)
  pw <- p$e[, i]
  keep <- pw > 0L
  e <- p$e[keep, , drop = FALSE]
  e[, i] <- e[, i] - 1L
  qp_norm(e, p$c[keep] * gmp::as.bigq(pw[keep]))
}

qp_vars_used <- function(p) {
  if (qp_is_zero(p)) return(integer())
  which(colSums(p$e) > 0L)
}

# Content: gcd of numerators over lcm of denominators (positive rational).
qp_content <- function(p) {
  if (qp_is_zero(p)) return(gmp::as.bigq(0))
  nums <- abs(gmp::numerator(p$c))
  dens <- gmp::denominator(p$c)
  g <- nums[1L]
  if (length(nums) > 1L) for (k in 2:length(nums)) {
    g <- gmp::gcd.bigz(g, nums[k])
    if (g == 1) break
  }
  l <- dens[1L]
  if (length(dens) > 1L) for (k in 2:length(dens)) l <- gmp::lcm.bigz(l, dens[k])
  gmp::as.bigq(g, l)
}

# Primitive integer polynomial with positive leading coefficient.
qp_primitive <- function(p) {
  if (qp_is_zero(p)) return(list(content = gmp::as.bigq(0), poly = p))
  ct <- qp_content(p)
  if (qp_lc(p) < 0) ct <- -ct
  list(content = ct, poly = qp(p$e, p$c / ct))
}

# Does monomial a divide monomial b (vectors)?
.mono_divides <- function(a, b) all(a <= b)

# Multivariate division of f by the list gs (sympy rings.div).
qp_div_multi <- function(f, gs, order = .ORDER$current) {
  n <- qp_nvars(f)
  qs <- lapply(gs, function(g) qp_zero(n))
  r_e <- list()
  r_c <- list()
  p <- f
  lms <- lapply(gs, qp_lm)
  lcs <- lapply(gs, qp_lc)
  while (!qp_is_zero(p)) {
    lmp <- p$e[1L, ]
    lcp <- p$c[1L]
    divided <- FALSE
    for (k in seq_along(gs)) {
      if (.mono_divides(lms[[k]], lmp)) {
        mono <- lmp - lms[[k]]
        cf <- lcp / lcs[[k]]
        qs[[k]] <- qp_add(qs[[k]], qp(matrix(mono, 1L), cf))
        p <- qp_sub(p, qp_mul_term(gs[[k]], mono, cf))
        divided <- TRUE
        break
      }
    }
    if (!divided) {
      r_e[[length(r_e) + 1L]] <- lmp
      r_c[[length(r_c) + 1L]] <- lcp
      p <- qp(p$e[-1L, , drop = FALSE], p$c[-1L])
    }
  }
  r <- if (length(r_e)) qp_norm(do.call(rbind, r_e), do.call(c, r_c), order) else qp_zero(n)
  list(q = qs, r = r)
}

# Remainder of f modulo the list gs (normal form); faster than full division.
qp_rem_multi <- function(f, gs) {
  if (!length(gs)) return(f)
  n <- qp_nvars(f)
  p <- f
  lms <- do.call(rbind, lapply(gs, qp_lm))
  lcs <- lapply(gs, qp_lc)
  rem_e <- list()
  rem_c <- list()
  while (!qp_is_zero(p)) {
    lmp <- p$e[1L, ]
    ok <- which(colSums(t(lms) <= lmp) == n)
    if (length(ok)) {
      k <- ok[1L]
      cf <- p$c[1L] / lcs[[k]]
      p <- qp_sub(p, qp_mul_term(gs[[k]], lmp - lms[k, ], cf))
    } else {
      rem_e[[length(rem_e) + 1L]] <- lmp
      rem_c[[length(rem_c) + 1L]] <- p$c[1L]
      p <- qp(p$e[-1L, , drop = FALSE], p$c[-1L])
    }
  }
  if (!length(rem_e)) return(qp_zero(n))
  qp_norm(do.call(rbind, rem_e), do.call(c, rem_c))
}

# Exact division f / g, or NULL when g does not divide f.
qp_exquo <- function(f, g) {
  if (qp_is_zero(g)) stop("polynomial division by zero")
  d <- qp_div_multi(f, list(g))
  if (!qp_is_zero(d$r)) return(NULL)
  d$q[[1L]]
}

# ---------------------------------------------------------------------------
# Conversion between expressions and qpolys.

# Split an expanded term into (rational coefficient, list of factors).
.term_factors <- function(t) {
  if (is_exact_number(t)) return(list(coef = t$value, facs = list()))
  if (t$type == "Mul") {
    if (is_exact_number(t$args[[1L]])) return(list(coef = t$args[[1L]]$value, facs = t$args[-1L]))
    return(list(coef = gmp::as.bigq(1), facs = t$args))
  }
  list(coef = gmp::as.bigq(1), facs = list(t))
}

# Base and nonnegative integer exponent of a factor, or NULL.
# decompose_power of _parallel_dict_from_expr_no_gens: x^(3/2) is
# sqrt(x)^3, exp(2*x) is exp(x)^2, and a negative power is one of 1/base
.gen_power <- function(f) {
  if (f$type == "Pow" && f$args[[2L]]$type == "Integer" && f$args[[2L]]$value > 0)
    return(list(f$args[[1L]], as.integer(as.character(f$args[[2L]]$value))))
  if (!(f$type == "Pow" || (f$type == "Function" && f$name == "exp"))) return(list(f, 1L))
  dp <- .decompose_power(f)
  b <- dp[[1L]]; e <- dp[[2L]]
  if (e < 0) { e <- -e; b <- .pow(b, S_NegativeOne) }
  list(b, as.integer(as.character(e)))
}

# Collect polynomial generators of expressions (sympy _parallel_dict_from_expr).
.find_gens <- function(exprs) {
  gens <- list()
  add_gen <- function(g) gens[[g$key]] <<- g
  for (ex in exprs) {
    terms <- if (ex$type == "Add") ex$args else list(ex)
    for (t in terms) {
      tf <- .term_factors(t)
      for (f in tf$facs) {
        gp <- .gen_power(f)
        if (is_number_type(gp[[1L]])) {
          if (gp[[1L]]$type == "Float") next
          if (is_exact_number(gp[[1L]]) && f$type != "Pow") next
        }
        add_gen(gp[[1L]])
      }
    }
  }
  gens <- unname(gens)
  .poly_sort_gens(gens)
}

# SymPy's _sort_gens: x, y, z first, then p to w, then a to o, then other
# names; a trailing number orders numerically.
.GENS_ORDER <- c(a = 301, b = 302, c = 303, d = 304, e = 305, f = 306, g = 307, h = 308, i = 309, j = 310,
                 k = 311, l = 312, m = 313, n = 314, o = 315, p = 216, q = 217, r = 218, s = 219, t = 220,
                 u = 221, v = 222, w = 223, x = 124, y = 125, z = 126)
.poly_sort_gens <- function(gens) {
  if (length(gens) < 2L) return(gens)
  nm <- vapply(gens, function(g) gsub("^", "**", .str(g), fixed = TRUE), "")
  m <- regmatches(nm, regexec("^(.*?)([0-9]*)$", nm))
  base <- vapply(m, `[`, "", 2L)
  idx <- vapply(m, function(z) if (nzchar(z[3L])) as.numeric(z[3L]) else 0, 0)
  ord <- ifelse(base %in% names(.GENS_ORDER), .GENS_ORDER[base], 1000)
  gens[order(ord, base, idx, method = "radix")]
}

.sort_gens <- function(gens) {
  if (length(gens) < 2L) return(gens)
  rank <- vapply(gens, function(g) if (is_symbol_like(g)) 0L
                 else if (length(free_symbols(g)) == 0L) 2L else 1L, 0L)
  nm <- vapply(gens, function(g) if (is_symbol_like(g)) g$name else .str(g), "")
  nat <- .natural_key(nm)
  gens[order(rank, nat, nm, method = "radix")]
}

# Sort key placing x2 before x10.
.natural_key <- function(nm) {
  m <- regmatches(nm, regexec("^(.*?)([0-9]+)$", nm))
  vapply(seq_along(nm), function(i) {
    mm <- m[[i]]
    if (length(mm)) sprintf("%s%012d", mm[2], as.numeric(mm[3])) else nm[i]
  }, "")
}

#' @noRd
expr_to_qpoly <- function(ex, gens, strict = TRUE) {
  ex <- expand(as_sym(ex))
  n <- length(gens)
  gk <- vapply(gens, .k, "")
  terms <- if (ex$type == "Add") ex$args else list(ex)
  if (length(terms) == 1L && is_zero_num(terms[[1L]])) return(qp_zero(n))
  e <- matrix(0L, length(terms), n)
  cs <- vector("list", length(terms))
  for (ti in seq_along(terms)) {
    t <- terms[[ti]]
    if (t$type == "Float") {
      cs[[ti]] <- gmp::as.bigq(t$value)
      next
    }
    tf <- .term_factors(t)
    coef <- tf$coef
    for (f in tf$facs) {
      if (f$type == "Float") {
        coef <- coef * gmp::as.bigq(f$value)
        next
      }
      gp <- .gen_power(f)
      j <- match(gp[[1L]]$key, gk)
      if (is.na(j)) {
        j2 <- match(f$key, gk)
        if (!is.na(j2)) {
          e[ti, j2] <- e[ti, j2] + 1L
          next
        }
        if (strict) return(NULL)
        stop("expression is not a polynomial in the given generators: ", .str(ex))
      }
      e[ti, j] <- e[ti, j] + gp[[2L]]
    }
    cs[[ti]] <- coef
  }
  qp_norm(e, do.call(c, cs))
}

#' @noRd
qpoly_to_expr <- function(p, gens) {
  if (qp_is_zero(p)) return(S_Zero)
  terms <- vector("list", length(p$c))
  for (i in seq_along(p$c)) {
    parts <- list(make_rational(p$c[i]))
    for (j in which(p$e[i, ] > 0L)) parts[[length(parts) + 1L]] <- .pow(gens[[j]], .int(p$e[i, j]))
    terms[[i]] <- .mul(parts)
  }
  .add(terms)
}

# Convert several expressions to qpolys over common generators.
.to_qpolys <- function(exprs, gens = NULL) {
  exprs <- lapply(exprs, function(z) expand(as_sym(z)))
  if (is.null(gens)) gens <- .find_gens(exprs)
  polys <- lapply(exprs, expr_to_qpoly, gens = gens, strict = TRUE)
  if (any(vapply(polys, is.null, TRUE))) {
    extra <- .find_gens(exprs)
    gk <- vapply(gens, .k, "")
    for (g in extra) if (!(g$key %in% gk)) gens[[length(gens) + 1L]] <- g
    polys <- lapply(exprs, expr_to_qpoly, gens = gens, strict = FALSE)
  }
  list(polys = polys, gens = gens)
}

# ---------------------------------------------------------------------------
# Dense univariate polynomials: bigq vectors, highest degree first.

up_trim <- function(a) {
  if (!length(a)) return(a)
  nz <- which(!(a == 0))
  if (!length(nz)) return(.bq0())
  a[nz[1L]:length(a)]
}

up_deg <- function(a) length(a) - 1L

up_add <- function(a, b) {
  la <- length(a)
  lb <- length(b)
  if (la < lb) a <- c(gmp::as.bigq(rep(0, lb - la)), a)
  if (lb < la) b <- c(gmp::as.bigq(rep(0, la - lb)), b)
  up_trim(a + b)
}

up_sub <- function(a, b) up_add(a, -b)

up_mul <- function(a, b) {
  la <- length(a)
  lb <- length(b)
  if (!la || !lb) return(.bq0())
  out <- gmp::as.bigq(rep(0, la + lb - 1L))
  if (la < lb) {
    tmp <- a
    a <- b
    b <- tmp
    lb <- length(b)
  }
  for (j in seq_len(lb)) {
    if (b[j] == 0) next
    idx <- j:(j + length(a) - 1L)
    out[idx] <- out[idx] + a * b[j]
  }
  up_trim(out)
}

up_divmod <- function(a, b) {
  b <- up_trim(b)
  if (!length(b)) stop("polynomial division by zero")
  a <- up_trim(a)
  db <- up_deg(b)
  if (up_deg(a) < db) return(list(q = .bq0(), r = a))
  q <- gmp::as.bigq(rep(0, up_deg(a) - db + 1L))
  r <- a
  lb <- b[1L]
  for (i in seq_along(q)) {
    cf <- r[i] / lb
    q[i] <- cf
    if (cf != 0) {
      idx <- i:(i + db)
      r[idx] <- r[idx] - b * cf
    }
  }
  rr <- if (db > 0) r[(length(r) - db + 1L):length(r)] else .bq0()
  list(q = up_trim(q), r = up_trim(rr))
}

up_monic <- function(a) if (!length(a)) a else a / a[1L]

up_gcd <- function(a, b) {
  a <- up_trim(a)
  b <- up_trim(b)
  while (length(b)) {
    r <- up_divmod(a, b)$r
    a <- b
    b <- r
  }
  up_monic(a)
}

up_diff <- function(a) {
  n <- up_deg(a)
  if (n <= 0) return(.bq0())
  up_trim(a[-length(a)] * gmp::as.bigq(n:1))
}

up_eval <- function(a, v) {
  v <- gmp::as.bigq(v)
  r <- gmp::as.bigq(0)
  for (k in seq_along(a)) r <- r * v + a[k]
  r
}

up_from_qpoly <- function(p, i = 1L) {
  if (qp_is_zero(p)) return(.bq0())
  d <- max(p$e[, i])
  out <- gmp::as.bigq(rep(0, d + 1L))
  for (k in seq_along(p$c)) {
    pos <- d - p$e[k, i] + 1L
    out[pos] <- out[pos] + p$c[k]
  }
  up_trim(out)
}

up_to_qpoly <- function(a, i = 1L, n = 1L) {
  a <- up_trim(a)
  if (!length(a)) return(qp_zero(n))
  d <- up_deg(a)
  e <- matrix(0L, d + 1L, n)
  e[, i] <- d:0
  qp_norm(e, a)
}

up_to_expr <- function(a, x) {
  a <- up_trim(a)
  if (!length(a)) return(S_Zero)
  d <- up_deg(a)
  terms <- lapply(seq_along(a), function(k) .mul(list(make_rational(a[k]), .pow(x, .int(d - k + 1L)))))
  .add(terms)
}

# Integer primitive part with positive leading coefficient.
up_primitive <- function(a) {
  a <- up_trim(a)
  if (!length(a)) return(list(content = gmp::as.bigq(0), poly = a))
  p <- qp(matrix(seq_along(a), ncol = 1L), a)
  ct <- qp_content(p)
  if (a[1L] < 0) ct <- -ct
  list(content = ct, poly = a / ct)
}

# Square-free decomposition (Yun) of a primitive polynomial over QQ.
# Returns list of list(factor, multiplicity).
up_sqf_list <- function(f) {
  f <- up_trim(f)
  if (up_deg(f) <= 0) return(list())
  out <- list()
  fp <- up_diff(f)
  g <- up_gcd(f, fp)
  if (up_deg(g) == 0) return(list(list(up_primitive(f)$poly, 1L)))
  w <- up_divmod(f, g)$q
  y <- up_divmod(fp, g)$q
  i <- 1L
  repeat {
    z <- up_sub(y, up_diff(w))
    if (!length(z)) {
      out[[length(out) + 1L]] <- list(up_primitive(w)$poly, i)
      break
    }
    h <- up_gcd(w, z)
    if (up_deg(h) > 0) out[[length(out) + 1L]] <- list(up_primitive(h)$poly, i)
    w <- up_divmod(w, h)$q
    y <- up_divmod(z, h)$q
    i <- i + 1L
    if (up_deg(w) <= 0) break
  }
  out
}
