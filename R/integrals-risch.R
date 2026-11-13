# The Risch algorithm for transcendental functions, after
# sympy/integrals/risch.py: DifferentialExtension, derivation, splitfactor,
# Hermite reduction, residue reduction and the primitive, hyperexponential
# and nonlinear cases, on the Poly layer of R/polys-rpoly.R.

.ri_error <- function(cls, msg) stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))
.nonelem <- function(msg = "") .ri_error("NonElementaryIntegralException", msg)
.ri_notimpl <- function(msg) .ri_error("NotImplementedError", msg)

.is_rational_function_in <- function(e, syms) {
  e <- as_sym(e)
  if (!length(syms)) return(TRUE)
  keys <- vapply(syms, `[[`, "", "key")
  rec <- function(z) {
    if (z$type %in% c("Symbol", "Dummy")) return(TRUE)
    if (is_number_type(z) || z$type == "Constant") return(TRUE)
    if (z$type %in% c("Add", "Mul")) return(all(vapply(z$args, rec, NA)))
    if (z$type == "Pow") {
      if (!any(vapply(keys, function(k) .contains_key(z$args[[1L]], k), NA)) &&
          !any(vapply(keys, function(k) .contains_key(z$args[[2L]], k), NA))) return(TRUE)
      return(z$args[[2L]]$type == "Integer" && rec(z$args[[1L]]))
    }
    !any(vapply(keys, function(k) .contains_key(z, k), NA))
  }
  rec(e)
}
.has_any_of <- function(e, syms) any(vapply(syms, function(s) .has1(e, s), NA))

# integer_powers
integer_powers <- function(exprs) {
  terms <- list()
  for (term in exprs) {
    found <- FALSE
    for (k in seq_along(terms)) {
      a <- cancel(.div(term, terms[[k]]$base))
      if (a$type %in% c("Integer", "Rational")) { terms[[k]]$list[[length(terms[[k]]$list) + 1L]] <- list(term, a); found <- TRUE; break }
    }
    if (!found) terms[[length(terms) + 1L]] <- list(base = term, list = list(list(term, S_One)))
  }
  out <- lapply(terms, function(tl) {
    cd <- Reduce(function(a, b) .dom_lcm(ZZ, a, b), lapply(tl$list, function(p) .as_numer_denom(p[[2L]])[[2L]]), S_One)
    newterm <- .div(tl$base, cd)
    list(newterm, lapply(tl$list, function(p) list(p[[1L]], .mul(list(p[[2L]], cd)))))
  })
  .key_sort(lapply(out, function(o) .sort_key(o[[1L]])), out)
}

# update_sets: keep the old members still among the atoms, add the new atoms
# satisfying func
.update_sets <- function(seq, atoms, func) {
  ak <- vapply(atoms, `[[`, "", "key"); sk <- vapply(seq, `[[`, "", "key")
  s <- atoms[ak %in% sk]
  new <- atoms[!(ak %in% sk)]
  c(s, Filter(func, new))
}
.uniq_syms <- function(l) l[!duplicated(vapply(l, `[[`, "", "key"))]
.atoms_fn <- function(e, name) .uniq_syms(Filter(function(z) z$type == "Function" && z$name == name, .all_subexprs(e)))
.atoms_pow <- function(e) .uniq_syms(Filter(function(z) z$type == "Pow", .all_subexprs(e)))
.all_subexprs <- function(e) { out <- list(); rec <- function(z) { out[[length(out) + 1L]] <<- z; for (a in z$args) rec(a) }; rec(e); out }

# ---------------------------------------------------------------------------
# DifferentialExtension

DifferentialExtension <- function(f = NULL, x = NULL, handle_first = "log", dummy = FALSE, extension = NULL, rewrite_complex = NULL) {
  DE <- new.env()
  class(DE) <- c("DifferentialExtension", "environment")
  if (!is.null(extension)) {
    if (is.null(extension$D)) stop("At least the key D must be included with the extension flag to DifferentialExtension.")
    for (nm in names(extension)) assign(nm, extension[[nm]], envir = DE)
    .de_auto_attrs(DE)
    return(DE)
  }
  if (is.null(f) || is.null(x)) stop("Either both f and x or a manual extension must be given.")
  DE$f <- as_sym(f); DE$x <- as_sym(x); DE$dummy <- dummy
  .de_reset(DE)
  exp_new <- TRUE; log_new <- TRUE
  if (is.null(rewrite_complex)) rewrite_complex <- any(vapply(.ri_atoms(DE$f), function(a) identical(a$key, S_I$key), NA))
  trigs <- c("sin", "cos", "cot", "tan", "sinh", "cosh", "coth", "tanh")
  atrigs <- c("asin", "acos", "acot", "atan")
  if (rewrite_complex) {
    DE$newf <- .map_bottom_up(DE$newf, function(z) {
      if (z$type == "Function" && z$name %in% trigs) return(rewrite(z, "exp"))
      if (z$type == "Function" && z$name %in% atrigs) return(rewrite(z, "log"))
      z
    })
    DE$newf <- cancel(DE$newf)
  } else if (any(vapply(Filter(function(z) z$type == "Function" && z$name %in% c(trigs, atrigs), .all_subexprs(DE$f)),
                        function(z) .has1(z, DE$x), NA))) {
    .ri_notimpl("Trigonometric and hyperbolic extensions are not supported (yet!).  Try rewriting in terms of exp and log, or using rewrite_complex=True.")
  }
  st <- new.env()
  st$exps <- list(); st$pows <- list(); st$numpows <- list(); st$sympows <- list(); st$logs <- list(); st$symlogs <- list()
  repeat {
    if (.is_rational_function_in(DE$newf, DE$T)) break
    if (!exp_new && !log_new) .ri_notimpl(paste0("Couldn't find an elementary transcendental extension for ", .str(f), ".  Try using a manual extension with the extension flag."))
    log_new <- .de_rewrite_exps_pows(DE, st, log_new)
    .de_rewrite_logs(DE, st)
    if (handle_first == "exp" || !log_new) {
      exp_new <- .de_exp_part(DE, st$exps)
      if (is.null(exp_new)) {
        DE$f <- DE$newf
        .de_reset(DE)
        exp_new <- TRUE
        next
      }
    }
    if (handle_first == "log" || !exp_new) log_new <- .de_log_part(DE, st$logs)
  }
  fr <- frac_in(DE$newf, DE$t)
  DE$fa <- fr[[1L]]; DE$fd <- fr[[2L]]
  .de_auto_attrs(DE)
  DE
}

.de_reset <- function(DE) {
  DE$t <- DE$x
  DE$T <- list(DE$x)
  DE$D <- list(.P(S_One, DE$x))
  DE$level <- -1L
  DE$exts <- list(NULL)
  DE$extargs <- list(NULL)
  DE$tcount <- 0L
  DE$backsubs <- list()
  DE$Tfuncs <- list()
  DE$newf <- DE$f
}
.de_next_t <- function(DE) {
  nm <- paste0("t", DE$tcount); DE$tcount <- DE$tcount + 1L
  if (isTRUE(DE$dummy)) Dummy(nm) else Symbol(nm)
}
.de_auto_attrs <- function(DE) {
  if (is.null(DE$T)) DE$T <- lapply(DE$D, function(d) d$x)
  if (is.null(DE$x)) DE$x <- DE$T[[1L]]
  DE$cases <- Map(get_case, DE$D, DE$T)
  DE$level <- -1L
  .de_set_level(DE)
}
.de_idx <- function(DE) length(DE$T) + DE$level + 1L
.de_set_level <- function(DE) {
  i <- .de_idx(DE)
  DE$t <- DE$T[[i]]; DE$d <- DE$D[[i]]; DE$case <- DE$cases[[i]]
}
.de_increment <- function(DE) {
  if (DE$level >= -1L) stop("The level of the differential extension cannot be incremented any further.")
  DE$level <- DE$level + 1L; .de_set_level(DE)
}
.de_decrement <- function(DE) {
  if (DE$level <= -length(DE$T)) stop("The level of the differential extension cannot be decremented any further.")
  DE$level <- DE$level - 1L; .de_set_level(DE)
}
# with DecrementLevel(DE): ...
.with_dec <- function(DE, expr) {
  .de_decrement(DE)
  on.exit(.de_increment(DE))
  force(expr)
}
.de_indices <- function(DE, ext) which(vapply(DE$exts, function(e) identical(e, ext), NA))
.de_backsubs_apply <- function(e, DE) {
  for (bs in DE$backsubs) e <- subs(e, bs[[1L]], bs[[2L]])
  e
}
# list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))
.de_s <- function(DE) {
  Tr <- rev(DE$T); Fr <- rev(lapply(DE$Tfuncs, function(lam) .lambda_call(lam, DE$x)))
  n <- min(length(Tr), length(Fr))
  lapply(seq_len(n), function(i) list(Tr[[i]], Fr[[i]]))
}
.lambda_call <- function(lam, v) subs(lam$args[[2L]], lam$args[[1L]]$args[[1L]], v)
.subs_seq <- function(e, pairs) { for (p in pairs) e <- subs(e, p[[1L]], p[[2L]]); e }

.de_rewrite_exps_pows <- function(DE, st, log_new) {
  Tl <- DE$T
  ratpows <- Filter(function(i) i$args[[1L]]$type == "Function" && i$args[[1L]]$name == "exp" && i$args[[2L]]$type %in% c("Integer", "Rational"),
                    .atoms_pow(DE$newf))
  repl <- lapply(ratpows, function(i) list(i, fcall("exp", list(.mul(list(i$args[[2L]], i$args[[1L]]$args[[1L]]))))))
  DE$backsubs <- c(DE$backsubs, lapply(repl, function(p) list(p[[2L]], p[[1L]])))
  for (p in repl) DE$newf <- xreplace(DE$newf, list(p))
  rf <- function(i, e) .is_rational_function_in(e, Tl) && .has_any_of(e, Tl)
  st$exps <- .update_sets(st$exps, .atoms_fn(DE$newf, "exp"), function(i) rf(i, i$args[[1L]]))
  st$pows <- .update_sets(st$pows, Filter(function(z) !(z$args[[1L]]$type == "Function" && FALSE), .atoms_pow(DE$newf)), function(i) rf(i, i$args[[2L]]))
  st$numpows <- .update_sets(st$numpows, st$pows, function(i) !.has_any_of(i$args[[1L]], Tl))
  nk <- vapply(st$numpows, `[[`, "", "key")
  st$sympows <- .update_sets(st$sympows, Filter(function(z) !(z$key %in% nk), st$pows),
                             function(i) .is_rational_function_in(i$args[[1L]], Tl) && i$args[[2L]]$type != "Integer")
  sk <- vapply(st$sympows, `[[`, "", "key")
  for (i in .bool_print_order(st$pows)) {
    old <- i
    new <- fcall("exp", list(.mul(list(i$args[[2L]], fcall("log", list(i$args[[1L]]))))))
    if (i$key %in% sk) {
      if (i$args[[2L]]$type %in% c("Integer", "Rational")) .ri_notimpl(paste0("Algebraic extensions are not supported (", .str(i), ")."))
      bfr <- frac_in(i$args[[1L]], DE$t)
      A <- is_deriv_k(bfr[[1L]], bfr[[2L]], DE)
      if (is.null(A)) {
        DE$newf <- xreplace(DE$newf, list(list(old, new)))
        DE$backsubs[[length(DE$backsubs) + 1L]] <- list(new, old)
        log_new <- .de_log_part(DE, list(fcall("log", list(i$args[[1L]]))))
        st$exps <- .update_sets(st$exps, .atoms_fn(DE$newf, "exp"), function(z) rf(z, z$args[[1L]]))
        next
      }
      newterm <- fcall("exp", list(.mul(list(i$args[[2L]], .add(list(fcall("log", list(A[[3L]])), A[[2L]]))))))
      DE$newf <- xreplace(DE$newf, list(list(i, newterm)))
    } else if (!(i$key %in% nk)) {
      next
    } else newterm <- new
    DE$backsubs[[length(DE$backsubs) + 1L]] <- list(new, old)
    DE$newf <- xreplace(DE$newf, list(list(old, newterm)))
    st$exps[[length(st$exps) + 1L]] <- newterm
  }
  log_new
}

.de_rewrite_logs <- function(DE, st) {
  Tl <- DE$T
  atoms <- .atoms_fn(DE$newf, "log")
  st$logs <- .update_sets(st$logs, atoms, function(i) .is_rational_function_in(i$args[[1L]], Tl) && .has_any_of(i$args[[1L]], Tl))
  st$symlogs <- .update_sets(st$symlogs, atoms, function(i) .has_any_of(i, Tl) && i$args[[1L]]$type == "Pow" &&
                               .is_rational_function_in(i$args[[1L]]$args[[1L]], Tl) && i$args[[1L]]$args[[2L]]$type != "Integer")
  for (i in .bool_print_order(st$symlogs)) {
    lbase <- fcall("log", list(i$args[[1L]]$args[[1L]]))
    st$logs[[length(st$logs) + 1L]] <- lbase
    new <- .mul(list(i$args[[1L]]$args[[2L]], lbase))
    DE$newf <- xreplace(DE$newf, list(list(i, new)))
    DE$backsubs[[length(DE$backsubs) + 1L]] <- list(new, i)
  }
  lg <- .uniq_syms(st$logs)
  st$logs <- .key_sort(lapply(lg, .sort_key), lg)
}

.de_exp_part <- function(DE, exps) {
  new_extension <- FALSE; restart <- FALSE
  expargs <- lapply(exps, function(i) i$args[[1L]])
  ip <- integer_powers(expargs)
  for (pr in ip) {
    arg <- pr[[1L]]; others <- pr[[2L]]
    others <- others[order(vapply(others, function(o) as.double(gmp::as.bigq(o[[2L]]$value)), 0))]
    afr <- frac_in(arg, DE$t)
    A <- is_log_deriv_k_t_radical(afr[[1L]], afr[[2L]], DE)
    if (!is.null(A)) {
      ans <- A[[1L]]; u <- A[[2L]]; n <- A[[3L]]; const <- A[[4L]]
      if (is_const_num(n, -1)) {
        n <- S_One; u <- .pow(u, S_NegativeOne); const <- .neg(const)
        ans <- lapply(ans, function(p) list(p[[1L]], .neg(p[[2L]])))
      }
      if (is_one_num(n)) {
        prod <- .mul(lapply(ans, function(p) .pow(p[[1L]], p[[2L]])))
        DE$newf <- xreplace(DE$newf, list(list(fcall("exp", list(arg)), .mul(list(fcall("exp", list(const)), prod)))))
        DE$newf <- xreplace(DE$newf, lapply(others, function(o) list(fcall("exp", list(.mul(list(o[[2L]], o[[1L]])))),
                                                                     .mul(list(fcall("exp", list(.mul(list(const, o[[2L]])))), prod)))))
        next
      } else if (!is_zero_num(const) || length(ans) > 1L) {
        rad <- .mul(lapply(ans, function(p) .pow(p[[1L]], .div(p[[2L]], n))))
        DE$newf <- xreplace(DE$newf, lapply(others, function(o) list(fcall("exp", list(.mul(list(o[[2L]], o[[1L]])))),
                                                                     .mul(list(fcall("exp", list(.mul(list(const, o[[2L]])))), rad)))))
        Tr <- rev(DE$T); Fr <- rev(lapply(DE$Tfuncs, function(l) .lambda_call(l, DE$x)))
        DE$newf <- xreplace(DE$newf, lapply(seq_len(min(length(Tr), length(Fr))), function(k) list(Tr[[k]], Fr[[k]])))
        restart <- TRUE
        break
      } else .ri_notimpl("Cannot integrate over algebraic extensions.")
    } else {
      afr <- frac_in(arg, DE$t)
      arga <- afr[[1L]]; argd <- afr[[2L]]
      darga <- P_sub(P_mul(argd, derivation(.P(P_expr(arga), DE$t), DE)), P_mul(arga, derivation(.P(P_expr(argd), DE$t), DE)))
      dargd <- P_pow(argd, 2L)
      cc <- P_cancel(darga, dargd, include = TRUE)
      darg <- .div(P_expr(cc[[1L]]), P_expr(cc[[2L]]))
      DE$t <- .de_next_t(DE)
      DE$T[[length(DE$T) + 1L]] <- DE$t
      DE$extargs[length(DE$extargs) + 1L] <- list(arg)
      DE$exts[length(DE$exts) + 1L] <- list("exp")
      DE$D[[length(DE$D) + 1L]] <- P_mul(.P(darg, DE$t), .P(DE$t, DE$t))
      iv <- if (isTRUE(DE$dummy)) Dummy("i") else Symbol("i")
      DE$Tfuncs[[length(DE$Tfuncs) + 1L]] <- Lambda(iv, fcall("exp", list(subs(arg, DE$x, iv))))
      DE$newf <- xreplace(DE$newf, lapply(others, function(o) list(fcall("exp", list(o[[1L]])), .pow(DE$t, o[[2L]]))))
      new_extension <- TRUE
    }
  }
  if (restart) return(NULL)
  new_extension
}

.de_log_part <- function(DE, logs) {
  new_extension <- FALSE
  logargs <- lapply(logs, function(i) i$args[[1L]])
  for (arg in .bool_print_order(logargs)) {
    afr <- frac_in(arg, DE$t)
    A <- is_deriv_k(afr[[1L]], afr[[2L]], DE)
    if (!is.null(A)) {
      newterm <- .add(list(fcall("log", list(A[[3L]])), A[[2L]]))
      DE$newf <- xreplace(DE$newf, list(list(fcall("log", list(arg)), newterm)))
      next
    }
    afr <- frac_in(arg, DE$t)
    arga <- afr[[1L]]; argd <- afr[[2L]]
    darga <- P_sub(P_mul(argd, derivation(.P(P_expr(arga), DE$t), DE)), P_mul(arga, derivation(.P(P_expr(argd), DE$t), DE)))
    dargd <- P_pow(argd, 2L)
    darg <- .div(P_expr(darga), P_expr(dargd))
    DE$t <- .de_next_t(DE)
    DE$T[[length(DE$T) + 1L]] <- DE$t
    DE$extargs[length(DE$extargs) + 1L] <- list(arg)
    DE$exts[length(DE$exts) + 1L] <- list("log")
    DE$D[[length(DE$D) + 1L]] <- .P(cancel(.div(darg, arg)), DE$t)
    iv <- if (isTRUE(DE$dummy)) Dummy("i") else Symbol("i")
    DE$Tfuncs[[length(DE$Tfuncs) + 1L]] <- Lambda(iv, fcall("log", list(subs(arg, DE$x, iv))))
    DE$newf <- xreplace(DE$newf, list(list(fcall("log", list(arg)), DE$t)))
    new_extension <- TRUE
  }
  new_extension
}

is_const_num <- function(e, v) is_number_type(e) && e$type %in% c("Integer", "Rational") && as.double(gmp::as.bigq(e$value)) == v

# ---------------------------------------------------------------------------
# basic routines

gcdex_diophantine <- function(a, b, c) {
  sg <- P_half_gcdex(a, b)
  s <- P_mul(sg[[1L]], P_exquo(c, sg[[2L]]))
  if (!P_is_zero(s) && P_deg(s) >= P_deg(b)) s <- P_div(s, b)[[2L]]
  t <- P_exquo(P_sub(c, P_mul(s, a)), b)
  list(s, t)
}

frac_in <- function(f, t, cancel = FALSE, field = FALSE) {
  t <- as_sym(t)
  if (inherits(f, "dpoly")) f <- P_expr(f)
  else if (is.list(f) && !is_sym(f)) f <- .div(P_expr(f[[1L]]), P_expr(f[[2L]]))
  nd <- .as_numer_denom(as_sym(f))
  fa <- P_or_null(nd[[1L]], t, field = field); fd <- P_or_null(nd[[2L]], t, field = field)
  if (is.null(fa) || is.null(fd)) stop("Could not turn ", .str(f), " into a fraction in ", .str(t), ".")
  if (cancel) { cc <- P_cancel(fa, fd, include = TRUE); fa <- cc[[1L]]; fd <- cc[[2L]] }
  list(fa, fd)
}

# as_poly_1t: p in K[t, 1/t] as an expression in t and z = 1/t
as_poly_1t <- function(p, t, z) {
  fr <- frac_in(p, t, cancel = TRUE)
  pa <- fr[[1L]]; pd <- fr[[2L]]
  if (!P_is_monomial(pd)) stop("not an element of K[t, 1/t]")
  qr <- P_div(pa, pd)
  ans <- P_expr(qr[[1L]])
  rem <- qr[[2L]]
  if (!P_is_zero(rem)) {
    r <- P_deg(pd) - P_deg(rem)
    # remainder.transform(1, t) * t^r: the coefficients reversed
    n <- P_deg(rem)
    zpart <- .add(lapply(0:n, function(k) .mul(list(P_nth(rem, k), .pow(z, .int(n - k + r))))))
    zpart <- .div(zpart, P_LC(pd))
    ans <- .add(list(ans, expand(zpart)))
  }
  ans
}

derivation <- function(p, DE, coefficientD = FALSE, basic = FALSE) {
  t <- DE$t
  r <- if (basic) S_Zero else .P(S_Zero, t)
  if (coefficientD) {
    if (DE$level <= -length(DE$T)) return(r)
    .de_decrement(DE)
    on.exit(.de_increment(DE))
  }
  n <- .de_idx(DE)
  pe <- if (inherits(p, "dpoly")) P_expr(p) else as_sym(p)
  for (k in seq_len(n)) {
    d <- DE$D[[k]]; v <- DE$T[[k]]
    pv <- if (basic) NULL else P_or_null(pe, v)
    dv <- if (is.null(pv)) .diff1(pe, v) else P_expr(P_diff(pv))
    if (basic) r <- .add(list(r, .mul(list(P_expr(d), dv))))
    else r <- P_add(r, .P(.mul(list(P_expr(d), dv)), t))
  }
  if (basic) r <- cancel(r)
  r
}

get_case <- function(d, t) {
  if (!.has1(P_expr(d), t)) return(if (P_is_one(d)) "base" else "primitive")
  if (P_is_zero(P_rem(d, .P(t, t)))) return("exp")
  if (P_is_zero(P_rem(d, .P(.add(list(S_One, .sq(t))), t)))) return("tan")
  if (P_deg(d) > 1) return("other_nonlinear")
  "other_linear"
}

# gcd of two expressions as polynomials in the generators gens (which may be
# reciprocals 1/v of symbols)
.ri_gens_gcd <- function(a, b, gens) {
  ws <- list(); vs <- list()
  for (g in gens) {
    if (g$type == "Pow" && is_const_num(g$args[[2L]], -1)) { w <- Dummy("w"); ws[[length(ws) + 1L]] <- list(g$args[[1L]], w); vs[[length(vs) + 1L]] <- w }
    else vs[[length(vs) + 1L]] <- g
  }
  tr <- function(e) {
    e <- expand(as_sym(e))
    for (pr in ws) e <- .map_bottom_up(e, function(z) {
      if (z$type == "Pow" && identical(z$args[[1L]]$key, pr[[1L]]$key) && z$args[[2L]]$type == "Integer" && z$args[[2L]]$value < 0)
        return(.pow(pr[[2L]], .neg(z$args[[2L]])))
      z
    })
    e
  }
  A <- tr(a); B <- tr(b)
  na <- .as_numer_denom(together(A)); nb <- .as_numer_denom(together(B))
  # the other generators of the domain: symbols and constant generators such
  # as exp(1/2) (construct_domain)
  others <- .poly_sort_gens(Filter(function(s) !any(vapply(vs, function(v) identical(v$key, s$key), NA)),
                                   .find_gens(lapply(c(na, nb), expand))))
  allg <- c(vs, others)
  field <- .has_any_of(na[[2L]], others) || .has_any_of(nb[[2L]], others)
  if (field) {
    # the ground domain is a fraction field: the monic gcd in the generators
    h <- .mpoly_gcd(expand(na[[1L]]), expand(nb[[1L]]), allg, FALSE)
    hp <- expr_to_qpoly(h, allg, strict = FALSE)
    nv <- length(vs)
    if (!qp_is_zero(hp) && nv) {
      lead <- hp$e[1L, seq_len(nv)]
      rows <- which(apply(hp$e[, seq_len(nv), drop = FALSE], 1L, function(r) all(r == lead)))
      lc <- qpoly_to_expr(qp(hp$e[rows, , drop = FALSE] * rep(c(rep(0L, nv), rep(1L, length(others))), each = length(rows)), hp$c[rows]), allg)
      h <- expand(cancel(.div(h, lc)))
    }
  } else {
    isint <- function(e) { p <- tryCatch(expr_to_qpoly(e, allg, strict = FALSE), error = function(err) NULL); !is.null(p) && all(gmp::denominator(p$c) == 1) }
    h <- .mpoly_gcd(A, B, allg, !(isint(A) && isint(B)))
  }
  for (pr in ws) h <- subs(h, pr[[2L]], .inv(pr[[1L]]))
  h
}

splitfactor <- function(p, DE, coefficientD = FALSE, z = NULL) {
  kinv <- lapply(DE$T[seq_len(.de_idx(DE) - 1L)], .inv)
  if (!is.null(z)) kinv <- c(kinv, list(z))
  One <- .P(S_One, DE$t, domain = p$K)
  Dp <- derivation(p, DE, coefficientD = coefficientD)
  if (P_is_zero(p)) return(list(p, One))
  if (!.has1(P_expr(p), DE$t)) {
    s <- .P(.ri_gens_gcd(P_expr(p), P_expr(Dp), kinv), DE$t)
    n <- P_exquo(p, s)
    return(list(n, s))
  }
  if (!P_is_zero(Dp)) {
    h <- P_to_field(P_gcd(p, Dp))
    g <- P_to_field(P_gcd(p, P_diff(p)))
    s <- P_exquo(h, g)
    if (P_deg(s) == 0) return(list(p, One))
    q <- splitfactor(P_exquo(p, s), DE, coefficientD = coefficientD)
    return(list(q[[1L]], P_mul(q[[2L]], s)))
  }
  list(p, One)
}

splitfactor_sqf <- function(p, DE, coefficientD = FALSE, z = NULL, basic = FALSE) {
  kkinv <- c(lapply(DE$T[seq_len(.de_idx(DE) - 1L)], .inv), DE$T[seq_len(.de_idx(DE) - 1L)])
  if (!is.null(z)) kkinv <- list(z)
  S <- list(); N <- list()
  if (P_is_zero(p)) return(list(list(list(p, 1L)), list()))
  for (pr in P_sqf_list_include(p)) {
    pi <- pr[[1L]]; i <- pr[[2L]]
    dpi <- derivation(pi, DE, coefficientD = coefficientD, basic = basic)
    Si <- .P(.ri_gens_gcd(P_expr(pi), if (inherits(dpi, "dpoly")) P_expr(dpi) else dpi, kkinv), DE$t)
    pi <- .P(P_expr(pi), DE$t)
    Ni <- P_exquo(pi, Si)
    if (!P_is_one(Si)) S[[length(S) + 1L]] <- list(Si, i)
    if (!P_is_one(Ni)) N[[length(N) + 1L]] <- list(Ni, i)
  }
  list(N, S)
}

canonical_representation <- function(a, d, DE) {
  l <- .P(.inv(P_LC(d)), DE$t)
  a <- P_mul(a, l); d <- P_mul(d, l)
  qr <- P_div(a, d)
  sp <- splitfactor(d, DE)
  dn <- sp[[1L]]; ds <- sp[[2L]]
  bc <- gcdex_diophantine(.P(dn, DE$t), .P(ds, DE$t), .P(qr[[2L]], DE$t))
  list(qr[[1L]], list(.P(bc[[1L]], DE$t), ds), list(.P(bc[[2L]], DE$t), dn))
}

hermite_reduce <- function(a, d, DE) {
  l <- .P(.inv(P_LC(d)), DE$t)
  a <- P_mul(a, l); d <- P_mul(d, l)
  cr <- canonical_representation(a, d, DE)
  fp <- cr[[1L]]; fs <- cr[[2L]]; fn <- cr[[3L]]
  a <- fn[[1L]]; d <- fn[[2L]]
  l <- .P(.inv(P_LC(d)), DE$t)
  a <- P_mul(a, l); d <- P_mul(d, l)
  ga <- .P(S_Zero, DE$t); gd <- .P(S_One, DE$t)
  dd <- derivation(d, DE)
  dm <- .P(P_gcd(P_to_field(d), P_to_field(dd)), DE$t)
  ds <- P_div(d, dm)[[1L]]
  while (P_deg(dm) > 0) {
    ddm <- derivation(dm, DE)
    dm2 <- P_gcd(P_to_field(dm), P_to_field(ddm))
    dms <- P_div(dm, dm2)[[1L]]
    ds_ddm <- P_mul(ds, ddm)
    ds_ddm_dm <- P_div(ds_ddm, dm)[[1L]]
    bc <- gcdex_diophantine(.P(P_neg(ds_ddm_dm), DE$t), .P(dms, DE$t), .P(a, DE$t))
    b <- .P(bc[[1L]], DE$t); c <- .P(bc[[2L]], DE$t)
    db <- .P(derivation(b, DE), DE$t)
    ds_dms <- P_div(ds, dms)[[1L]]
    a <- P_sub(.P(c, DE$t), .P(P_mul(db, ds_dms), DE$t))
    ga <- P_add(P_mul(ga, dm), P_mul(b, gd))
    gd <- P_mul(gd, dm)
    gg <- P_cancel(ga, gd, include = TRUE); ga <- gg[[1L]]; gd <- gg[[2L]]
    dm <- dm2
  }
  qr <- P_div(a, ds)
  gg <- P_cancel(ga, gd, include = TRUE); ga <- gg[[1L]]; gd <- gg[[2L]]
  rd <- P_cancel(qr[[2L]], ds, include = TRUE)
  rra <- P_add(P_add(P_mul(qr[[1L]], fs[[2L]]), P_mul(fp, fs[[2L]])), fs[[1L]])
  rrd <- fs[[2L]]
  rr <- P_cancel(rra, rrd, include = TRUE)
  list(list(ga, gd), list(rd[[1L]], rd[[2L]]), list(rr[[1L]], rr[[2L]]))
}

polynomial_reduce <- function(p, DE) {
  q <- .P(S_Zero, DE$t)
  while (P_deg(p) >= P_deg(DE$d)) {
    m <- P_deg(p) - P_deg(DE$d) + 1L
    q0 <- P_mul(.P(.pow(DE$t, .int(m)), DE$t), .P(.div(P_LC(.P(P_expr(p), DE$t)), .mul(list(.int(m), P_LC(DE$d)))), DE$t))
    q <- P_add(q, q0)
    p <- P_sub(p, derivation(q0, DE))
  }
  list(q, p)
}

recognize_log_derivative <- function(a, d, DE, z = NULL) {
  if (is.null(z)) z <- Dummy("z")
  cc <- P_cancel(a, d, include = TRUE); a <- cc[[1L]]; d <- cc[[2L]]
  a <- P_div(a, d)[[2L]]
  pz <- .P(z, DE$t)
  Dd <- derivation(d, DE)
  q <- P_sub(a, P_mul(pz, Dd))
  r <- .P(P_resultant(d, q, includePRS = TRUE)[[1L]], z)
  sp <- splitfactor_sqf(r, DE, coefficientD = TRUE, z = z)
  if (any(vapply(sp[[1L]], function(s) P_deg(.P(P_expr(s[[1L]]), z)) > 0, NA))) return(FALSE)
  for (s in sp[[2L]]) {
    sz <- .P(P_expr(s[[1L]]), z)
    rs <- roots(P_expr(sz), z)
    if (sum(attr(rs, "multiplicity")) != P_deg(sz) || !all(vapply(rs, function(j) j$type == "Integer", NA))) return(FALSE)
  }
  TRUE
}

residue_reduce <- function(a, d, DE, z = NULL, invert = TRUE) {
  if (is.null(z)) z <- Dummy("z")
  cc <- P_cancel(a, d, include = TRUE); a <- cc[[1L]]; d <- cc[[2L]]
  lc <- P_LC(d)
  a <- P_mul_ground(P_to_field(a), .inv(lc)); d <- P_mul_ground(P_to_field(d), .inv(lc))
  kkinv <- c(lapply(DE$T[seq_len(.de_idx(DE) - 1L)], .inv), DE$T[seq_len(.de_idx(DE) - 1L)])
  if (P_is_zero(a)) return(list(list(), TRUE))
  a <- P_div(a, d)[[2L]]
  pz <- .P(z, DE$t)
  Dd <- derivation(d, DE)
  q <- P_sub(a, P_mul(pz, Dd))
  rr <- if (P_deg(Dd) <= P_deg(d)) P_resultant(d, q, includePRS = TRUE) else P_resultant(q, d, includePRS = TRUE)
  r <- rr[[1L]]; R <- rr[[2L]]
  Rmap <- list()
  for (i in R) Rmap[[as.character(P_deg(i))]] <- i
  H <- list()
  r <- .P(r, z)
  sp <- splitfactor_sqf(r, DE, coefficientD = TRUE, z = z)
  Np <- sp[[1L]]; Sp <- sp[[2L]]
  for (si in Sp) {
    s <- si[[1L]]; i <- si[[2L]]
    if (i == P_deg(d)) {
      s <- P_monic(.P(P_expr(s), z))
      H[[length(H) + 1L]] <- list(s, d)
    } else {
      h <- Rmap[[as.character(i)]]
      if (is.null(h)) next
      h_lc <- .P(P_LC(.P(h, DE$t)), DE$t, field = TRUE)
      for (aj in P_sqf_list_include(h_lc, all = TRUE)) {
        g <- .ri_gens_gcd(P_expr(aj[[1L]]), P_expr(P_pow(s, aj[[2L]])), c(list(z), kkinv))
        h <- P_exquo(.P(h, DE$t, field = TRUE), .P(g, DE$t))
      }
      s <- P_monic(.P(P_expr(s), z))
      if (invert) {
        h_lc <- .P(P_LC(.P(h, DE$t)), DE$t, field = TRUE)
        inv <- P_invert(.P(P_expr(h_lc), z, field = TRUE), s)
        coeffs <- list(S_One)
        for (cf in P_coeffs(h)[-1L]) {
          L <- P_rem(P_mul(inv, .P(cf, z)), s)
          coeffs[[length(coeffs) + 1L]] <- P_expr(L)
        }
        ms <- P_monoms(h)
        d0 <- list(); for (k in seq_along(ms)) d0[[as.character(ms[k])]] <- coeffs[[k]]
        h <- P_from_dict(d0, DE$t)
      }
      H[[length(H) + 1L]] <- list(s, h)
    }
  }
  b <- !any(vapply(Np, function(pr) { ce <- cancel(P_expr(pr[[1L]])); .has1(ce, DE$t) || .has1(ce, z) }, NA))
  list(H, b)
}

residue_reduce_to_basic <- function(H, DE, z) {
  i <- Dummy("i")
  s <- .de_s(DE)
  .add(lapply(H, function(a) RootSum(P_expr(a[[1L]]), Lambda(i, .mul(list(i, .subs_seq(subs(fcall("log", list(P_expr(a[[2L]]))), z, i), s)))), z)))
}

residue_reduce_derivation <- function(H, DE, z) {
  i <- Dummy("i")
  .add(lapply(H, function(a) RootSum(P_expr(a[[1L]]), Lambda(i, .mul(list(i, .div(subs(P_expr(derivation(a[[2L]], DE)), z, i), subs(P_expr(a[[2L]]), z, i))))), z)))
}

integrate_primitive_polynomial <- function(p, DE) {
  Zero <- .P(S_Zero, DE$t)
  q <- .P(S_Zero, DE$t)
  if (!.has1(P_expr(p), DE$t)) return(list(Zero, p, TRUE))
  repeat {
    if (!.has1(P_expr(p), DE$t)) return(list(q, p, TRUE))
    Dt <- frac_in(DE$d, DE$T[[.de_idx(DE) - 1L]])
    res <- .with_dec(DE, {
      a <- P_LC(p)
      af <- frac_in(a, DE$t)
      tryCatch({
        rv <- limited_integrate(af[[1L]], af[[2L]], list(list(Dt[[1L]], Dt[[2L]])), DE)
        if (is.null(rv)) .nonelem()
        rv
      }, NonElementaryIntegralException = function(e) NULL)
    })
    if (is.null(res)) return(list(q, p, FALSE))
    ba <- res[[1L]][[1L]]; bd <- res[[1L]][[2L]]; cc <- res[[2L]]
    m <- P_deg(p)
    q0 <- P_add(P_mul(.P(cc[[1L]], DE$t), .P(.div(.pow(DE$t, .int(m + 1L)), .int(m + 1L)), DE$t)),
                P_mul(.P(.div(P_expr(ba), P_expr(bd)), DE$t), .P(.pow(DE$t, .int(m)), DE$t)))
    p <- P_sub(p, derivation(q0, DE))
    q <- P_add(q, q0)
  }
}

integrate_primitive <- function(a, d, DE, z = NULL) {
  if (is.null(z)) z <- Dummy("z")
  s <- .de_s(DE)
  hr <- hermite_reduce(a, d, DE)
  g1 <- hr[[1L]]; h <- hr[[2L]]; r <- hr[[3L]]
  rr <- residue_reduce(h[[1L]], h[[2L]], DE, z = z)
  g2 <- rr[[1L]]; b <- rr[[2L]]
  if (!b) {
    i <- cancel(.sub(.sub(.div(P_expr(a), P_expr(d)),
                          .div(P_expr(P_sub(P_mul(g1[[2L]], derivation(g1[[1L]], DE)), P_mul(g1[[1L]], derivation(g1[[2L]], DE)))), P_expr(P_pow(g1[[2L]], 2L)))),
                     residue_reduce_derivation(g2, DE, z)))
    i <- NonElementaryIntegral(.subs_seq(cancel(i), s), DE$x)
    return(list(.add(list(.subs_seq(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), s), residue_reduce_to_basic(g2, DE, z))), i, b))
  }
  p <- cancel(.add(list(.div(P_expr(h[[1L]]), P_expr(h[[2L]])), .neg(residue_reduce_derivation(g2, DE, z)), .div(P_expr(r[[1L]]), P_expr(r[[2L]])))))
  p <- .P(p, DE$t)
  ip <- integrate_primitive_polynomial(p, DE)
  q <- ip[[1L]]; i <- ip[[2L]]; b <- ip[[3L]]
  ret <- .add(list(.subs_seq(.add(list(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), P_expr(q))), s), residue_reduce_to_basic(g2, DE, z)))
  i <- if (!b) NonElementaryIntegral(.subs_seq(cancel(P_expr(i)), s), DE$x) else cancel(P_expr(i))
  list(ret, i, b)
}

# p is an expression in t and z = 1/t
integrate_hyperexponential_polynomial <- function(p, DE, z) {
  t1 <- DE$t
  dtt <- P_exquo(DE$d, .P(DE$t, DE$t))
  qa <- .P(S_Zero, DE$t); qd <- .P(S_One, DE$t)
  b <- TRUE
  if (is_zero_num(p)) return(list(qa, qd, b))
  pz <- .P(p, z); pt <- .P(p, t1)
  degz <- .ri_deg_in(p, z); degt <- .ri_deg_in(p, t1)
  .de_decrement(DE)
  on.exit(.de_increment(DE))
  for (i in (-degz):degt) {
    if (i == 0L) next
    a <- if (i < 0L) .ri_nth_in(p, z, -i) else .ri_nth_in(p, t1, i)
    af <- frac_in(a, DE$t, field = TRUE)
    cc <- P_cancel(af[[1L]], af[[2L]], include = TRUE)
    aa <- cc[[1L]]; ad <- cc[[2L]]
    iDt <- P_mul(.P(.int(i), t1), dtt)
    idf <- frac_in(iDt, DE$t, field = TRUE)
    rv <- tryCatch(rischDE(idf[[1L]], idf[[2L]], .P(P_expr(aa), DE$t), .P(P_expr(ad), DE$t), DE), NonElementaryIntegralException = function(e) NULL)
    if (is.null(rv)) { b <- FALSE; next }
    vv <- frac_in(list(rv[[1L]], rv[[2L]]), t1, cancel = TRUE)
    va <- vv[[1L]]; vd <- vv[[2L]]
    ti <- if (i > 0L) .P(.pow(t1, .int(i)), t1) else .P(.pow(z, .int(-i)), z)
    qa <- .P(.add(list(.mul(list(P_expr(qa), P_expr(vd))), .mul(list(P_expr(va), P_expr(ti), P_expr(qd))))), t1)
    qd <- P_mul(qd, vd)
  }
  list(qa, qd, b)
}
.ri_deg_in <- function(e, v) { e <- expand(e); if (!.has1(e, v)) 0L else degree(e, v) }
.ri_nth_in <- function(e, v, k) { cl <- .coeff_list(expand(e), v); n <- length(cl) - 1L; if (k > n) S_Zero else cl[[n - k + 1L]] }

integrate_hyperexponential <- function(a, d, DE, z = NULL, conds = "piecewise") {
  if (is.null(z)) z <- Dummy("z")
  s <- c(list(list(z, .inv(DE$t))), .de_s(DE))
  hr <- hermite_reduce(a, d, DE)
  g1 <- hr[[1L]]; h <- hr[[2L]]; r <- hr[[3L]]
  rr <- residue_reduce(h[[1L]], h[[2L]], DE, z = z)
  g2 <- rr[[1L]]; b <- rr[[2L]]
  if (!b) {
    i <- cancel(.sub(.sub(.div(P_expr(a), P_expr(d)),
                          .div(P_expr(P_sub(P_mul(g1[[2L]], derivation(g1[[1L]], DE)), P_mul(g1[[1L]], derivation(g1[[2L]], DE)))), P_expr(P_pow(g1[[2L]], 2L)))),
                     residue_reduce_derivation(g2, DE, z)))
    i <- NonElementaryIntegral(cancel(.subs_seq(i, s)), DE$x)
    return(list(.add(list(.subs_seq(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), s), residue_reduce_to_basic(g2, DE, z))), i, b))
  }
  p <- cancel(.add(list(.div(P_expr(h[[1L]]), P_expr(h[[2L]])), .neg(residue_reduce_derivation(g2, DE, z)), .div(P_expr(r[[1L]]), P_expr(r[[2L]])))))
  pp <- as_poly_1t(p, DE$t, z)
  hp <- integrate_hyperexponential_polynomial(pp, DE, z)
  qa <- hp[[1L]]; qd <- hp[[2L]]; b <- hp[[3L]]
  i <- .ri_nth_in(.ri_nth_in(pp, z, 0L), DE$t, 0L)
  ret <- .add(list(.subs_seq(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), s), residue_reduce_to_basic(g2, DE, z)))
  qas <- .subs_seq(P_expr(qa), s); qds <- .subs_seq(P_expr(qd), s)
  if (conds == "piecewise" && !.has1(qds, DE$x)) {
    ret <- .add(list(ret, Piecewise(list(.div(qas, qds), Ne(qds, 0)),
                                    list(integrate(.subs_seq(subs(.sub(p, i), DE$t, S_One), s), DE$x), S_true))))
  } else ret <- .add(list(ret, .div(qas, qds)))
  if (!b) {
    i <- .sub(p, .div(P_expr(P_sub(P_mul(qd, derivation(qa, DE)), P_mul(qa, derivation(qd, DE)))), P_expr(P_pow(qd, 2L))))
    i <- NonElementaryIntegral(.subs_seq(cancel(i), s), DE$x)
  }
  list(ret, i, b)
}

integrate_hypertangent_polynomial <- function(p, DE) {
  qr <- polynomial_reduce(p, DE)
  a <- P_exquo(DE$d, .P(.add(list(.sq(DE$t), S_One)), DE$t))
  c <- .P(.div(P_nth(qr[[2L]], 1L), .mul(list(S_Two, P_expr(a)))), DE$t)
  list(qr[[1L]], c)
}

integrate_nonlinear_no_specials <- function(a, d, DE, z = NULL) {
  if (is.null(z)) z <- Dummy("z")
  s <- .de_s(DE)
  hr <- hermite_reduce(a, d, DE)
  g1 <- hr[[1L]]; h <- hr[[2L]]; r <- hr[[3L]]
  rr <- residue_reduce(h[[1L]], h[[2L]], DE, z = z)
  g2 <- rr[[1L]]; b <- rr[[2L]]
  if (!b) return(list(.add(list(.subs_seq(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), s), residue_reduce_to_basic(g2, DE, z))), b))
  p <- .P(cancel(.add(list(.div(P_expr(h[[1L]]), P_expr(h[[2L]])), .neg(residue_reduce_derivation(g2, DE, z)), .div(P_expr(r[[1L]]), P_expr(r[[2L]]))))), DE$t)
  q12 <- polynomial_reduce(p, DE)
  b <- !.has1(P_expr(q12[[2L]]), DE$t)
  ret <- .add(list(.subs_seq(cancel(.add(list(.div(P_expr(g1[[1L]]), P_expr(g1[[2L]])), P_expr(q12[[1L]])))), s), residue_reduce_to_basic(g2, DE, z)))
  list(ret, b)
}

# NonElementaryIntegral: an Integral that prints and behaves as one
NonElementaryIntegral <- function(f, ...) {
  r <- Integral(f, ...)
  new_sym("Integral", r$args, extra = "NonElementary", key = paste0("NonElementary", r$key))
}

#' The Risch algorithm
#'
#' A port of SymPy's `risch_integrate()` (sympy/integrals/risch.py, rde.py
#' and prde.py) for elementary transcendental functions: exponentials and
#' logarithms over a tower built by `DifferentialExtension()`. Returns the
#' antiderivative; a part proven nonelementary is an unevaluated integral.
#'
#' @param f The integrand.
#' @param x The variable.
#' @param extension A `DifferentialExtension` to use.
#' @param handle_first `"log"` or `"exp"`: which extensions are built first.
#' @param separate_integral Return `list(elementary, nonelementary)`.
#' @param rewrite_complex Rewrite trigonometric functions with `exp`.
#' @param conds `"piecewise"` for conditional results.
#' @return An expression or a list.
#' @examples
#' x <- Symbol("x")
#' risch_integrate(exp(x^2)*x, x)
#' @export
risch_integrate <- function(f, x, extension = NULL, handle_first = "log", separate_integral = FALSE,
                            rewrite_complex = NULL, conds = "piecewise") {
  f <- as_sym(f); x <- as_sym(x)
  DE <- if (!is.null(extension)) extension else DifferentialExtension(f, x, handle_first = handle_first, dummy = TRUE, rewrite_complex = rewrite_complex)
  fa <- DE$fa; fd <- DE$fd
  result <- S_Zero
  for (case in rev(unlist(DE$cases))) {
    if (!.has1(P_expr(fa), DE$t) && !.has1(P_expr(fd), DE$t) && case != "base") {
      .de_decrement(DE)
      fr <- frac_in(list(fa, fd), DE$t)
      fa <- fr[[1L]]; fd <- fr[[2L]]
      next
    }
    cc <- P_cancel(fa, fd, include = TRUE); fa <- cc[[1L]]; fd <- cc[[2L]]
    if (case == "exp") { r <- integrate_hyperexponential(fa, fd, DE, conds = conds); ans <- r[[1L]]; i <- r[[2L]]; b <- r[[3L]] }
    else if (case == "primitive") { r <- integrate_primitive(fa, fd, DE); ans <- r[[1L]]; i <- r[[2L]]; b <- r[[3L]] }
    else if (case == "base") { ans <- integrate(.div(P_expr(fa), P_expr(fd)), DE$x, risch = FALSE); b <- FALSE; i <- S_Zero }
    else .ri_notimpl("Only exponential and logarithmic extensions are currently supported.")
    result <- .add(list(result, ans))
    if (b) {
      .de_decrement(DE)
      fr <- frac_in(i, DE$t)
      fa <- fr[[1L]]; fd <- fr[[2L]]
    } else {
      result <- .de_backsubs_apply(result, DE)
      if (!is_zero_num(i)) i <- NonElementaryIntegral(.de_backsubs_apply(i$args[[1L]], DE), i$args[[2L]])
      if (!separate_integral) return(.add(list(result, i)))
      if (identical(i$extra, "NonElementary")) return(list(result, i))
      return(list(result, S_Zero))
    }
  }
}
