# Canonical constructors for Add, Mul and Pow, following sympy/core/add.py,
# sympy/core/mul.py and sympy/core/power.py.

.nc_state <- new.env()
.nc_state$any <- FALSE

.is_comm <- function(e) {
  if (!.nc_state$any) return(TRUE)
  !identical(ask_fact(e, "commutative"), FALSE)
}

# Split a term into numeric coefficient and the rest: 3*x*y -> (3, x*y).
as_coeff_mul_split <- function(a) {
  if (a$type == "Mul" && is_number_type(a$args[[1L]])) {
    rest <- a$args[-1L]
    list(a$args[[1L]],
         if (length(rest) == 1L) rest[[1L]] else new_sym("Mul", rest))
  } else if (is_number_type(a)) {
    list(a, S_One)
  } else {
    list(S_One, a)
  }
}

# Rebuild c*t for canonical t without a numeric coefficient.
.coeff_times <- function(cf, t) {
  if (t$type == "Integer" && t$value == 1) return(cf)
  if (cf$type != "Float" && num_is_one(cf)) return(t)
  if (num_is_zero(cf)) return(S_Zero)
  if (t$type == "Mul") return(new_sym("Mul", c(list(cf), t$args)))
  if (t$type == "Add") return(.mul(list(cf, t)))
  new_sym("Mul", list(cf, t))
}

.sort_by_key <- function(lst) {
  if (length(lst) < 2L) return(lst)
  lst[order(vapply(lst, .k, ""), method = "radix")]
}

# ---------------------------------------------------------------------------
# Add

# AssocOp.__new__ (sympy/core/operations.py:74)
.rel_in_assoc <- function(cls) stop(.sympy_condition("TypeError", paste0("Relational cannot be used in ", cls)))

.add <- function(args) {
  n <- length(args)
  if (n == 0L) return(S_Zero)
  if (n == 1L) return(args[[1L]])
  types <- vapply(args, `[[`, "", "type")
  if (any(types %in% .REL_TYPES)) .rel_in_assoc("Add")
  isadd <- types == "Add"
  if (any(isadd)) {
    # Flatten, reusing the type and key vectors cached on Add nodes.
    parts_args <- vector("list", n)
    parts_types <- vector("list", n)
    parts_keys <- vector("list", n)
    for (i in seq_len(n)) {
      a <- args[[i]]
      if (isadd[i]) {
        parts_args[[i]] <- a$args
        parts_types[[i]] <- a$extra$types
        parts_keys[[i]] <- a$extra$keys
      } else {
        parts_args[[i]] <- list(a)
        parts_types[[i]] <- types[i]
        parts_keys[[i]] <- a$key
      }
    }
    args <- do.call(c, parts_args)
    types <- unlist(parts_types, use.names = FALSE)
    keys <- unlist(parts_keys, use.names = FALSE)
  } else {
    keys <- vapply(args, .k, "")
  }
  if (any(types == "Order")) return(.add_with_order(args, types))
  if (any(types == "AccumulationBounds")) {
    # bounds absorb the numbers of the sum, as Add.flatten does
    fo <- .ab_fold(args, "add")
    if (!is.null(fo)) {
      if (!length(fo$rest)) return(fo$acc)
      if (!.is_ab(fo$acc)) return(.add(c(fo$rest, list(fo$acc))))
      args <- c(fo$rest, list(fo$acc))
      types <- vapply(args, `[[`, "", "type")
      keys <- vapply(args, .k, "")
    }
  }
  isnum <- types == "Integer" | types == "Rational" | types == "Float"
  coeff <- S_Zero
  if (any(isnum)) for (a in args[isnum]) coeff <- num_add(coeff, a)
  infs <- character()
  isconst <- types == "Constant"
  if (any(isconst)) {
    cn <- vapply(args[isconst], `[[`, "", "name")
    if (any(cn == "nan")) return(S_NaN)
    special <- cn %in% .INF_NAMES
    infs <- cn[special]
    isnum[which(isconst)[special]] <- TRUE
  }
  terms <- args[!isnum]
  ttypes <- types[!isnum]
  tkeys <- keys[!isnum]
  inf_term <- NULL
  if (length(infs)) {
    u <- unique(infs)
    if (length(u) > 1L || (u == "zoo" && length(infs) > 1L)) return(S_NaN)
    inf_term <- new_constant(u)
    coeff <- S_Zero
  }
  m <- length(terms)
  if (m > 1L) {
    rk <- tkeys
    cs <- NULL
    ts <- NULL
    ismul <- which(ttypes == "Mul")
    if (length(ismul)) {
      lead <- ismul[vapply(terms[ismul], function(a) {
        t1 <- a$args[[1L]]$type
        t1 == "Integer" || t1 == "Rational" || t1 == "Float"
      }, TRUE)]
      if (length(lead)) {
        cs <- vector("list", m)
        ts <- terms
        for (i in lead) {
          sp <- as_coeff_mul_split(terms[[i]])
          cs[[i]] <- sp[[1L]]
          ts[[i]] <- sp[[2L]]
          rk[i] <- sp[[2L]]$key
        }
      }
    }
    if (anyDuplicated(rk)) {
      if (is.null(cs)) {
        cs <- vector("list", m)
        ts <- terms
      }
      u <- unique(rk)
      idx <- match(rk, u)
      dup <- which(tabulate(idx, length(u)) > 1L)
      keepi <- !(idx %in% dup)
      newterms <- list()
      for (g in dup) {
        w <- which(idx == g)
        cf <- NULL
        for (j in w) {
          cj <- if (is.null(cs[[j]])) S_One else cs[[j]]
          cf <- if (is.null(cf)) cj else num_add(cf, cj)
        }
        if (num_is_zero(cf)) next
        newterms[[length(newterms) + 1L]] <- .coeff_times(cf, ts[[w[1L]]])
      }
      ntypes <- vapply(newterms, `[[`, "", "type")
      if (length(newterms) && any(ntypes %in% c("Add", .NUMBER_TYPES))) {
        extra <- c(terms[keepi], newterms)
        if (!is.null(inf_term)) extra <- c(extra, list(inf_term))
        return(.add(c(list(coeff), extra)))
      }
      terms <- c(terms[keepi], newterms)
      ttypes <- c(ttypes[keepi], ntypes)
      tkeys <- c(tkeys[keepi], vapply(newterms, .k, ""))
    }
  }
  # oo absorbs terms known to be real or nonnegative, -oo those real or
  # nonpositive, and zoo finite terms of known reality (Add.flatten).
  if (!is.null(inf_term) && length(terms)) {
    drop <- vapply(terms, function(t) {
      f <- .facts(t)
      switch(inf_term$name,
             oo = isTRUE(f[["nonnegative"]]) || isTRUE(f[["real"]]),
             "-oo" = isTRUE(f[["nonpositive"]]) || isTRUE(f[["real"]]),
             zoo = isTRUE(f[["finite"]]) && (!is.na(f[["extended_real"]]) || !is.na(f[["real"]])),
             FALSE)
    }, TRUE)
    terms <- terms[!drop]
    ttypes <- ttypes[!drop]
    tkeys <- tkeys[!drop]
  }
  if (length(terms) > 1L) {
    o <- order(tkeys, method = "radix")
    terms <- terms[o]
    ttypes <- ttypes[o]
    tkeys <- tkeys[o]
  }
  head <- if (!is.null(inf_term)) inf_term
          else if (!num_is_zero(coeff) || (coeff$type == "Float" && !length(terms))) coeff
          else NULL
  if (!is.null(head)) {
    terms <- c(list(head), terms)
    ttypes <- c(head$type, ttypes)
    tkeys <- c(head$key, tkeys)
  }
  if (!length(terms)) return(S_Zero)
  if (length(terms) == 1L) return(terms[[1L]])
  new_sym("Add", terms, extra = list(types = ttypes, keys = tkeys),
          key = paste0("Add(", paste(tkeys, collapse = ","), ")"))
}

# ---------------------------------------------------------------------------
# Mul

.base_exp <- function(f) {
  if (f$type == "Pow") return(list(f$args[[1L]], f$args[[2L]]))
  if (f$type == "Function" && f$name == "exp") return(list(S_E, f$args[[1L]]))
  # ExpBase.as_base_exp: exp_polar(a) = exp_polar(1)^a
  # (exp_polar(0) is special: it is its own base)
  if (f$type == "Function" && f$name == "exp_polar" && !is_one_num(f$args[[1L]]) && !is_zero_num(f$args[[1L]])) return(list(.EXP_POLAR_1, f$args[[1L]]))
  list(f, S_One)
}
.EXP_POLAR_1 <- new_sym("Function", list(new_sym("Integer", value = gmp::as.bigq(1), key = "Z1")), name = "exp_polar",
                        key = "exp_polar(Z1)")

.mul <- function(args, depth = 0L) {
  n <- length(args)
  if (n == 0L) return(S_One)
  if (n == 1L) return(args[[1L]])
  types <- vapply(args, `[[`, "", "type")
  if (any(types %in% .REL_TYPES)) .rel_in_assoc("Mul")
  if (any(types == "Mul")) {
    args <- do.call(c, lapply(args, function(a)
      if (a$type == "Mul") a$args else list(a)))
    types <- vapply(args, `[[`, "", "type")
  }
  if (any(types == "Order")) return(.mul_with_order(args, types))
  if (any(types == "AccumulationBounds")) {
    # bounds absorb the numbers of the product, as Mul.flatten does
    fo <- .ab_fold(args, "mul")
    if (!is.null(fo)) {
      if (!length(fo$rest)) return(fo$acc)
      if (!.is_ab(fo$acc)) return(.mul(c(fo$rest, list(fo$acc)), depth))
      args <- c(fo$rest, list(fo$acc))
      types <- vapply(args, `[[`, "", "type")
    }
  }
  coeff <- S_One
  infs <- character()
  cfac <- list()
  ncfac <- list()
  for (i in seq_along(args)) {
    a <- args[[i]]
    t <- types[i]
    if (t == "Integer" || t == "Rational" || t == "Float") {
      coeff <- num_mul(coeff, a)
    } else if (t == "Constant" && a$name %in% c(.INF_NAMES, "nan")) {
      if (a$name == "nan") return(S_NaN)
      infs <- c(infs, a$name)
    } else if (.is_comm(a)) {
      cfac[[length(cfac) + 1L]] <- if (t == "Pow" || t == "Function") .reeval_exp(a) else a
    } else {
      ncfac[[length(ncfac) + 1L]] <- a
    }
  }
  inf_head <- NULL
  if (length(infs)) {
    if (num_is_zero(coeff)) return(S_NaN)
    if ("zoo" %in% infs) {
      inf_head <- S_ComplexInfinity
    } else {
      s <- num_sign(coeff) * prod(ifelse(infs == "oo", 1, -1))
      inf_head <- if (s > 0) S_Infinity else S_NegativeInfinity
    }
    coeff <- S_One
  } else if (num_is_zero(coeff)) {
    return(S_Zero)
  }

  if (length(cfac) > 1L || (length(cfac) == 1L && depth == 0L &&
                            FALSE)) {
    m <- length(cfac)
    bs <- vector("list", m)
    ecs <- vector("list", m)
    ets <- vector("list", m)
    for (i in seq_len(m)) {
      be <- .base_exp(cfac[[i]])
      bs[[i]] <- be[[1L]]
      sp <- as_coeff_mul_split(be[[2L]])
      ecs[[i]] <- sp[[1L]]
      ets[[i]] <- sp[[2L]]
    }
    gk <- paste(vapply(bs, .k, ""), vapply(ets, .k, ""), sep = "|")
    if (anyDuplicated(gk)) {
      u <- unique(gk)
      idx <- match(gk, u)
      newf <- vector("list", length(u))
      for (g in seq_along(u)) {
        w <- which(idx == g)
        if (length(w) == 1L) {
          newf[[g]] <- cfac[[w]]
          next
        }
        ec <- ecs[[w[1L]]]
        for (j in w[-1L]) ec <- num_add(ec, ecs[[j]])
        newf[[g]] <- .pow(bs[[w[1L]]], .coeff_times(ec, ets[[w[1L]]]))
      }
      cfac <- newf
    }
    cfac <- .combine_num_exp(cfac)
    cfac <- .combine_rational_radicals(cfac)
    ft <- vapply(cfac, `[[`, "", "type")
    if (depth < 8L && (any(ft %in% c("Mul", .NUMBER_TYPES)) ||
                       any(vapply(cfac, is_inf_const, TRUE)) ||
                       any(vapply(cfac, is_nan, TRUE)))) {
      return(.mul(c(list(coeff), if (!is.null(inf_head)) list(inf_head),
                    cfac, ncfac), depth + 1L))
    }
  }

  if (length(ncfac) > 1L) {
    m <- .merge_nc(ncfac)
    ncfac <- m$nc
    # a combined power that commutes joins the commutative factors
    if (length(m$c)) return(.mul(c(list(coeff), if (!is.null(inf_head)) list(inf_head), cfac, m$c, ncfac), depth + 1L))
  }

  if (is.null(inf_head) && !length(ncfac) && length(cfac) == 1L &&
      cfac[[1L]]$type == "Add" && !(coeff$type != "Float" && num_is_one(coeff))) {
    a <- cfac[[1L]]
    return(.add(lapply(a$args, function(t) .mul(list(coeff, t)))))
  }

  if (!is.null(inf_head) && inf_head$name %in% c("oo", "-oo")) {
    # Mul.flatten's _handle_for_oo: factors of known sign are absorbed into oo
    sgn <- if (inf_head$name == "oo") 1L else -1L
    absorb <- function(fs) Filter(function(t) {
      if (isTRUE(ask_fact(t, "positive"))) return(FALSE)
      if (isTRUE(ask_fact(t, "negative"))) { sgn <<- -sgn; return(FALSE) }
      TRUE
    }, fs)
    cfac <- absorb(cfac); ncfac <- absorb(ncfac)
    inf_head <- if (sgn > 0L) S_Infinity else S_NegativeInfinity
  }
  if (!is.null(inf_head) && inf_head$name == "zoo") {
    # Mul.flatten: zoo absorbs the nonzero factors of known reality
    keep <- function(fs) Filter(function(t) !(identical(ask_fact(t, "zero"), FALSE) && !is.na(ask_fact(t, "extended_real"))), fs)
    cfac <- keep(cfac); ncfac <- keep(ncfac)
  }
  cfac <- .sort_by_key(cfac)
  head <- if (!is.null(inf_head)) list(inf_head)
          else if (coeff$type == "Float" || !num_is_one(coeff)) list(coeff)
          else list()
  out <- c(head, cfac, ncfac)
  if (!length(out)) return(S_One)
  if (length(out) == 1L) return(out[[1L]])
  # constructor postprocessors (Basic._constructor_postprocessor_mapping)
  if (length(ncfac) && !is.null(.arith_hooks$mul_post)) return(.arith_hooks$mul_post(out))
  new_sym("Mul", out)
}

# Mul.flatten rebuilds each factor as Pow(b, c*t) from as_coeff_Mul() of its
# exponent, so an unevaluated exponent 2*(1 - a) becomes 2 - 2*a.
.reeval_exp <- function(f) {
  if (f$type == "Pow") x <- f$args[[2L]]
  else if (f$name == "exp") x <- f$args[[1L]]
  else return(f)
  if (x$type != "Mul" || length(x$args) != 2L || !is_number_type(x$args[[1L]]) || x$args[[2L]]$type != "Add") return(f)
  nx <- .mul(x$args)
  if (identical(nx$key, x$key)) return(f)
  if (f$type == "Pow") .pow(f$args[[1L]], nx) else fcall("exp", list(nx))
}

# Mul.flatten's num_exp: powers of numbers with the same non-rational
# exponent combine, 2^x*3^x -> 6^x (a negative base only with an integer
# exponent).
.combine_num_exp <- function(fs) {
  ok <- vapply(fs, function(f) f$type == "Pow" && is_number_type(f$args[[1L]]) &&
                 !(f$args[[2L]]$type %in% c("Integer", "Rational")) &&
                 (num_sign(f$args[[1L]]) > 0 || isTRUE(ask_fact(f$args[[2L]], "integer"))), TRUE)
  if (sum(ok) < 2L) return(fs)
  ek <- vapply(fs[ok], function(f) f$args[[2L]]$key, "")
  if (!anyDuplicated(ek)) return(fs)
  ne <- fs[ok]
  out <- fs[!ok]
  for (u in unique(ek)) {
    w <- which(ek == u)
    b <- ne[[w[1L]]]$args[[1L]]
    for (j in w[-1L]) b <- num_mul(b, ne[[j]]$args[[1L]])
    out[[length(out) + 1L]] <- if (length(w) == 1L) ne[[w]] else .pow(b, ne[[w[1L]]]$args[[2L]])
  }
  out
}

# sqrt(2)*sqrt(3) -> sqrt(6): positive rational bases sharing an exponent.
.combine_rational_radicals <- function(fs) {
  if (length(fs) < 2L) return(fs)
  fs <- .combine_rr_positive(fs)
  # Mul.flatten: (-1)^r joins a positive number raised to the same r
  neg <- which(vapply(fs, function(f) f$type == "Pow" && is_negone_num(f$args[[1L]]) && f$args[[2L]]$type == "Rational", TRUE))
  if (length(neg) == 1L) {
    r <- fs[[neg]]$args[[2L]]
    j <- which(vapply(fs, function(f) f$type == "Pow" && is_exact_number(f$args[[1L]]) && f$args[[1L]]$value > 0 &&
                        identical(f$args[[2L]]$key, r$key), TRUE))
    if (length(j)) {
      j <- j[1L]
      fs[[j]] <- new_sym("Pow", list(make_rational(-fs[[j]]$args[[1L]]$value), r))
      fs <- fs[-neg]
    }
  }
  fs
}
.combine_rr_positive <- function(fs) {
  is_rr <- vapply(fs, function(f) f$type == "Pow" &&
                    is_exact_number(f$args[[1L]]) && f$args[[1L]]$value > 0 &&
                    f$args[[2L]]$type == "Rational", TRUE)
  if (sum(is_rr) < 2L) return(fs)
  ek <- vapply(fs[is_rr], function(f) f$args[[2L]]$key, "")
  if (!anyDuplicated(ek)) return(fs)
  rr <- fs[is_rr]
  out <- fs[!is_rr]
  for (u in unique(ek)) {
    w <- which(ek == u)
    if (length(w) == 1L) {
      out[[length(out) + 1L]] <- rr[[w]]
      next
    }
    prodv <- rr[[w[1L]]]$args[[1L]]$value
    for (j in w[-1L]) prodv <- prodv * rr[[j]]$args[[1L]]$value
    out[[length(out) + 1L]] <- .pow(make_rational(prodv), rr[[w[1L]]]$args[[2L]])
  }
  out
}

.merge_nc <- function(fs) {
  out <- list()
  cpart <- list()
  todo <- fs
  while (length(todo)) {
    f <- todo[[1L]]
    todo <- todo[-1L]
    if (!length(out)) { out <- list(f); next }
    last <- out[[length(out)]]
    b1 <- .base_exp(last)
    b2 <- .base_exp(f)
    # Mul.flatten: noncommutative powers combine unless the new exponent is
    # a sum (a^x*a^y stays, a^2*a^3 = a^5); a commutative result leaves the
    # noncommutative part, any other is combined again with the factor before
    ne <- .add(list(b1[[2L]], b2[[2L]]))
    if (identical(b1[[1L]]$key, b2[[1L]]$key) && ne$type != "Add") {
      out[[length(out)]] <- NULL
      p <- .pow(b1[[1L]], ne)
      if (.is_comm(p)) cpart[[length(cpart) + 1L]] <- p else todo <- c(list(p), todo)
    } else {
      out[[length(out) + 1L]] <- f
    }
  }
  list(nc = out, c = cpart)
}

# ---------------------------------------------------------------------------
# Pow

.pow <- function(b, e) {
  if (is_sym(b) && b$type == "Order") { r <- .order_pow(b, as_sym(e)); if (!is.null(r)) return(r) }
  if (!is.null(.arith_hooks$pow_pre)) { r <- .arith_hooks$pow_pre(b, e); if (!is.null(r)) return(r) }
  bt <- b$type
  et <- e$type
  if (et == "Integer" || et == "Rational" || et == "Float") {
    if (num_is_zero(e)) return(S_One)
    if (et == "Integer" && e$value == 1) return(b)
  }
  if (is_nan(b) || is_nan(e)) return(S_NaN)
  # Dimension._eval_power (physics-units-dimensions.R)
  if (bt == "Dimension") return(.dim_eval_power(b, e))
  if (bt == "AccumulationBounds") return(.ab_pow(b, e))
  if (bt == "Integer" && b$value == -1 && !(et %in% c("Integer", "Rational", "Float"))) {
    # Pow.__new__: (-1)^e is 1 for even e and -1 for odd e
    if (isTRUE(ask_fact(e, "even"))) return(S_One)
    if (isTRUE(ask_fact(e, "odd"))) return(S_NegativeOne)
  }
  if (et == "AccumulationBounds" && !(bt == "Constant" && b$name == "E")) {
    r <- .ab_rpow(b, e)
    return(if (is.null(r)) new_sym("Pow", list(b, e)) else r)
  }
  # Pow.__new__: a negative number to an integer symbol sheds its sign
  # for even and odd exponents
  if (et %in% c("Symbol", "Dummy") && (bt %in% c("Integer", "Rational", "Float") || (bt == "Mul" && !length(free_symbols(b)))) &&
      isTRUE(ask_fact(e, "integer")) && could_extract_minus_sign(b)) {
    if (isTRUE(ask_fact(e, "even"))) { b <- .neg(b); bt <- b$type }
    else if (isTRUE(ask_fact(e, "odd"))) return(.neg(.pow(.neg(b), e)))
  }
  # Integer._eval_power and Rational._eval_power: (-2)^k is 2^k for even k
  if (bt %in% c("Integer", "Rational") && !(et %in% c("Integer", "Rational", "Float")) &&
      num_sign(b) < 0 && isTRUE(ask_fact(e, "even"))) { b <- .neg(b); bt <- b$type }
  if (bt == "Constant" && b$name == "E") return(fcall("exp", list(e)))
  if (bt == "Integer" && b$value == 1) {
    if (is_inf_const(e)) return(S_NaN)
    return(S_One)
  }
  # Pow.__new__: recognize base as E, b^(c*f/log(b)) -> E^(c*f)
  if (length(e$args) && !(bt == "Function" && b$name == "exp_polar") && .has_fun(e, "log")) {
    cm <- as_coeff_mul_split(factor_terms(e, sign = FALSE))
    nd <- fraction(cm[[2L]])
    den <- nd[[2L]]
    if (den$type == "Function" && den$name == "log" && length(den$args) == 1L && identical(den$args[[1L]]$key, b$key))
      return(fcall("exp", list(.mul(list(cm[[1L]], nd[[1L]])))))
    if (den$type == "Add") {
      s <- fcall("sign", list(fcall("im", list(b))))
      if (is_number_type(s) && !is_zero_num(s) &&
          identical(den$key, .add(list(fcall("log", list(.neg(factor_terms(b, sign = FALSE)))), .mul(list(s, S_I, S_Pi))))$key))
        return(fcall("exp", list(.mul(list(cm[[1L]], nd[[1L]])))))
    }
  }
  b_num <- bt == "Integer" || bt == "Rational" || bt == "Float"
  e_num <- et == "Integer" || et == "Rational" || et == "Float"
  if (b_num && e_num) return(.pow_num(b, e))
  if (is_inf_const(e)) {
    r <- .pow_to_inf(b, e)
    if (!is.null(r)) return(r)
  }
  if (is_inf_const(b)) {
    r <- .pow_inf_base(b, e)
    if (!is.null(r)) return(r)
  }
  if (bt == "Integer" && b$value == 0) {
    if (isTRUE(ask_fact(e, "positive"))) return(S_Zero)
    if (isTRUE(ask_fact(e, "negative"))) return(S_ComplexInfinity)
  }
  if (bt == "Constant" && b$name == "I" && et == "Integer") {
    k <- .int_mod(e, 4L)
    return(switch(k + 1L, S_One, S_I, S_NegativeOne,
                  new_sym("Mul", list(S_NegativeOne, S_I))))
  }
  # ImaginaryUnit._eval_power for a rational exponent: i, r = divmod(e, 2),
  # I**r unevaluated, negated for odd i
  if (bt == "Constant" && b$name == "I" && et == "Rational") {
    q <- gmp::as.bigq(e$value)
    i <- gmp::as.bigz(floor(q / 2))
    r <- q - 2 * gmp::as.bigq(i)
    rv <- new_sym("Pow", list(S_I, make_rational(r)))
    if (i %% 2 != 0) return(new_sym("Mul", list(S_NegativeOne, rv)))
    return(rv)
  }
  # KroneckerDelta._eval_power
  if (bt == "Function" && b$name == "KroneckerDelta") {
    if (isTRUE(ask_fact(e, "positive"))) return(b)
    if (isTRUE(ask_fact(e, "negative")) && !is_negone_num(e)) return(.inv(b))
  }
  # Piecewise._eval_power: the power goes into every branch
  if (bt == "Piecewise") return(do.call(Piecewise, lapply(b$args, function(pc) list(.pow(pc$args[[1L]], e), pc$args[[2L]]))))
  if (bt == "Add" && e_num) {
    r <- .add_eval_power(b, e)
    if (!is.null(r)) return(r)
  }
  if (bt == "Pow") {
    r <- .pow_eval_power(b$args[[1L]], b$args[[2L]], e)
    if (!is.null(r)) return(r)
  }
  # a polar base: (exp_polar(a))^e = exp_polar(a*e)
  if (bt == "Function" && b$name == "exp_polar") return(fcall("exp_polar", list(.mul(list(b$args[[1L]], e)))))
  if (bt == "Function" && b$name == "exp") {
    # ExpBase._eval_power through Pow._eval_power of E^a: an integer
    # exponent, a real a, or a nonreal a whose imaginary part keeps
    # floor(1/2 - im(a)/(2*pi)) = 0
    a <- b$args[[1L]]
    ar <- ask_fact(a, "extended_real")
    if (et == "Integer" || isTRUE(ask_fact(e, "integer")) || isTRUE(ar)) {
      return(fcall("exp", list(.mul(list(a, e)))))
    }
    if (isFALSE(ar)) {
      ia <- fcall("im", list(a))
      if (!length(free_symbols(ia))) {
        v <- tryCatch(Re(.neval(ia, list())), error = function(err) NA_real_)
        if (!is.na(v) && floor(0.5 - v / (2 * pi)) == 0) return(fcall("exp", list(.mul(list(a, e)))))
      }
    }
  }
  if (bt == "Mul") {
    fs <- b$args
    comm <- vapply(fs, .is_comm, TRUE)
    if (et == "Integer") {
      if (all(comm)) return(.mul(lapply(fs, function(f) .pow(f, e))))
      if (any(comm)) {
        return(.mul(c(lapply(fs[comm], function(f) .pow(f, e)),
                      list(new_sym("Pow", list(.mul(fs[!comm]), e))))))
      }
    } else if (all(comm)) {
      # Mul._eval_power: sqrt(2*I) = 1 + I, then for a rational or float
      # exponent Pow._eval_expand_power_base; other exponents keep (a*b)^e
      if (et == "Rational" && gmp::denominator(e$value) == 2) {
        r <- .mul_pow_imag_half(b, e)
        if (!is.null(r)) return(r)
      }
      if (et %in% c("Rational", "Float")) return(.expand_power_base_rat(fs, e))
    }
  }
  new_sym("Pow", list(b, e))
}

# Pow._eval_power: (b^e)^expt = s*b^(e*expt) when the branch factor s is known
.pow_eval_power <- function(b, e, expt) {
  if (is_nan(b)) return(NULL)
  num_lt <- function(x, k) {
    if (length(free_symbols(x))) return(FALSE)
    d <- tryCatch(.mp_num_sign(.sub(x, .int(k))), error = function(err) NA)
    isTRUE(d < 0)
  }
  half <- function(x) {
    if (x$type == "Rational" && gmp::denominator(x$value) == 2) return(TRUE)
    nd <- .as_numer_denom(x)
    isTRUE(ask_fact(nd[[1L]], "integer")) && identical(nd[[2L]]$key, S_Two$key)
  }
  # s when it is real and numerically sign(s)
  real_sign <- function(s) {
    if (!isTRUE(ask_fact(s, "extended_real"))) return(NULL)
    sg <- fcall("sign", list(s))
    d <- tryCatch(.neval(.sub(sg, s), list()), error = function(err) NA)
    if (!is.na(d) && Mod(d) < 1e-2) sg else NULL
  }
  s <- NULL
  if (isTRUE(ask_fact(expt, "integer"))) s <- S_One
  else if (.is_polar(b)) s <- S_One
  else {
    er <- ask_fact(e, "extended_real")
    if (isTRUE(er)) {
      if (is_negone_num(e)) {
        if (half(expt)) {
          if (isTRUE(ask_fact(b, "negative"))) return(.mul(list(.pow(S_NegativeOne, expt), .pow(.neg(b), .mul(list(e, expt))))))
          if (isFALSE(ask_fact(b, "negative"))) return(.pow(b, .neg(expt)))
        }
      } else if (isTRUE(ask_fact(e, "even"))) {
        if (isTRUE(ask_fact(b, "extended_real"))) b <- fcall("Abs", list(b))
        if (isTRUE(ask_fact(b, "imaginary"))) b <- .mul(list(fcall("Abs", list(fcall("im", list(b)))), S_I))
      }
      if (num_lt(fcall("Abs", list(e)), 1L) || is_one_num(e)) s <- S_One
      else if (isTRUE(ask_fact(b, "extended_nonnegative"))) s <- S_One
      else if (isTRUE(ask_fact(fcall("re", list(b)), "extended_nonnegative")) && num_lt(fcall("Abs", list(e)), 2L)) s <- S_One
      else if (half(expt)) {
        fl <- fcall("floor", list(.sub(S_Half, .div(.mul(list(e, fcall("arg", list(b)))), .mul(list(S_Two, S_Pi))))))
        s <- real_sign(fcall("exp", list(.mul(list(S_Two, S_Pi, S_I, expt, fl)))))
      }
    } else if (isFALSE(er)) {
      fl <- fcall("floor", list(.sub(S_Half, .div(fcall("im", list(.mul(list(e, fcall("log", list(b)))))), .mul(list(S_Two, S_Pi))))))
      s <- real_sign(fcall("exp", list(.mul(list(S_Two, S_I, S_Pi, expt, fl)))))
    }
  }
  if (is.null(s)) return(NULL)
  .mul(list(s, .pow(b, .mul(list(e, expt)))))
}

# Mul._eval_power for an imaginary base and a half-integer exponent:
# (c*I)^(p/2) with |c/2| a rational square.
.mul_pow_imag_half <- function(b, e) {
  if (!isTRUE(ask_fact(b, "imaginary"))) return(NULL)
  a <- .imag_coeff(b)
  if (is.null(a) || !(a$type %in% c("Integer", "Rational"))) return(NULL)
  h <- abs(gmp::as.bigq(a$value) / 2)
  n <- .iroot(gmp::numerator(h), 2); d <- .iroot(gmp::denominator(h), 2)
  if (!n$exact || !d$exact) return(NULL)
  p <- .int(as.character(gmp::numerator(e$value)))
  r <- make_rational(gmp::as.bigq(n$root, d$root))
  sg <- if (a$value > 0) S_I else .neg(S_I)
  f1 <- .pow(r, p); f2 <- .pow(.add(list(S_One, sg)), p)
  # _unevaluated_Mul: 2*(1 + I) is not distributed
  if (is_one_num(f1)) f2 else if (f2$type == "Add") new_sym("Mul", list(f1, f2)) else .mul(list(f1, f2))
}

# Pow._eval_expand_power_base (force = FALSE) for commutative factors and a
# rational or float exponent that is not an integer.
.expand_power_base_rat <- function(fs, e) {
  other <- list(); nonneg <- list(); neg <- list()
  for (f in fs) {
    # polar factors split out of powers as nonnegative ones do
    if (.is_polar(f)) { nonneg[[length(nonneg) + 1L]] <- f; next }
    if (isFALSE(ask_fact(f, "real"))) { other[[length(other) + 1L]] <- f; next }
    nn <- ask_fact(f, "nonnegative")
    if (isTRUE(nn)) nonneg[[length(nonneg) + 1L]] <- f
    else if (isFALSE(nn)) neg[[length(neg) + 1L]] <- f
    else other[[length(other) + 1L]] <- f
  }
  if (length(neg) > 1L) {
    o <- S_One
    if (!length(other) && is_number_type(neg[[1L]])) { o <- neg[[1L]]; neg <- neg[-1L] }
    if (length(neg) %% 2L == 1L) o <- .neg(o)
    for (n in neg) nonneg[[length(nonneg) + 1L]] <- .neg(n)
    if (!is_one_num(o)) other[[length(other) + 1L]] <- o
  } else if (length(neg) && length(other)) {
    if (is_number_type(neg[[1L]]) && !is_negone_num(neg[[1L]])) {
      other[[length(other) + 1L]] <- S_NegativeOne
      nonneg[[length(nonneg) + 1L]] <- .neg(neg[[1L]])
    } else other <- c(other, neg)
  } else other <- c(other, neg)
  parts <- lapply(nonneg, function(f) .pow(f, e))
  if (length(other)) {
    om <- .mul(other)
    parts[[length(parts) + 1L]] <- if (is_number_type(om)) .pow(om, e) else new_sym("Pow", list(om, e))
  }
  .mul(parts)
}

# Remainder of a sym Integer modulo a small positive integer.
.int_mod <- function(e, m) {
  as.integer(as.character(gmp::mod.bigz(gmp::numerator(e$value), m)))
}

.pow_to_inf <- function(b, e) {
  if (!is_number_type(b)) return(NULL)
  v <- num_value_dbl(b)
  # exact bases near 1 must not round to 1
  if (b$type != "Float" && v %in% c(-1, 1)) {
    q <- b$value
    v <- if (q == 1) 1 else if (q == -1) -1 else if (q > 1) 2 else if (q < -1) -2 else if (q > 0) 0.5 else -0.5
  }
  if (e$name == "oo") {
    if (v > 1) return(S_Infinity)
    if (v > -1 && v < 1) return(S_Zero)
    if (v == 1 || v == -1) return(S_NaN)
    return(S_ComplexInfinity)
  }
  if (e$name == "-oo") {
    if (v > 1 || v < -1) return(S_Zero)
    if (v > 0 && v < 1) return(S_Infinity)
    if (v == 0) return(S_ComplexInfinity)
    return(S_NaN)
  }
  S_NaN
}

.pow_inf_base <- function(b, e) {
  if (!is_number_type(e)) return(NULL)
  s <- num_sign(e)
  if (s < 0) return(S_Zero)
  if (b$name == "oo") return(S_Infinity)
  if (b$name == "zoo") return(S_ComplexInfinity)
  if (e$type == "Integer") {
    return(if (.int_mod(e, 2L) == 0L) S_Infinity else S_NegativeInfinity)
  }
  NULL
}

.pow_num <- function(b, e) {
  if (b$type == "Float" || e$type == "Float") {
    # Float powers at the Floats' precision, as mpf_pow and mpc_pow_mpf
    p <- .floats_prec(list(b, e))
    bzero <- if (b$type == "Float") .float_mp(b)$m == 0 else is_zero_num(b)
    if (bzero && num_sign(e) < 0) return(S_ComplexInfinity)
    r <- tryCatch(.mpeval(new_sym("Pow", list(b, e)), p + 10, list()), error = function(err) NULL)
    if (!is.null(r)) return(.mpc_sym(.mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p)), p))
    bd <- num_value_dbl(b)
    ed <- num_value_dbl(e)
    if (bd == 0 && ed < 0) return(S_ComplexInfinity)
    if (bd >= 0 || ed == round(ed)) return(Float(bd^ed))
    r <- abs(bd)^ed
    th <- pi * ed
    return(.add(list(Float(r * cos(th)), .mul(list(Float(r * sin(th)), S_I)))))
  }
  if (e$type == "Integer") {
    if (b$value == 0) {
      if (e$value < 0) return(S_ComplexInfinity)
      return(S_Zero)
    }
    if (b$value == 1) return(S_One)
    if (b$value == -1) return(.pow_negone(e$value))
    if (abs(e$value) > 1e7) stop("exponent too large for exact power")
    return(make_rational(.bigq_pow_int(b$value, as.numeric(e$value))))
  }
  if (b$value == 0) {
    if (e$value > 0) return(S_Zero)
    return(S_ComplexInfinity)
  }
  if (b$value < 0) {
    # Integer and Rational _eval_power for a negative base
    ev <- e$value
    if (b$value == -1) return(.pow_negone(ev))
    if (ev < 0) return(.mul(list(.pow_negone(ev), .pow_num(make_rational(1 / (-b$value)), make_rational(-ev)))))
    if (gmp::denominator(b$value) != 1)
      return(.mul(list(.pow_num(make_rational(gmp::as.bigq(gmp::numerator(b$value))), e),
                       .pow_num(make_rational(1 / gmp::as.bigq(gmp::denominator(b$value))), e))))
    if (ev == gmp::as.bigq(1, 2)) return(.mul(list(S_I, .pow_num(make_rational(-b$value), e))))
    n <- gmp::as.bigz(gmp::numerator(-b$value))
    a <- .int_pow_rat(n, ev)
    # nothing extracted: (-2)^(1/3) stays whole
    if (is.null(a)) return(new_sym("Pow", list(b, e)))
    return(.mul(list(.pow_negone(ev), a)))
  }
  .root_rational(b$value, e$value)
}

.pow_negone <- function(r) {
  # NegativeOne._eval_power: I^p for halves, else (-1)^i * (-1)^(r - i), i = floor(r)
  r <- gmp::as.bigq(r)
  p <- gmp::numerator(r); q <- gmp::denominator(r)
  if (q == 1) return(if (gmp::mod.bigz(p, 2) == 0) S_One else S_NegativeOne)
  if (q == 2) return(.pow(S_I, .int(as.character(p))))
  i <- .bigq_floor(r)
  base <- new_sym("Pow", list(S_NegativeOne, make_rational(r - gmp::as.bigq(i))))
  if (gmp::mod.bigz(i, 2) != 0) .mul(list(S_NegativeOne, base)) else base
}

# Factor an integer n >= 1 into list(primes, exponents); primes as bigz.
.factor_small <- function(n) {
  n <- gmp::as.bigz(n)
  if (n == 1) return(list(p = gmp::as.bigz(integer()), e = integer()))
  if (n < gmp::as.bigz("1000000000000000000000000000000")) {
    f <- gmp::factorize(n)
    tab <- table(as.character(f))
    ps <- gmp::as.bigz(names(tab))
    o <- order(as.double(ps))
    return(list(p = ps[o], e = as.integer(tab)[o], rest = gmp::as.bigz(1)))
  }
  # ponytail: very large integers only get trial division by small primes;
  # the cofactor is kept whole. Use full factorization if needed.
  ps <- c()
  es <- c()
  rest <- n
  for (p in c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)) {
    k <- 0L
    while (gmp::mod.bigz(rest, p) == 0) {
      rest <- gmp::divq.bigz(rest, p)
      k <- k + 1L
    }
    if (k) {
      ps <- c(ps, p)
      es <- c(es, k)
    }
  }
  list(p = gmp::as.bigz(ps), e = as.integer(es), rest = rest)
}

# n^(p/q) for positive integer n, as list(coef = bigq, base = bigz,
# frac = bigq in [0, 1)) with n^(p/q) = coef * base^frac.
# Floor of the k-th root of a nonnegative bigz, and whether it is exact.
.iroot <- function(n, k) {
  n <- gmp::as.bigz(n)
  if (n < 2) return(list(root = n, exact = TRUE))
  bits <- nchar(as.character(n, b = 2))
  x <- gmp::as.bigz(2)^(ceiling(bits / k))
  repeat {
    y <- gmp::divq.bigz((k - 1) * x + gmp::divq.bigz(n, x^(k - 1)), k)
    if (y >= x) break
    x <- y
  }
  list(root = x, exact = x^k == n)
}

# n^(p/q) for an integer n >= 2 and rational p/q > 0 with q > 1, after
# Integer._eval_power: perfect roots, then the factors of n up to 2^15
# (or the base of a perfect power) give integer parts and radicals.
# NULL when nothing simplifies.
.int_pow_rat <- function(n, ex) {
  n <- gmp::as.bigz(n)
  p <- gmp::numerator(ex); q <- gmp::denominator(ex)
  big <- function(b, k) if (as.double(k) * .bitlen(gmp::as.bigz(b)) > 1e7) stop("exponent too large for exact power") else b^k
  # an exact q-th root needs q <= log2(n)
  if (q <= .bitlen(n)) {
    rr <- .iroot(n, as.integer(as.character(q)))
    if (rr$exact) return(.int(big(rr$root, p)))
  }
  pp <- .perfect_power(n)
  fd <- if (!isFALSE(pp)) stats::setNames(as.numeric(as.character(pp[[2L]])), as.character(pp[[1L]]))
        else .factorint(n, 2^15)
  out_int <- gmp::as.bigz(1); out_rad <- list(); sqr <- list()
  for (pr in names(fd)) {
    pz <- gmp::as.bigz(pr)
    E <- gmp::as.bigz(fd[[pr]]) * p
    de <- gmp::divq.bigz(E, q); dm <- gmp::mod.bigz(E, q)
    if (de > 0) out_int <- out_int * big(pz, de)
    if (dm > 0) {
      g <- gmp::gcd.bigz(dm, q)
      if (g != 1) out_rad[[length(out_rad) + 1L]] <- new_sym("Pow", list(.int(pz), make_rational(gmp::as.bigq(dm / g, q / g))))
      else sqr[[pr]] <- dm
    }
  }
  sg <- gmp::as.bigz(0)
  for (ex2 in sqr) { sg <- if (sg == 0) ex2 else gmp::gcd.bigz(sg, ex2); if (sg == 1) break }
  sqr_int <- gmp::as.bigz(1)
  for (pr in names(sqr)) sqr_int <- sqr_int * big(gmp::as.bigz(pr), gmp::divq.bigz(sqr[[pr]], sg))
  if (sqr_int == n && out_int == 1 && !length(out_rad)) return(NULL)
  last <- if (sqr_int == 1 || sg == 0) S_One else .pow(.int(sqr_int), make_rational(gmp::as.bigq(sg, q)))
  .mul(c(list(.int(out_int)), out_rad, list(last)))
}

# A positive rational to a positive non-integer rational power, after
# Rational._eval_power and Integer._eval_power.
.root_rational <- function(qv, ex) {
  # (3/4)^(-r) = (4/3)^r
  if (ex < 0) return(.root_rational(1 / gmp::as.bigq(qv), -ex))
  num <- gmp::numerator(qv); den <- gmp::denominator(qv)
  ipow <- function(n, r) {
    if (n == 1) return(S_One)
    v <- .int_pow_rat(n, r)
    if (is.null(v)) new_sym("Pow", list(.int(n), make_rational(r))) else v
  }
  if (den == 1) return(ipow(num, ex))
  p <- gmp::numerator(ex); q <- gmp::denominator(ex)
  ip <- gmp::divq.bigz(p, q)
  if (ip > 0) { ip <- ip + 1; rem <- ip * q - p } else { ip <- gmp::as.bigz(1); rem <- q - p }
  .mul(list(ipow(num, ex), ipow(den, gmp::as.bigq(rem, q)), make_rational(gmp::as.bigq(1, den^ip))))
}

# Add._eval_power: a + b*I with an infinite part, (r + i*I)^(p/2) for a
# Pythagorean triple, and (r + i*I)^-1
.add_eval_power <- function(b, e) {
  if (length(b$args) == 2L && any(vapply(b$args, function(a) isTRUE(ask_fact(a, "infinite")), TRUE))) {
    if (isFALSE(ask_fact(e, "zero")) && !(is_number_type(e) && is_one_num(e))) {
      a1 <- b$args[[1L]]; b1 <- b$args[[2L]]
      if (!is_zero_num(coeff(a1, S_I))) { t <- a1; a1 <- b1; b1 <- t }
      ico <- coeff(b1, S_I)
      if (!is_zero_num(ico) && isTRUE(ask_fact(ico, "extended_real")) && isTRUE(ask_fact(a1, "extended_real"))) {
        if (isTRUE(ask_fact(e, "extended_negative"))) return(S_Zero)
        if (isTRUE(ask_fact(e, "extended_positive"))) return(S_ComplexInfinity)
      }
    }
    return(NULL)
  }
  et <- e$type
  if (et %in% c("Integer", "Rational") && !length(free_symbols(b))) {
    h <- if (is_number_type(b$args[[1L]])) b$args[[1L]] else NULL
    if (!is.null(h) && length(b$args) == 2L) {
      sp <- as_coeff_mul_split(b$args[[2L]])
      if (is_const(sp[[2L]], "I") && is_number_type(sp[[1L]])) {
        r <- h; i <- sp[[1L]]
        q <- gmp::denominator(e$value); p <- gmp::numerator(e$value)
        if (q == 2) {
          D <- .sqrt(.add(list(.pow(r, S_Two), .pow(i, S_Two))))
          if (.is_rat(D)) {
            root <- .pow(.sqrt(factor_terms(.div(.sub(D, r), S_Two))), .int(p))
            si <- if (num_sign(i) > 0) S_One else S_NegativeOne
            return(.mul(list(root, expand_multinomial(.pow(.add(list(.div(.add(list(D, r)), .abs_num(i)), .mul(list(si, S_I)))), .int(p))))))
          }
        } else if (q == 1 && p == -1) {
          return(.umul(.sub(r, .mul(list(i, S_I))), .inv(.add(list(.pow(r, S_Two), .pow(i, S_Two))))))
        }
      }
    }
  }
  NULL
}
