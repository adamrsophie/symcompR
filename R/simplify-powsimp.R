# powsimp and powdenest after sympy/simplify/powsimp.py.

#' Combine powers
#'
#' `powsimp()` combines powers with the same base (`x^a*x^b = x^(a+b)`) and,
#' when valid or with `force = TRUE`, powers with the same exponent
#' (`x^a*y^a = (x*y)^a`), as SymPy's `powsimp()`. `powdenest()` rewrites
#' `(x^a)^b` as `x^(a*b)` when valid. `logcombine()` merges sums of
#' logarithms.
#'
#' @param e An expression.
#' @param deep Also simplify inside function arguments, bases and
#'   exponents.
#' @param combine `"all"`, `"exp"` or `"base"`.
#' @param force Ignore assumptions.
#' @param measure The complexity measure (unused by the combination itself,
#'   it orders the radical bases).
#' @param polar Treat the expression as one on the Riemann surface of the
#'   logarithm (`powdenest()`).
#' @return An expression.
#' @examples
#' x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b")
#' powsimp(x^a * x^b)
#' logcombine(log(x) + log(Symbol("y")), force = TRUE)
#' @export
powsimp <- function(e, deep = FALSE, combine = "all", force = FALSE, measure = count_ops) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) powsimp(z, deep, combine, force, measure)))
  .powsimp(as_sym(e), deep, combine, force, measure)
}

.ps_y <- Dummy("y")

# Expr.as_base_exp
.abe <- function(z) {
  if (z$type == "Function" && z$name == "exp") return(list(S_E, z$args[[1L]]))
  if (z$type == "Function" && z$name == "exp_polar" && !is_zero_num(z$args[[1L]])) return(list(fcall("exp_polar", list(S_One)), z$args[[1L]]))
  # Pow.as_base_exp is its args (the reference SymPy no longer turns
  # (1/q)^e into q^(-e))
  if (z$type == "Pow") return(list(z$args[[1L]], z$args[[2L]]))
  list(z, S_One)
}

.pw <- function(b, e) if (is_const(b, "E")) fcall("exp", list(e)) else .pow(b, e)

.fa <- function(z, nm) isTRUE(ask_fact(z, nm))

.powsimp <- function(expr, deep, combine, force, measure) {
  recurse <- function(arg, deep_ = deep, combine_ = combine) .powsimp(arg, deep_, combine_, force, measure)
  if (!is_sym(expr) || !length(expr$args) || expr$type == "MatrixSymbol") return(expr)
  if (expr$type == "Function" && expr$name == "exp_polar" && (is_zero_num(expr$args[[1L]]) || is_one_num(expr$args[[1L]]))) return(expr)
  y_in <- expr$type == "Mul" && any(vapply(expr$args, function(a) identical(a$key, .ps_y$key), TRUE))
  if (deep || expr$type == "Add" || (expr$type == "Mul" && !y_in)) {
    args <- lapply(expr$args, recurse)
    if (!identical(vapply(args, .k, ""), vapply(expr$args, .k, ""))) expr <- sym_rebuild(expr, args)
  }
  if (expr$type == "Pow") return(.div(recurse(.mul(list(expr, .ps_y)), deep_ = FALSE), .ps_y))
  if (expr$type != "Mul") return(expr)
  if (combine %in% c("exp", "all")) {
    cb <- list(); ce <- list()
    nc_part <- list(); coeff <- S_One
    for (term in .sympy_args(expr)) {
      if (term$type %in% c("Integer", "Rational")) { coeff <- .mul(list(coeff, term)); next }
      if (term$type == "Pow") term <- .denest_pow(term)
      if (.is_comm(term)) {
        be <- .abe(term); b <- be[[1L]]; e <- be[[2L]]
        if (deep) { b <- recurse(b); e <- recurse(e) }
        if (b$type == "Pow" || (b$type == "Function" && b$name == "exp")) { b <- .pw(b, e); e <- S_One }
        k <- b$key
        if (is.null(cb[[k]])) { cb[[k]] <- b; ce[[k]] <- list(e) } else ce[[k]][[length(ce[[k]]) + 1L]] <- e
      } else {
        if (length(nc_part)) {
          b1 <- .abe(nc_part[[length(nc_part)]]); b2 <- .abe(term)
          if (identical(b1[[1L]]$key, b2[[1L]]$key) && .is_comm(b1[[2L]]) && .is_comm(b2[[2L]])) {
            nc_part[[length(nc_part)]] <- .pow(b1[[1L]], .add(list(b1[[2L]], b2[[2L]])))
            next
          }
        }
        nc_part[[length(nc_part) + 1L]] <- term
      }
    }
    # a rational base absorbs its powers from the coefficient
    for (b in ordered(unname(cb))) {
      k <- b$key; e <- ce[[k]]
      if (!is_zero_num(b) && b$type %in% c("Integer", "Rational") && !all(vapply(e, is_number_type, TRUE)) &&
          !is_one_num(coeff) && !(is_one_num(b) || identical(b$key, S_NegativeOne$key))) {
        m <- multiplicity(.abs_num(b), .abs_num(coeff))
        if (m != 0) {
          e[[length(e) + 1L]] <- .int(m)
          coeff <- .div(coeff, .pow(b, .int(m)))
        }
      }
      ce[[k]] <- .add(e)
    }
    if (!is_one_num(coeff)) {
      k <- coeff$key
      if (!is.null(cb[[k]])) ce[[k]] <- .add(list(ce[[k]], S_One)) else { cb[[k]] <- coeff; ce[[k]] <- S_One }
    }
    # b^e with 1/b^e
    skip <- character()
    ce0 <- ce
    for (k in names(cb)) {
      if (k %in% skip || is.null(cb[[k]])) next
      b <- cb[[k]]; e <- ce0[[k]]
      if (.fa(b, "positive") || .fa(b, "polar")) {
        binv <- .inv(b)
        if (b$type == "Float" && .float_is_one(b)) { ce[[k]] <- S_One; next }
        ki <- binv$key
        if (!identical(ki, k) && !is.null(cb[[ki]])) {
          if (is_one_num(.as_numer_denom(b)[[1L]])) {
            cb[[k]] <- NULL; ce[[ki]] <- .sub(ce[[ki]], e); ce[[k]] <- NULL
          } else {
            skip <- c(skip, ki)
            e2 <- ce[[ki]]; cb[[ki]] <- NULL; ce[[ki]] <- NULL
            ce[[k]] <- .sub(ce[[k]], e2)
          }
        }
      }
    }
    # b^e with (-b)^e
    ce0 <- ce
    for (k in names(cb)) {
      b <- cb[[k]]
      if (is.null(b) || !(b$type %in% c("Symbol", "Add"))) next
      nb <- .neg(b); kn <- nb$key
      if (is.null(cb[[kn]])) next
      e <- ce0[[k]]
      if (!is.na(.tv(ask_fact(b, "positive"))) || .fa(e, "integer")) {
        if (.fa(e, "integer") || .fa(b, "negative")) {
          ce[[kn]] <- .add(list(ce[[kn]], e)); cb[[k]] <- NULL; ce[[k]] <- NULL
        } else {
          e <- ce[[kn]]; cb[[kn]] <- NULL; ce[[kn]] <- NULL
          ce[[k]] <- .add(list(ce[[k]], e))
        }
        kk <- S_NegativeOne$key
        if (!is.null(cb[[kk]])) ce[[kk]] <- .add(list(ce[[kk]], e)) else { cb[[kk]] <- S_NegativeOne; ce[[kk]] <- e }
      }
    }
    cp <- list()
    for (k in names(cb)) if (!is_zero_num(ce[[k]])) cp[[length(cp) + 1L]] <- list(cb[[k]], ce[[k]])
    cp <- .powsimp_radicals(cp, measure)
    newexpr <- .mul(lapply(cp, function(p) .pw(p[[1L]], p[[2L]])))
    if (combine == "exp") return(.mul(list(newexpr, .mul(nc_part))))
    return(.mul(list(recurse(.mul(nc_part), combine_ = "base"), recurse(newexpr, combine_ = "base"))))
  }
  if (combine != "base") stop("combine must be one of 'all', 'exp', 'base'")
  .powsimp_base(expr, deep, force, recurse)
}

.tv <- function(x) if (is.null(x) || is.na(x)) NA else x
.abs_num <- function(z) if (num_sign(z) < 0) .neg(z) else z
.float_is_one <- function(b) isTRUE(tryCatch(.sort_num(b) == 1, error = function(e) FALSE))

# bkey of powsimp: the base with the denominator of a rational exponent,
# and the numerator
.ps_bkey <- function(b, e = NULL) {
  if (is.null(e)) { be <- .abe(b); return(.ps_bkey(be[[1L]], be[[2L]])) }
  if (e$type == "Integer") return(list(list(b, S_One), e))
  if (e$type == "Rational") return(list(list(b, .int(gmp::denominator(e$value))), .int(gmp::numerator(e$value))))
  sp <- as_coeff_mul_split(e); c <- sp[[1L]]; m <- sp[[2L]]
  if (c$type %in% c("Integer", "Rational") && !is_one_num(c)) {
    cq <- .int(gmp::denominator(c$value)); cp <- .int(gmp::numerator(c$value))
    if (.fa(m, "integer")) return(list(list(b, cq), .mul(list(m, cp))))
    return(list(list(.pw(b, m), cq), cp))
  }
  list(list(.pw(b, e), S_One), S_One)
}
.bk <- function(bq) paste0(bq[[1L]]$key, "|", bq[[2L]]$key)
.bigint <- function(z) gmp::as.bigz(z$value)

# radicals with the same denominator joined: sqrt(x)*sqrt(y) = sqrt(x*y)
.powsimp_radicals <- function(cp, measure) {
  cbq <- list(); cbe <- list(); bases <- list()
  for (p in cp) {
    r <- .ps_bkey(p[[1L]], p[[2L]]); bq <- r[[1L]]; k <- .bk(bq)
    if (!is.null(cbq[[k]])) cbe[[k]] <- .add(list(cbe[[k]], r[[2L]])) else { cbq[[k]] <- bq; cbe[[k]] <- r[[2L]] }
    if (!is_one_num(bq[[2L]]) && bq[[1L]]$type == "Mul") bases[[length(bases) + 1L]] <- bq
  }
  update <- function(bq) {
    k <- .bk(bq); v <- cbe[[k]]
    if (v$type != "Integer") return(invisible())
    q <- .bigint(bq[[2L]]); vv <- .bigint(v)
    newe <- vv %/% q; r <- vv %% q
    if (r != 0) return(invisible())
    cbq[[k]] <<- NULL; cbe[[k]] <<- NULL
    if (newe != 0) {
      pm <- .pw(bq[[1L]], .int(newe))
      for (m in (if (pm$type == "Mul") pm$args else list(pm))) {
        rr <- .ps_bkey(m); b2 <- rr[[1L]]; k2 <- .bk(b2)
        if (is.null(cbq[[k2]])) { cbq[[k2]] <<- b2; cbe[[k2]] <<- S_Zero }
        cbe[[k2]] <<- .add(list(cbe[[k2]], rr[[2L]]))
        if (!is_one_num(b2[[2L]])) bases[[length(bases) + 1L]] <<- b2
      }
    }
    invisible()
  }
  if (length(bases)) {
    bases <- bases[!duplicated(vapply(bases, .bk, ""))]
    bases <- .key_sort(lapply(bases, function(bq) list(.sort_key(bq[[1L]]), .sort_key(bq[[2L]]))), bases)
    ms <- vapply(bases, function(bq) as.numeric(measure(bq[[1L]])) + as.numeric(measure(bq[[2L]])), 0)
    bases <- bases[order(-ms)]
  }
  done <- list()
  ratq <- function(x) .ps_bkey(x)[[1L]][[2L]]
  lt <- function(a, b) { d <- .sub(a, b); if (!is_number_type(d)) NA else num_sign(d) < 0 }
  for (base in bases) {
    kb <- .bk(base)
    if (is.null(cbq[[kb]])) next
    b <- base[[1L]]; exponent <- base[[2L]]
    last <- FALSE; qlcm <- gmp::as.bigz(1)
    repeat {
      bstart <- b; qstart <- qlcm
      bb <- list(); ee <- list()
      for (bi in (if (b$type == "Mul") b$args else list(b))) {
        rr <- .ps_bkey(bi); bib <- rr[[1L]]; bie <- rr[[2L]]; kk <- .bk(bib)
        if (is.null(cbq[[kk]]) || !isFALSE(lt(cbe[[kk]], bie))) { ee <- list(); bb <- list(); break }
        ee[[length(ee) + 1L]] <- list(bie, cbe[[kk]]); bb[[length(bb) + 1L]] <- bib
      }
      if (length(ee) && all(vapply(ee, function(p) p[[1L]]$type == "Integer" && p[[2L]]$type == "Integer", TRUE))) {
        fd <- function(p) .bigint(p[[2L]]) %/% .bigint(p[[1L]])
        min1 <- fd(ee[[1L]]); ok <- TRUE
        for (i in seq_along(ee)[-1L]) { rat <- fd(ee[[i]]); if (rat < 1) { ok <- FALSE; break }; min1 <- min(min1, rat) }
        if (ok) {
          for (i in seq_along(bb)) {
            kk <- .bk(bb[[i]])
            cbe[[kk]] <- .sub(cbe[[kk]], .mul(list(.int(min1), ee[[i]][[1L]])))
            update(bb[[i]])
          }
          cbe[[kb]] <- .add(list(cbe[[kb]], .mul(list(.int(min1), .int(qstart), exponent))))
        }
      }
      if (last || length(cbq) == 1L || all(vapply(cbq, function(bq) is_one_num(bq[[2L]]), TRUE))) break
      qs <- vapply((if (bstart$type == "Mul") bstart$args else list(bstart)), function(bi) as.numeric(.bigint(ratq(bi))), 0)
      ql <- Reduce(function(a, b) a * b / gmp::gcd.bigz(a, b), lapply(qs, gmp::as.bigz))
      if (ql == 1) break
      b <- .pw(bstart, .int(ql))
      qlcm <- ql * qstart
      if (all(vapply((if (b$type == "Mul") b$args else list(b)), function(bi) is_one_num(ratq(bi)), TRUE))) last <- TRUE
    }
    if (is.null(cbq[[kb]])) next
    done[[length(done) + 1L]] <- list(base[[1L]], .mul(list(cbe[[kb]], make_rational(gmp::as.bigq(1, .bigint(base[[2L]]))))))
    cbq[[kb]] <- NULL; cbe[[kb]] <- NULL
  }
  for (k in names(cbq)) {
    b <- cbq[[k]][[1L]]; q <- cbq[[k]][[2L]]
    if ((b$type == "Pow" || (b$type == "Function" && b$name == "exp")) && !is_one_num(q) && !(.abe(b)[[2L]]$type %in% c("Integer", "Rational"))) {
      be <- .abe(b); b <- .pw(be[[1L]], .div(be[[2L]], q))
    } else if (!is_one_num(q)) b <- .pow(b, .inv(q))
    done[[length(done) + 1L]] <- list(b, cbe[[k]])
  }
  done
}

.coeff_isneg <- function(a) {
  if (a$type == "Mul") a <- a$args[[1L]]
  is_number_type(a) && isTRUE(ask_fact(a, "extended_negative"))
}

.powsimp_base <- function(expr, deep, force, recurse) {
  cp <- list(); nc_part <- list()
  for (term in expr$args) if (.is_comm(term)) cp[[length(cp) + 1L]] <- .abe(term) else nc_part[[length(nc_part) + 1L]] <- term
  for (i in seq_along(cp)) {
    b <- cp[[i]][[1L]]; e <- cp[[i]][[2L]]
    if (!(all(vapply(.as_numer_denom(b), .fa, TRUE, nm = "nonnegative")) || .fa(e, "integer") || force || .fa(b, "polar"))) next
    sp <- as_coeff_mul_split(e)
    if (sp[[1L]]$type %in% c("Integer", "Rational") && !is_one_num(sp[[1L]]) && !is_one_num(sp[[2L]]))
      cp[[i]] <- list(.pw(b, sp[[1L]]), sp[[2L]])
  }
  ek <- list(); eb <- list()
  for (p in cp) {
    b <- p[[1L]]; e <- p[[2L]]
    if (deep) e <- recurse(e)
    if (e$type == "Add" && (.fa(b, "positive") || .fa(e, "integer"))) {
      e <- factor_terms(e)
      if (.coeff_isneg(e)) { e <- .neg(e); b <- .inv(b) }
    }
    k <- e$key
    if (is.null(ek[[k]])) { ek[[k]] <- e; eb[[k]] <- list(b) } else eb[[k]][[length(eb[[k]]) + 1L]] <- b
  }
  pk <- list(); pe <- list()
  addp <- function(b, e) { k <- b$key; if (is.null(pk[[k]])) { pk[[k]] <<- b; pe[[k]] <<- list(e) } else pe[[k]][[length(pe[[k]]) + 1L]] <<- e }
  nterms <- function(z) if (z$type == "Add") sum(vapply(z$args, nterms, 0)) else if (z$type == "Mul") prod(vapply(z$args, nterms, 0)) else 1
  for (k in names(ek)) {
    e <- ek[[k]]; bases <- eb[[k]]
    if (length(bases) == 1L) new_base <- bases[[1L]]
    else if (.fa(e, "integer") || force) new_base <- .mul(bases)
    else {
      unk <- list(); nonneg <- list(); neg <- list()
      for (bi in bases) {
        if (.fa(bi, "negative")) neg[[length(neg) + 1L]] <- bi
        else if (.fa(bi, "nonnegative") || .fa(bi, "polar")) nonneg[[length(nonneg) + 1L]] <- bi
        else unk[[length(unk) + 1L]] <- bi
      }
      if ((length(unk) == 1L && !length(neg)) || (length(neg) == 1L && !length(unk))) {
        nonneg <- c(nonneg, unk, neg); unk <- list(); neg <- list()
      } else if (length(neg)) {
        nd <- .as_numer_denom(e)
        israt <- e$type %in% c("Integer", "Rational") || (.fa(nd[[1L]], "integer") && .fa(nd[[2L]], "integer"))
        if (israt) { neg <- lapply(neg, .neg); unk <- c(unk, rep(list(S_NegativeOne), length(neg))) }
        else { unk <- c(unk, neg); neg <- list() }
      }
      for (b in unk) addp(b, e)
      new_base <- .mul(c(nonneg, neg))
      xnb <- expand_mul(new_base, deep = FALSE)
      if ((if (xnb$type == "Add") length(xnb$args) else 1L) < nterms(new_base)) new_base <- factor_terms(xnb)
    }
    addp(new_base, e)
  }
  c_part <- list()
  for (k in names(pk)) for (ei in pe[[k]]) c_part[[length(c_part) + 1L]] <- .pw(pk[[k]], ei)
  .mul(c(c_part, nc_part))
}

#' @rdname powsimp
#' @export
powdenest <- function(e, force = FALSE, polar = FALSE) {
  eq <- as_sym(e)
  if (force) {
    den <- function(b, ex) {
      if (!(b$type == "Pow" || (b$type == "Function" && b$name == "exp"))) return(list(.tv(ask_fact(b, "positive")), new_sym("Pow", list(b, ex))))
      be <- .abe(b)
      den(be[[1L]], .mul(list(be[[2L]], ex)))
    }
    reps <- list()
    for (p in c(atoms(eq, "Pow"), Filter(function(f) f$name == "exp", atoms(eq, "Function")))) {
      pb <- .abe(p)
      if (pb[[1L]]$type == "Pow" || (pb[[1L]]$type == "Function" && pb[[1L]]$name == "exp")) {
        r <- den(pb[[1L]], pb[[2L]])
        if (!isFALSE(r[[1L]])) reps[[length(reps) + 1L]] <- list(p, r[[2L]])
      }
    }
    if (length(reps)) eq <- subs(eq, reps)
    ps <- posify(eq)
    return(xreplace(powdenest(ps$expr, force = FALSE, polar = polar), ps$reps))
  }
  if (polar) {
    pr <- polarify(eq)
    return(unpolarify(powdenest(unpolarify(pr[[1L]], exponents_only = TRUE)), pr[[2L]]))
  }
  new <- powsimp(eq)
  .xreplace_td(new, function(m) m$type == "Pow" || (m$type == "Function" && m$name == "exp"), .denest_pow)
}

# xreplace with a Transform: top down, a replaced node is not entered
.xreplace_td <- function(e, filter, f) {
  if (!is_sym(e)) return(e)
  if (filter(e)) {
    r <- f(e)
    if (!identical(r$key, e$key)) return(r)
  }
  if (!length(e$args)) return(e)
  args <- lapply(e$args, .xreplace_td, filter = filter, f = f)
  if (identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e else sym_rebuild(e, args)
}

.denest_pow <- function(eq) {
  be <- .abe(eq); b <- be[[1L]]; e <- be[[2L]]
  if (is_const(b, "E") && e$type == "Mul") {
    lg <- list(); other <- list()
    for (ei in e$args) {
      if (any(vapply(if (ei$type == "Add") ei$args else list(ei), function(ai) ai$type == "Function" && ai$name == "log", TRUE)))
        lg[[length(lg) + 1L]] <- ei else other[[length(other) + 1L]] <- ei
    }
    return(.pow(fcall("exp", list(logcombine(.mul(lg)))), .mul(other)))
  }
  bbe <- .abe(b)[[2L]]
  if (is_one_num(bbe) && !(b$type == "Mul" || (b$type == "Rational") || .fa(b, "positive"))) return(eq)
  polars <- list(); nonpolars <- list()
  for (bb in (if (b$type == "Mul") b$args else list(b))) {
    if (.fa(bb, "polar")) polars[[length(polars) + 1L]] <- .abe(bb) else nonpolars[[length(nonpolars) + 1L]] <- bb
  }
  if (length(polars) == 1L && polars[[1L]][[1L]]$type != "Mul")
    return(.mul(list(.pw(polars[[1L]][[1L]], .mul(list(polars[[1L]][[2L]], e))), powdenest(.pow(.mul(nonpolars), e)))))
  if (length(polars))
    return(.mul(c(lapply(polars, function(p) powdenest(.pow(p[[1L]], .mul(list(p[[2L]], e))))), list(powdenest(.pow(.mul(nonpolars), e))))))
  if (b$type == "Integer") {
    logb <- expand_log(fcall("log", list(b)))
    if (logb$type == "Mul") {
      c <- logb$args[[1L]]; logb <- logb$args[[2L]]
      return(.pow(logb$args[[1L]], .mul(list(e, c))))
    }
  }
  bf <- if (b$type == "Mul") b$args else list(b)
  if (b$type != "Mul" || any(vapply(bf, function(s) !length(s$args), TRUE))) return(eq)
  nc_gcd <- function(aa, bb) {
    a <- as_coeff_mul_split(aa); b2 <- as_coeff_mul_split(bb)
    cc <- .as_numer_denom(gcd(a[[1L]], b2[[1L]]))[[1L]]
    fa_ <- if (a[[2L]]$type == "Mul") a[[2L]]$args else list(a[[2L]])
    fb_ <- if (b2[[2L]]$type == "Mul") b2[[2L]]$args else list(b2[[2L]])
    common <- Filter(function(z) z$key %in% vapply(fb_, .k, ""), fa_)
    .keep_coeff(cc, .mul(common))
  }
  glogb <- expand_log(fcall("log", list(b)))
  if (glogb$type == "Add") {
    args <- glogb$args
    g <- Reduce(nc_gcd, args)
    if (!is_one_num(g)) {
      sp <- as_coeff_mul_split(g)
      glogb <- .keep_coeff(sp[[1L]], .mul(list(sp[[2L]], .add(lapply(args, function(a) .div(a, g))))))
    }
  }
  if ((glogb$type == "Function" && glogb$name == "log") || glogb$type != "Mul") {
    a0 <- glogb$args[[1L]]
    if (length(glogb$args) && (a0$type == "Pow" || (a0$type == "Function" && a0$name == "exp"))) {
      glogb <- .denest_pow(a0)
      gb <- .abe(glogb)
      ab <- tryCatch(.sort_num(fcall("Abs", list(gb[[2L]]))), error = function(err) NA)
      if (isTRUE(ab < 1)) return(.pw(gb[[1L]], .mul(list(gb[[2L]], e))))
    }
    return(eq)
  }
  add <- list(); other <- list()
  for (a in glogb$args) if (a$type == "Add") add[[length(add) + 1L]] <- a else other[[length(other) + 1L]] <- a
  .pow(fcall("exp", list(logcombine(.mul(add)))), .mul(list(e, .mul(other))))
}

#' @rdname powsimp
#' @export
logcombine <- function(e, force = FALSE) {
  expr <- as_sym(e)
  islog <- function(z) z$type == "Function" && z$name == "log" && length(z$args) == 1L
  f <- function(rv) {
    if (!(rv$type %in% c("Add", "Mul"))) return(rv)
    gooda <- function(a) !identical(a$key, S_NegativeOne$key) &&
      (isTRUE(ask_fact(a, "extended_real")) || (force && !isFALSE(ask_fact(a, "extended_real"))))
    goodlog <- function(l) { a <- l$args[[1L]]; isTRUE(ask_fact(a, "positive")) || (force && !isFALSE(ask_fact(a, "nonpositive"))) }
    other <- list(); logs <- list()
    # log1: key (a tuple of factors) -> list of list(coefficients, log)
    l1k <- list(); l1v <- list()
    add1 <- function(kfac, item) {
      k <- paste0("t:", paste(vapply(kfac, .k, ""), collapse = "\u0001"))
      if (is.null(l1k[[k]])) { l1k[[k]] <<- kfac; l1v[[k]] <<- list() }
      l1v[[k]][[length(l1v[[k]]) + 1L]] <<- item
    }
    for (a in .aargs(rv)) {
      if (islog(a) && goodlog(a)) add1(list(), list(list(), a))
      else if (a$type != "Mul") other[[length(other) + 1L]] <- a
      else {
        ot <- list(); co <- list(); lo <- list()
        for (ai in a$args) {
          if (.is_rat(ai) && num_sign(ai) < 0) { ot[[length(ot) + 1L]] <- S_NegativeOne; co[[length(co) + 1L]] <- .neg(ai) }
          else if (islog(ai) && goodlog(ai)) lo[[length(lo) + 1L]] <- ai
          else if (gooda(ai)) co[[length(co) + 1L]] <- ai
          else ot[[length(ot) + 1L]] <- ai
        }
        if (length(lo) > 1L) logs[[length(logs) + 1L]] <- list(ot, co, lo)
        else if (length(lo)) add1(ot, list(co, lo[[1L]]))
        else other[[length(other) + 1L]] <- a
      }
    }
    if (length(other) == 1L && islog(other[[1L]])) { add1(list(), list(list(), other[[1L]])); other <- list() }
    if (!length(logs) && !length(l1k)) return(rv)
    for (t in logs) {
      o <- t[[1L]]; ee <- t[[2L]]; l <- ordered(t[[3L]])
      e2 <- fcall("log", list(.pow(l[[1L]]$args[[1L]], .mul(ee))))
      for (li in l[-1L]) e2 <- fcall("log", list(.pow(li$args[[1L]], e2)))
      c <- .mul(o)
      if (islog(e2)) add1(list(c), list(list(), e2)) else other[[length(other) + 1L]] <- .mul(list(c, e2))
    }
    # the coefficient keys become expressions
    ek <- list(); ev <- list()
    for (k in names(l1k)) {
      kk <- .mul(l1k[[k]])
      val <- fnode("log", list(logcombine(.mul(lapply(l1v[[k]], function(cl) .pow(cl[[2L]]$args[[1L]], .mul(cl[[1L]])))), force = force)))
      ek[[kk$key]] <- kk; ev[[kk$key]] <- val
    }
    for (kk in ordered(unname(ek))) {
      k <- kk$key
      if (is.null(ek[[k]])) next
      nk <- .neg(kk)$key
      if (!is.null(ek[[nk]])) {
        num <- kk; den <- ek[[nk]]
        if (count_ops(num) > count_ops(den)) { t <- num; num <- den; den <- t }
        lnum <- ev[[num$key]]; lden <- ev[[den$key]]
        ek[[num$key]] <- NULL; ev[[num$key]] <- NULL; ek[[den$key]] <- NULL; ev[[den$key]] <- NULL
        other[[length(other) + 1L]] <- .mul(list(num, fnode("log", list(.div(lnum$args[[1L]], lden$args[[1L]])))))
      } else {
        other[[length(other) + 1L]] <- .mul(list(kk, ev[[k]]))
        ek[[k]] <- NULL; ev[[k]] <- NULL
      }
    }
    .add(other)
  }
  if (expr$type == "Mul") {
    um <- .umul(expr$args)
    cr <- as_coeff_mul_split(um); coefficient <- cr[[1L]]; rest <- cr[[2L]]
    if (coefficient$type == "Rational") {
      factors <- .margs(rest)
      for (i in seq_along(factors)) {
        term <- factors[[i]]
        if (term$type != "Add") next
        cd <- .as_coeff_dict(term)
        oneK <- S_One$key
        constant <- if (!is.null(cd$k[[oneK]])) cd$v[[oneK]] else S_Zero
        ks <- setdiff(names(cd$k), oneK)
        if (length(ks) && all(vapply(ks, function(k) { kk <- cd$k[[k]]; islog(kk) && !length(free_symbols(kk$args[[1L]])) && isTRUE(ask_fact(kk$args[[1L]], "positive")) }, TRUE))) {
          factors[[i]] <- .add(c(list(.mul(list(coefficient, constant))), lapply(ks, function(k) .keep_coeff(.mul(list(coefficient, cd$v[[k]])), cd$k[[k]]))))
          expr <- .umul(factors)
          break
        }
      }
    }
  }
  .bottom_up_args(expr, f)
}
