# Order terms, after sympy/series/order.py, with the Order handling of
# Add.flatten and Mul.flatten, Expr.getO/removeO/getn and the leading term
# machinery of Expr.as_leading_term and Add.extract_leading_order.
# An Order is new_sym("Order", list(expr, Tuple(x1, p1), ...)).

.so_vars <- function(o) lapply(o$args[-1L], function(t) t$args[[1L]])
.so_points <- function(o) lapply(o$args[-1L], function(t) t$args[[2L]])
.so_is_order <- function(e) is_sym(e) && e$type == "Order"
.so_is_inf_point <- function(p) is_const(p, "oo") || identical(p$key, .mul(list(S_Infinity, S_I))$key)
.so_is_ninf_point <- function(p) is_const(p, "-oo") || identical(p$key, .mul(list(S_NegativeInfinity, S_I))$key)

# Order.__new__
.so_new <- function(expr, vars = NULL, points = NULL) {
  expr <- as_sym(expr)
  if (is.null(vars)) {
    if (.so_is_order(expr)) { vars <- .so_vars(expr); points <- .so_points(expr) }
    else { vars <- free_symbols(expr); points <- rep(list(S_Zero), length(vars)) }
  } else if (is.null(points)) points <- rep(list(S_Zero), length(vars))
  vars <- lapply(vars, as_sym); points <- lapply(points, as_sym)
  if (!all(vapply(vars, function(v) v$type %in% c("Symbol", "Dummy", "Wild"), NA))) stop(.sympy_condition("TypeError", paste0("Variables are not symbols, got [", paste(vapply(vars, .str, ""), collapse = ", "), "]")))
  if (anyDuplicated(vapply(vars, .k, ""))) stop(.sympy_condition("ValueError", paste0("Variables are supposed to be unique symbols, got [", paste(vapply(vars, .str, ""), collapse = ", "), "]")))
  if (.so_is_order(expr)) {
    evp <- stats::setNames(.so_points(expr), vapply(.so_vars(expr), .k, ""))
    evv <- stats::setNames(.so_vars(expr), names(evp))
    nvp <- evp; nvv <- evv
    for (i in seq_along(vars)) {
      k <- vars[[i]]$key
      if (!is.null(nvp[[k]])) { if (!identical(nvp[[k]]$key, points[[i]]$key)) .notimpl_error("Mixing Order at different points is not supported.") }
      else { nvp[[k]] <- points[[i]]; nvv[[k]] <- vars[[i]] }
    }
    if (setequal(names(evp), names(nvp))) return(expr)
    vars <- unname(nvv); points <- unname(nvp[names(nvv)])
  }
  if (is_nan(expr)) return(S_NaN)
  if (any(vapply(points, function(p) any(vapply(vars, function(v) .has1(p, v), NA)), NA))) stop("Got a point that contains a variable.", call. = FALSE)
  if (length(vars)) {
    if (any(vapply(points, function(p) !identical(p$key, points[[1L]]$key), NA))) .notimpl_error("Multivariable orders at different points are not supported.")
    p0 <- points[[1L]]
    if (.so_is_inf_point(p0) || .so_is_ninf_point(p0)) {
      sg <- if (.so_is_inf_point(p0)) S_One else S_NegativeOne
      ds <- lapply(vars, function(v) Dummy())
      s <- lapply(seq_along(vars), function(i) list(vars[[i]], .mul(list(sg, .inv(ds[[i]])))))
      # rs = {sg/v: sg/k} with v = sg/d: the dummy d goes back to sg/k
      rs <- lapply(seq_along(vars), function(i) list(ds[[i]], .mul(list(sg, .inv(vars[[i]])))))
      ps <- rep(list(S_Zero), length(points))
      args <- ds
    } else if (!is_zero_num(p0)) {
      ds <- lapply(vars, function(v) Dummy())
      s <- lapply(seq_along(vars), function(i) list(vars[[i]], .add(list(ds[[i]], p0))))
      rs <- lapply(seq_along(vars), function(i) list(ds[[i]], .sub(vars[[i]], p0)))
      ps <- rep(list(S_Zero), length(points))
      args <- ds
    } else { s <- list(); rs <- list(); ps <- points; args <- vars }
    if (length(s)) expr <- subs(expr, s)
    if (expr$type == "Add") expr <- sym_factor(expr)
    if (length(vars) > 1L) expr <- expand(expr)
    old_key <- ""
    while (!identical(old_key, expr$key)) {
      old_key <- expr$key
      if (expr$type == "Add") {
        lst <- .so_extract_leading_order(expr, args)
        expr <- .add(lapply(lst, function(p) p[[2L]]$args[[1L]]))
      } else if (!is_zero_num(expr)) {
        lt <- tryCatch(.so_alt_multi(expr, args), PoleError = function(err) NULL)
        if (!is.null(lt)) expr <- lt
        else if (!(expr$type == "Function" || all(vapply(expr$args, function(a) a$type == "Function", NA)))) {
          pts <- lapply(seq_along(args), function(i) list(args[[i]], ps[[i]]))
          orders <- lapply(expr$args, function(a) {
            lta <- tryCatch(.so_alt_multi(a, args), PoleError = function(err) a)
            if (any(vapply(args, function(v) identical(v$key, lta$key), NA))) .so_new(lta, lapply(pts, `[[`, 1L), lapply(pts, `[[`, 2L)) else .so_new(lta)
          })
          if (expr$type == "Add") {
            ne <- .so_new(.add(orders), lapply(pts, `[[`, 1L), lapply(pts, `[[`, 2L))
            if (ne$type == "Add") ne <- .so_new(.add(lapply(ne$args, function(a) a$args[[1L]])), lapply(pts, `[[`, 1L), lapply(pts, `[[`, 2L))
            expr <- ne$args[[1L]]
          } else if (expr$type == "Mul") expr <- .mul(lapply(orders, function(a) a$args[[1L]]))
          else if (expr$type == "Pow") expr <- fcall("exp", list(.mul(list(expr$args[[2L]], fcall("log", list(expr$args[[1L]]))))))
        }
        if (isTRUE(ask_fact(expr, "zero"))) expr <- S_Zero
        else expr <- .mul(.so_indep(expr, args)[[2L]])
        expr <- expand_power_base(expr)
        expr <- expand_log(expr)
        if (length(args) == 1L) {
          x <- args[[1L]]
          margs <- .mi_make_args_mul(.mul(.so_indep(expr, list(x))[[2L]]))
          for (i in seq_along(margs)) {
            t <- margs[[i]]
            if (t$type == "Pow") {
              b <- t$args[[1L]]; q <- t$args[[2L]]
              isx <- function(z) identical(z$key, x$key) || identical(z$key, .neg(x)$key)
              if (isx(b) && isTRUE(ask_fact(q, "real")) && !.has1(q, x)) margs[[i]] <- .pow(x, q)
              else if (b$type == "Pow" && !.has1(b$args[[2L]], x)) {
                bb <- b$args[[1L]]; r <- b$args[[2L]]
                if (isx(bb) && isTRUE(ask_fact(r, "real"))) margs[[i]] <- .pow(x, .mul(list(r, q)))
              } else if (b$type == "Mul" && identical(b$args[[1L]]$key, S_NegativeOne$key)) {
                b <- .neg(b)
                if (b$type == "Pow" && !.has1(b$args[[2L]], x)) {
                  bb <- b$args[[1L]]; r <- b$args[[2L]]
                  if (isx(bb) && isTRUE(ask_fact(r, "real"))) margs[[i]] <- .pow(x, .mul(list(r, q)))
                }
              }
            }
          }
          expr <- .mul(margs)
        }
      }
    }
    if (length(rs)) expr <- subs(expr, rs)
  }
  if (.so_is_order(expr)) expr <- expr$args[[1L]]
  if (!any(vapply(vars, function(v) .has1(expr, v), NA)) && !isTRUE(ask_fact(expr, "zero"))) expr <- S_One
  if (length(vars)) {
    o <- order(vapply(vars, function(v) paste(.sort_key(v), collapse = "\r"), ""))
    vars <- vars[o]; points <- points[o]
  }
  new_sym("Order", c(list(expr), lapply(seq_along(vars), function(i) Tuple(vars[[i]], points[[i]]))))
}

# expr.as_independent(*syms, as_Add=False): list(independent, dependent factors)
.so_indep <- function(e, syms) {
  fs <- .mi_make_args_mul(e)
  dep <- vapply(fs, function(f) any(vapply(syms, function(s) .has1(f, s), NA)), NA)
  list(.mul(fs[!dep]), fs[dep])
}

.so_alt_multi <- function(e, syms) { for (x in syms) e <- as_leading_term(e, x); e }

# Add.extract_leading_order
.so_extract_leading_order <- function(e, syms, points = NULL) {
  if (is.null(points)) points <- rep(list(S_Zero), length(syms))
  lst <- list()
  for (f in e$args) {
    of <- .so_new(f, syms, points)
    drop <- FALSE
    for (p in lst) if (isTRUE(.so_contains(p[[2L]], of)) && !identical(p[[2L]]$key, of$key)) { drop <- TRUE; break }
    if (drop) next
    nl <- list(list(f, of))
    for (p in lst) {
      if (isTRUE(.so_contains(of, p[[2L]])) && !identical(p[[2L]]$key, of$key)) next
      nl[[length(nl) + 1L]] <- p
    }
    lst <- nl
  }
  lst
}

# Order.contains: TRUE, FALSE or NA
.so_contains <- function(o, expr) {
  expr <- as_sym(expr)
  if (isTRUE(ask_fact(expr, "zero"))) return(TRUE)
  if (is_nan(expr)) return(FALSE)
  pts <- .so_points(o)
  point <- if (length(pts)) pts[[1L]] else S_Zero
  ovars <- .so_vars(o)
  oexpr <- o$args[[1L]]
  if (.so_is_order(expr)) {
    if (any(vapply(.so_points(expr), function(p) !identical(p$key, point$key), NA)) || any(vapply(pts, function(p) !identical(p$key, point$key), NA))) return(NA)
    eexpr <- expr$args[[1L]]
    if (identical(eexpr$key, oexpr$key)) {
      ok <- vapply(o$args[-1L], .k, "")
      return(all(vapply(expr$args[-1L], function(t) t$key %in% ok, NA)))
    }
    if (eexpr$type == "Add") {
      r <- vapply(eexpr$args, function(a) .so_contains(o, a), NA)
      return(all(r))
    }
    if (oexpr$type == "Add" && is_zero_num(point)) {
      for (a in oexpr$args) if (isTRUE(.so_contains(new_sym("Order", c(list(a), o$args[-1L])), expr))) return(TRUE)
      return(FALSE)
    }
    evars <- .so_vars(expr)
    ek <- vapply(evars, .k, ""); ovk <- vapply(ovars, .k, "")
    common <- if (length(ovars) && length(evars)) ovars[ovk %in% ek] else if (length(ovars)) ovars else evars
    if (!length(common)) return(NA)
    if (oexpr$type == "Pow" && length(ovars) == 1L && identical(ovk, ek)) {
      symbol <- ovars[[1L]]
      other <- .mul(.so_indep(eexpr, list(symbol))[[2L]])
      if (other$type == "Pow" && identical(other$args[[1L]]$key, symbol$key) && identical(oexpr$args[[1L]]$key, symbol$key)) {
        rv <- NA
        if (is_zero_num(point)) rv <- ask_fact(.sub(oexpr$args[[2L]], other$args[[2L]]), "nonpositive")
        if (isTRUE(ask_fact(point, "infinite"))) rv <- ask_fact(.sub(oexpr$args[[2L]], other$args[[2L]]), "nonnegative")
        if (!is.na(rv)) return(rv)
      }
    }
    ratio <- powsimp(.div(oexpr, eexpr), deep = TRUE, combine = "exp")
    r <- NULL
    for (s in common) {
      l <- tryCatch(.so_limit_noheur(ratio, s, point), error = function(err) NULL)
      lv <- if (is.null(l) || l$type == "Limit") NA else !is_zero_num(l)
      if (is.null(r)) r <- lv
      else if (!identical(r, lv)) return(NA)
    }
    return(r)
  }
  if (oexpr$type == "Pow" && length(ovars) == 1L) {
    symbol <- ovars[[1L]]
    other <- .mul(.so_indep(expr, list(symbol))[[2L]])
    if (other$type == "Pow" && identical(other$args[[1L]]$key, symbol$key) && identical(oexpr$args[[1L]]$key, symbol$key)) {
      rv <- NA
      if (is_zero_num(point)) rv <- ask_fact(.sub(oexpr$args[[2L]], other$args[[2L]]), "nonpositive")
      if (isTRUE(ask_fact(point, "infinite"))) rv <- ask_fact(.sub(oexpr$args[[2L]], other$args[[2L]]), "nonnegative")
      if (!is.na(rv)) return(rv)
    }
  }
  .so_contains(o, .so_new(expr, ovars, pts))
}

# Limit(ratio, s, point).doit(heuristics=False); the limit engine stands in
# until Limit.doit is ported
.so_limit_noheur <- function(e, s, p) limit(e, s, p)

# Order.as_expr_variables
.so_expr_variables <- function(o, order_symbols) {
  if (is.null(order_symbols)) return(list(o$args[[1L]], o$args[-1L]))
  osp <- lapply(order_symbols, function(t) t$args[[2L]])
  if (!all(vapply(osp, function(p) identical(p$key, osp[[1L]]$key), NA)) && !all(vapply(.so_points(o), function(p) identical(p$key, .so_points(o)[[1L]]$key), NA)))
    .notimpl_error("Order at points other than 0 or oo not supported.")
  if (length(order_symbols) && length(o$args) > 1L && !identical(osp[[1L]]$key, .so_points(o)[[1L]]$key))
    .notimpl_error("Multiplying Order at different points is not supported.")
  d <- stats::setNames(order_symbols, vapply(order_symbols, function(t) t$args[[1L]]$key, ""))
  for (t in o$args[-1L]) if (is.null(d[[t$args[[1L]]$key]])) d[[t$args[[1L]]$key]] <- t
  d <- unname(d)
  d <- d[order(vapply(d, function(t) paste(.sort_key(t$args[[1L]]), collapse = "\r"), ""))]
  list(o$args[[1L]], d)
}

# Expr.getO, removeO and getn
.so_getO <- function(e) {
  if (.so_is_order(e)) return(e)
  if (e$type == "Add") {
    os <- Filter(.so_is_order, e$args)
    if (length(os)) return(if (length(os) == 1L) os[[1L]] else new_sym("Add", os))
  }
  NULL
}
.so_removeO <- function(e) {
  if (.so_is_order(e)) return(S_Zero)
  if (e$type == "Add") return(.add(Filter(Negate(.so_is_order), e$args)))
  e
}
.so_getn <- function(e) {
  o <- .so_getO(e)
  if (is.null(o)) return(NULL)
  if (.so_is_order(o)) {
    oe <- o$args[[1L]]
    if (is_one_num(oe)) return(S_Zero)
    if (oe$type %in% c("Symbol", "Dummy")) return(S_One)
    if (oe$type == "Pow") return(oe$args[[2L]])
    if (oe$type == "Mul") for (oi in oe$args) {
      if (oi$type %in% c("Symbol", "Dummy")) return(S_One)
      if (oi$type == "Pow") {
        syms <- Filter(function(z) z$type %in% c("Symbol", "Dummy"), atoms(oi))
        if (length(syms) == 1L) {
          oi2 <- subs(oi, syms[[1L]], Dummy("x", positive = TRUE))
          if (oi2$type == "Pow" && oi2$args[[1L]]$type %in% c("Symbol", "Dummy") && oi2$args[[2L]]$type %in% c("Integer", "Rational")) return(abs(oi2$args[[2L]]))
        }
      }
    }
  }
  .notimpl_error(paste("not sure of order of", .str(o)))
}

# ---------------------------------------------------------------------------
# Expr.as_leading_term and its _eval_as_leading_term methods

.pole_error <- function(msg) stop(structure(class = c("PoleError", "error", "condition"), list(message = msg, call = NULL)))

#' @rdname series
#' @export
as_leading_term <- function(e, x, logx = NULL, cdir = 0) {
  e <- as_sym(e); x <- as_sym(x); cdir <- as_sym(cdir)
  if (!(x$type %in% c("Symbol", "Dummy", "Wild"))) stop("expecting a Symbol but got ", .str(x), call. = FALSE)
  if (!(x$key %in% names(.free_sym_list(e)))) return(e)
  obj <- .so_eval_alt(e, x, logx, cdir)
  if (is.null(obj)) .notimpl_error(paste0("as_leading_term(", .str(e), ", ", .str(x), ")"))
  powsimp(obj, deep = TRUE, combine = "exp")
}


.so_eval_alt <- function(e, x, logx, cdir) {
  switch(e$type,
    Add = .so_alt_add(e, x, logx, cdir),
    Mul = .mul(lapply(.sympy_args(e), function(t) as_leading_term(t, x, logx, cdir))),
    Pow = .so_alt_pow(e, x, logx, cdir),
    Order = e,
    Function = {
      h <- .so_alt_hooks[[e$name]]
      if (!is.null(h)) h(e, x, logx, cdir) else .so_alt_function(e, x, logx, cdir)
    },
    {
      h <- .so_alt_hooks[[e$type]]
      if (!is.null(h)) h(e, x, logx, cdir) else e
    })
}

# Function._eval_as_leading_term
.so_alt_function <- function(e, x, logx, cdir) {
  args <- lapply(e$args, function(a) as_leading_term(a, x, logx))
  o <- .so_new(S_One, list(x))
  if (any(vapply(args, function(a) .has1(a, x) && isTRUE(.so_contains(o, a)), NA)))
    .notimpl_error(paste(e$name, "has no _eval_as_leading_term routine"))
  e
}

# Add._eval_as_leading_term
.so_alt_add <- function(e, x, logx, cdir) {
  o <- .so_getO(e)
  if (is.null(o)) o <- new_sym("Order", list(S_Zero))
  old <- .so_removeO(e)
  if (.has_type(old, "Piecewise")) old <- piecewise_fold(old)
  if (any(vapply(e$args, function(a) a$type == "Function" && a$name == "log", NA)))
    old <- expand(old, deep = TRUE, log = TRUE, mul = FALSE, power_exp = FALSE, power_base = FALSE, multinomial = FALSE, basic = FALSE, force = FALSE, factor = FALSE)
  expr <- expand_mul(old)
  if (expr$type != "Add") return(as_leading_term(expr, x, logx, cdir))
  infinite <- Filter(function(t) isTRUE(ask_fact(t, "infinite")), .sympy_args(expr))
  .logx <- if (is.null(logx)) Dummy("logx") else logx
  lts <- lapply(.sympy_args(expr), function(t) as_leading_term(t, x, .logx, cdir))
  mn <- NULL; new_expr <- S_Zero
  for (term in lts) {
    od <- .so_new(term, list(x))
    if (is.null(mn) || !isTRUE(.so_contains(mn, od))) {
      # `order not in min`: a None from contains is a TypeError in SymPy
      if (!is.null(mn) && is.na(.so_contains(mn, od))) return(expr)
      mn <- od; new_expr <- term
    } else {
      cm <- .so_contains(od, mn)
      if (is.na(cm)) return(expr)
      if (cm) new_expr <- .add(list(new_expr, term))
    }
  }
  if (is.null(logx)) new_expr <- subs(new_expr, .logx, fcall("log", list(x)))
  iz <- ask_fact(new_expr, "zero")
  if (is.na(iz)) {
    new_expr <- cancel(trigsimp(new_expr))
    iz <- ask_fact(new_expr, "zero")
  }
  if (isTRUE(iz)) {
    n0 <- tryCatch(.so_getn(.add(list(S_Zero, mn))), NotImplementedError = function(err) S_One)
    if (is.null(n0) || .has_type(n0, "Symbol")) n0 <- S_One
    res <- .so_new(S_One)
    incr <- S_One
    while (.so_is_order(res)) {
      res <- trigsimp(powsimp(cancel(.ns_eval(old, x, .add(list(n0, incr)), logx, cdir))))
      incr <- .mul(list(incr, S_Two))
    }
    return(as_leading_term(res, x, logx, cdir))
  }
  if (is_nan(new_expr)) return(.add(c(infinite, list(o))))
  new_expr
}

# Pow._eval_as_leading_term
.so_alt_pow <- function(e, x, logx, cdir) {
  b <- e$args[[1L]]; ex <- e$args[[2L]]
  if (.has1(ex, x)) return(as_leading_term(fcall("exp", list(.mul(list(ex, fcall("log", list(b)))))), x, logx, cdir))
  f <- tryCatch(as_leading_term(b, x, logx, cdir), PoleError = function(err) NULL)
  if (is.null(f)) return(e)
  if (!isTRUE(ask_fact(ex, "integer")) && isTRUE(ask_fact(f, "negative")) && !.has1(f, x)) {
    ndir <- .so_dir(.sub(b, f), x, cdir)
    imd <- ask_fact(fcall("im", list(ndir)), "negative")
    if (isTRUE(imd)) return(.mul(list(.pow(f, ex), .pow(S_NegativeOne, .mul(list(.int(-2), ex))))))
    if (isTRUE(ask_fact(fcall("im", list(ndir)), "zero"))) {
      ll <- .so_eval_alt(fcall("log", list(b)), x, logx, cdir)
      if (identical(ask_fact(ll, "infinite"), FALSE)) return(fcall("exp", list(.mul(list(ex, ll)))))
    }
  }
  .pow(f, ex)
}

# Expr.dir: the direction in which the expression leaves its value at
# x = 0, from the first nonzero derivative there
.so_dir <- function(e, x, cdir) {
  if (isTRUE(ask_fact(e, "zero"))) return(S_Zero)
  minexp <- S_Zero
  arg <- e
  coeff <- S_Zero
  while (!is_zero_num(arg)) {
    minexp <- .add(list(minexp, S_One))
    arg <- .diff1(arg, x)
    coeff <- subs(arg, x, S_Zero)
    if (is_nan(coeff)) coeff <- limit(arg, x, S_Zero)
    if (identical(coeff$key, S_ComplexInfinity$key)) {
      coeff <- tryCatch({
        lt <- leadterm(arg, x)
        c0 <- if (is.list(lt) && !is_sym(lt)) lt[[1L]] else lt
        if (.has1(c0, fcall("log", list(x)))) stop("log in leadterm")
        c0
      }, error = function(err) limit(arg, x, S_Zero))
    }
    if (!is_zero_num(coeff)) break
  }
  .mul(list(coeff, .pow(as_sym(cdir), minexp)))
}
