# Limits, after sympy/series/limits.py: Limit.doit tries the leading term
# of the series expansion, then the Gruntz algorithm (R/series-gruntz.R),
# then heuristics.

#' Limits
#'
#' `limit(e, x, x0, dir)` computes the limit of `e` as `x` approaches `x0`,
#' as SymPy's `limit()` (`Limit(e, x, x0, dir).doit(deep = FALSE)`). `dir`
#' is `"+"` (default), `"-"` or `"+-"` (two sided). `x0` may be `oo` or
#' `-oo`. `Limit()` builds an unevaluated limit. A limit that cannot be
#' computed is returned unevaluated.
#'
#' @param e An expression.
#' @param x The variable.
#' @param x0 The limit point.
#' @param dir Direction.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' limit(sin(x) / x, x, 0)
#' limit((1 + 1/x)^x, x, oo)
#' limit(x^2 * exp(-x), x, oo)
#' @export
limit <- function(e, x, x0, dir = "+") .limit_doit(Limit(e, x, x0, dir), deep = FALSE)

#' @rdname limit
#' @export
Limit <- function(e, x, x0, dir = "+") {
  e <- as_sym(e); x <- as_sym(x); x0 <- as_sym(x0)
  # Limit.__new__: oo is approached from below, -oo from above
  if (is_const(x0, "oo") || identical(x0$key, .mul(list(S_I, S_Infinity))$key)) dir <- "-"
  else if (is_const(x0, "-oo") || identical(x0$key, .mul(list(S_I, S_NegativeInfinity))$key)) dir <- "+"
  if (.has1(x0, x)) .notimpl_error(paste0("Limits approaching a variable point are not supported (", .str(x), " -> ", .str(x0), ")"))
  if (!(dir %in% c("+", "-", "+-"))) stop("direction must be one of '+', '-' or '+-', not ", dir, call. = FALSE)
  new_sym("Limit", list(e, x, x0, new_sym("LimitDir", name = dir, key = paste0("LD", dir))))
}
register_str_printer("Limit", function(e) paste0("Limit(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", ",
                                                 .str(e$args[[3L]]), ", dir='", e$args[[4L]]$name, "')"))
register_bound_vars("Limit", function(e) e$args[[2L]]$key)
register_doit("Limit", function(e) .limit_doit(e))

.lim_is_limit <- function(v) is_sym(v) && v$type == "Limit"
.lim_illegal <- function(e) any(vapply(list(S_NaN, S_Infinity, S_NegativeInfinity, S_ComplexInfinity), function(c) .contains_key(e, c$key), NA))

# Limit.doit
.limit_doit <- function(self, deep = TRUE) {
  e <- self$args[[1L]]; z <- self$args[[2L]]; z0 <- self$args[[3L]]; dir <- self$args[[4L]]$name
  if (dir == "+-") {
    r <- limit(e, z, z0, "+")
    l <- limit(e, z, z0, "-")
    if (.lim_is_limit(r) && .lim_is_limit(l) && identical(r$args[[1L]]$key, l$args[[1L]]$key)) return(self)
    if (identical(r$key, l$key)) return(l)
    if (isTRUE(ask_fact(r, "infinite")) && isTRUE(ask_fact(l, "infinite"))) return(S_ComplexInfinity)
    .gr_value_error(paste0("The limit does not exist since left hand limit = ", .str(l), " and right hand limit = ", .str(r)))
  }
  if (identical(z0$key, S_ComplexInfinity$key)) .notimpl_error("Limits at complex infinity are not implemented")
  if (isTRUE(ask_fact(z0, "infinite"))) {
    cdir <- sign(z0)
    cdir <- .div(cdir, Abs(cdir))
    e <- subs(e, z, .mul(list(cdir, z)))
    dir <- "-"
    z0 <- S_Infinity
  }
  if (deep) { e <- doit(e); z <- doit(z); z0 <- doit(z0) }
  if (identical(e$key, z$key)) return(z0)
  if (!.has1(e, z)) return(e)
  if (is_nan(z0)) return(S_NaN)
  if (.lim_illegal(e)) return(self)
  if (e$type == "Order") return(sym_rebuild(e, c(list(limit(e$args[[1L]], z, z0)), e$args[-1L])))
  cdir <- if (dir == "+") S_One else if (dir == "-") S_NegativeOne else S_Zero

  set_signs <- function(expr) {
    if (!length(expr$args) || expr$type %in% c("LimitDir")) return(expr)
    newargs <- lapply(expr$args, set_signs)
    if (!identical(vapply(newargs, .k, ""), vapply(expr$args, .k, ""))) expr <- sym_rebuild(expr, newargs)
    fl <- if (expr$type == "Function") expr$name else ""
    if (fl %in% c("Abs", "sign", "arg")) {
      sig <- tryCatch({
        s <- limit(expr$args[[1L]], z, z0, dir)
        if (isTRUE(ask_fact(s, "zero"))) s <- limit(.inv(expr$args[[1L]]), z, z0, dir)
        s
      }, NotImplementedError = function(err) NULL)
      if (is.null(sig)) return(expr)
      if (isTRUE(ask_fact(sig, "extended_real"))) {
        if (.is_true(Lt(sig, 0))) return(switch(fl, Abs = .neg(expr$args[[1L]]), sign = S_NegativeOne, arg = S_Pi))
        if (.is_true(Gt(sig, 0))) return(switch(fl, Abs = expr$args[[1L]], sign = S_One, arg = S_Zero))
      }
    }
    expr
  }
  if (.has_type(e, "Float")) e <- nsimplify(e)
  e <- set_signs(e)

  if (isTRUE(.so_is_meromorphic(e, z, z0))) {
    if (is_const(z0, "oo")) { newe <- subs(e, z, .inv(z)); cdir <- .neg(cdir) }
    else newe <- subs(e, z, .add(list(z, z0)))
    lt <- tryCatch(.so_leadterm(newe, z, NULL, cdir), ValueError = function(err) NULL)
    if (!is.null(lt)) {
      coeff <- lt[[1L]]; ex <- lt[[2L]]
      if (.is_true(Gt(ex, 0))) return(S_Zero)
      if (is_zero_num(ex)) return(coeff)
      odd <- is_number_type(ex) && (as.integer(trunc(num_value_dbl(ex))) %% 2L != 0L)
      if (is_one_num(cdir) || !odd) return(.mul(list(S_Infinity, sign(coeff))))
      if (identical(cdir$key, S_NegativeOne$key)) return(.mul(list(S_NegativeInfinity, sign(coeff))))
      return(S_ComplexInfinity)
    }
  }
  if (is_const(z0, "oo")) {
    if (e$type == "Mul") e <- factor_terms(e)
    asm <- list()
    for (k in c("positive", "negative", "real")) { v <- ask_fact(z, k); if (!is.na(v)) asm[[k]] <- v }
    dummy <- do.call(Dummy, c(list("z"), asm))
    newe <- subs(e, z, .inv(dummy))
    cdir <- .neg(cdir)
    newz <- dummy
  } else {
    newe <- subs(e, z, .add(list(z, z0)))
    newz <- z
  }
  lt <- tryCatch(.so_leadterm(newe, newz, NULL, cdir), ValueError = function(err) NULL,
                 NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
  if (is.null(lt)) {
    e <- powsimp(e)
    if (e$type == "Pow") {
      r <- .lim_pow_heuristics(e, z, z0, dir)
      if (!is.null(r)) return(r)
    }
    r <- tryCatch({
      coeff <- as_leading_term(newe, newz, NULL, cdir)
      if (!identical(coeff$key, newe$key) && (.has_fun(coeff, "exp") || .contains_key(coeff, S_E$key)))
        gruntz(coeff, newz, S_Zero, if (isTRUE(ask_fact(.se_re(cdir), "negative"))) "-" else "+")
      else NULL
    }, ValueError = function(err) NULL, NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
    if (!is.null(r)) return(r)
  } else {
    coeff <- lt[[1L]]; ex <- lt[[2L]]
    if (coeff$type == "AccumulationBounds" && is_zero_num(ex)) return(coeff)
    if (.lim_illegal(coeff)) return(self)
    if (!.has1(coeff, newz)) {
      if (isTRUE(ask_fact(ex, "positive"))) return(S_Zero)
      if (is_zero_num(ex)) return(coeff)
      if (isTRUE(ask_fact(ex, "negative"))) {
        if (is_one_num(cdir)) return(.mul(list(S_Infinity, sign(coeff))))
        if (identical(cdir$key, S_NegativeOne$key))
          return(.mul(list(S_NegativeInfinity, sign(coeff), .pow(S_NegativeOne, .add(list(S_One, ex))))))
        return(S_ComplexInfinity)
      }
      .notimpl_error(paste("Not sure of sign of", .str(ex)))
    }
  }
  if (isTRUE(ask_fact(z0, "extended_nonnegative")))
    e <- .map_bottom_up(e, function(t) if (t$type == "Function" && t$name == "factorial") fcall("gamma", list(.add(list(t$args[[1L]], S_One)))) else t)
  r <- tryCatch({
    v <- gruntz(e, z, z0, dir)
    if (is_nan(v)) .pole_error("")
    v
  }, PoleError = function(err) NULL, ValueError = function(err) NULL)
  if (is.null(r)) {
    r <- .lim_heuristics(e, z, z0, dir)
    if (is.null(r)) return(self)
  }
  r
}
.has_fun <- function(e, name) {
  if (e$type == "Function" && e$name == name) return(TRUE)
  any(vapply(e$args, .has_fun, NA, name = name))
}

# Limit.pow_heuristics
.lim_pow_heuristics <- function(e, z, z0, dir) {
  b1 <- e$args[[1L]]; e1 <- e$args[[2L]]
  if (!.has1(b1, z)) return(exp(limit(.mul(list(e1, log(b1))), z, z0, dir)))
  ex_lim <- limit(e1, z, z0, dir)
  base_lim <- limit(b1, z, z0, dir)
  if (is_one_num(base_lim) && (is_const(ex_lim, "oo") || is_const(ex_lim, "-oo")))
    return(exp(limit(.mul(list(e1, .sub(b1, S_One))), z, z0, dir)))
  if (is_const(base_lim, "-oo") && is_const(ex_lim, "oo")) return(S_ComplexInfinity)
  NULL
}

# heuristics
.lim_heuristics <- function(e, z, z0, dir) {
  rv <- NULL
  if (is_const(z0, "oo")) {
    rv <- limit(subs(e, z, .inv(z)), z, S_Zero, "+")
    if (.lim_is_limit(rv)) return(NULL)
  } else if (e$type %in% c("Mul", "Add", "Pow") || (e$type == "Function" && !is.null(.FUNS[[e$name]]))) {
    r <- list()
    for (a in e$args) {
      l <- limit(a, z, z0, dir)
      if (.contains_key(l, S_Infinity$key) && is.na(ask_fact(l, "finite"))) {
        if (e$type == "Add") {
          m <- factor_terms(e)
          if (m$type != "Mul") m <- together(m)
          if (m$type != "Mul") m <- sym_factor(e)
          if (m$type == "Mul") return(.lim_heuristics(m, z, z0, dir))
          return(NULL)
        }
        return(NULL)
      }
      if (.lim_is_limit(l) || is_nan(l)) return(NULL)
      r[[length(r) + 1L]] <- l
    }
    if (length(r)) {
      rv <- sym_rebuild(e, r)
      if (is_nan(rv) && e$type == "Mul" && any(vapply(r, function(v) v$type == "AccumulationBounds", NA))) {
        isab <- vapply(r, function(v) v$type == "AccumulationBounds", NA)
        if (any(!isab)) {
          e3 <- simplify(.mul(e$args[!isab]))
          rv <- .mul(c(list(limit(e3, z, z0, dir)), r[isab]))
        }
      }
      if (is_nan(rv)) {
        rat_e <- tryCatch(ratsimp(e), PolynomialError = function(err) NULL)
        if (is.null(rat_e) || is_nan(rat_e) || identical(rat_e$key, e$key)) return(NULL)
        return(limit(rat_e, z, z0, dir))
      }
    }
  }
  rv
}
