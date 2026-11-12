# Order terms and series expansion: the public API over R/series-order.R
# and R/series-nseries.R, after sympy/series/order.py and Expr.series in
# sympy/core/expr.py.

#' Order terms
#'
#' `O(expr, ...)` (alias `Order()`) is SymPy's big O term. The extra
#' arguments are the variables (each a symbol, or `list(x, x0)` for an
#' expansion point other than 0); by default the free symbols of `expr` at
#' 0. Added to an expression it absorbs the terms it contains;
#' `removeO()` drops it and `getO()` returns it (or `NULL`).
#'
#' @param expr An expression.
#' @param ... Variables, or `list(variable, point)` pairs.
#' @param e An expression with an order term.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' O(x^3)
#' x + x^2 + O(x^2)
#' O(1/x, list(x, Inf))
#' @export
O <- function(expr, ...) {
  args <- list(...)
  if (!length(args)) return(.so_new(as_sym(expr)))
  vars <- list(); pts <- list()
  for (a in args) {
    if (is.null(a)) next
    if (is_sym(a) && a$type == "Tuple") { vars <- c(vars, list(a$args[[1L]])); pts <- c(pts, list(a$args[[2L]])) }
    else if (is.list(a) && !is_sym(a)) { vars <- c(vars, list(as_sym(a[[1L]]))); pts <- c(pts, list(as_sym(a[[2L]]))) }
    else { vars <- c(vars, list(as_sym(a))); pts <- c(pts, list(S_Zero)) }
  }
  if (!length(vars)) return(.so_new(as_sym(expr)))
  .so_new(as_sym(expr), vars, pts)
}

#' @rdname O
#' @export
Order <- O

#' @rdname O
#' @export
removeO <- function(e) .so_removeO(as_sym(e))

#' @rdname O
#' @export
getO <- function(e) .so_getO(as_sym(e))

register_str_printer("Order", function(e) {
  vars <- .so_vars(e); pts <- .so_points(e)
  if (!length(vars) || all(vapply(pts, is_zero_num, NA))) {
    if (length(vars) <= 1L) return(paste0("O(", .str(e$args[[1L]]), ")"))
    return(paste0("O(", paste(vapply(c(list(e$args[[1L]]), vars), .str, ""), collapse = ", "), ")"))
  }
  paste0("O(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
})
register_rebuilder("Order", function(e, args) .so_new(args[[1L]], lapply(args[-1L], function(t) t$args[[1L]]), lapply(args[-1L], function(t) t$args[[2L]])))
register_diff_hook("Order", function(e, x) {
  d <- .diff1(e$args[[1L]], x)
  if (is_zero_num(d)) return(e)
  .so_new(d, .so_vars(e), .so_points(e))
})

# Order._eval_subs
register_subs_hook("Order", function(e, old, new) {
  vars <- .so_vars(e); pts <- .so_points(e)
  i <- match(old$key, vapply(vars, .k, ""))
  if (is.na(i)) return(NULL)
  newexpr <- subs(e$args[[1L]], old, new)
  newvars <- vars; newpt <- pts
  if (new$type %in% c("Symbol", "Dummy")) newvars[[i]] <- new
  else {
    syms <- free_symbols(new)
    if (length(syms) == 1L || old$key %in% vapply(syms, .k, "")) {
      var <- if (old$key %in% vapply(syms, .k, "")) vars[[i]] else syms[[1L]]
      no <- .so_getO(new)
      if (!is.null(no) && identical(limit(no$args[[1L]], var, .so_points(no)[[1L]])$key, pts[[i]]$key))
        return(.so_new(newexpr, list(var), list(.so_points(no)[[1L]])))
      point <- subs(new, var, pts[[i]])
      if (!identical(point$key, pts[[i]]$key)) {
        d <- Dummy()
        sol <- solveset(.sub(old, subs(new, var, d)), d)
        sv <- .so_first_solution(sol)
        point <- limit(sv, old, pts[[i]])
      }
      newvars[[i]] <- var; newpt[[i]] <- point
    } else if (!(old$key %in% vapply(syms, .k, ""))) {
      newvars <- newvars[-i]; newpt <- newpt[-i]
      if (!length(syms) && identical(new$key, pts[[i]]$key)) {
        newvars <- c(newvars, syms); newpt <- c(newpt, rep(list(S_Zero), length(syms)))
      }
    } else return(NULL)
  }
  .so_new(newexpr, newvars, newpt)
})

# the first element of iter(sol), as dict(zip((d,), sol)) takes it
.so_first_solution <- function(sol) {
  t <- sol$type
  if (t == "FiniteSet" && length(sol$args)) return(sol$args[[1L]])
  if (t == "Complement" && sol$args[[1L]]$type == "FiniteSet") {
    rest <- Filter(function(z) !isTRUE(.is_true(contains(sol$args[[2L]], z))), sol$args[[1L]]$args)
    if (length(rest)) return(rest[[1L]])
  }
  if (t == "ImageSet") {
    lam <- sol$args[[1L]]; base <- sol$args[[2L]]
    v <- lam$args[[1L]]$args
    first <- switch(base$name %||% "", Integers = S_Zero, Naturals0 = S_Zero, Naturals = S_One, NULL)
    if (length(v) == 1L && !is.null(first)) return(subs(lam$args[[2L]], v[[1L]], first))
  }
  if (t == "Union") return(.so_first_solution(sol$args[[1L]]))
  stop("could not invert the substitution", call. = FALSE)
}

# Add.flatten with Order terms: an order absorbs the terms (and orders) it
# contains; the orders come last
.add_with_order <- function(args, types) {
  ords <- list()
  for (o in args[types == "Order"]) {
    if (isTRUE(ask_fact(o$args[[1L]], "zero"))) next
    if (any(vapply(ords, function(o1) isTRUE(.so_contains(o1, o)), NA))) next
    ords <- c(list(o), Filter(function(o1) !isTRUE(.so_contains(o, o1)), ords))
  }
  rest <- .add(args[types != "Order"])
  terms <- if (is_zero_num(rest)) list() else if (rest$type == "Add") rest$args else list(rest)
  kept <- Filter(function(t) !any(vapply(ords, function(o) isTRUE(.so_contains(o, t)), NA)), terms)
  if (!length(ords)) return(.add(kept))
  if (!length(kept) && length(ords) == 1L) return(ords[[1L]])
  all_terms <- c(kept, ords)
  new_sym("Add", all_terms, extra = list(types = vapply(all_terms, `[[`, "", "type"), keys = vapply(all_terms, .k, "")))
}

# Mul.flatten with Order factors: Order(product, *order_symbols)
.mul_with_order <- function(args, types) {
  order_symbols <- NULL
  exprs <- list()
  for (a in args) {
    if (a$type == "Order") {
      r <- .so_expr_variables(a, order_symbols)
      exprs[[length(exprs) + 1L]] <- r[[1L]]; order_symbols <- r[[2L]]
    } else exprs[[length(exprs) + 1L]] <- a
  }
  prod <- .mul(exprs)
  if (!length(order_symbols)) return(.so_new(prod))
  .so_new(prod, lapply(order_symbols, function(t) t$args[[1L]]), lapply(order_symbols, function(t) t$args[[2L]]))
}

# Order._eval_power
.order_pow <- function(b, e) {
  if (is_number_type(e) && isTRUE(ask_fact(e, "nonnegative"))) return(.so_new(.pow(b$args[[1L]], e), .so_vars(b), .so_points(b)))
  if (e$type == "Order" && is_one_num(e$args[[1L]])) return(b)
  NULL
}

# A Python generator: nxt() returns the next value, or NULL once exhausted.
.py_generator <- function(nxt) structure(list(nxt = nxt), class = "sym_generator")
.gen_map <- function(g, f) .py_generator(function() { v <- g$nxt(); if (is.null(v)) NULL else f(v) })
.gen_list <- function(l) { i <- 0L; .py_generator(function() { i <<- i + 1L; if (i > length(l)) NULL else l[[i]] }) }
# yield_lseries(self.removeO()._eval_lseries(x, logx=logx, cdir=cdir)) (sympy/core/expr.py:3190)
.lseries_gen <- function(e, x, logx, cdir) {
  k <- 0L; buf <- list()
  .py_generator(function() {
    k <<- k + 1L
    # Known ceiling: each call recomputes the first k terms (O(k^2) nseries calls);
    # keep the _eval_lseries state between calls if long expansions matter.
    if (length(buf) < k) buf <<- .so_lseries(e, x, k, logx, cdir, nmax = Inf)
    if (length(buf) < k) NULL else buf[[k]]
  })
}

# Expr.series; n = NULL gives the lazy series (a sym_generator)
.ns_series_driver <- function(e, x, x0, n, dir, logx, cdir) {
  e <- as_sym(e)
  lazy <- is.null(n)
  if (!(x$key %in% names(.free_sym_list(e)))) return(if (lazy) .gen_list(list(e)) else e)
  if (!(dir %in% c("+", "-"))) stop(.sympy_condition("ValueError", "Dir must be '+' or '-'"))
  if (!lazy) {
    n <- as.integer(n)
    if (n < 0L) stop(.sympy_condition("ValueError", "Number of terms should be nonnegative"))
  }
  x0 <- as_sym(x0); cdir <- as_sym(cdir)
  if (!is_zero_num(cdir)) {
    if (isTRUE(ask_fact(cdir, "real"))) dir <- if (isTRUE(ask_fact(cdir, "positive"))) "+" else "-"
    else dir <- if (isTRUE(ask_fact(fcall("im", list(cdir)), "positive"))) "+" else "-"
  } else if (!is_zero_num(x0) && isTRUE(ask_fact(x0, "infinite"))) cdir <- simplify(fcall("sign", list(x0)))
  else if (dir == "+") cdir <- S_One
  else if (dir == "-") cdir <- S_NegativeOne
  cdir <- .div(cdir, Abs(cdir))
  if (!is_zero_num(x0) && isTRUE(ask_fact(x0, "infinite"))) {
    if (lazy) return(.gen_map(.ns_series_driver(subs(e, x, .div(cdir, x)), x, S_Zero, n, "+", NULL, 1), function(si) subs(si, x, .div(cdir, x))))
    s <- tryCatch(.ns_series_driver(subs(e, x, .div(cdir, x)), x, S_Zero, n, "+", NULL, 1), PoleError = function(err) NULL)
    if (is.null(s)) return(subs(aseries(subs(e, x, .mul(list(cdir, x))), x, n), x, .mul(list(cdir, x))))
    return(subs(s, x, .div(cdir, x)))
  }
  if (!is_zero_num(x0) || !is_one_num(cdir)) {
    s <- .ns_series_driver(subs(e, x, .add(list(x0, .mul(list(cdir, x))))), x, S_Zero, n, "+", logx, 1)
    if (lazy) return(.gen_map(s, function(si) subs(si, x, .sub(.div(x, cdir), .div(x0, cdir)))))
    return(subs(s, x, .sub(.div(x, cdir), .div(x0, cdir))))
  }
  if ((is.na(ask_fact(x, "positive")) && is.na(ask_fact(x, "negative"))) || !(x$type %in% c("Symbol", "Dummy"))) {
    xpos <- Dummy("x", positive = TRUE)
    rv <- .ns_series_driver(subs(e, x, xpos), xpos, x0, n, dir, logx, cdir)
    if (lazy) return(.gen_map(rv, function(si) subs(si, xpos, x)))
    return(subs(rv, xpos, x))
  }
  if (lazy) return(.lseries_gen(.so_removeO(e), x, logx, cdir))
  s1 <- .ns_eval(e, x, .int(n), logx, cdir)
  o <- .so_getO(s1)
  if (!is.null(o)) {
    ngot <- .so_getn(s1)
    if (.is_true(Gt(ngot, n))) {
      if (n != 0L) s1 <- .add(list(s1, subs(o, x, .pow(x, .div(.int(n), ngot)))))
      else s1 <- .add(list(s1, .so_new(S_One, list(x))))
    } else if (.is_true(Lt(ngot, n))) {
      done <- FALSE
      for (more in 1:8) {
        s1 <- .ns_eval(e, x, .int(n + more), logx, cdir)
        newn <- .so_getn(s1)
        if (!identical(newn$key, ngot$key)) {
          ndo <- n + as.integer(ceiling(num_value_dbl(.div(.mul(list(.sub(.int(n), ngot), .int(more))), .sub(newn, ngot)))))
          s1 <- .ns_eval(e, x, .int(ndo), logx, cdir)
          while (.is_true(Lt(.so_getn(s1), n))) { s1 <- .ns_eval(e, x, .int(ndo), logx, cdir); ndo <- ndo + 1L }
          done <- TRUE
          break
        }
      }
      if (!done) stop(.sympy_condition("ValueError", paste0("Could not calculate ", n, " terms for ", .str(e))))
      s1 <- .add(list(s1, .so_new(.pow(x, .int(n)), list(x))))
    }
    o <- .so_getO(s1)
    s1 <- .so_removeO(s1)
  } else if (.has_type(s1, "Order")) return(s1)
  else {
    o <- .so_new(.pow(x, .int(n)), list(x))
    s1done <- doit(s1)
    r <- tryCatch(identical(.so_removeO(.add(list(s1done, o)))$key, s1done$key), NotImplementedError = function(err) NA)
    if (is.na(r)) return(s1)
    if (r) o <- NULL
  }
  out <- tryCatch(collect(s1, x), NotImplementedError = function(err) s1)
  if (is.null(o)) out else .add(list(out, o))
}

#' Series expansion
#'
#' `series(e, x, x0, n)` is SymPy's `Expr.series`: the expansion of `e`
#' about `x = x0` up to (but excluding) `(x - x0)^n`, with an [O()] term
#' (none when the expansion is exact). `dir = "-"` expands from the left;
#' `x0` may be `oo`. With `n = NULL` it returns SymPy's lazy series, a
#' generator of the terms one at a time: [take()] its first terms.
#' `nseries()` is `Expr.nseries` and `leadterm()` returns
#' the coefficient and exponent of the leading term at 0;
#' `as_leading_term()` the leading term itself.
#'
#' @param e An expression.
#' @param x The variable (default: the only free symbol).
#' @param x0 The expansion point.
#' @param n Number of terms (order of truncation), or `NULL` for the lazy
#'   series.
#' @param dir `"+"` or `"-"`.
#' @param logx Optional symbol standing for `log(x)`.
#' @param cdir The direction as a complex number (overrides `dir`).
#' @return An expression (a generator for `n = NULL`); for `leadterm()` a
#'   list `(coeff, exponent)`.
#' @examples
#' x <- Symbol("x")
#' series(sin(x), x, 0, 8)
#' series(exp(x) / x, x, 0, 3)
#' series(log(1 + x), x)
#' take(series(cos(x), x, 1, n = NULL), 3)
#' @export
series <- function(e, x = NULL, x0 = 0, n = 6, dir = "+", logx = NULL, cdir = 0) {
  e <- as_sym(e)
  if (is.null(x)) {
    fs <- free_symbols(e)
    if (!length(fs)) return(e)
    if (length(fs) > 1L) stop(.sympy_condition("ValueError", "x must be given for multivariate functions."))
    x <- fs[[1L]]
  }
  .ns_series_driver(e, as_sym(x), x0, n, dir, logx, cdir)
}

#' @rdname series
#' @export
nseries <- function(e, x = NULL, x0 = 0, n = 6, dir = "+", logx = NULL, cdir = 0) {
  e <- as_sym(e)
  if (!is.null(x) && !(as_sym(x)$key %in% names(.free_sym_list(e)))) return(e)
  if (is.null(x) || !is_zero_num(as_sym(x0)) || dir != "+") return(series(e, x, x0, n, dir, cdir = cdir))
  .ns_eval(e, as_sym(x), .int(n), logx, cdir)
}

#' @rdname series
#' @export
leadterm <- function(e, x, logx = NULL, cdir = 0) {
  r <- .so_leadterm(as_sym(e), as_sym(x), logx, cdir)
  list(coeff = r[[1L]], exponent = r[[2L]])
}

#' @rdname series
#' @param bound,hir Arguments of SymPy's `aseries()`: the recursion bound for
#'   exponential terms and whether to return the hierarchical form.
#' @export
aseries <- function(e, x, n = 6, bound = 0, hir = FALSE) {
  # Expr.aseries: the asymptotic expansion at x = oo through the mrv set
  e <- as_sym(e); x <- as_sym(x)
  if (is.na(ask_fact(x, "positive")) && is.na(ask_fact(x, "negative"))) {
    xpos <- Dummy("x", positive = TRUE)
    return(subs(aseries(subs(e, x, xpos), xpos, n, bound, hir), xpos, x))
  }
  r <- tryCatch(.mrv(e, x), PoleError = function(err) NULL)
  if (is.null(r)) return(e)
  om <- r$s; exps <- r$e
  if (x$key %in% .ss_keys(om)) {
    s <- subs(aseries(subs(e, x, exp(x)), x, n, bound, hir), x, log(x))
    if (!is.null(.so_getO(s))) return(.add(list(s, .so_new(.pow(x, .int(-n)), list(x), list(S_Infinity)))))
    return(s)
  }
  k <- Dummy("k", positive = TRUE)
  rw <- .gr_rewrite(exps, om, x, k)
  func <- rw$f; logw <- rw$logw
  if (e$key %in% .ss_keys(om)) {
    if (bound <= 0) return(e)
    s <- aseries(.gr_exp_arg(e), x, n, bound = bound)
    s <- sym_rebuild(s, lapply(s$args, .so_removeO))
    res <- tryCatch(exp(subs(as_leading_term(subs(s, x, .inv(x)), x), x, .inv(x))), PoleError = function(err) e)
    func <- .div(exp(.sub(e$args[[1L]], res$args[[1L]])), k)
    logw <- log(.inv(res))
  }
  s <- expand_mul(series(func, k, 0, n))
  if (hir) return(subs(s, k, exp(logw)))
  o <- .so_getO(s)
  rest <- .so_removeO(s)
  terms <- if (rest$type == "Add") rest$args else list(rest)
  ex <- vapply(terms, function(t) as.integer(trunc(num_value_dbl(.so_as_coeff_exponent(t, k)[[2L]]))), 0L)
  terms <- terms[order(ex)]
  s <- S_Zero; has_ord <- FALSE
  for (t in terms) {
    ce <- .so_as_coeff_exponent(t, k)
    if (.has1(ce[[1L]], x)) {
      snew <- aseries(ce[[1L]], x, n, bound = bound - 1)
      if (has_ord && !is.null(.so_getO(snew))) break
      if (!is.null(.so_getO(snew))) has_ord <- TRUE
      s <- .add(list(s, .mul(list(snew, .pow(k, ce[[2L]])))))
    } else s <- .add(list(s, t))
  }
  if (is.null(o) || has_ord) return(subs(s, k, exp(logw)))
  subs(.add(list(s, o)), k, exp(logw))
}
