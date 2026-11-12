# collect and rcollect, after sympy/simplify/radsimp.py.

#' Collect terms
#'
#' `collect(e, syms)` groups the additive terms of `e` by powers of the
#' symbols (or expressions) in `syms`, as SymPy's `collect()`; it does not
#' expand `e` first. With `evaluate = FALSE` it returns a named list mapping
#' the collected powers (by their printed form) to their coefficients.
#' `rcollect(e, ...)` collects recursively inside the arguments.
#'
#' @param e An expression.
#' @param syms A symbol, expression or list of them.
#' @param func Optional function applied to each coefficient.
#' @param evaluate Return the collected expression (default) or a list.
#' @param exact Collect only exact powers of the patterns.
#' @param distribute_order_term Distribute an order term over the coefficients.
#' @param ... Symbols for `rcollect()`.
#' @return An expression or a named list.
#' @examples
#' x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b")
#' collect(a*x^2 + b*x^2 + a*x - b*x + 1, x)
#' @export
collect <- function(e, syms, func = NULL, evaluate = TRUE, exact = FALSE, distribute_order_term = TRUE) {
  expr <- as_sym(e)
  syms <- if (is_sym(syms)) list(syms) else lapply(syms, as_sym)
  issym <- function(z) z$type %in% c("Symbol", "Dummy") || .neg(z)$type %in% c("Symbol", "Dummy") || .has_type(z, "Wild")
  nonsyms <- Filter(Negate(issym), syms)
  if (length(nonsyms)) {
    reps <- lapply(nonsyms, function(z) list(z, Dummy()))
    rk <- vapply(reps, function(r) r[[1L]]$key, "")
    syms <- lapply(syms, function(z) { i <- match(z$key, rk); if (is.na(i)) z else reps[[i]][[2L]] })
    rv <- collect(subs(expr, reps), syms, func = func, evaluate = evaluate, exact = exact, distribute_order_term = distribute_order_term)
    urep <- lapply(reps, function(r) list(r[[2L]], r[[1L]]))
    if (is_sym(rv)) return(xreplace(rv, urep))
    out <- lapply(rv, function(v) xreplace(v, urep))
    names(out) <- vapply(attr(rv, "keys"), function(k) .str(xreplace(k, urep)), "")
    attr(out, "keys") <- lapply(attr(rv, "keys"), function(k) xreplace(k, urep))
    return(out)
  }
  if (is.null(exact)) {
    ss <- list()
    for (i in .mi_make_args_add(expr)) {
      if (!any(vapply(syms, function(s) .has1(i, s), NA)) || any(vapply(syms, function(s) identical(s$key, i$key), NA))) next
      if (i$type != "Mul") ss[[length(ss) + 1L]] <- i
      else {
        dep <- Filter(function(f) any(vapply(syms, function(s) .has1(f, s), NA)), .sympy_args(i))
        g <- .mul(dep)
        if (!any(vapply(syms, function(s) identical(s$key, g$key), NA))) ss[[length(ss) + 1L]] <- g
      }
    }
    ss <- .uniq_syms(ss)
    simple <- all(vapply(ss, function(i) i$type == "Pow" && any(vapply(syms, function(s) identical(s$key, i$args[[1L]]$key), NA)), NA))
    syms <- c(syms, .mi_sort_ordered(ss))
    if (!simple) return(collect(expr, syms, func = func, evaluate = evaluate, exact = FALSE, distribute_order_term = distribute_order_term))
  }
  make_expression <- function(terms) {
    product <- list()
    for (tm in terms) {
      term <- tm[[1L]]; rat <- tm[[2L]]; sym <- tm[[3L]]; deriv <- tm[[4L]]
      if (!is.null(deriv)) for (k in seq_len(.int_val(deriv[[2L]]))) term <- Derivative(term, deriv[[1L]])
      product[[length(product) + 1L]] <- if (is.null(sym)) (if (is_one_num(rat)) term else .pow(term, rat)) else .pow(term, .mul(list(rat, sym)))
    }
    .mul(product)
  }
  parse_derivative <- function(d) {
    dvars <- unlist(lapply(d$args[-1L], function(t) rep(list(t$args[[1L]]), .int_val(t$args[[2L]]))), recursive = FALSE)
    ex <- d$args[[1L]]; sym <- dvars[[1L]]; order <- 1L
    for (s in dvars[-1L]) if (identical(s$key, sym$key)) order <- order + 1L else .notimpl_error("Improve MV Derivative support in collect")
    while (ex$type == "Derivative") {
      ev <- unlist(lapply(ex$args[-1L], function(t) rep(list(t$args[[1L]]), .int_val(t$args[[2L]]))), recursive = FALSE)
      s0 <- ev[[1L]]
      if (any(vapply(ev, function(s) !identical(s$key, s0$key), NA))) .notimpl_error("Improve MV Derivative support in collect")
      if (identical(s0$key, sym$key)) { order <- order + length(ev); ex <- ex$args[[1L]] } else break
    }
    list(ex, list(sym, .int(order)))
  }
  parse_term <- function(ex) {
    rat <- S_One; symx <- NULL; sexpr <- ex; deriv <- NULL
    if (ex$type == "Pow") {
      if (ex$args[[1L]]$type == "Derivative") { pd <- parse_derivative(ex$args[[1L]]); sexpr <- pd[[1L]]; deriv <- pd[[2L]] }
      else sexpr <- ex$args[[1L]]
      p <- ex$args[[2L]]
      if (is_number_type(p)) rat <- p
      else {
        cm <- as_coeff_mul_split(p)
        if (is_number_type(cm[[1L]]) && !is_one_num(cm[[1L]])) { rat <- cm[[1L]]; symx <- cm[[2L]] }
        else symx <- p
      }
    } else if (ex$type == "Function" && ex$name == "exp") {
      a <- ex$args[[1L]]
      if (a$type %in% c("Integer", "Rational")) { sexpr <- S_E; rat <- a }
      else if (a$type == "Mul") {
        cm <- .mi_as_coeff_Mul(a)
        if (cm[[1L]]$type %in% c("Integer", "Rational")) { sexpr <- fcall("exp", list(cm[[2L]])); rat <- cm[[1L]] }
        else sexpr <- fcall("exp", list(a))
      }
    } else if (ex$type == "Derivative") { pd <- parse_derivative(ex); sexpr <- pd[[1L]]; deriv <- pd[[2L]] }
    list(sexpr, rat, symx, deriv)
  }
  matches <- function(a, b) !is.null(tryCatch(sym_match(a, b), error = function(err) NULL))
  parse_expression <- function(terms, pattern) {
    pattern <- .mi_make_args_mul(pattern)
    if (length(terms) < length(pattern)) return(NULL)
    pattern <- lapply(pattern, parse_term)
    elems <- list(); has_deriv <- FALSE
    for (pe in pattern) {
      if (is_number_type(pe[[1L]]) && is_one_num(pe[[2L]]) && is.null(pe[[3L]])) next
      found <- FALSE
      for (j in seq_along(terms)) {
        tj <- terms[[j]]
        if (is.null(tj)) next
        if (!is.null(tj[[4L]])) has_deriv <- TRUE
        symok <- (is.null(tj[[3L]]) && is.null(pe[[3L]])) || (!is.null(tj[[3L]]) && !is.null(pe[[3L]]) && (identical(tj[[3L]]$key, pe[[3L]]$key) || matches(tj[[3L]], pe[[3L]])))
        if (matches(tj[[1L]], pe[[1L]]) && symok) {
          if (!isFALSE(exact) && (!identical(pe[[2L]]$key, tj[[2L]]$key) || !identical(.mr_key(pe[[4L]]), .mr_key(tj[[4L]])))) next
          elems[[length(elems) + 1L]] <- tj
          terms[j] <- list(NULL)
          found <- TRUE
          break
        }
      }
      if (!found) return(NULL)
    }
    list(Filter(Negate(is.null), terms), elems, has_deriv)
  }
  if (evaluate) {
    if (expr$type == "Add") {
      o <- getO(expr)
      parts <- Filter(function(a) is.null(o) || !identical(a$key, o$key), expr$args)
      expr <- .add(c(lapply(parts, function(a) collect(a, syms, func, TRUE, exact, distribute_order_term)), if (is.null(o)) list() else list(o)))
    } else if (expr$type == "Mul") {
      return(.mul(lapply(expr$args, function(tm) collect(tm, syms, func, TRUE, exact, distribute_order_term))))
    } else if (expr$type == "Pow") {
      return(.pow(collect(expr$args[[1L]], syms, func, TRUE, exact, distribute_order_term), expr$args[[2L]]))
    }
  }
  syms <- lapply(syms, function(i) expand_power_base(i, deep = FALSE))
  order_term <- NULL
  if (distribute_order_term) {
    order_term <- getO(expr)
    if (!is.null(order_term)) {
      if (any(vapply(syms, function(s) .has1(order_term, s), NA))) order_term <- NULL
      else expr <- removeO(expr)
    }
  }
  summa <- lapply(.mi_make_args_add(expr), function(i) expand_power_base(i, deep = FALSE))
  ckeys <- list(); cvals <- list(); disliked <- S_Zero
  for (product in summa) {
    fs <- .mi_make_args_mul(product)
    cm <- vapply(fs, .is_comm, NA)
    args <- c(.mi_sort_ordered(fs[cm]), fs[!cm])
    terms <- lapply(args, parse_term)
    hit <- FALSE
    for (symbol in syms) {
      if (symbol$type == "Derivative") terms <- rev(terms)
      result <- parse_expression(terms, symbol)
      if (!is.null(result)) {
        if (!.is_comm(symbol)) stop("Can not collect noncommutative symbol", call. = FALSE)
        terms <- result[[1L]]; elems <- result[[2L]]
        index <- if (!result[[3L]]) .mul(lapply(elems, function(el) .pow(el[[1L]], if (is.null(el[[3L]])) el[[2L]] else .mul(list(el[[2L]], el[[3L]])))))
                 else make_expression(elems)
        tv <- expand_power_base(make_expression(terms), deep = FALSE)
        index <- expand_power_base(index, deep = FALSE)
        k <- index$key
        if (is.null(ckeys[[k]])) { ckeys[[k]] <- index; cvals[[k]] <- list() }
        cvals[[k]][[length(cvals[[k]]) + 1L]] <- tv
        hit <- TRUE
        break
      }
    }
    if (!hit) disliked <- .add(list(disliked, product))
  }
  collected <- lapply(cvals, .add)
  if (!is_zero_num(disliked)) {
    k1 <- S_One$key
    if (is.null(ckeys[[k1]])) ckeys[[k1]] <- S_One
    collected[[k1]] <- disliked
  }
  if (!is.null(order_term)) collected <- lapply(collected, function(v) .add(list(v, order_term)))
  if (!is.null(func)) collected <- lapply(collected, func)
  if (evaluate) return(.add(lapply(names(collected), function(k) .mul(list(ckeys[[k]], collected[[k]])))))
  out <- unname(collected)
  names(out) <- vapply(names(collected), function(k) .str(ckeys[[k]]), "")
  attr(out, "keys") <- unname(ckeys[names(collected)])
  out
}

#' @rdname collect
#' @export
rcollect <- function(e, ...) {
  e <- as_sym(e)
  vars <- lapply(list(...), as_sym)
  if (!length(e$args) || !any(vapply(vars, function(v) .has1(e, v), NA))) return(e)
  e <- sym_rebuild(e, lapply(e$args, function(a) do.call(rcollect, c(list(a), vars))))
  if (e$type == "Add") return(collect(e, vars))
  e
}
