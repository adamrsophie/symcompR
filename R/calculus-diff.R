# Differentiation, after sympy/core/function.py (Derivative).

#' Tuples
#'
#' A `Tuple` is an ordered container of expressions, used for derivative
#' variable counts, integration limits and solution sets.
#'
#' @param ... Expressions.
#' @return A `sym` Tuple.
#' @export
Tuple <- function(...) new_sym("Tuple", lapply(list(...), as_sym))

register_str_printer("Tuple", function(e) {
  inner <- paste(vapply(e$args, .str, ""), collapse = ", ")
  if (length(e$args) == 1L) paste0("(", inner, ",)") else paste0("(", inner, ")")
})

# A derivative count: an R integer, or a sym when it is symbolic.
.dcount <- function(c) {
  if (!is_sym(c)) return(as.integer(c))
  if (c$type == "Integer") as.integer(as.numeric(c$value)) else c
}
.dc_sym <- function(c) if (is_sym(c)) c else .int(c)
.dc_add <- function(a, b) if (is_sym(a) || is_sym(b)) .dcount(.add(list(.dc_sym(a), .dc_sym(b)))) else a + b

.diff_spec <- function(vars) {
  # Normalize diff() variable arguments into list(list(var, count), ...).
  out <- list()
  i <- 1L
  while (i <= length(vars)) {
    v <- vars[[i]]
    if (is.list(v) && !inherits(v, "sym")) {
      out[[length(out) + 1L]] <- list(as_sym(v[[1L]]), .dcount(as_sym(v[[2L]])))
      i <- i + 1L
      next
    }
    if (inherits(v, "sym") && v$type == "Tuple") {
      out[[length(out) + 1L]] <- list(v$args[[1L]], .dcount(v$args[[2L]]))
      i <- i + 1L
      next
    }
    v <- as_sym(v)
    cnt <- 1L
    if (i < length(vars)) {
      nx <- vars[[i + 1L]]
      if (is.numeric(nx) || (inherits(nx, "sym") && nx$type == "Integer")) {
        cnt <- as.integer(as.numeric(nx))
        i <- i + 1L
      }
    }
    out[[length(out) + 1L]] <- list(v, cnt)
    i <- i + 1L
  }
  out
}

#' Differentiate an expression
#'
#' `diff(f, x)` differentiates with respect to `x`; `diff(f, x, 3)` gives the
#' third derivative and `diff(f, x, y)` a mixed partial derivative. Without
#' variables, `f` must have exactly one free symbol. Derivatives of unknown
#' functions stay as [Derivative()] objects. This is a method for the base R
#' generic [base::diff()].
#'
#' @param x A `sym` expression.
#' @param ... Variables, optionally each followed by a count.
#' @param evaluate If `FALSE`, return an unevaluated [Derivative()].
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' diff(sin(x) * exp(x), x)
#' diff(x^3 * y^2, x, 2, y)
#' f <- Function("f")
#' diff(f(x)^2, x)
#' @export
diff.sym <- function(x, ..., evaluate = TRUE) {
  f <- x
  vars <- list(...)
  if (!length(vars)) {
    fs <- free_symbols(f)
    if (length(fs) != 1L) stop("specify differentiation variables")
    vars <- list(fs[[1L]])
  }
  spec <- .diff_spec(vars)
  if (!evaluate) return(.derivative_node(f, spec, sort = FALSE))
  .diff_spec_eval(f, spec)
}

# Derivative(f, *spec, evaluate=True): differentiate pair by pair; a pair
# with no closed form leaves it and the pairs after it unevaluated.
.diff_spec_eval <- function(f, spec) {
  for (i in seq_along(spec)) {
    r <- .diff_count(f, spec[[i]][[1L]], spec[[i]][[2L]])
    if (is.null(r)) return(.derivative_node(f, spec[i:length(spec)]))
    f <- r
    if (is_zero_num(f)) return(f)
  }
  f
}

# _eval_derivative_n_times for one variable: a loop for an integer count;
# for a symbolic count zero when the count is positive and v is absent, the
# AtomicExpr Piecewise for an atom, else NULL.
.diff_count <- function(f, v, cnt) {
  if (!is_sym(cnt)) {
    for (k in seq_len(cnt)) {
      f <- .diff1(f, v)
      if (is_zero_num(f)) break
    }
    return(f)
  }
  if (is_symbol_like(v) && isTRUE(ask_fact(cnt, "positive")) && is.null(.free_sym_list(f)[[v$key]])) return(S_Zero)
  if (length(f$args)) return(NULL)
  if (identical(f$key, v$key))
    return(Piecewise(list(f, Eq(cnt, 0)), list(S_One, Eq(cnt, 1)), list(S_Zero, S_true)))
  Piecewise(list(f, Eq(cnt, 0)), list(S_Zero, S_true))
}

.derivative_node <- function(f, spec, sort = TRUE) {
  if (!length(spec)) return(f)
  if (f$type == "Derivative") {
    inner <- f$args[[1L]]
    old <- lapply(f$args[-1L], function(t) list(t$args[[1L]], .dcount(t$args[[2L]])))
    spec <- c(old, spec)
    f <- inner
  }
  if (sort) {
    keys <- vapply(spec, function(s) s[[1L]]$key, "")
    u <- unique(keys)
    merged <- lapply(u, function(k) {
      w <- which(keys == k)
      list(spec[[w[1L]]][[1L]], Reduce(.dc_add, lapply(spec[w], `[[`, 2L)))
    })
    if (all(vapply(merged, function(s) is_symbol_like(s[[1L]]), TRUE))) {
      o <- order(vapply(merged, function(s) s[[1L]]$key, ""), method = "radix")
      merged <- merged[o]
    }
  } else {
    # Derivative.__new__: zero counts drop and only contiguous identical
    # variables merge, in the order given
    merged <- list()
    for (sp in spec) {
      if (identical(sp[[2L]], 0L)) next
      m <- length(merged)
      if (m && identical(merged[[m]][[1L]]$key, sp[[1L]]$key)) {
        c <- .dc_add(merged[[m]][[2L]], sp[[2L]])
        if (identical(c, 0L)) merged[[m]] <- NULL else merged[[m]][[2L]] <- c
      } else merged[[m + 1L]] <- sp
    }
    if (!length(merged)) return(f)
  }
  new_sym("Derivative", c(list(f), lapply(merged, function(s)
    new_sym("Tuple", list(s[[1L]], .dc_sym(s[[2L]]))))))
}

#' Unevaluated derivative
#'
#' `Derivative(f, x, ...)` represents a derivative without computing it.
#' Use [doit()] to evaluate.
#'
#' @param f Expression.
#' @param ... Variables and counts, as for [diff()].
#' @return A `sym` Derivative.
#' @export
Derivative <- function(f, ...) {
  f <- as_sym(f)
  vars <- list(...)
  if (!length(vars)) {
    # Derivative(expr): the single free symbol is the variable
    fs <- free_symbols(f)
    if (length(fs) != 1L) {
      if (.expr_is_number(f)) return(S_Zero)
      if (!length(fs)) stop("Since there are no variables in the expression, the variable(s) of differentiation must be supplied to differentiate ", .str(f), call. = FALSE)
      stop("Since there is more than one variable in the expression, the variable(s) of differentiation must be supplied to differentiate ", .str(f), call. = FALSE)
    }
    vars <- list(fs[[1L]])
  }
  .derivative_node(f, .diff_spec(vars), sort = FALSE)
}

register_str_printer("Derivative", function(e) {
  parts <- vapply(e$args[-1L], function(t) {
    n <- t$args[[2L]]
    if (is_one_num(n)) .str(t$args[[1L]]) else paste0("(", .str(t$args[[1L]]), ", ", .str(n), ")")
  }, "")
  paste0("Derivative(", .str(e$args[[1L]]), ", ", paste(parts, collapse = ", "), ")")
})

register_bound_vars("Derivative", function(e) character())

.diff1 <- function(e, x) {
  if (identical(e$key, x$key)) return(S_One)
  t <- e$type
  if (!length(e$args)) return(S_Zero)
  if (is_symbol_like(x) && is.null(.free_sym_list(e)[[x$key]])) return(S_Zero)
  if (!is_symbol_like(x) && !.contains_key(e, x$key)) return(S_Zero)
  switch(t,
         Add = .add(lapply(e$args, .diff1, x)),
         Mul = {
           fs <- e$args
           terms <- list()
           for (i in seq_along(fs)) {
             d <- .diff1(fs[[i]], x)
             if (is_zero_num(d)) next
             fi <- fs
             fi[[i]] <- d
             terms[[length(terms) + 1L]] <- .mul(fi)
           }
           .add(terms)
         },
         Pow = {
           b <- e$args[[1L]]
           p <- e$args[[2L]]
           db <- .diff1(b, x)
           dp <- .diff1(p, x)
           if (is_zero_num(dp)) {
             # Pow._eval_derivative: self*(dexp*log(b) + db*p/b)
             .mul(list(e, .div(.mul(list(db, p)), b)))
           } else if (is_zero_num(db)) {
             .mul(list(e, fcall("log", list(b)), dp))
           } else {
             .mul(list(e, .add(list(.mul(list(dp, fcall("log", list(b)))),
                                   .div(.mul(list(p, db)), b)))))
           }
         },
         Function = .diff_function(e, x),
         Derivative = {
           inner <- e$args[[1L]]
           if (!is_symbol_like(x)) {
             # d/dv Derivative(f, t) = Derivative(df/dv, t), as in SymPy.
             r <- .diff1(inner, x)
             for (tv in e$args[-1L])
               for (i in seq_len(as.integer(as.numeric(tv$args[[2L]])))) r <- .diff1(r, tv$args[[1L]])
             r
           } else if (.is_undefined_application(inner)) {
             .derivative_node(e, list(list(x, 1L)))
           } else {
             ev <- doit(e)
             # an unevaluated derivative only gains one more order
             if (ev$type == "Derivative" && identical(ev$key, e$key)) .derivative_node(e, list(list(x, 1L)))
             else .diff1(ev, x)
           }
         },
         Equality = Eq(.diff1(e$args[[1L]], x), .diff1(e$args[[2L]], x)),
         {
           h <- .diff_hooks[[t]]
           if (!is.null(h)) h(e, x) else .derivative_node(e, list(list(x, 1L)))
         })
}

.is_undefined_application <- function(e) {
  e$type == "Function" && is.null(.FUNS[[e$name]])
}

# The Dummy SymPy's Function.fdiff uses for argument i equal to A: the same
# (i, A) always gives the same Dummy, so equal chain rule terms cancel.
.fdiff_dummy <- function(i, A) {
  nm <- paste0("xi_", i)
  new_sym("Dummy", name = nm, extra = list(facts = .symbol_facts(logical()), id = 0L),
          key = paste0("Dfd", i, "|", A$key, ":", nchar(nm), ":", nm))
}

# Chain rule for an undefined function, as SymPy's Function._eval_derivative
# and fdiff: differentiate with respect to an argument when it is a symbol
# found in no other argument, else through Subs(Derivative(f(xi), xi), xi, A).
.diff_undefined <- function(e, x) {
  args <- e$args
  terms <- list()
  for (i in seq_along(args)) {
    A <- args[[i]]
    da <- .diff1(A, x)
    if (is_zero_num(da)) next
    wrt <- is_symbol_like(A) || .is_undefined_application(A) || A$type == "Derivative"
    direct <- wrt && (length(args) == 1L || !is_symbol_like(A) ||
                      !any(vapply(args[-i], function(v) .has1(v, A), TRUE)))
    df <- if (direct) .derivative_node(e, list(list(A, 1L))) else {
      D <- .fdiff_dummy(i, A)
      a2 <- args
      a2[[i]] <- D
      Subs(.derivative_node(sym_rebuild(e, a2), list(list(D, 1L))), D, A)
    }
    terms[[length(terms) + 1L]] <- .mul(list(df, da))
  }
  .add(terms)
}

.diff_function <- function(e, x) {
  d <- .FUNS[[e$name]]
  args <- e$args
  if (is.null(d)) return(.diff_undefined(e, x))
  # functions with Tuple arguments differentiate themselves
  if (!is.null(d$diffx)) return(d$diffx(e, x))
  terms <- list()
  for (i in seq_along(args)) {
    da <- .diff1(args[[i]], x)
    if (is_zero_num(da)) next
    pd <- if (!is.null(d$diff)) d$diff(args, i) else NULL
    if (is.null(pd)) {
      if (is_symbol_like(args[[i]]) && identical(args[[i]]$key, x$key)) {
        pd <- .derivative_node(e, list(list(x, 1L)))
        terms[[length(terms) + 1L]] <- pd
        next
      }
      return(.derivative_node(e, list(list(x, 1L))))
    }
    terms[[length(terms) + 1L]] <- .mul(list(pd, da))
  }
  .add(terms)
}

# ---------------------------------------------------------------------------
# doit: evaluate unevaluated objects.


#' Evaluate unevaluated objects
#'
#' Evaluates [Derivative()], `Integral`, `Sum`, `Product` and `Limit`
#' objects inside an expression, like SymPy's `doit()`.
#'
#' @param e A `sym` expression.
#' @param deep Evaluate nested objects too.
#' @return A `sym` expression.
#' @export
doit <- function(e, deep = TRUE) {
  if (inherits(e, "TensExpr")) return(.tn_doit(e, deep))
  if (inherits(e, "symmatrix")) return(.mat_map(e, doit))
  e <- as_sym(e)
  if (!length(e$args)) return(e)
  # Piecewise.doit: a piece whose doit gives back this Piecewise is kept
  if (e$type == "Piecewise") {
    if (!deep) return(e)
    return(.piecewise(lapply(.pw_pairs(e), function(p) {
      ne <- doit(p[[1L]])
      list(if (identical(ne$key, e$key)) p[[1L]] else ne, doit(p[[2L]]))
    })))
  }
  if (deep) {
    # an explicit matrix result stays an argument as an ImmutableMatrix node
    args <- lapply(e$args, function(a) { r <- doit(a); if (inherits(r, "symmatrix")) .as_block(r) else r })
    if (!identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e <- sym_rebuild(e, args)
  }
  t <- e$type
  if (t == "Derivative")
    return(.diff_spec_eval(e$args[[1L]], lapply(e$args[-1L], function(tp) list(tp$args[[1L]], .dcount(tp$args[[2L]])))))
  h <- .doit_hooks[[t]]
  if (!is.null(h)) return(if ("deep" %in% names(formals(h))) h(e, deep = deep) else h(e))
  e
}

register_rebuilder("Derivative", function(e, args) {
  spec <- lapply(args[-1L], function(t) list(t$args[[1L]], .dcount(t$args[[2L]])))
  .derivative_node(args[[1L]], spec, sort = FALSE)
})

register_neval("Tuple", function(e, env) stop("cannot evaluate a Tuple numerically"))

#' Jacobian, Hessian, gradient and Wronskian
#'
#' @param f An expression or list of expressions.
#' @param vars A list of variables. For `gradient()` without `vars`, the
#'   gradient of a [CoordSys3D()] scalar or vector field (`sympy.vector`).
#' @param doit Evaluate the derivatives (coordinate system fields only).
#' @return A symbolic matrix (see [Matrix()]).
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' jacobian(list(x^2 * y, 5 * x + sin(y)), list(x, y))
#' hessian(x^3 + x * y^2, list(x, y))
#' @export
jacobian <- function(f, vars) {
  if (inherits(f, "symmatrix")) f <- f$data
  if (!is.list(f) || inherits(f, "sym")) f <- list(f)
  f <- lapply(f, as_sym)
  vars <- lapply(if (is.list(vars) && !inherits(vars, "sym")) vars else list(vars), as_sym)
  rows <- lapply(f, function(fi) lapply(vars, function(v) .diff1(fi, v)))
  Matrix(rows)
}

#' @rdname jacobian
#' @export
hessian <- function(f, vars) {
  f <- as_sym(f)
  vars <- lapply(if (is.list(vars) && !inherits(vars, "sym")) vars else list(vars), as_sym)
  rows <- lapply(vars, function(v) {
    dv <- .diff1(f, v)
    lapply(vars, function(w) .diff1(dv, w))
  })
  Matrix(rows)
}

#' @rdname jacobian
#' @export
gradient <- function(f, vars, doit = TRUE) {
  if (missing(vars)) return(.vec_gradient(f, doit = doit))
  f <- as_sym(f)
  vars <- lapply(if (is.list(vars) && !inherits(vars, "sym")) vars else list(vars), as_sym)
  Matrix(lapply(vars, function(v) list(.diff1(f, v))))
}

#' @rdname jacobian
#' @param x The variable for the Wronskian.
#' @export
wronskian <- function(f, x) {
  f <- lapply(f, as_sym)
  x <- as_sym(x)
  n <- length(f)
  rows <- list()
  cur <- f
  for (i in seq_len(n)) {
    rows[[i]] <- cur
    cur <- lapply(cur, .diff1, x)
  }
  det(Matrix(rows))
}

#' Implicit differentiation
#'
#' Returns `dy/dx` (or the `n`-th derivative) for `eq = 0` defining `y` as a
#' function of `x`, like SymPy's `idiff()`.
#'
#' @param eq An expression equal to zero.
#' @param y The dependent symbol.
#' @param x The independent symbol.
#' @param n Order of the derivative.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' idiff(x^2 + y^2 - 4, y, x)
#' @export
idiff <- function(eq, y, x, n = 1) {
  eq <- as_sym(eq)
  if (is_relational(eq)) eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  y <- as_sym(y)
  x <- as_sym(x)
  fy <- .diff1(eq, y)
  fx <- .diff1(eq, x)
  dydx <- .neg(.div(fx, fy))
  if (n == 1) return(simplify(dydx))
  cur <- dydx
  for (k in seq_len(n - 1)) {
    cur <- .add(list(.diff1(cur, x), .mul(list(.diff1(cur, y), dydx))))
  }
  simplify(cur)
}
