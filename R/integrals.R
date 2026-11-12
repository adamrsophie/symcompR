# integrate() and Integral, after sympy/integrals/integrals.py; the
# evaluation is Integral.doit (R/integrals-doit.R) with _eval_integral
# (R/integrals-pipeline.R).

#' @noRd
integrate <- function(f, ..., meijerg = NULL, conds = "piecewise", risch = NULL, heurisch = NULL, manual = NULL) {
  if (is.function(f)) return(stats::integrate(f, ...))
  if (is.character(f) && length(f) == 1L && ...length() >= 2L && is.numeric(..1) && is.numeric(..2))
    return(stats::integrate(f, ...))
  if (inherits(f, "symmatrix")) {
    a <- list(...)
    return(.mat_map(f, function(z) do.call(integrate, c(list(z), a, list(meijerg = meijerg, conds = conds, risch = risch, heurisch = heurisch, manual = manual)))))
  }
  f <- as_sym(f)
  # integrate() of a vector keeps the integrals of its components (Integral._eval_Integral)
  if (.is_bd(f)) return(.bd_integral(f, ...))
  hints <- list(meijerg = meijerg, risch = risch, manual = manual, heurisch = heurisch, conds = conds)
  I <- .integral_new(f, list(...))
  if (!is_sym(I)) return(I)
  if (I$type == "Integral") return(.integral_doit(I, hints, deep = FALSE))
  if (I$type == "Equality") return(sym_rebuild(I, lapply(I$args, function(z) if (z$type == "Integral") .integral_doit(z, hints, deep = FALSE) else z)))
  I
}

# Integral(f, *limits): AddWithLimits._common_new (denesting, folding the
# Piecewise parts that contain an integration variable)
.integral_new <- function(f, specs) {
  f <- as_sym(f)
  if (f$type == "Equality") return(Eq(.integral_new(f$args[[1L]], specs), .integral_new(f$args[[2L]], specs)))
  if (is_nan(f)) return(S_NaN)
  if (length(specs)) lims <- lapply(specs, .lim_tuple)
  else {
    fs <- free_symbols(f)
    if (length(fs) != 1L) stop("specify dummy variables for ", .str(f), call. = FALSE)
    lims <- list(Tuple(fs[[1L]]))
  }
  while (f$type == "Integral") { lims <- c(f$args[-1L], lims); f <- f$args[[1L]] }
  vars <- lapply(lims, function(l) l$args[[1L]])
  pws <- .uniq_syms(Filter(function(z) z$type == "Piecewise", .all_subexprs(f)))
  reps <- list()
  for (p in pws) if (!any(vapply(vars, function(v) .has1(p, v), NA))) reps[[length(reps) + 1L]] <- list(p, Dummy())
  if (length(reps)) f <- xreplace(f, reps)
  f <- piecewise_fold(f)
  if (length(reps)) f <- xreplace(f, lapply(reps, function(r) list(r[[2L]], r[[1L]])))
  new_sym("Integral", c(list(f), lims))
}

#' @rdname sym_integrate
#' @export
Integral <- function(f, ...) {
  f <- as_sym(f)
  if (.is_bd(f)) return(.bd_integral(f, ...))
  specs <- list(...)
  # Integral(f) of a single free symbol integrates with respect to it
  if (!length(specs) && length(free_symbols(f)) == 1L) specs <- free_symbols(f)
  lims <- lapply(specs, .lim_tuple)
  new_sym("Integral", c(list(f), lims))
}

.lim_tuple <- function(sp) {
  if (is.list(sp) && !inherits(sp, "sym")) return(do.call(Tuple, lapply(sp, as_sym)))
  sp <- as_sym(sp)
  if (sp$type == "Tuple") return(sp)
  Tuple(sp)
}

register_bound_vars("Integral", function(e) {
  vapply(e$args[-1L], function(t) if (length(t$args) >= 2L) t$args[[1L]]$key else "", "")
})
register_str_printer("Integral", function(e) {
  lims <- vapply(e$args[-1L], function(t) if (length(t$args) == 1L) .str(t$args[[1L]]) else .str(t), "")
  paste0("Integral(", .str(e$args[[1L]]), ", ", paste(lims, collapse = ", "), ")")
})
register_doit("Integral", function(e) .integral_doit(e))
register_diff_hook("Integral", function(e, x) {
  f <- e$args[[1L]]
  lims <- e$args[-1L]
  if (length(lims) == 1L && length(lims[[1L]]$args) == 3L) {
    v <- lims[[1L]]$args[[1L]]; a <- lims[[1L]]$args[[2L]]; b <- lims[[1L]]$args[[3L]]
    if (identical(v$key, x$key)) {
      # the integration variable is a dummy: rename it before Leibniz' rule
      d <- Dummy(v$name)
      f <- xreplace(f, list(list(v, d)))
      v <- d
      lims[[1L]] <- Tuple(d, a, b)
    }
    out <- list(.mul(list(subs(f, v, b), .diff1(b, x))), .neg(.mul(list(subs(f, v, a), .diff1(a, x)))))
    if (.has1(f, x)) out[[3L]] <- new_sym("Integral", list(.diff1(f, x), lims[[1L]]))
    return(.add(out))
  }
  if (length(lims) == 1L && length(lims[[1L]]$args) == 1L && identical(lims[[1L]]$args[[1L]]$key, x$key)) return(f)
  .derivative_node(e, list(list(x, 1L)))
})

# the antiderivative by SymPy's _eval_integral, NULL when none is found
# ponytail: a nesting guard of 20 integrate() calls stands in for Python's
# recursion limit; raise it if a legitimate integral nests deeper
.antideriv_sympy <- function(f, x, hints = list()) {
  depth <- .mi_state$int_depth %||% 0L
  if (depth > 20L) return(NULL)
  .mi_state$int_depth <- depth + 1L
  on.exit(.mi_state$int_depth <- depth)
  do.call(.sympy_eval_integral, c(list(f, x), hints))
}

# The real set where a boolean condition in x holds.
.cond_set <- function(c, x) {
  if (c$type == "BooleanAtom") return(if (.is_true(c)) Reals else S_EmptySet)
  switch(c$type,
    And = Reduce(Intersection, lapply(c$args, .cond_set, x = x)),
    Or = do.call(Union, lapply(c$args, .cond_set, x = x)),
    Not = Complement(Reals, .cond_set(c$args[[1L]], x)),
    solve_univariate_inequality(c, x, relational = FALSE))
}

.real_singularities <- function(f, x, real = TRUE) {
  out <- list()
  dens <- list()
  walk <- function(z) {
    if (z$type == "Pow" && .is_neg_exp(z$args[[2L]]) && .has1(z$args[[1L]], x)) dens[[length(dens) + 1L]] <<- z$args[[1L]]
    if (z$type == "Function" && z$name == "log" && .has1(z$args[[1L]], x)) dens[[length(dens) + 1L]] <<- z$args[[1L]]
    for (a in z$args) walk(a)
  }
  walk(f)
  for (d in dens) {
    if (!is_polynomial(d, x)) next
    rs <- tryCatch(roots(d, x), error = function(err) list())
    for (r in rs) if (!real || isTRUE(ask_fact(r, "real"))) out[[length(out) + 1L]] <- r
  }
  out[!duplicated(vapply(out, .k, ""))]
}

# ---------------------------------------------------------------------------
# Antiderivatives

# Linear argument a*x + b of a function, or NULL.
.linear_arg <- function(g, x) {
  if (!is_polynomial(g, x)) return(NULL)
  cl <- tryCatch(.coeff_list(g, x), error = function(err) NULL)
  if (is.null(cl) || length(cl) != 2L) return(NULL)
  list(a = cl[[1L]], b = cl[[2L]])
}

#' Line integrals
#'
#' `line_integrate(f, curve, vars)` integrates `f` along a parametric curve
#' given as `list(param, a, b, x(param), y(param), ...)`.
#'
#' @param f An expression.
#' @param curve A list `list(param, a, b, coord1, coord2, ...)`.
#' @param vars The coordinate variables.
#' @return An expression.
#' @examples
#' t <- Symbol("t")
#' line_integrate(1, list(t, 0, 1, t, t), list(Symbol("x"), Symbol("y")))
#' @export
line_integrate <- function(f, curve, vars) {
  f <- as_sym(f)
  p <- as_sym(curve[[1L]])
  a <- curve[[2L]]; b <- curve[[3L]]
  coords <- lapply(curve[-(1:3)], as_sym)
  vars <- lapply(vars, as_sym)
  ds <- .sqrt(.add(lapply(coords, function(cc) .sq(.diff1(cc, p)))))
  g <- subs(f, lapply(seq_along(vars), function(i) list(vars[[i]], coords[[i]])))
  integrate(.mul(list(g, ds)), list(p, a, b))
}

