# Approximations for numerical evaluation, after
# sympy/codegen/approximations.py.

# The bounds of SumApprox and SeriesApprox: a list of pairs
# list(expr, list(lo, hi)) (the R form of a Python dict)
.cg_bounds <- function(bounds) {
  if (inherits(bounds, "cg_dict")) return(bounds)
  .cgd(lapply(bounds, function(pr) list(as_sym(pr[[1L]]), lapply(pr[[2L]], as_sym))))
}
.cg_dbl <- function(z) as.numeric(N(as_sym(z)))

#' Approximations of sums and of functions
#'
#' `SumApprox(bounds, reltol)` is an optimization (see [sym_optimize()])
#' dropping the terms of a sum that are small next to the others, given
#' bounds on the symbols, as SymPy's `SumApprox`. `SeriesApprox(bounds,
#' reltol)` replaces functions of one bounded symbol by their truncated
#' series when that is accurate enough, as SymPy's `SeriesApprox`.
#'
#' @param bounds A list of pairs `list(symbol, list(lo, hi))`.
#' @param reltol The relative tolerance.
#' @param max_order The largest order of the series.
#' @param n_point_checks The number of points (even) where the series is
#'   checked.
#' @param cost_function,priority As in [ReplaceOptim()].
#' @return An optimization (an R function of an expression).
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' sa <- SumApprox(list(list(x, list(0, 0.001)), list(y, list(100, 1000))), reltol = 0.01)
#' sym_optimize(x^2 + y^2 + 1, list(sa))
#' @export
SumApprox <- function(bounds, reltol, cost_function = NULL, priority = 1) {
  st <- new.env(parent = emptyenv())
  st$bounds <- .cg_bounds(bounds)
  rt <- .cg_dbl(reltol)
  inb <- function(t) .cgd_has(st$bounds, t)
  value <- function(add) {
    for (term in add$args) {
      fs <- free_symbols(term)
      if (.expr_is_number(term) || inb(term) || length(fs) != 1L) next
      fs <- fs[[1L]]
      if (!inb(fs)) next
      b <- .cgd_get(st$bounds, fs)
      intrvl <- Interval(b[[1L]], b[[2L]])
      at <- function(v) subs(term, fs, v)
      if (isTRUE(is_increasing(term, intrvl, fs))) {
        st$bounds <- .cgd_merge(st$bounds, .cgd(list(list(term, list(at(b[[1L]]), at(b[[2L]]))))))
      } else if (isTRUE(is_decreasing(term, intrvl, fs))) {
        st$bounds <- .cgd_merge(st$bounds, .cgd(list(list(term, list(at(b[[2L]]), at(b[[1L]]))))))
      } else return(add)
    }
    if (!all(vapply(add$args, function(t) .expr_is_number(t) || inb(t), TRUE))) return(add)
    bnds <- lapply(add$args, function(t) if (.expr_is_number(t)) list(t, t) else .cgd_get(st$bounds, t))
    largest <- 0
    for (b in bnds) {
      lo <- .cg_dbl(b[[1L]]); hi <- .cg_dbl(b[[2L]])
      if (lo <= 0 && 0 <= hi) next
      largest <- max(largest, min(abs(lo), abs(hi)))
    }
    keep <- vapply(bnds, function(b) max(abs(.cg_dbl(b[[1L]])), abs(.cg_dbl(b[[2L]]))) >= largest * rt, TRUE)
    .add(add$args[keep])
  }
  .cg_optim_obj(function(expr) .cg_replace(sym_factor(as_sym(expr)), function(e) e$type == "Add", value),
                cost_function, priority, cls = "SumApprox", bounds = st, reltol = reltol)
}

#' @rdname SumApprox
#' @export
SeriesApprox <- function(bounds, reltol, max_order = 4, n_point_checks = 4, cost_function = NULL, priority = 1) {
  bounds <- .cg_bounds(bounds)
  if (n_point_checks %% 2 == 1)
    .cg_err("ValueError", "Checking the solution at expansion point is not helpful")
  rt <- .cg_dbl(reltol)
  prec <- ceiling(-log10(rt))
  value <- function(fexpr) {
    fs <- free_symbols(fexpr)
    if (length(fs) != 1L) return(fexpr)
    symb <- fs[[1L]]
    if (!.cgd_has(bounds, symb)) return(fexpr)
    b <- .cgd_get(bounds, symb)
    lo <- b[[1L]]; hi <- b[[2L]]
    x0 <- (lo + hi) / 2L
    cheapest <- NULL
    for (n in seq(max_order + 1L, 1L)) {
      fseri <- removeO(series(fexpr, symb, x0 = x0, n = n))
      n_ok <- TRUE
      for (idx in seq_len(n_point_checks) - 1L) {
        x <- lo + idx * (hi - lo) / (n_point_checks - 1L)
        val <- xreplace(fseri, list(list(symb, x)))
        ref <- xreplace(fexpr, list(list(symb, x)))
        if (abs(as.numeric(N(1L - val / ref, prec))) > rt) { n_ok <- FALSE; break }
      }
      if (n_ok) cheapest <- fseri else break
    }
    if (is.null(cheapest)) fexpr else cheapest
  }
  query <- function(e) .cgp_is_function(e) && length(e$args) == 1L
  .cg_optim_obj(function(expr) .cg_replace(sym_factor(as_sym(expr)), query, value),
                cost_function, priority, cls = "SeriesApprox", bounds = bounds, reltol = reltol,
                max_order = max_order, n_point_checks = n_point_checks)
}
