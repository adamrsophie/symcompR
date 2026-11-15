# Plotting with base graphics. The probability code that used to live here
# is in the stats-*.R files, a port of sympy.stats.


#' Plot symbolic expressions with base graphics
#'
#' `plot(expr, range)` draws a function of one variable, with `range` given
#' as `list(x, a, b)` (default `-10` to `10`). `plot_parametric()`,
#' `plot3d()` (surface), `plot_implicit()` (zero contour), `plot_contour()`,
#' `plot3d_parametric_line()` and `plot3d_parametric_surface()` are the
#' analogues of SymPy's plotting functions. Extra arguments go to the base
#' graphics functions.
#'
#' @param x A `sym` expression, for `plot()`.
#' @param range `list(var, a, b)`.
#' @param ... Graphical parameters.
#' @param fx,fy Coordinate expressions.
#' @param e An expression.
#' @param range1,range2 Ranges for two variables.
#' @param n Number of sample points.
#' @return The sampled data, invisibly.
#' @examples
#' x <- Symbol("x")
#' if (interactive()) plot(sin(x) / x, list(x, -10, 10))
#' @export
plot.sym <- function(x, range = NULL, ..., n = 400) {
  if (is.null(range)) {
    fs <- free_symbols(x)
    if (length(fs) != 1L) stop("plot() needs a range for multivariate expressions")
    range <- list(fs[[1L]], -10, 10)
  }
  v <- as_sym(range[[1L]])
  xs <- seq(as.numeric(N(as_sym(range[[2L]]))), as.numeric(N(as_sym(range[[3L]]))), length.out = n)
  ys <- .eval_grid(x, list(v), list(xs))
  graphics::plot(xs, ys, type = "l", xlab = .str(v), ylab = .str(x), ...)
  invisible(data.frame(x = xs, y = ys))
}

.eval_grid <- function(e, vars, vals) {
  out <- rep(NA_real_, length(vals[[1L]]))
  f <- tryCatch(do.call(lambdify, list(vars, e)), error = function(err) NULL)
  if (!is.null(f)) {
    r <- tryCatch(suppressWarnings(do.call(f, vals)), error = function(err) NULL)
    if (!is.null(r)) {
      r <- rep_len(r, length(out))
      return(if (is.complex(r)) ifelse(abs(Im(r)) < 1e-12, Re(r), NA_real_) else as.numeric(r))
    }
  }
  nms <- vapply(vars, function(s) s$name, "")
  for (i in seq_along(out)) {
    env <- stats::setNames(lapply(vals, `[`, i), nms)
    v <- tryCatch(.neval(e, env), error = function(err) NA_real_)
    out[i] <- if (is.complex(v)) { if (abs(Im(v)) < 1e-12) Re(v) else NA_real_ } else as.numeric(v)
  }
  out
}

#' @rdname plot.sym
#' @export
plot_parametric <- function(fx, fy, range, ..., n = 400) {
  v <- as_sym(range[[1L]])
  ts <- seq(as.numeric(N(as_sym(range[[2L]]))), as.numeric(N(as_sym(range[[3L]]))), length.out = n)
  xs <- .eval_grid(as_sym(fx), list(v), list(ts))
  ys <- .eval_grid(as_sym(fy), list(v), list(ts))
  graphics::plot(xs, ys, type = "l", xlab = .str(as_sym(fx)), ylab = .str(as_sym(fy)), ...)
  invisible(data.frame(t = ts, x = xs, y = ys))
}

.grid2 <- function(e, r1, r2, n) {
  v1 <- as_sym(r1[[1L]]); v2 <- as_sym(r2[[1L]])
  xs <- seq(as.numeric(N(as_sym(r1[[2L]]))), as.numeric(N(as_sym(r1[[3L]]))), length.out = n)
  ys <- seq(as.numeric(N(as_sym(r2[[2L]]))), as.numeric(N(as_sym(r2[[3L]]))), length.out = n)
  g <- expand.grid(xs, ys)
  z <- matrix(.eval_grid(as_sym(e), list(v1, v2), list(g[[1L]], g[[2L]])), n, n)
  list(x = xs, y = ys, z = z)
}

#' @rdname plot.sym
#' @export
plot3d <- function(e, range1, range2, ..., n = 50) {
  g <- .grid2(e, range1, range2, n)
  graphics::persp(g$x, g$y, g$z, xlab = .str(as_sym(range1[[1L]])), ylab = .str(as_sym(range2[[1L]])),
                  zlab = "z", theta = 30, phi = 25, ...)
  invisible(g)
}

#' @rdname plot.sym
#' @export
plot_implicit <- function(e, range1, range2, ..., n = 200) {
  e <- .eq_to_expr(e)
  g <- .grid2(e, range1, range2, n)
  graphics::contour(g$x, g$y, g$z, levels = 0, drawlabels = FALSE,
                    xlab = .str(as_sym(range1[[1L]])), ylab = .str(as_sym(range2[[1L]])), ...)
  invisible(g)
}

#' @rdname plot.sym
#' @param levels Contour levels; by default chosen by [graphics::contour()].
#' @export
plot_contour <- function(e, range1, range2, ..., levels = NULL, n = 100) {
  g <- .grid2(e, range1, range2, n)
  args <- list(g$x, g$y, g$z, xlab = .str(as_sym(range1[[1L]])), ylab = .str(as_sym(range2[[1L]])), ...)
  if (!is.null(levels)) args$levels <- levels
  do.call(graphics::contour, args)
  invisible(g)
}

# an empty perspective box around the points, for drawing 3D lines and facets
.persp_box <- function(xs, ys, zs, labs, ...) {
  rg <- function(v) { r <- range(v, finite = TRUE); if (diff(r) == 0) r + c(-1, 1) else r }
  graphics::persp(rg(xs), rg(ys), matrix(NA_real_, 2, 2), zlim = rg(zs), xlab = labs[1L], ylab = labs[2L],
                  zlab = labs[3L], theta = 30, phi = 25, ...)
}

#' @rdname plot.sym
#' @param fz The third coordinate expression.
#' @export
plot3d_parametric_line <- function(fx, fy, fz, range, ..., n = 300) {
  v <- as_sym(range[[1L]])
  ts <- seq(as.numeric(N(as_sym(range[[2L]]))), as.numeric(N(as_sym(range[[3L]]))), length.out = n)
  xyz <- lapply(list(fx, fy, fz), function(f) .eval_grid(as_sym(f), list(v), list(ts)))
  pm <- .persp_box(xyz[[1L]], xyz[[2L]], xyz[[3L]], vapply(list(fx, fy, fz), function(f) .str(as_sym(f)), ""), ...)
  graphics::lines(grDevices::trans3d(xyz[[1L]], xyz[[2L]], xyz[[3L]], pm))
  invisible(data.frame(t = ts, x = xyz[[1L]], y = xyz[[2L]], z = xyz[[3L]]))
}

#' @rdname plot.sym
#' @export
plot3d_parametric_surface <- function(fx, fy, fz, range1, range2, ..., n = 30) {
  gs <- lapply(list(fx, fy, fz), function(f) .grid2(f, range1, range2, n))
  X <- gs[[1L]]$z; Y <- gs[[2L]]$z; Z <- gs[[3L]]$z
  pm <- .persp_box(X, Y, Z, vapply(list(fx, fy, fz), function(f) .str(as_sym(f)), ""), ...)
  # facets drawn back to front would need depth sorting; a wire mesh needs none
  for (i in seq_len(n)) {
    graphics::lines(grDevices::trans3d(X[i, ], Y[i, ], Z[i, ], pm))
    graphics::lines(grDevices::trans3d(X[, i], Y[, i], Z[, i], pm))
  }
  invisible(list(u = gs[[1L]]$x, v = gs[[1L]]$y, x = X, y = Y, z = Z))
}

#' @rdname plot.sym
#' @param nrows,ncols The grid shape.
#' @details `PlotGrid(nrows, ncols, ...)` draws several plots in one
#'   figure; each argument in `...` is a function of no arguments that draws
#'   one plot, as SymPy's `PlotGrid` takes plot objects.
#' @export
PlotGrid <- function(nrows, ncols, ...) {
  old <- graphics::par(mfrow = c(nrows, ncols))
  on.exit(graphics::par(old))
  invisible(lapply(list(...), function(f) if (is.function(f)) f() else f))
}
