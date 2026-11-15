# Power series rings in a named variable, after sympy/polys/series/ring.py
# (PowerSeriesRingRing, PowerSeriesRingField and PowerSeriesElement). The
# arithmetic is the port of ringpython.py in polys-series.R; this file adds
# the ring of SymPy's power_series_ring(), whose elements carry their ring
# and support the arithmetic operators, and the conversions to and from
# expressions.

#' Power series rings
#'
#' A port of SymPy's `sympy.polys.series`: `power_series_ring(symbol, K,
#' prec)` returns the ring of power series in one variable over `ZZ` or
#' `QQ` truncated at `O(x^prec)`, together with its generator, as
#' `list(R, x)`. Series are added, multiplied, divided and raised to integer
#' powers with the usual operators; an exact polynomial has no `O()` term
#' until it grows past the ring precision. The ring `R` holds the methods
#' of SymPy's ring as functions: `R$compose(s, t)`, `R$inverse(s)`,
#' `R$reversion(s)`, `R$differentiate(s)`, `R$truncate(s, n)`,
#' `R$square(s)`, `R$from_list(coeffs, prec)`, `R$from_expr(e)`,
#' `R$to_list(s)`, `R$to_expr(s)` and `R$order_term()`; over `QQ` also
#' `R$integrate()`, `R$sqrt()`, `R$log()`, `R$log1p()`, `R$exp()`,
#' `R$expm1()`, `R$sin()`, `R$cos()`, `R$tan()`, `R$sinh()`, `R$cosh()`,
#' `R$tanh()`, `R$asin()`, `R$asinh()`, `R$atan()`, `R$atanh()` and
#' `R$hypot()`. `PowerSeriesRingRing()` and `PowerSeriesRingField()` build
#' the rings over `ZZ` and `QQ` directly.
#'
#' @param symbol The variable name, or a symbol.
#' @param K The ground domain, `ZZ` or `QQ`.
#' @param prec The ring precision.
#' @param domain The ground domain.
#' @return `power_series_ring()` a list of the ring and its generator; the
#'   other constructors a ring.
#' @examples
#' r <- power_series_ring("x", QQ, 6)
#' R <- r[[1]]; x <- r[[2]]
#' (1 + x + x^2)^3
#' R$sin(x)
#' R$compose(R$exp(x) - 1, x + x^2)
#' R$reversion(R$sin(x))
#' R$to_expr(R$tan(x))
#' @name power_series_ring
#' @export
power_series_ring <- function(symbol, K, prec = 6) {
  kind <- attr(K, "dom")$kind
  R <- if (identical(kind, "ZZ")) PowerSeriesRingRing(K, symbol, prec)
       else if (identical(kind, "QQ")) PowerSeriesRingField(K, symbol, prec)
       else stop("ValueError: Unsupported ground domain: ", format(K), call. = FALSE)
  list(R, R$gen)
}

#' @rdname power_series_ring
#' @export
PowerSeriesRingRing <- function(domain = ZZ, symbol = "x", prec = 6) {
  if (!identical(attr(domain, "dom")$kind, "ZZ")) stop("ValueError: Unsupported ground domain: ", format(domain), call. = FALSE)
  .ps_ring(domain, symbol, prec, field = FALSE)
}
#' @rdname power_series_ring
#' @export
PowerSeriesRingField <- function(domain = QQ, symbol = "x", prec = 6) {
  if (!identical(attr(domain, "dom")$kind, "QQ")) stop("ValueError: Unsupported ground domain: ", format(domain), call. = FALSE)
  .ps_ring(domain, symbol, prec, field = TRUE)
}

.ps_ring <- function(domain, symbol, prec, field) {
  if (!(is.character(symbol) || is_sym(symbol))) stop("GeneratorsError: Symbol must be a string or a SymPy expression.", call. = FALSE)
  sym <- if (is.character(symbol)) Symbol(symbol) else symbol
  P <- .us_ring(if (field) "QQ" else "ZZ", prec)
  R <- new.env()
  R$domain <- domain; R$symbol <- sym; R$prec <- as.integer(prec); R$field <- field; R$P <- P
  el <- function(s) structure(list(s = s, ring = R), class = "PowerSeriesElement")
  conv <- function(v) {
    q <- if (is_sym(v)) gmp::as.bigq(sstr(v)) else gmp::as.bigq(v)
    if (!field) {
      if (gmp::denominator(q) != 1) stop("CoercionFailed: cannot convert ", .de_str(q), " to ZZ", call. = FALSE)
      return(gmp::as.bigz(q))
    }
    q
  }
  R$el <- el; R$conv <- conv
  R$from_list <- function(coeffs, prec = NULL) el(P$from_list(lapply(coeffs, conv), prec))
  R$from_ground <- function(a) R$from_list(list(a))
  R$one <- el(P$one); R$zero <- el(P$zero); R$gen <- el(P$gen)
  R$order_term <- function() el(P$from_list(list(), R$prec))
  R$truncate <- function(s, prec) el(P$truncate(s$s, prec))
  R$to_list <- function(s) lapply(P$to_list(s$s), .ps_tosym)
  R$to_expr <- function(s) {
    cf <- P$to_list(s$s); terms <- list()
    for (i in seq_along(cf)) if (cf[[i]] != 0) terms[[length(terms) + 1L]] <- .ps_tosym(cf[[i]]) * sym^(i - 1L)
    if (!is.null(P$series_prec(s$s))) terms[[length(terms) + 1L]] <- O(sym^P$series_prec(s$s))
    if (!length(terms)) return(S_Zero)
    Reduce(`+`, terms)
  }
  R$from_expr <- function(e) {
    e <- as_sym(e); prec <- NULL
    if (.has_type(e, "Order")) {
      ex <- getO(e)$args[[1L]]
      b <- if (ex$type == "Pow") ex$args[[1L]] else ex
      pw <- if (ex$type == "Pow") as.integer(ex$args[[2L]]$value) else 1L
      if (!identical(b$key, sym$key)) stop("ValueError: Order contains different symbol than ring.", call. = FALSE)
      prec <- min(pw, R$prec)
      e <- removeO(e)
    }
    rebuild <- function(e) {
      if (e$type == "Symbol") {
        if (identical(e$key, sym$key)) return(R$gen)
        stop("ValueError: Expr contains different symbol than ring.", call. = FALSE)
      }
      if (e$type == "Add") return(Reduce(`+`, lapply(e$args, rebuild)))
      if (e$type == "Mul") return(Reduce(`*`, lapply(e$args, rebuild)))
      if (e$type == "Pow") {
        if (e$args[[2L]]$type != "Integer") stop("ValueError: Invalid exponent in power series", call. = FALSE)
        return(rebuild(e$args[[1L]])^as.integer(e$args[[2L]]$value))
      }
      if (!e$type %in% c("Integer", "Rational")) stop("CoercionFailed: cannot convert ", sstr(e), " to ", format(domain), call. = FALSE)
      R$from_ground(e)
    }
    s <- rebuild(e)
    if (is.null(prec)) s else el(P$from_list(P$to_list(s$s), prec))
  }
  for (nm in c("square", "inverse", "reversion", "differentiate")) local({
    f <- nm; R[[f]] <- function(s) el(do.call(`$`, list(P, f))(s$s))
  })
  R$compose <- function(s, t) el(P$compose(s$s, t$s))
  if (field) {
    for (nm in c("integrate", "sqrt", "log", "log1p", "exp", "expm1", "sin", "cos", "tan", "sinh", "cosh", "tanh",
                 "asin", "asinh", "atan", "atanh")) local({
      f <- nm; R[[f]] <- function(s) el(do.call(`$`, list(P, f))(s$s))
    })
    R$hypot <- function(s, t) el(P$hypot(s$s, t$s))
  }
  structure(R, class = "PowerSeriesRing")
}

.ps_tosym <- function(q) { q <- gmp::as.bigq(q); Rational(as.character(gmp::numerator(q)), as.character(gmp::denominator(q))) }

#' @export
format.PowerSeriesRing <- function(x, ...)
  sprintf("Power Series Ring in %s over %s of size %d", sstr(x$symbol), format(x$domain), x$prec)
#' @export
print.PowerSeriesRing <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
#' @export
`$.PowerSeriesElement` <- function(x, name) {
  v <- unclass(x); P <- v$ring$P
  switch(name, prec = P$series_prec(v$s), series = v$s, is_ground = P$is_ground(v$s),
         constant_coefficient = .ps_tosym(P$constant_coefficient(v$s)),
         as_expr = v$ring$to_expr(x), removeO = v$ring$from_list(v$ring$to_list(x)),
         v[[name]])
}
# series_pprint with ^ for powers, as the R printers write them
#' @export
format.PowerSeriesElement <- function(x, ...) gsub("**", "^", x$ring$P$pretty(x$s, symbol = sstr(x$ring$symbol)), fixed = TRUE)
#' @export
print.PowerSeriesElement <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

#' @export
chooseOpsMethod.PowerSeriesElement <- function(x, y, mx, my, cl, reverse) TRUE
#' @export
Ops.PowerSeriesElement <- function(e1, e2) {
  isel <- function(v) inherits(v, "PowerSeriesElement")
  R <- if (isel(e1)) e1$ring else e2$ring
  P <- R$P; el <- R$el
  if (missing(e2)) return(switch(.Generic, "-" = el(P$negative(e1$s)), "+" = el(P$positive(e1$s)),
                                 stop("unsupported unary operator ", .Generic, call. = FALSE)))
  if (.Generic == "^") {
    n <- as.integer(if (is_sym(e2)) e2$value else e2)
    return(if (n == 0L) R$one else if (n == 1L) e1 else el(P$pow_int(e1$s, n)))
  }
  if (isel(e1) && isel(e2) && !identical(e1$ring, e2$ring)) stop("ValueError: Cannot combine power series with different rings.", call. = FALSE)
  if (.Generic %in% c("==", "!=")) {
    eq <- isel(e1) && isel(e2) && P$equal_repr(e1$s, e2$s)
    return(if (.Generic == "==") eq else !eq)
  }
  a <- if (isel(e1)) e1$s else NULL; b <- if (isel(e2)) e2$s else NULL
  switch(.Generic,
    "+" = el(if (is.null(b)) P$add_ground(a, R$conv(e2)) else if (is.null(a)) P$add_ground(b, R$conv(e1)) else P$add(a, b)),
    "-" = el(if (is.null(b)) P$subtract_ground(a, R$conv(e2)) else if (is.null(a)) P$rsubtract_ground(b, R$conv(e1)) else P$subtract(a, b)),
    "*" = el(if (is.null(b)) P$multiply_ground(a, R$conv(e2)) else if (is.null(a)) P$multiply_ground(b, R$conv(e1)) else P$multiply(a, b)),
    "/" = el(P$divide(if (is.null(a)) R$from_ground(e1)$s else a, if (is.null(b)) R$from_ground(e2)$s else b)),
    stop("operator ", .Generic, " is not defined for power series", call. = FALSE))
}
