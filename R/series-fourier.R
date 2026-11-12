# Fourier series, after sympy/series/fourier.py.

.fourier_limits <- function(f, limits) {
  find_x <- function() {
    fs <- free_symbols(f)
    if (length(fs) == 1L) return(fs[[1L]])
    if (!length(fs)) return(Dummy("k"))
    stop(" specify dummy variables for ", .str(f), ". If the function contains more than one free symbol, ",
         "a dummy variable should be supplied explicitly e.g. FourierSeries(m*n**2, (n, -pi, pi))", call. = FALSE)
  }
  if (is.null(limits)) return(list(find_x(), .neg(S_Pi), S_Pi))
  limits <- lapply(if (is_sym(limits) && limits$type == "Tuple") limits$args else limits, as_sym)
  lim <- if (length(limits) == 3L) limits else if (length(limits) == 2L) c(list(find_x()), limits) else NULL
  if (is.null(lim) || !is_symbol_like(lim[[1L]])) stop("Invalid limits given: ", deparse(limits), call. = FALSE)
  if (is_inf_const(lim[[2L]]) || is_inf_const(lim[[3L]]))
    stop("Both the start and end value should be bounded", call. = FALSE)
  lim
}

# A trigonometric polynomial with period 2*L as list(a0, an, bn), where an
# and bn are lists of list(k, coefficient); NULL when f is not one. The
# terms go through exponentials: c*exp(I*(k*pi*x/L + b)).
.fourier_finite <- function(f, x, L) {
  trig <- Filter(function(z) .has1(z, x), atoms(f, "Function"))
  if (!all(vapply(trig, function(z) z$name %in% c("sin", "cos"), TRUE))) return(NULL)
  e <- expand(rewrite(f, "exp"))
  terms <- if (e$type == "Add") e$args else list(e)
  C <- list()
  w <- .div(.mul(list(S_I, S_Pi)), L)
  for (t in terms) {
    facs <- if (t$type == "Mul") t$args else list(t)
    dep <- vapply(facs, function(z) .has1(z, x), TRUE)
    if (sum(dep) > 1L) return(NULL)
    if (!any(dep)) {
      C[["0"]] <- .add(list(if (is.null(C[["0"]])) S_Zero else C[["0"]], t))
      next
    }
    ex <- facs[[which(dep)]]
    if (!(ex$type == "Function" && ex$name == "exp")) return(NULL)
    arg <- expand(ex$args[[1L]])
    cl <- tryCatch(.coeff_list(arg, x), error = function(e) NULL)
    if (length(cl) != 2L) return(NULL)
    k <- simplify(.div(cl[[1L]], w))
    if (k$type != "Integer" || is_zero_num(k)) return(NULL)
    c <- .mul(c(facs[!dep], list(fcall("exp", list(cl[[2L]])))))
    kk <- .str(k)
    C[[kk]] <- .add(list(if (is.null(C[[kk]])) S_Zero else C[[kk]], c))
  }
  tidy <- function(v) {
    v <- simplify(expand(v))
    if (.has_function(v, "exp")) v <- simplify(rewrite(v, "cos"))
    v
  }
  get <- function(k) { v <- C[[as.character(k)]]; if (is.null(v)) S_Zero else v }
  # harmonic indices stay exact integers
  ks <- unique(abs(gmp::as.bigz(setdiff(names(C), "0"))))
  ks <- if (length(ks)) sort(ks) else ks
  an <- list(); bn <- list()
  for (i in seq_along(ks)) {
    k <- ks[i]
    a <- tidy(.add(list(get(k), get(-k))))
    b <- tidy(.mul(list(S_I, .sub(get(k), get(-k)))))
    if (!is_zero_num(a)) an[[length(an) + 1L]] <- list(k, a)
    if (!is_zero_num(b)) bn[[length(bn) + 1L]] <- list(k, b)
  }
  list(tidy(get(0)), an, bn)
}

.fourier_L <- function(lim) .div(fcall("Abs", list(.sub(lim[[3L]], lim[[2L]]))), S_Two)

#' Fourier series
#'
#' `fourier_series(f, limits, finite)` is the Fourier series of `f` over
#' `limits = list(x, a, b)` (default `list(x, -pi, pi)`). Trigonometric
#' polynomials give a `FiniteFourierSeries`. Use [truncate()] for the first
#' `n` nonzero terms, `sigma_approximation()` for the Lanczos smoothed sum,
#' and `shift()`, `shiftx()`, `scale(s, k)`, `scalex()` to transform the series
#' as SymPy does; `+` and `-` add series with the same period.
#'
#' @param f An expression.
#' @param limits `list(x, a, b)` or `list(a, b)`.
#' @param finite Detect trigonometric polynomials.
#' @param s,x,con A Fourier series (`con` for [truncate()], `x` for [scale()]).
#' @param n Number of terms; `NULL` in `truncate()` returns all terms of a
#'   finite series.
#' @param k The shift or scale amount.
#' @param center The factor for `scale(s, k)` (the second argument of the
#'   [base::scale()] generic).
#' @param scale Ignored; part of the [base::scale()] generic.
#' @param ... Unused.
#' @return A Fourier series or an expression.
#' @examples
#' x <- Symbol("x")
#' s <- fourier_series(x^2, list(x, -pi, pi))
#' truncate(s, 4)
#' fourier_series(sin(x)^2)
#' @export
fourier_series <- function(f, limits = NULL, finite = TRUE) {
  f <- as_sym(f)
  lim <- .fourier_limits(f, limits)
  x <- lim[[1L]]
  if (!.has1(f, x)) return(f)
  L <- .fourier_L(lim)
  if (finite) {
    fin <- .fourier_finite(f, x, L)
    if (!is.null(fin)) return(.finite_fourier(f, lim, fin))
  }
  n <- Dummy("n")
  center <- .div(.add(list(lim[[2L]], lim[[3L]])), S_Two)
  per <- .sub(lim[[3L]], lim[[2L]])
  cos_t <- fcall("cos", list(.div(.mul(list(S_Two, n, S_Pi, x)), per)))
  sin_t <- fcall("sin", list(.div(.mul(list(S_Two, n, S_Pi, x)), per)))
  lims <- list(x, lim[[2L]], lim[[3L]])
  # the integral with symbolic n, with SymPy's Piecewise for the integer
  # indices where the generic formula has a pole
  coef <- function(trig) {
    g <- .mul(list(f, trig))
    I <- do.call(sym_integrate, list(g, lims))
    den <- tryCatch(fraction(together(I))[[2L]], error = function(e) S_One)
    if (!.has1(den, n) || !is_polynomial(den, n)) return(I)
    rts <- tryCatch(roots(den, n), error = function(e) list())
    # integer roots and parameter-dependent ones (SymPy keeps both)
    ks <- Filter(function(r) r$type == "Integer" || (length(free_symbols(r)) && !isFALSE(ask_fact(r, "integer"))), rts)
    if (!length(ks)) return(I)
    num <- vapply(ks, function(r) r$type == "Integer", TRUE)
    ks <- c(ks[num][order(vapply(ks[num], function(r) as.numeric(.Z(r)), 0))], ks[!num])
    pairs <- lapply(ks, function(k) list(do.call(sym_integrate, list(subs(g, n, k), lims)), Eq(n, k)))
    do.call(Piecewise, c(pairs, list(list(I, S_true))))
  }
  cos_seq <- function() {
    an <- SeqFormula(.div(.mul(list(S_Two, cos_t, coef(cos_t))), per),
                     list(n, S_One, S_Infinity))
    a0 <- simplify(.div(do.call(sym_integrate, list(f, lims)), per))
    list(a0, an)
  }
  sin_seq <- function()
    SeqFormula(.div(.mul(list(S_Two, sin_t, coef(sin_t))), per),
               list(n, S_One, S_Infinity))
  zero_seq <- function() SeqFormula(S_Zero, list(S_One, S_Infinity))
  if (is_zero_num(simplify(center))) {
    nf <- subs(f, x, .neg(x))
    if (identical(expand(f)$key, expand(nf)$key)) {
      ca <- cos_seq()
      return(.fourier_node(f, lim, ca[[1L]], ca[[2L]], zero_seq()))
    }
    if (identical(expand(f)$key, expand(.neg(nf))$key))
      return(.fourier_node(f, lim, S_Zero, zero_seq(), sin_seq()))
  }
  ca <- cos_seq()
  .fourier_node(f, lim, ca[[1L]], ca[[2L]], sin_seq())
}

.fourier_node <- function(f, lim, a0, an, bn)
  new_sym("FourierSeries", list(f, do.call(Tuple, lim), Tuple(a0, an, bn)))

.dict_tuple <- function(d) do.call(Tuple, lapply(d, function(p) Tuple(.zs(gmp::as.bigz(p[[1L]])), p[[2L]])))
.finite_fourier <- function(f, lim, fin)
  new_sym("FiniteFourierSeries", list(f, do.call(Tuple, lim), Tuple(fin[[1L]], .dict_tuple(fin[[2L]]), .dict_tuple(fin[[3L]]))))

.is_fourier <- function(s) is_sym(s) && s$type %in% c("FourierSeries", "FiniteFourierSeries")

# Term pt of the series (pt = 0 is a0).
.fourier_term <- function(s, pt) {
  ex <- s$args[[3L]]$args
  if (identical(as.character(pt), "0")) return(ex[[1L]])
  x <- s$args[[2L]]$args[[1L]]
  lim <- s$args[[2L]]$args
  if (s$type == "FourierSeries") {
    # resonant indices are resolved by the Piecewise built in fourier_series()
    return(.add(list(.seq_coeff(ex[[2L]], pt), .seq_coeff(ex[[3L]], pt))))
  }
  L <- .fourier_L(lim)
  pk <- as.character(pt)
  lookup <- function(d) {
    for (p in d$args) if (.str(p$args[[1L]]) == pk) return(p$args[[2L]])
    S_Zero
  }
  th <- .div(.mul(list(.zs(gmp::as.bigz(pt)), S_Pi, x)), L)
  .add(list(.mul(list(lookup(ex[[2L]]), fcall("cos", list(th)))), .mul(list(lookup(ex[[3L]]), fcall("sin", list(th))))))
}

# Harmonic indices of a finite series, sorted, as strings.
.fourier_keys <- function(s) {
  ks <- unlist(lapply(s$args[[3L]]$args[2:3], function(d) vapply(d$args, function(p) .str(p$args[[1L]]), "")))
  if (!length(ks)) return(character())
  ks <- unique(ks)
  ks[order(gmp::as.bigz(ks))]
}

#' @rdname fourier_series
#' @export
truncate.sym <- function(con, n, ...) {
  # One method serves both series classes: FormalPowerSeries truncates to
  # 6 terms by default, as FormalPowerSeries.truncate; Fourier series to 3.
  if (con$type == "FormalPowerSeries" || .ffps_is(con)) return(fps_truncate(con, if (missing(n)) 6 else n))
  if (missing(n)) n <- 3
  s <- con
  if (!.is_fourier(s)) stop("truncate() needs a Fourier or formal power series", call. = FALSE)
  finite <- s$type == "FiniteFourierSeries"
  if (!finite && is.null(n)) stop("n = NULL needs a finite series", call. = FALSE)
  terms <- list()
  # a finite series visits only its harmonics, which may be huge
  pts <- if (finite) c("0", .fourier_keys(s)) else NULL
  i <- 0L
  repeat {
    if (!is.null(n) && length(terms) >= n) break
    i <- i + 1L
    if (finite && i > length(pts)) break
    pt <- if (finite) pts[i] else i - 1L
    t <- .fourier_term(s, pt)
    if (!is_zero_num(t)) terms[[length(terms) + 1L]] <- t
  }
  .add(terms)
}

#' @rdname fourier_series
#' @export
sigma_approximation <- function(s, n = 3) {
  terms <- list()
  for (i in seq_len(n) - 1L) {
    t <- .fourier_term(s, i)
    if (!is_zero_num(t)) terms[[length(terms) + 1L]] <- .mul(list(fcall("sinc", list(.div(.mul(list(S_Pi, .int(i))), .int(n)))), t))
  }
  .add(terms)
}

.fourier_check_indep <- function(s, k) {
  k <- as_sym(k)
  x <- s$args[[2L]]$args[[1L]]
  if (.has1(k, x)) stop("'", .str(k), "' should be independent of ", .str(x), call. = FALSE)
  k
}

# Rebuild a finite series from a transformed truncated expression.
.fourier_refinite <- function(s, f, e) {
  lim <- s$args[[2L]]$args
  x <- lim[[1L]]; L <- .fourier_L(lim)
  # terms that are no longer harmonics of the period join a0, as SymPy
  e <- expand(e)
  parts <- if (e$type == "Add") e$args else list(e)
  ok <- vapply(parts, function(t) !is.null(.fourier_finite(t, x, L)), TRUE)
  fin <- if (any(ok)) .fourier_finite(.add(parts[ok]), x, L) else list(S_Zero, list(), list())
  fin[[1L]] <- .add(c(list(fin[[1L]]), parts[!ok]))
  .finite_fourier(f, lim, fin)
}

#' @rdname fourier_series
#' @export
shift <- function(s, k) {
  k <- .fourier_check_indep(s, k)
  ex <- s$args[[3L]]$args
  if (s$type == "FiniteFourierSeries")
    return(.fourier_refinite(s, .add(list(s$args[[1L]], k)), .add(list(truncate(s, NULL), k))))
  .fourier_node(.add(list(s$args[[1L]], k)), s$args[[2L]]$args, .add(list(ex[[1L]], k)), ex[[2L]], ex[[3L]])
}

#' @rdname fourier_series
#' @export
shiftx <- function(s, k) {
  k <- .fourier_check_indep(s, k)
  x <- s$args[[2L]]$args[[1L]]
  ex <- s$args[[3L]]$args
  xk <- .add(list(x, k))
  if (s$type == "FiniteFourierSeries")
    return(.fourier_refinite(s, subs(s$args[[1L]], x, xk), subs(truncate(s, NULL), x, xk)))
  .fourier_node(subs(s$args[[1L]], x, xk), s$args[[2L]]$args, ex[[1L]], subs(ex[[2L]], x, xk), subs(ex[[3L]], x, xk))
}

#' @rdname fourier_series
#' @export
scale.sym <- function(x, center = TRUE, scale = TRUE) {
  s <- x
  if (!.is_fourier(s)) stop("scale() needs a Fourier series", call. = FALSE)
  # base::scale() has no `...`: the factor is the second argument
  k <- .fourier_check_indep(s, if (isTRUE(center)) 1 else center)
  ex <- s$args[[3L]]$args
  if (s$type == "FiniteFourierSeries")
    return(.fourier_refinite(s, .mul(list(s$args[[1L]], k)), .mul(list(truncate(s, NULL), k))))
  .fourier_node(.mul(list(s$args[[1L]], k)), s$args[[2L]]$args, .mul(list(ex[[1L]], k)),
                .seq_coeff_mul(ex[[2L]], k), .seq_coeff_mul(ex[[3L]], k))
}

#' @rdname fourier_series
#' @export
scalex <- function(s, k) {
  k <- .fourier_check_indep(s, k)
  x <- s$args[[2L]]$args[[1L]]
  ex <- s$args[[3L]]$args
  xk <- .mul(list(x, k))
  if (s$type == "FiniteFourierSeries")
    return(.fourier_refinite(s, subs(s$args[[1L]], x, xk), subs(truncate(s, NULL), x, xk)))
  .fourier_node(subs(s$args[[1L]], x, xk), s$args[[2L]]$args, ex[[1L]], subs(ex[[2L]], x, xk), subs(ex[[3L]], x, xk))
}

.fourier_ops <- function(gen, e1, e2) {
  if (gen == "-") e2 <- if (.is_fourier(e2)) scale(e2, -1) else .neg(as_sym(e2))
  if (!(gen %in% c("+", "-"))) return(NULL)
  if (!(.is_fourier(e1) && .is_fourier(e2))) return(.add(list(as_sym(e1), as_sym(e2))))
  l1 <- e1$args[[2L]]$args; l2 <- e2$args[[2L]]$args
  if (!identical(l1[[2L]]$key, l2[[2L]]$key) || !identical(l1[[3L]]$key, l2[[3L]]$key))
    stop("Both the series should have same periods", call. = FALSE)
  x <- l1[[1L]]
  fn <- .add(list(e1$args[[1L]], subs(e2$args[[1L]], l2[[1L]], x)))
  if (!.has1(fn, x)) return(fn)
  if (e1$type == "FiniteFourierSeries" || e2$type == "FiniteFourierSeries") return(fourier_series(fn, l1))
  a <- e1$args[[3L]]$args; b <- e2$args[[3L]]$args
  .fourier_node(fn, l1, .add(list(a[[1L]], b[[1L]])), SeqAdd(a[[2L]], b[[2L]]), SeqAdd(a[[3L]], b[[3L]]))
}

register_bound_vars("FourierSeries", function(e) e$args[[2L]]$args[[1L]]$key)
register_bound_vars("FiniteFourierSeries", function(e) e$args[[2L]]$args[[1L]]$key)
register_str_printer("FourierSeries", function(e)
  paste0("FourierSeries(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", ", .str(e$args[[3L]]), ")"))
register_str_printer("FiniteFourierSeries", function(e) {
  ex <- e$args[[3L]]$args
  d <- function(t) paste0("{", paste(vapply(t$args, function(p) paste0(.str(p$args[[1L]]), ": ", .str(p$args[[2L]])), ""),
                                     collapse = ", "), "}")
  paste0("FiniteFourierSeries(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", (", .str(ex[[1L]]), ", ",
         d(ex[[2L]]), ", ", d(ex[[3L]]), "))")
})
for (.t in c("FourierSeries", "FiniteFourierSeries"))
  register_latex_printer(.t, function(e) paste0(.latex(truncate(e)), " + \\ldots"))
rm(.t)
