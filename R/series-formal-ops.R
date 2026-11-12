# FormalPowerSeries arithmetic, after sympy/series/formal.py: derivative,
# integration, sums and constant multiples, product, composition and
# inverse through incomplete Bell polynomials.

# bell._bell_incomplete_poly: the incomplete Bell polynomial B_{n,k}
.bell_incomplete <- function(n, k, syms) {
  if (n == 0L && k == 0L) return(S_One)
  if (n == 0L || k == 0L) return(S_Zero)
  s <- S_Zero; a <- S_One
  for (m in seq_len(n - k + 1L)) {
    s <- .add(list(s, .mul(list(a, .bell_incomplete(n - m, k - 1L, syms), syms[[m]]))))
    a <- .div(.mul(list(a, .int(n - m))), .int(m))
  }
  expand_mul(s)
}

# FormalPowerSeries methods
.fps_new <- function(f, x, x0, dir, ak, xk, ind) new_sym("FormalPowerSeries", list(f, x, x0, dir, Tuple(ak, xk, ind)))
.fps_k <- function(seq) .seq_var(seq)
.fps_ak_seq <- function(p) { k <- .fps_k(p$ak); SeqFormula(p$ak$args[[1L]], list(k, 1, oo)) }
.fps_fact_seq <- function(p) { k <- .fps_k(p$ak); SeqFormula(fcall("factorial", list(k)), list(k, 1, oo)) }
.fps_bell_coeff_seq <- function(p) SeqMul(.fps_ak_seq(p), .fps_fact_seq(p))
.fps_sign_seq <- function(p) { k <- .fps_k(p$ak); SeqPer(list(-1, 1), list(k, 1, oo)) }
.fps_zero_coeff <- function(s) .fps_eval_term(s, S_Zero)
.fps_seq_head <- function(sq, n) lapply(seq_len(n) - 1L, function(i) .seq_coeff(sq, .seq_ith_point(sq, i)))

# FormalPowerSeries._eval_derivative
.fps_derivative <- function(s, x) {
  p <- .fps_parts(s)
  f <- diff(p$f, x); ind <- diff(p$ind, x)
  pow_xk <- .fps_get_pow_x(p$xk$args[[1L]], p$x)
  ak <- p$ak; k <- .fps_k(ak); b <- .seq_bounds(ak)
  form <- ak$args[[1L]]
  if (.has1(form, x)) {
    pieces <- lapply(form$args, function(ec) {
      temp <- .add(lapply(.fps_add_args(ec$args[[1L]]), function(t) .mul(list(t, .add(list(pow_xk, .fps_get_pow_x(t, p$x)))))))
      list(temp, ec$args[[2L]])
    })
    ak <- SeqFormula(subs(.piecewise(pieces), k, .add(list(k, S_One))), list(k, .sub(b[[1L]], S_One), b[[2L]]))
  } else ak <- SeqFormula(subs(.mul(list(form, pow_xk)), k, .add(list(k, S_One))), list(k, .sub(b[[1L]], S_One), b[[2L]]))
  .fps_new(f, p$x, p$x0, p$dir, ak, p$xk, ind)
}
register_diff_hook("FormalPowerSeries", function(e, x) if (identical(x$key, e$args[[2L]]$key)) .fps_derivative(e, x) else S_Zero)

#' @rdname fps
#' @export
fps_integrate <- function(s, x = NULL) {
  p <- .fps_parts(s)
  x <- if (is.null(x)) p$x else as_sym(x)
  f <- integrate(p$f, x)
  ind <- integrate(p$ind, x)
  ind <- .add(list(ind, limit(.sub(f, ind), x, S_Zero)))
  pow_xk <- .fps_get_pow_x(p$xk$args[[1L]], p$x)
  ak <- p$ak; k <- .fps_k(ak); b <- .seq_bounds(ak)
  form <- ak$args[[1L]]
  if (.has1(form, x)) {
    pieces <- lapply(form$args, function(ec) {
      temp <- .add(lapply(.fps_add_args(ec$args[[1L]]), function(t) .div(t, .add(list(pow_xk, .fps_get_pow_x(t, p$x), S_One)))))
      list(temp, ec$args[[2L]])
    })
    ak <- SeqFormula(subs(.piecewise(pieces), k, .sub(k, S_One)), list(k, .add(list(b[[1L]], S_One)), b[[2L]]))
  } else ak <- SeqFormula(subs(.div(form, .add(list(pow_xk, S_One))), k, .sub(k, S_One)), list(k, .add(list(b[[1L]], S_One)), b[[2L]]))
  .fps_new(f, p$x, p$x0, p$dir, ak, p$xk, ind)
}

# FormalPowerSeries.__add__, __neg__, __mul__ (by a constant)
.fps_add <- function(a, b) {
  pa <- .fps_parts(a)
  if (is_sym(b) && b$type == "FormalPowerSeries") {
    pb <- .fps_parts(b)
    if (!identical(pa$dir$key, pb$dir$key)) stop(.sympy_condition("ValueError", "Both series should be calculated from the same direction."))
    if (!identical(pa$x0$key, pb$x0$key)) stop(.sympy_condition("ValueError", "Both series should be calculated about the same point."))
    f <- .add(list(pa$f, subs(pb$f, pb$x, pa$x)))
    if (!.has1(f, pa$x)) return(f)
    ak <- SeqAdd(pa$ak, pb$ak)
    sa <- .seq_bounds(pa$ak)[[1L]]; sb <- .seq_bounds(pb$ak)[[1L]]
    if (.is_true(Gt(sa, sb))) { sq <- pb$ak; s <- sb; e <- sa } else { sq <- pa$ak; s <- sa; e <- sb }
    si <- as.integer(num_value_dbl(s)); ei <- as.integer(num_value_dbl(e))
    save <- if (ei > si) .add(lapply(seq_len(ei - si) - 1L, function(i) .mul(list(.seq_coeff(sq, .int(si + i)), .seq_coeff(pa$xk, .int(si + i)))))) else S_Zero
    return(.fps_new(f, pa$x, pa$x0, pa$dir, ak, pa$xk, .add(list(pa$ind, pb$ind, save))))
  }
  b <- as_sym(b)
  if (!.has1(b, pa$x)) return(.fps_new(.add(list(pa$f, b)), pa$x, pa$x0, pa$dir, pa$ak, pa$xk, .add(list(pa$ind, b))))
  new_sym("Add", list(a, b), extra = list(types = c(a$type, b$type), keys = c(a$key, b$key)))
}
.fps_neg <- function(a) {
  p <- .fps_parts(a)
  .fps_new(.neg(p$f), p$x, p$x0, p$dir, .seq_coeff_mul(p$ak, S_NegativeOne), p$xk, .neg(p$ind))
}
.fps_scale <- function(a, c) {
  p <- .fps_parts(a)
  if (.has1(c, p$x)) return(NULL)
  .fps_new(.mul(list(p$f, c)), p$x, p$x0, p$dir, .seq_coeff_mul(p$ak, c), p$xk, .mul(list(p$ind, c)))
}
.fps_is <- function(e) is_sym(e) && e$type == "FormalPowerSeries"
.ops_fps <- function(gen, e1, e2) {
  a <- if (is_sym(e1)) e1 else as_sym(e1); b <- if (is_sym(e2)) e2 else as_sym(e2)
  if (!(.fps_is(a) || .fps_is(b))) return(NULL)
  switch(gen,
    "+" = if (.fps_is(a)) .fps_add(a, b) else .fps_add(b, a),
    "-" = if (.fps_is(a)) .fps_add(a, if (.fps_is(b)) .fps_neg(b) else .neg(b)) else .fps_add(.fps_neg(b), a),
    "*" = if (.fps_is(a) && !.fps_is(b)) .fps_scale(a, b) else if (.fps_is(b) && !.fps_is(a)) .fps_scale(b, a) else NULL,
    NULL)
}

# FormalPowerSeries.coeff_bell
.fps_coeff_bell <- function(s, n) {
  p <- .fps_parts(s)
  bc <- .fps_seq_head(.fps_bell_coeff_seq(p), n)
  inner <- lapply(seq_len(n), function(j) .bell_incomplete(n, j, bc[seq_len(n - j + 1L)]))
  SeqPer(inner, list(Dummy("k"), 1, oo))
}
.fps_check_pair <- function(a, b) {
  pa <- .fps_parts(a); pb <- .fps_parts(b)
  if (!identical(pa$dir$key, pb$dir$key)) stop(.sympy_condition("ValueError", "Both series should be calculated from the same direction."))
  if (!identical(pa$x0$key, pb$x0$key)) stop(.sympy_condition("ValueError", "Both series should be calculated about the same point."))
  if (!identical(pa$x$key, pb$x$key)) stop(.sympy_condition("ValueError", "Both series should have the same symbol."))
}

# FiniteFormalPowerSeries and its subclasses: lazy products, compositions
# and inverses whose terms come from truncate().
.ffps_types <- c("FormalPowerSeriesProduct", "FormalPowerSeriesCompose", "FormalPowerSeriesInverse")
.ffps_is <- function(e) is_sym(e) && e$type %in% .ffps_types
for (.t in .ffps_types) local({
  nm <- .t
  register_str_printer(nm, function(e) paste0(nm, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
})
rm(.t)

#' @rdname fps
#' @param other A second `FormalPowerSeries`.
#' @export
fps_product <- function(s, other, n = 6) {
  if (!.fps_is(other)) stop(.sympy_condition("ValueError", "Both series should be an instance of FormalPowerSeries class."))
  .fps_check_pair(s, other)
  new_sym("FormalPowerSeriesProduct", list(s, other))
}

#' @rdname fps
#' @export
fps_compose <- function(s, other, n = 6) {
  if (!.fps_is(other)) stop(.sympy_condition("ValueError", "Both series should be an instance of FormalPowerSeries class."))
  .fps_check_pair(s, other)
  # other._eval_term(0).as_coeff_mul(other.x)[0]: the factors free of x
  z <- .fps_zero_coeff(other); gx <- .fps_parts(other)$x
  cf <- if (!.has1(z, gx)) z else if (z$type == "Mul") .mul(Filter(function(a) !.has1(a, gx), z$args)) else S_One
  if (!is_zero_num(cf))
    stop(.sympy_condition("ValueError", "The formal power series of the inner function should not have any constant coefficient term."))
  new_sym("FormalPowerSeriesCompose", list(s, other))
}

#' @rdname fps
#' @export
fps_inverse <- function(s, n = 6) {
  if (isTRUE(ask_fact(.fps_zero_coeff(s), "zero")))
    stop(.sympy_condition("ValueError", "Constant coefficient should exist for an inverse of a formal power series to exist."))
  new_sym("FormalPowerSeriesInverse", list(s))
}

# FiniteFormalPowerSeries._eval_terms of each subclass
.ffps_terms <- function(e, n) {
  pf <- .fps_parts(e$args[[1L]])
  switch(e$type,
    FormalPowerSeriesProduct = {
      pg <- .fps_parts(e$args[[2L]])
      c1 <- SeqFormula(pf$ak$args[[1L]], list(.fps_k(pf$ak), 0, oo))
      c2 <- SeqFormula(pg$ak$args[[1L]], list(.fps_k(pg$ak), 0, oo))
      aks <- convolution(.fps_seq_head(c1, n), .fps_seq_head(c2, n))
      .add(lapply(seq_len(n) - 1L, function(i) .mul(list(as_sym(aks[[i + 1L]]), .seq_coeff(pf$xk, .int(i))))))
    },
    FormalPowerSeriesCompose = {
      f <- e$args[[1L]]; g <- e$args[[2L]]
      terms <- list(.fps_zero_coeff(f))
      for (i in seq_len(n - 1L)) {
        sq <- SeqMul(.fps_bell_coeff_seq(pf), .fps_coeff_bell(g, i))
        terms[[length(terms) + 1L]] <- .mul(list(.div(.add(.fps_seq_head(sq, i)), fcall("factorial", list(.int(i)))), .seq_coeff(pf$xk, .int(i))))
      }
      .add(terms)
    },
    FormalPowerSeriesInverse = {
      f <- e$args[[1L]]
      k <- .fps_k(pf$xk)
      inv_seq <- SeqFormula(.pow(.fps_zero_coeff(f), .neg(.add(list(k, S_One)))), list(k, 1, oo))
      aux <- SeqMul(SeqMul(.fps_sign_seq(pf), .fps_fact_seq(pf)), inv_seq)
      terms <- list(.fps_zero_coeff(f))
      for (i in seq_len(n - 1L)) {
        sq <- SeqMul(aux, .fps_coeff_bell(f, i))
        terms[[length(terms) + 1L]] <- .mul(list(.div(.add(.fps_seq_head(sq, i)), fcall("factorial", list(.int(i)))), .seq_coeff(pf$xk, .int(i))))
      }
      .add(terms)
    })
}
.ffps_truncate <- function(e, n = 6) {
  pf <- .fps_parts(e$args[[1L]])
  n <- as.integer(n)
  .add(list(.ffps_terms(e, n), .so_new(.seq_coeff(pf$xk, .int(n)), list(pf$x), list(pf$x0))))
}
