# Complex regions, after ComplexRegion in sympy/sets/fancysets.py and its
# intersection and union handlers: subsets of the complex plane given as
# products of real sets, in rectangular (x + y*I) or polar form.

#' Regions of the complex plane
#'
#' `ComplexRegion(sets, polar = FALSE)` is the set of `x + y*I` with
#' `(x, y)` in `sets`, a [ProductSet()] of two real sets or a union of
#' them, as SymPy's `CartesianComplexRegion`. With `polar = TRUE` the pairs
#' are `(r, theta)` and the points `r*(cos(theta) + I*sin(theta))`, with
#' theta normalized into `[0, 2*pi)` (`PolarComplexRegion`).
#' `complex_region_from_real(s)` embeds a real set (SymPy's `from_real`);
#' `a_interval()` and `b_interval()` are the unions of the first and second
#' factors, and `psets()` lists the products. [contains()],
#' [Intersection()] and [Union()] accept regions.
#'
#' @param sets A product of two real sets, or a union of products.
#' @param polar Polar form.
#' @param s A real set.
#' @param region A complex region.
#' @param theta A set of angles.
#' @return A set.
#' @examples
#' a <- Interval(2, 3); b <- Interval(4, 6)
#' C1 <- ComplexRegion(ProductSet(a, b))
#' contains(C1, Rational(5, 2) + Rational(9, 2)*sym_I)
#' ComplexRegion(ProductSet(Interval(0, 1), Interval(0, 2*sym_pi)), polar = TRUE)
#' Intersection(ComplexRegion(ProductSet(Interval(1, 3), Interval(-1, 1))), Reals)
#' @name ComplexRegion
NULL

#' @rdname ComplexRegion
#' @export
ComplexRegion <- function(sets, polar = FALSE) {
  sets <- as_sym(sets)
  if (!isTRUE(polar)) {
    if (sets$type == "ProductSet" && length(sets$args) == 2L && all(vapply(sets$args, .reals_like, TRUE)))
      return(Complexes)
    # CartesianComplexRegion.__new__: a product of two finite sets is the
    # finite set of the points x + I*y
    if (length(sets$args) == 2L && all(vapply(sets$args, function(a) a$type == "FiniteSet", TRUE))) {
      pts <- list()
      for (a in sets$args[[1L]]$args) for (b in sets$args[[2L]]$args) pts[[length(pts) + 1L]] <- .add(list(a, .mul(list(S_I, b))))
      return(FiniteSet(pts))
    }
    return(new_sym("CartesianComplexRegion", list(sets)))
  }
  ps <- if (sets$type == "ProductSet") list(sets) else sets$args
  ps <- lapply(ps, function(v) ProductSet(v$args[[1L]], normalize_theta_set(v$args[[2L]])))
  new_sym("PolarComplexRegion", list(do.call(Union, ps)))
}

.is_cregion <- function(z) is_sym(z) && z$type %in% c("CartesianComplexRegion", "PolarComplexRegion")
.cr_polar <- function(z) z$type == "PolarComplexRegion"

#' @rdname ComplexRegion
#' @export
psets <- function(region) { s <- region$args[[1L]]; if (s$type == "ProductSet") list(s) else s$args }
#' @rdname ComplexRegion
#' @export
a_interval <- function(region) do.call(Union, lapply(psets(region), function(p) p$args[[1L]]))
#' @rdname ComplexRegion
#' @export
b_interval <- function(region) do.call(Union, lapply(psets(region), function(p) p$args[[2L]]))

#' @rdname ComplexRegion
#' @export
complex_region_from_real <- function(s) {
  s <- as_sym(s)
  if (!isTRUE(is_subset(s, Reals))) stop("sets must be a subset of the real line", call. = FALSE)
  new_sym("CartesianComplexRegion", list(ProductSet(s, FiniteSet(S_Zero))))
}

# SymPy's _pi_coeff for a numeric angle: c with angle = c*pi, reduced mod 2.
.theta_coeff <- function(v) {
  v <- as_sym(v)
  if (is_zero_num(v)) return(gmp::as.bigq(0))
  c <- .div(v, S_Pi)
  if (!(c$type %in% c("Integer", "Rational"))) return(NULL)
  q <- gmp::as.bigq(c$value)
  q - 2 * .bigq_floor(q / 2)
}

#' @rdname ComplexRegion
#' @export
normalize_theta_set <- function(theta) {
  theta <- as_sym(theta)
  tp <- .mul(list(S_Two, S_Pi))
  if (theta$type == "Interval") {
    i <- .iv(theta)
    len <- .sub(i$b, i$a)
    if (.is_true(Ge(len, tp))) {
      if (identical(simplify(.sub(len, tp))$key, S_Zero$key) && i$lo && i$ro) {
        k <- .theta_coeff(i$a)
        kp <- .mul(list(make_rational(k), S_Pi))
        return(Union(Interval(S_Zero, kp, FALSE, TRUE), Interval(kp, tp, TRUE, TRUE)))
      }
      return(Interval(S_Zero, tp, FALSE, TRUE))
    }
    ks <- .theta_coeff(i$a); ke <- .theta_coeff(i$b)
    if (is.null(ks) || is.null(ke)) stop("Normalizing theta without pi as coefficient is not yet implemented", call. = FALSE)
    ns <- .mul(list(make_rational(ks), S_Pi)); ne <- .mul(list(make_rational(ke), S_Pi))
    if (ks > ke) return(Union(Interval(S_Zero, ne, FALSE, i$ro), Interval(ns, tp, i$lo, TRUE)))
    return(Interval(ns, ne, i$lo, i$ro))
  }
  if (theta$type == "FiniteSet") {
    return(FiniteSet(lapply(theta$args, function(el) {
      k <- .theta_coeff(el)
      if (is.null(k)) stop("Normalizing theta without pi as coefficient, is not Implemented.", call. = FALSE)
      .mul(list(make_rational(k), S_Pi))
    })))
  }
  if (theta$type == "Union") return(do.call(Union, lapply(theta$args, normalize_theta_set)))
  if (isTRUE(is_subset(theta, Reals))) stop("Normalizing theta of this type is not implemented", call. = FALSE)
  stop(.str(theta), " is not a real set", call. = FALSE)
}

# fuzzy_or of fuzzy_and over the products: TRUE, FALSE or NA
.cr_fuzzy <- function(ps, test) {
  vals <- vapply(ps, function(p) {
    a <- test(p)
    if (all(a %in% TRUE)) TRUE else if (any(a %in% FALSE)) FALSE else NA
  }, NA)
  if (any(vals %in% TRUE)) TRUE else if (all(vals %in% FALSE)) FALSE else NA
}
.cr_in <- function(S, v) { r <- .rel_to_r(.contains(S, v)); if (is_sym(r)) NA else r }

.cr_contains <- function(R, other) {
  ps <- psets(R)
  tup <- other$type == "Tuple"
  if (tup && length(other$args) != 2L) stop("expecting Tuple of length 2", call. = FALSE)
  v <- if (!.cr_polar(R)) {
    ri <- if (tup) other$args else .re_im(other)
    .cr_fuzzy(ps, function(p) c(.cr_in(p$args[[1L]], ri[[1L]]), .cr_in(p$args[[2L]], ri[[2L]])))
  } else if (!tup && is_zero_num(other)) {
    .cr_fuzzy(ps, function(p) .cr_in(p$args[[1L]], S_Zero))
  } else {
    rt <- if (tup) other$args else list(fcall("Abs", list(other)), fcall("arg", list(other)))
    th <- rt[[2L]]
    if (!(isTRUE(ask_fact(th, "real")) && !length(free_symbols(th)))) NA
    else {
      tp <- .mul(list(S_Two, S_Pi))
      th <- .sub(th, .mul(list(tp, fcall("floor", list(.div(th, tp))))))
      .cr_fuzzy(ps, function(p) c(.cr_in(p$args[[1L]], rt[[1L]]), .cr_in(p$args[[2L]], th)))
    }
  }
  if (is.na(v)) new_sym("Contains", list(other, R)) else if (v) S_true else S_false
}

.cr_intersect <- function(A, B) {
  if (.is_cregion(B)) {
    if (!.cr_polar(A) && !.cr_polar(B)) return(ComplexRegion(Intersection(A$args[[1L]], B$args[[1L]])))
    if (.cr_polar(A) && .cr_polar(B)) {
      t1 <- b_interval(A); t2 <- b_interval(B)
      th <- Intersection(t1, t2)
      tp <- .mul(list(S_Two, S_Pi))
      if ((isTRUE(.cr_in(t1, tp)) && isTRUE(.cr_in(t2, S_Zero))) || (isTRUE(.cr_in(t2, tp)) && isTRUE(.cr_in(t1, S_Zero))))
        th <- Union(th, FiniteSet(S_Zero))
      return(ComplexRegion(ProductSet(Intersection(a_interval(A), a_interval(B)), th), polar = TRUE))
    }
    return(NULL)
  }
  if (!isTRUE(is_subset(B, Reals))) return(NULL)
  parts <- list()
  for (p in psets(A)) {
    if (isTRUE(.cr_in(p$args[[2L]], S_Zero))) parts[[length(parts) + 1L]] <- p$args[[1L]]
    if (.cr_polar(A)) {
      if (isTRUE(.cr_in(p$args[[2L]], S_Pi))) {
        xx <- Dummy("x", real = TRUE)
        parts[[length(parts) + 1L]] <- ImageSet(Lambda(xx, .neg(xx)), p$args[[1L]])
      }
      if (isTRUE(.cr_in(p$args[[1L]], S_Zero))) parts[[length(parts) + 1L]] <- FiniteSet(S_Zero)
    }
  }
  Intersection(do.call(Union, parts), B)
}

.cr_union <- function(A, B) {
  if (!.is_cregion(B) && isTRUE(is_subset(B, Reals))) B <- complex_region_from_real(B)
  if (!.is_cregion(B)) return(NULL)
  if (!.cr_polar(A) && !.cr_polar(B)) return(ComplexRegion(Union(A$args[[1L]], B$args[[1L]])))
  if (.cr_polar(A) && .cr_polar(B)) return(ComplexRegion(Union(A$args[[1L]], B$args[[1L]]), polar = TRUE))
  NULL
}

for (.t in c("CartesianComplexRegion", "PolarComplexRegion"))
  register_str_printer(.t, local({ tt <- .t; function(e) paste0(tt, "(", .str(e$args[[1L]]), ")") }))
rm(.t)
.latex_set_product <- function(s) {
  if (s$type == "ProductSet") return(paste(vapply(s$args, .latex, ""), collapse = " \\times "))
  if (s$type == "Union") return(paste(vapply(s$args, .latex_set_product, ""), collapse = " \\cup "))
  .latex(s)
}
register_latex_printer("CartesianComplexRegion", function(e)
  paste0("\\left\\{x + y i\\; \\middle|\\; x, y \\in ", .latex_set_product(e$args[[1L]]), " \\right\\}"))
register_latex_printer("PolarComplexRegion", function(e)
  paste0("\\left\\{r \\left(i \\sin{\\left(\\theta \\right)} + \\cos{\\left(\\theta \\right)}\\right)\\; \\middle|\\; r, \\theta \\in ",
         .latex_set_product(e$args[[1L]]), " \\right\\}"))
