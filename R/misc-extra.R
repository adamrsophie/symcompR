# B-splines (sympy/functions/special/bsplines.py) and small utilities:
# polarify, prod, reshape, parallel_poly_from_expr, sqf_norm, gff, inverse
# monomial orderings and the carmichael helpers.

.bs_ivl <- function(cond, x) {
  a <- cond$args[[1L]]; b <- cond$args[[2L]]
  lts <- function(r) if (r$type %in% c("GreaterThan", "StrictGreaterThan")) r$args[[2L]] else r$args[[1L]]
  gts <- function(r) if (r$type %in% c("GreaterThan", "StrictGreaterThan")) r$args[[1L]] else r$args[[2L]]
  if (identical(lts(a)$key, x$key)) { tmp <- a; a <- b; b <- tmp }
  list(lts(a), gts(b))
}

.pw_expand <- function(p) if (p$type == "Piecewise") .pw_build(lapply(.pw_pairs(p), function(pr) list(expand(pr[[1L]]), pr[[2L]]))) else expand(p)

.bs_add <- function(c, b1, d, b2, x) {
  if (is_zero_num(b1) || is_zero_num(c)) return(.pw_expand(piecewise_fold(.mul(list(d, b2)))))
  if (is_zero_num(b2) || is_zero_num(d)) return(.pw_expand(piecewise_fold(.mul(list(c, b1)))))
  p1 <- .pw_pairs(piecewise_fold(.mul(list(c, b1))))
  p2 <- .pw_pairs(piecewise_fold(.mul(list(d, b2))))
  p2args <- p2[-length(p2)]
  out <- list()
  for (arg in p1[-length(p1)]) {
    expr <- arg[[1L]]; cond <- arg[[2L]]
    lower <- .bs_ivl(cond, x)[[1L]]
    for (i in seq_along(p2args)) {
      c2 <- p2args[[i]][[2L]]
      iv2 <- .bs_ivl(c2, x)
      if (identical(c2$key, cond$key)) {
        expr <- .add(list(expr, p2args[[i]][[1L]]))
        p2args[[i]] <- NULL
        break
      } else if (.is_true(Lt(iv2[[1L]], lower)) && .is_true(Le(iv2[[2L]], lower))) {
        out[[length(out) + 1L]] <- p2args[[i]]
        p2args[[i]] <- NULL
        break
      }
    }
    out[[length(out) + 1L]] <- list(expr, cond)
  }
  out <- c(out, p2args, list(list(S_Zero, S_true)))
  .pw_expand(.pw_build(out))
}

.bspline <- function(d, knots, n, x) {
  nint <- length(knots) - 1L
  if (n + d + 1L > nint) stop("n + d + 1 must not exceed len(knots) - 1", call. = FALSE)
  if (d == 0L) return(.pw_build(list(list(S_One, .contains(Interval(knots[[n + 1L]], knots[[n + 2L]]), x)), list(S_Zero, S_true))))
  den <- .sub(knots[[n + d + 2L]], knots[[n + 2L]])
  if (!is_zero_num(den)) { B <- .div(.sub(knots[[n + d + 2L]], x), den); b2 <- .bspline(d - 1L, knots, n + 1L, x) }
  else { B <- S_Zero; b2 <- S_Zero }
  den <- .sub(knots[[n + d + 1L]], knots[[n + 1L]])
  if (!is_zero_num(den)) { A <- .div(.sub(x, knots[[n + 1L]]), den); b1 <- .bspline(d - 1L, knots, n, x) }
  else { A <- S_Zero; b1 <- S_Zero }
  .bs_add(A, b1, B, b2, x)
}

#' B-splines
#'
#' `bspline_basis(d, knots, n, x)` is the `n`-th B-spline of degree `d` on
#' the knots (zero based `n`, as in SymPy), a Piecewise in `x`.
#' `bspline_basis_set(d, knots, x)` lists all of them, and
#' `interpolating_spline(d, x, X, Y)` is the spline of degree `d` through
#' the points `(X, Y)`.
#'
#' @param d The degree.
#' @param knots The knots.
#' @param n The index, from 0.
#' @param x The variable.
#' @param X,Y Coordinates of the interpolation points.
#' @return A Piecewise, or a list of them.
#' @examples
#' x <- Symbol("x")
#' bspline_basis(1, c(0, 1, 2), 0, x)
#' interpolating_spline(3, x, c(-2, 0, 1, 3, 4), c(4, 2, 1, 1, 3))
#' @name bsplines
NULL

#' @rdname bsplines
#' @export
bspline_basis <- function(d, knots, n, x) {
  knots <- lapply(as.list(knots), as_sym)
  if (d < 0) stop("degree must be non-negative", call. = FALSE)
  # built on an assumption-free dummy, then x goes in, as SymPy does
  u <- Dummy("x")
  xreplace(.bspline(as.integer(d), knots, as.integer(n), u), list(list(u, as_sym(x))))
}

#' @rdname bsplines
#' @export
bspline_basis_set <- function(d, knots, x) {
  lapply(seq_len(length(knots) - d - 1L) - 1L, function(i) bspline_basis(d, knots, i, x))
}

#' @rdname bsplines
#' @export
interpolating_spline <- function(d, x, X, Y) {
  x <- as_sym(x)
  if (d < 1 || d != round(d)) stop("Spline degree must be a positive integer, not ", d, ".", call. = FALSE)
  if (length(X) != length(Y)) stop("Number of X and Y coordinates must be the same.", call. = FALSE)
  if (length(X) < d + 1) stop("Degree must be less than the number of control points.", call. = FALSE)
  X <- lapply(as.list(X), as_sym); Y <- lapply(as.list(Y), as_sym)
  if (length(X) > 1L && !all(vapply(seq_len(length(X) - 1L), function(i) .is_true(Lt(X[[i]], X[[i + 1L]])), TRUE)))
    stop("The x-coordinates must be strictly increasing.", call. = FALSE)
  m <- length(X)
  inner <- if (d %% 2 == 1) { j <- (d + 1) %/% 2; X[seq_len(m - 2L * j) + j] }
           else { j <- d %/% 2; lapply(seq_len(m - 2L * j - 1L), function(i) .div(.add(list(X[[i + j]], X[[i + j + 1L]])), S_Two)) }
  knots <- c(rep(X[1L], d + 1L), inner, rep(X[m], d + 1L))
  basis <- bspline_basis_set(d, knots, x)
  cs <- lapply(seq_len(m) - 1L, function(i) Dummy(paste0("c", i)))
  eqs <- lapply(seq_len(m), function(i) .sub(.add(Map(function(c, b) .mul(list(c, subs(b, x, X[[i]]))), cs, basis)), Y[[i]]))
  sol <- .solve_linear_system(eqs, cs)
  coef <- lapply(cs, function(c) sol[[c$name]])
  conds <- list()
  for (b in basis) for (pr in .pw_pairs(b)) if (!.is_true(pr[[2L]])) conds[[pr[[2L]]$key]] <- pr[[2L]]
  conds <- unname(conds)
  lo <- vapply(conds, function(c) num_value_dbl(N(.bs_ivl(c, x)[[1L]])), 0)
  conds <- conds[order(lo)]
  pieces <- lapply(conds, function(c) {
    val <- .add(Map(function(k, b) {
      e <- S_Zero
      for (pr in .pw_pairs(b)) if (identical(pr[[2L]]$key, c$key)) e <- pr[[1L]]
      .mul(list(k, e))
    }, coef, basis))
    list(expand(val), c)
  })
  .pw_build(pieces)
}

#' Small helpers
#'
#' `polarify(e, subs, lift)` replaces the symbols of `e` by polar dummies
#' and returns `list(expression, mapping)` (or the expression alone when
#' `subs = FALSE`). `sym_prod(l)` is the product of a list (SymPy's
#' `prod()`). `sym_reshape(seq, how)` (SymPy's `reshape()`, renamed so it does not mask [stats::reshape()]) groups a sequence by the counts in `how`.
#' `parallel_poly_from_expr(exprs)` builds polynomials over common
#' generators. `sqf_norm(f, extension)` returns `list(s, g, r)` with
#' `g = f(x - s*a)` for the algebraic number `a` and its square free norm
#' `r`. `gff()` is not implemented, as in SymPy. `ilex()`, `igrlex()` and
#' `igrevlex()` are the inverse monomial order keys. `carmichael` collects
#' SymPy's Carmichael helpers.
#'
#' @param e,f An expression.
#' @param subs Replace symbols by polar dummies.
#' @param lift Wrap numbers in `polar_lift()`.
#' @param l,seq A list.
#' @param how A vector of counts, possibly nested.
#' @param exprs Expressions.
#' @param extension An algebraic number.
#' @param m An exponent vector.
#' @param ... Generators.
#' @return Varies.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' polarify(exp(x)*y)
#' sym_reshape(1:6, 2)
#' sqf_norm(x^2 - 2, sqrt(Integer(3)))
#' @name misc-helpers
NULL

#' @rdname misc-helpers
#' @export
polarify <- function(e, subs = TRUE, lift = FALSE) {
  if (lift) subs <- FALSE
  e <- .polarify(as_sym(e), lift)
  if (!subs) return(e)
  syms <- free_symbols(e)
  ds <- lapply(syms, function(s) Dummy(s$name, polar = TRUE))
  # the mapping back, as substitution pairs list(dummy, symbol)
  list(xreplace(e, Map(list, syms, ds)), Map(list, ds, syms))
}

# SymPy's _polarify: numbers (and symbols when lifting) become polar_lift();
# the terms of a sum are left alone, but a lifted sum is lifted as a whole.
.polarify <- function(e, lift, pause = FALSE) {
  if (isTRUE(ask_fact(e, "polar"))) return(e)
  lp <- function(z) fcall("polar_lift", list(z))
  if (!pause && !length(free_symbols(e)) && is_sym(e) && !(e$type %in% c("BooleanAtom", "Tuple"))) return(lp(e))
  if (e$type %in% c("Symbol", "Dummy") && !pause && lift) return(lp(e))
  if (!length(e$args)) return(e)
  if (e$type == "Add") {
    r <- .add(lapply(e$args, .polarify, lift = lift, pause = TRUE))
    return(if (lift) lp(r) else r)
  }
  if (e$type == "Function") return(sym_rebuild(e, lapply(e$args, .polarify, lift = lift, pause = FALSE)))
  if (e$type == "Integral") {
    # the integration variables stay; the integrand and the limits are lifted
    lims <- lapply(e$args[-1L], function(t) do.call(Tuple, c(list(.polarify(t$args[[1L]], FALSE, pause)),
                                                          lapply(t$args[-1L], .polarify, lift = lift, pause = pause))))
    return(new_sym("Integral", c(list(.polarify(e$args[[1L]], lift, pause)), lims)))
  }
  sym_rebuild(e, lapply(e$args, .polarify, lift = lift, pause = pause))
}

#' @rdname misc-helpers
#' @export
sym_prod <- function(l) .mul(lapply(as.list(l), as_sym))

#' @rdname misc-helpers
#' @export
sym_reshape <- function(seq, how) {
  seq <- as.list(seq)
  m <- sum(unlist(how))
  if (any(unlist(how) < 0) || m <= 0 || length(seq) %% m != 0)
    stop("template must contain nonnegative counts with a positive sum that divides the sequence length", call. = FALSE)
  pos <- 0L
  # one group filled from seq; nested templates give nested lists, possibly empty
  fill <- function(tmpl) {
    grp <- list()
    for (hi in as.list(tmpl)) {
      if (is.list(hi)) grp[[length(grp) + 1L]] <- fill(hi)
      else { grp <- c(grp, seq[pos + seq_len(hi)]); pos <<- pos + hi }
    }
    grp
  }
  lapply(seq_len(length(seq) %/% m), function(k) fill(how))
}

#' @rdname misc-helpers
#' @export
parallel_poly_from_expr <- function(exprs, ...) {
  exprs <- lapply(exprs, as_sym)
  gens <- .gens_of(exprs, list(...))
  list(lapply(exprs, function(z) do.call(Poly, c(list(z), gens))), list(gens = gens))
}

#' @rdname misc-helpers
#' @export
sqf_norm <- function(f, extension) .ext_sqf_norm(f, extension)

# Rational coordinates of an expanded radical expression: monomial key ->
# list(monomial, coefficient).
.field_coords <- function(z) {
  out <- list()
  for (t in (function(e) if (e$type == "Add") e$args else list(e))(expand(z))) {
    sp <- as_coeff_mul_split(t)
    if (!(sp[[1L]]$type %in% c("Integer", "Rational"))) return(NULL)
    k <- sp[[2L]]$key
    out[[k]] <- if (is.null(out[[k]])) list(sp[[2L]], gmp::as.bigq(sp[[1L]]$value)) else list(sp[[2L]], out[[k]][[2L]] + gmp::as.bigq(sp[[1L]]$value))
  }
  out
}

# c as a polynomial p(y) of degree < n with c = p(a), by exact Gaussian
# elimination over the radical monomials of the powers of a; NULL if none.
.to_field_poly <- function(c, a, y, n) {
  # algebraic quotients are rationalized first: 1/(sqrt(2) + sqrt(3))
  c <- expand(c)
  if (.has_type(c, "Pow") && any(vapply(atoms(c, "Pow"), function(p) is_number_type(p$args[[2L]]) && num_sign(p$args[[2L]]) < 0 && p$args[[1L]]$type == "Add", TRUE)))
    c <- expand(radsimp(c))
  if (!length(free_symbols(c)) && is_exact_number(c)) return(c)
  vs <- lapply(seq_len(n) - 1L, function(i) .field_coords(.pow(a, .int(i))))
  vc <- .field_coords(c)
  if (is.null(vc) || any(vapply(vs, is.null, TRUE))) return(NULL)
  keys <- unique(c(unlist(lapply(vs, names)), names(vc)))
  getq <- function(v, k) if (is.null(v[[k]])) gmp::as.bigq(0) else v[[k]][[2L]]
  M <- lapply(keys, function(k) c(lapply(vs, getq, k = k), list(getq(vc, k))))
  rows <- length(M); r <- 1L; piv <- integer()
  for (j in seq_len(n)) {
    pr <- which(vapply(seq_len(rows), function(i) i >= r && M[[i]][[j]] != 0, TRUE))
    if (!length(pr)) next
    t <- M[[r]]; M[[r]] <- M[[pr[1L]]]; M[[pr[1L]]] <- t
    pv <- M[[r]][[j]]
    M[[r]] <- lapply(M[[r]], function(v) v / pv)
    for (i in seq_len(rows)) if (i != r && M[[i]][[j]] != 0) {
      f <- M[[i]][[j]]
      M[[i]] <- Map(function(u, w) u - f * w, M[[i]], M[[r]])
    }
    piv <- c(piv, j); r <- r + 1L
    if (r > rows) break
  }
  # consistency: no row 0 = nonzero
  for (i in seq_len(rows)) if (all(vapply(M[[i]][seq_len(n)], function(v) v == 0, TRUE)) && M[[i]][[n + 1L]] != 0) return(NULL)
  q <- rep(list(gmp::as.bigq(0)), n)
  for (k in seq_along(piv)) q[[piv[k]]] <- M[[k]][[n + 1L]]
  .add(lapply(seq_len(n), function(i) .mul(list(make_rational(q[[i]]), .pow(y, .int(i - 1L))))))
}

#' @rdname misc-helpers
#' @export
gff <- function(f, ...) stop("symbolic falling factorial", call. = FALSE)

#' @rdname misc-helpers
#' @export
ilex <- function(m) -m
#' @rdname misc-helpers
#' @export
igrlex <- function(m) c(-sum(m), -m)
#' @rdname misc-helpers
#' @export
igrevlex <- function(m) c(-sum(m), rev(m))

#' @rdname misc-helpers
#' @export
carmichael <- list(is_carmichael = function(n) is_carmichael(n),
                   find_carmichael_numbers_in_range = function(x, y) find_carmichael_numbers_in_range(x, y),
                   find_first_n_carmichaels = function(n) find_first_n_carmichaels(n))

#' Symbolic dictionaries
#'
#' `Dict(keys, values)` is SymPy's `Dict`, an immutable mapping between
#' expressions stored as `(key, value)` tuples. `keys` may also be a list of
#' `list(key, value)` pairs. A repeated key keeps its last value.
#' `dict_get()` looks a key up, and `dict_keys()`, `dict_values()` and
#' `dict_items()` list the contents in printing order.
#'
#' @param keys Keys, or a list of key and value pairs.
#' @param values Values, one per key.
#' @param d A `Dict`.
#' @param key A key.
#' @param default Returned when the key is absent.
#' @return A `Dict`, an expression, or a list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' d <- Dict(list(x, y), list(1, 2))
#' d
#' dict_get(d, x)
#' subs(Dict(list(x), list(y^2)), y, 3)
#' @export
Dict <- function(keys, values = NULL) {
  if (is.null(values)) {
    pairs <- as.list(keys)
    keys <- lapply(pairs, `[[`, 1L); values <- lapply(pairs, `[[`, 2L)
  }
  keys <- lapply(as.list(keys), as_sym); values <- lapply(as.list(values), as_sym)
  if (length(keys) != length(values)) stop("keys and values must have the same length", call. = FALSE)
  keep <- !duplicated(vapply(keys, .k, ""), fromLast = TRUE)
  items <- Map(function(k, v) Tuple(k, v), keys[keep], values[keep])
  new_sym("Dict", if (length(items)) .bool_print_order(items) else list())
}

.dict_sorted <- function(d) {
  ks <- lapply(d$args, function(t) t$args[[1L]])
  .key_sort(lapply(ks, .sort_key), d$args)
}

#' @rdname Dict
#' @export
dict_get <- function(d, key, default = NULL) {
  key <- as_sym(key)
  for (t in d$args) if (identical(t$args[[1L]]$key, key$key)) return(t$args[[2L]])
  default
}
#' @rdname Dict
#' @export
dict_keys <- function(d) lapply(.dict_sorted(d), function(t) t$args[[1L]])
#' @rdname Dict
#' @export
dict_values <- function(d) lapply(.dict_sorted(d), function(t) t$args[[2L]])
#' @rdname Dict
#' @export
dict_items <- function(d) .dict_sorted(d)

register_rebuilder("Dict", function(e, args) Dict(lapply(args, function(t) t$args)))
register_str_printer("Dict", function(e)
  paste0("{", paste(vapply(.dict_sorted(e), function(t) paste0(.str(t$args[[1L]]), ": ", .str(t$args[[2L]])), ""),
                    collapse = ", "), "}"))
register_latex_printer("Dict", function(e)
  paste0("\\left\\{ ", paste(vapply(.dict_sorted(e), function(t) paste0(.latex(t$args[[1L]]), " : ", .latex(t$args[[2L]])), ""),
                             collapse = ", \\  "), "\\right\\}"))
