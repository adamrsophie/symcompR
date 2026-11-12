# Accumulation bounds, after sympy/calculus/accumulationbounds.py: the
# closed real interval of values an expression takes, with interval
# arithmetic. The core Add, Mul and Pow fold them with numbers (see
# .ab_add_args and .ab_mul_args), and sin, cos, tan, cot, exp and log
# evaluate on them.

#' Accumulation bounds
#'
#' `AccumBounds(min, max)` (also `AccumulationBounds()`) is SymPy's closed
#' real interval `<min, max>` of the values an expression can take, for
#' example `sin(oo)` or `limit(sin(x), x, oo)`. Arithmetic with numbers and
#' other bounds is interval arithmetic: `AccumBounds(1, 2) + 3` is
#' `AccumBounds(4, 5)` and `AccumBounds(-1, 2)^2` is `AccumBounds(0, 4)`.
#' Comparisons decide when the bounds do, `contains()` tests membership and
#' `intersection()` intersects with bounds or a finite set.
#'
#' @param min,max The bounds; equal bounds give the value itself.
#' @return An `AccumulationBounds` expression, or a number.
#' @examples
#' AccumBounds(1, 2) + 3
#' AccumBounds(-1, 2) * AccumBounds(3, 4)
#' 1/AccumBounds(0, 2)
#' sin(AccumBounds(0, pi/2))
#' x <- Symbol("x")
#' limit(sin(x), x, Inf)
#' @export
AccumBounds <- function(min, max) {
  lo <- as_sym(min); hi <- as_sym(max)
  if (isFALSE(.ab_xreal(lo)) || isFALSE(.ab_xreal(hi))) stop("Only real AccumulationBounds are supported", call. = FALSE)
  if (identical(lo$key, hi$key)) return(hi)
  bad <- if (!length(free_symbols(lo)) && !length(free_symbols(hi))) .is_true(Lt(hi, lo))
         else isTRUE(ask_fact(.sub(hi, lo), "negative"))
  if (bad) stop("Lower limit should be smaller than upper limit", call. = FALSE)
  new_sym("AccumulationBounds", list(lo, hi))
}

#' @rdname AccumBounds
#' @export
AccumulationBounds <- AccumBounds

.is_ab <- function(z) is_sym(z) && z$type == "AccumulationBounds"
.ab_xreal <- function(z) {
  if (is_inf_const(z) && !is_const(z, "zoo")) return(TRUE)
  v <- ask_fact(z, "extended_real")
  if (!is.na(v)) v else ask_fact(z, "real")
}
# extended-real sign facts, with the infinities included; a constant's
# sign is decided numerically when its facts do not settle it
.ab_f <- function(z, f) {
  if (is_const(z, "oo")) return(f %in% c("positive", "nonnegative", "nonzero"))
  if (is_const(z, "-oo")) return(f %in% c("negative", "nonpositive", "nonzero"))
  v <- ask_fact(z, f)
  if (!is.na(v)) return(isTRUE(v))
  if (!length(free_symbols(z)) && f %in% c("positive", "negative", "nonnegative", "nonpositive")) {
    r <- switch(f, positive = Gt(z, S_Zero), negative = Lt(z, S_Zero), nonnegative = Ge(z, S_Zero), nonpositive = Le(z, S_Zero))
    return(.is_true(r))
  }
  FALSE
}
.ab_zero <- function(z) is_zero_num(z) || (z$type == "Float" && z$value == 0) || isTRUE(ask_fact(z, "zero"))
.ab_oo <- function() new_sym("AccumulationBounds", list(S_NegativeInfinity, S_Infinity))
.ab_lo <- function(z) z$args[[1L]]
.ab_hi <- function(z) z$args[[2L]]

# AccumBounds.__add__ (b is bounds or an expression)
.ab_add <- function(a, b) {
  if (.is_ab(b)) return(AccumBounds(.add(list(.ab_lo(a), .ab_lo(b))), .add(list(.ab_hi(a), .ab_hi(b)))))
  if ((is_const(b, "oo") && is_const(.ab_lo(a), "-oo")) || (is_const(b, "-oo") && is_const(.ab_hi(a), "oo"))) return(.ab_oo())
  if (isTRUE(.ab_xreal(b))) {
    if (is_const(.ab_lo(a), "-oo") && is_const(.ab_hi(a), "oo")) return(.ab_oo())
    if (is_const(.ab_lo(a), "-oo")) return(AccumBounds(S_NegativeInfinity, .add(list(.ab_hi(a), b))))
    if (is_const(.ab_hi(a), "oo")) return(AccumBounds(.add(list(.ab_lo(a), b)), S_Infinity))
    return(AccumBounds(.add(list(.ab_lo(a), b)), .add(list(.ab_hi(a), b))))
  }
  NULL
}

# AccumBounds.__mul__
.ab_mul <- function(a, b) {
  if (identical(a$key, .ab_oo()$key)) return(a)
  if (.is_ab(b)) {
    if (identical(b$key, .ab_oo()$key)) return(b)
    v <- list()
    for (e in a$args) {
      vi <- .ab_mul(b, e)
      if (is.null(vi)) return(NULL)
      v <- c(v, if (.is_ab(vi)) vi$args else list(vi))
    }
    return(AccumBounds(do.call(Min, v), do.call(Max, v)))
  }
  if (is_const(b, "oo")) {
    if (.ab_zero(.ab_lo(a))) return(AccumBounds(S_Zero, S_Infinity))
    if (.ab_zero(.ab_hi(a))) return(AccumBounds(S_NegativeInfinity, S_Zero))
  }
  if (is_const(b, "-oo")) {
    if (.ab_zero(.ab_lo(a))) return(AccumBounds(S_NegativeInfinity, S_Zero))
    if (.ab_zero(.ab_hi(a))) return(AccumBounds(S_Zero, S_Infinity))
  }
  if (isTRUE(.ab_xreal(b))) {
    if (.ab_zero(b)) {
      if (is_const(.ab_hi(a), "oo")) return(AccumBounds(S_Zero, S_Infinity))
      if (is_const(.ab_lo(a), "-oo")) return(AccumBounds(S_NegativeInfinity, S_Zero))
      return(S_Zero)
    }
    if (.ab_f(b, "positive")) return(AccumBounds(.mul(list(.ab_lo(a), b)), .mul(list(.ab_hi(a), b))))
    if (.ab_f(b, "negative")) return(AccumBounds(.mul(list(.ab_hi(a), b)), .mul(list(.ab_lo(a), b))))
  }
  NULL
}

# AccumBounds.__rtruediv__: other / a
.ab_rdiv <- function(a, other) {
  lo <- .ab_lo(a); hi <- .ab_hi(a)
  if (!isTRUE(.ab_xreal(other))) return(NULL)
  if (.ab_zero(other)) return(S_Zero)
  if (.ab_f(lo, "nonpositive") && .ab_f(hi, "nonnegative")) {
    if (.ab_zero(lo)) {
      if (.ab_f(other, "positive")) return(AccumBounds(.mul(list(other, .pow(hi, S_NegativeOne))), S_Infinity))
      if (.ab_f(other, "negative")) return(AccumBounds(S_NegativeInfinity, .mul(list(other, .pow(hi, S_NegativeOne)))))
    }
    if (.ab_zero(hi)) {
      if (.ab_f(other, "positive")) return(AccumBounds(S_NegativeInfinity, .mul(list(other, .pow(lo, S_NegativeOne)))))
      if (.ab_f(other, "negative")) return(AccumBounds(.mul(list(other, .pow(lo, S_NegativeOne))), S_Infinity))
    }
    return(.ab_oo())
  }
  q1 <- .mul(list(other, .pow(lo, S_NegativeOne))); q2 <- .mul(list(other, .pow(hi, S_NegativeOne)))
  AccumBounds(Min(q1, q2), Max(q1, q2))
}

# AccumBounds.__truediv__ (a / b with a bounds) and __rtruediv__ (b bounds)
.ab_div <- function(a, b) {
  if (!.is_ab(a)) return(.ab_rdiv(b, a))
  lo <- .ab_lo(a); hi <- .ab_hi(a)
  if (.is_ab(b)) {
    blo <- .ab_lo(b); bhi <- .ab_hi(b)
    if (.ab_f(blo, "positive") || .ab_f(bhi, "negative"))
      return(.ab_mul(a, AccumBounds(.pow(bhi, S_NegativeOne), .pow(blo, S_NegativeOne))))
    if (.ab_f(lo, "nonpositive") && .ab_f(hi, "nonnegative") && .ab_f(blo, "nonpositive") && .ab_f(bhi, "nonnegative")) {
      if (.ab_zero(lo) && .ab_zero(blo)) return(AccumBounds(S_Zero, S_Infinity))
      if (.ab_zero(hi) && .ab_zero(blo)) return(AccumBounds(S_NegativeInfinity, S_Zero))
      return(.ab_oo())
    }
    if (.ab_f(hi, "negative")) {
      if (.ab_f(blo, "negative")) {
        if (.ab_zero(bhi)) return(AccumBounds(.div(hi, blo), S_Infinity))
        if (.ab_f(bhi, "positive")) return(.ab_oo())
      }
      if (.ab_zero(blo) && .ab_f(bhi, "positive")) return(AccumBounds(S_NegativeInfinity, .div(hi, bhi)))
    }
    if (.ab_f(lo, "positive")) {
      if (.ab_f(blo, "negative")) {
        if (.ab_zero(bhi)) return(AccumBounds(S_NegativeInfinity, .div(lo, blo)))
        if (.ab_f(bhi, "positive")) return(.ab_oo())
      }
      if (.ab_zero(blo) && .ab_f(bhi, "positive")) return(AccumBounds(.div(lo, bhi), S_Infinity))
    }
  } else if (isTRUE(.ab_xreal(b))) {
    if (is_const(b, "oo") || is_const(b, "-oo")) {
      if (identical(a$key, .ab_oo()$key)) return(.ab_oo())
      if (is_const(hi, "oo")) return(AccumBounds(Min(S_Zero, b), Max(S_Zero, b)))
      if (is_const(lo, "-oo")) return(AccumBounds(Min(S_Zero, .neg(b)), Max(S_Zero, .neg(b))))
    }
    if (.ab_f(b, "positive")) return(AccumBounds(.div(lo, b), .div(hi, b)))
    if (.ab_f(b, "negative")) return(AccumBounds(.div(hi, b), .div(lo, b)))
  }
  NULL
}

# AccumBounds.__pow__
.ab_pow <- function(a, e) {
  lo <- .ab_lo(a); hi <- .ab_hi(a)
  raw <- function() new_sym("Pow", list(a, e))
  gt <- function(u, v) .is_true(Gt(u, v)); lt <- function(u, v) .is_true(Lt(u, v))
  if (is_const(e, "oo")) {
    if (.ab_f(lo, "nonnegative")) {
      if (lt(hi, S_One)) return(S_Zero)
      if (gt(lo, S_One)) return(S_Infinity)
      return(AccumBounds(S_Zero, S_Infinity))
    }
    if (.ab_f(hi, "negative")) {
      if (gt(lo, S_NegativeOne)) return(S_Zero)
      if (lt(hi, S_NegativeOne)) return(S_ComplexInfinity)
      return(S_NaN)
    }
    if (gt(lo, S_NegativeOne)) return(if (lt(hi, S_One)) S_Zero else AccumBounds(S_Zero, S_Infinity))
    return(.ab_oo())
  }
  if (is_const(e, "-oo")) return(.pow(.ab_rdiv(a, S_One), S_Infinity))
  if (.ab_f(lo, "nonnegative") && .ab_f(e, "nonnegative"))
    return(AccumBounds(.pow(lo, e), .pow(hi, e)))
  if (.ab_zero(e)) return(S_One)
  if (e$type == "Integer" || isTRUE(ask_fact(e, "integer"))) {
    p1 <- .pow(lo, e); p2 <- .pow(hi, e)
    if (.ab_f(lo, "positive") || .ab_f(hi, "negative")) return(AccumBounds(Min(p1, p2), Max(p1, p2)))
    if (e$type == "Integer") {
      neg <- .ab_f(e, "negative")
      if (.Z(e) %% 2 == 0) {
        if (neg) {
          if (.ab_zero(lo)) return(AccumBounds(p2, S_Infinity))
          if (.ab_zero(hi)) return(AccumBounds(p1, S_Infinity))
          return(.pow(.ab_rdiv(a, S_One), .neg(e)))
        }
        return(AccumBounds(S_Zero, Max(p1, p2)))
      }
      if (neg) {
        if (.ab_zero(lo)) return(AccumBounds(p2, S_Infinity))
        if (.ab_zero(hi)) return(AccumBounds(S_NegativeInfinity, p1))
        return(.pow(.ab_rdiv(a, S_One), .neg(e)))
      }
      return(AccumBounds(p1, p2))
    }
  }
  if (e$type == "Rational" && .ab_f(lo, "nonnegative")) {
    num <- .zs(gmp::numerator(e$value)); den <- .zs(gmp::denominator(e$value))
    if (is_one_num(num)) return(AccumBounds(.pow(lo, .pow(den, S_NegativeOne)), .pow(hi, .pow(den, S_NegativeOne))))
    return(.pow(.pow(a, num), .pow(den, S_NegativeOne)))
  }
  if (.is_ab(e) && (.ab_f(lo, "positive") || (.ab_f(lo, "nonnegative") && .ab_f(.ab_lo(e), "nonnegative")))) {
    p <- lapply(e$args, function(i) .pow(a, i))
    if (!any(vapply(p, function(z) z$type == "Pow", TRUE))) {
      v <- unlist(lapply(p, function(z) if (.is_ab(z)) z$args else list(z)), recursive = FALSE)
      return(AccumBounds(do.call(Min, v), do.call(Max, v)))
    }
  }
  raw()
}

# AccumBounds.__rpow__: b^a for a base b that is not bounds
.ab_rpow <- function(b, a) {
  if (!(isTRUE(ask_fact(b, "real")) && isTRUE(ask_fact(b, "nonnegative")) &&
        isTRUE(ask_fact(.sub(.ab_hi(a), .ab_lo(a)), "positive")))) return(NULL)
  if (is_one_num(b)) return(S_One)
  if (isTRUE(ask_fact(b, "positive"))) {
    u <- .pow(b, .ab_lo(a)); v <- .pow(b, .ab_hi(a))
    if (.is_true(Lt(v, u))) { t <- u; u <- v; v <- t }
    return(AccumBounds(u, v))
  }
  if (.ab_zero(b)) {
    if (.ab_zero(.ab_lo(a))) return(AccumBounds(S_Zero, S_One))
    if (.ab_f(.ab_lo(a), "positive")) return(S_Zero)
  }
  NULL
}

# Fold bounds with the terms of a sum or the factors of a product that
# interval arithmetic accepts (real terms; factors of known sign), as the
# operators of SymPy's AccumBounds do; the others stay as terms beside it.
.ab_fold <- function(args, op) {
  isab <- vapply(args, .is_ab, TRUE)
  acc <- args[[which(isab)[1L]]]
  rest <- list()
  for (z in c(args[isab][-1L], args[!isab])) {
    if (!.is_ab(acc)) {
      # the bounds collapsed to a value: the rest folds through the core
      acc <- if (op == "add") .add(list(acc, z)) else .mul(list(acc, z))
      next
    }
    # Mul.flatten and Add.flatten fold only Numbers (and other bounds) into
    # the bounds; other terms stay beside them
    if (!(.is_ab(z) || is_number_type(z) || is_inf_const(z) || is_nan(z))) { rest[[length(rest) + 1L]] <- z; next }
    r <- if (op == "add") .ab_add(acc, z) else .ab_mul(acc, z)
    if (is.null(r)) rest[[length(rest) + 1L]] <- z else acc <- r
  }
  list(acc = acc, rest = rest)
}

# Interval comparison: each side is [lo, hi], a single value for a non-bound.
.ab_rel <- function(type, a, b) {
  lo <- function(z) if (.is_ab(z)) .ab_lo(z) else z
  hi <- function(z) if (.is_ab(z)) .ab_hi(z) else z
  le <- function() if (.is_true(Le(hi(a), lo(b)))) TRUE else if (.is_true(Gt(lo(a), hi(b)))) FALSE else NA
  ge <- function() if (.is_true(Ge(lo(a), hi(b)))) TRUE else if (.is_true(Lt(hi(a), lo(b)))) FALSE else NA
  v <- switch(type, LessThan = le(), GreaterThan = ge(), StrictLessThan = !ge(), StrictGreaterThan = !le())
  if (is.na(v)) NULL else if (v) S_true else S_false
}

# sin, cos, tan, cot, exp and log on bounds, and the trigonometric
# functions at the infinities.
.AB_FUNS <- c("sin", "cos", "tan", "cot", "exp", "log")
.ab_fcall <- function(name, a) {
  if (!.is_ab(a)) {
    if (name %in% c("sin", "cos")) return(AccumBounds(S_NegativeOne, S_One))
    if (name %in% c("tan", "cot")) return(.ab_oo())
    return(NULL)
  }
  lo <- .ab_lo(a); hi <- .ab_hi(a)
  hits <- function(l, h, pts) {
    r <- vapply(pts, function(p) { v <- And(Le(l, p), Ge(h, p)); if (.is_true(v)) 1L else if (.is_false(v)) 0L else NA_integer_ }, 0L)
    if (anyNA(r)) NA else any(r == 1L)
  }
  shifted <- function(per) {
    d <- fcall("floor", list(.div(lo, per)))
    l <- if (is_const(lo, "-oo")) lo else .sub(lo, .mul(list(d, per)))
    h <- if (is_const(hi, "oo")) hi else .sub(hi, .mul(list(d, per)))
    list(l, h)
  }
  switch(name,
    exp = AccumBounds(fcall("exp", list(lo)), fcall("exp", list(hi))),
    log = if (.ab_f(lo, "positive")) AccumBounds(fcall("log", list(lo)), fcall("log", list(hi)))
          else if (.ab_zero(lo)) AccumBounds(S_NegativeInfinity, fcall("log", list(hi))) else S_NaN,
    sin = {
      s <- shifted(.mul(list(S_Two, S_Pi)))
      top <- hits(s[[1L]], s[[2L]], list(.div(S_Pi, S_Two), .mul(list(Rational(5, 2), S_Pi))))
      bot <- hits(s[[1L]], s[[2L]], list(.mul(list(Rational(3, 2), S_Pi)), .mul(list(Rational(7, 2), S_Pi))))
      bot2 <- hits(s[[1L]], s[[2L]], list(.mul(list(Rational(3, 2), S_Pi)), .mul(list(.int(4), S_Pi))))
      if (is.na(top) || is.na(bot) || is.na(bot2)) return(NULL)
      sl <- fcall("sin", list(s[[1L]])); sh <- fcall("sin", list(s[[2L]]))
      if (top && bot) AccumBounds(S_NegativeOne, S_One)
      else if (top) AccumBounds(Min(sl, sh), S_One)
      else if (bot2) AccumBounds(S_NegativeOne, Max(sl, sh))
      else AccumBounds(Min(sl, sh), Max(sl, sh))
    },
    cos = fcall("sin", list(.ab_add(a, .div(S_Pi, S_Two)))),
    tan = {
      s <- shifted(S_Pi)
      h <- hits(s[[1L]], s[[2L]], list(.div(S_Pi, S_Two), .mul(list(Rational(3, 2), S_Pi))))
      if (is.na(h)) return(NULL)
      if (h) .ab_oo() else AccumBounds(fcall("tan", list(s[[1L]])), fcall("tan", list(s[[2L]])))
    },
    cot = .neg(fcall("tan", list(.ab_add(a, .div(S_Pi, S_Two))))),
    NULL)
}

.ab_intersection <- function(a, b) {
  if (!.is_ab(a)) stop("Input must be AccumulationBounds", call. = FALSE)
  if (is_sym(b) && b$type == "FiniteSet") return(FiniteSet(Filter(function(i) isTRUE(.ab_contains(a, i)), b$args)))
  if (!.is_ab(b)) stop("Input must be AccumulationBounds or FiniteSet object", call. = FALSE)
  if (.is_true(Lt(.ab_hi(a), .ab_lo(b))) || .is_true(Gt(.ab_lo(a), .ab_hi(b)))) return(S_EmptySet)
  if (.is_true(Le(.ab_lo(a), .ab_lo(b)))) {
    if (.is_true(Le(.ab_hi(a), .ab_hi(b)))) return(AccumBounds(.ab_lo(b), .ab_hi(a)))
    if (.is_true(Gt(.ab_hi(a), .ab_hi(b)))) return(b)
  }
  if (.is_true(Le(.ab_lo(b), .ab_lo(a)))) {
    if (.is_true(Le(.ab_hi(b), .ab_hi(a)))) return(AccumBounds(.ab_lo(a), .ab_hi(b)))
    if (.is_true(Gt(.ab_hi(b), .ab_hi(a)))) return(a)
  }
  stop("cannot decide the intersection", call. = FALSE)
}

#' @rdname AccumBounds
#' @param a,b Bounds; `ab_union()` is SymPy's `AccumBounds.union()` (base R
#'   has `union()`), and [intersection()] takes bounds or a finite set.
#' @export
ab_union <- function(a, b) {
  if (!.is_ab(a) || !.is_ab(b)) stop("Input must be AccumulationBounds object", call. = FALSE)
  AccumBounds(Min(.ab_lo(a), .ab_lo(b)), Max(.ab_hi(a), .ab_hi(b)))
}

.ab_contains <- function(a, x) {
  x <- as_sym(x)
  if (is_const(x, "oo") || is_const(x, "-oo")) return(is_const(.ab_lo(a), "-oo") || is_const(.ab_hi(a), "oo"))
  v <- And(Le(.ab_lo(a), x), Ge(.ab_hi(a), x))
  if (!(v$type == "BooleanAtom")) stop("input failed to evaluate", call. = FALSE)
  .is_true(v)
}

register_str_printer("AccumulationBounds", function(e) paste0("AccumBounds(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))
register_latex_printer("AccumulationBounds", function(e) paste0("\\left\\langle ", .latex(e$args[[1L]]), ", ", .latex(e$args[[2L]]), "\\right\\rangle"))

register_rebuilder("AccumulationBounds", function(e, args) AccumBounds(args[[1L]], args[[2L]]))
