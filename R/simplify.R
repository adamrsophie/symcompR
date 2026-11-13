# Simplification, after sympy/simplify.

.measure <- function(e) .cops(e)

.shortest <- function(cands) {
  cands <- Filter(Negate(is.null), cands)
  ops <- vapply(cands, .measure, 0L)
  len <- vapply(cands, function(z) nchar(z$key), 0L)
  cands[[order(ops, len)[1L]]]
}


# ---------------------------------------------------------------------------
# Powers and logarithms


# ---------------------------------------------------------------------------
# Radicals


#' @rdname radsimp
#' @export
sqrtdenest <- function(e) {
  e <- as_sym(e)
  .map_bottom_up(e, function(z) {
    if (!(z$type == "Pow" && is_const_num_exp(z, gmp::as.bigq(1, 2)))) return(z)
    inner <- z$args[[1L]]
    if (inner$type != "Add" || length(inner$args) != 2L) return(z)
    a <- Filter(is_exact_number, inner$args)
    rad <- Filter(Negate(is_exact_number), inner$args)
    if (length(a) != 1L || length(rad) != 1L) return(z)
    a <- a[[1L]]
    sp <- as_coeff_mul_split(rad[[1L]])
    b <- sp[[1L]]
    r <- sp[[2L]]
    if (!(r$type == "Pow" && is_exact_number(r$args[[1L]]) && is_const_num_exp(r, gmp::as.bigq(1, 2)))) return(z)
    c <- r$args[[1L]]
    d2 <- .sub(.sq(a), .mul(list(.sq(b), c)))
    if (!is_exact_number(d2) || d2$value < 0 || a$value <= 0) return(z)
    d <- .sqrt(d2)
    if (!is_exact_number(d)) return(z)
    p1 <- .sqrt(.div(.add(list(a, d)), S_Two))
    p2 <- .sqrt(.div(.sub(a, d), S_Two))
    if (num_sign(b) < 0) p2 <- .neg(p2)
    .add(list(p1, p2))
  })
}

is_const_num_exp <- function(z, q) z$args[[2L]]$type == "Rational" && z$args[[2L]]$value == q

# ---------------------------------------------------------------------------
# Collection and rational simplification


#' Separate variables
#'
#' `separatevars(e)` factors `e` into a product of factors that each depend
#' on fewer symbols, when possible.
#'
#' @param e An expression.
#' @param syms Optional symbols; with `dict = TRUE` returns a named list.
#' @param dict Return a list of factors per symbol.
#' @return An expression or a list.
#' @export
separatevars <- function(e, syms = NULL, dict = FALSE) {
  e <- as_sym(e)
  f <- .separatevars(e)
  if (!dict) return(f)
  syms <- if (is.null(syms)) free_symbols(e) else .as_gens(list(syms))
  facs <- if (f$type == "Mul") f$args else list(f)
  out <- list(coeff = S_One)
  for (s in syms) out[[s$name]] <- S_One
  for (fc in facs) {
    dep <- Filter(function(s) .has1(fc, s), syms)
    if (length(dep) > 1L) return(NULL)
    if (!length(dep)) out$coeff <- .mul(list(out$coeff, fc))
    else out[[dep[[1L]]$name]] <- .mul(list(out[[dep[[1L]]$name]], fc))
  }
  out
}

# SymPy _separatevars: keep a Mul and separate its factors, separate the
# base of a power and expand it, else factor.
.separatevars <- function(e) {
  if (length(free_symbols(e)) < 2L) return(e)
  if (e$type == "Mul") return(.mul(lapply(e$args, .separatevars)))
  if (e$type == "Pow" && !is_const(e$args[[1L]], "E"))
    e <- expand_power_base(.pow(.separatevars(e$args[[1L]]), e$args[[2L]]))
  factor_sym(e)
}

#' Find simple exact forms of numbers
#'
#' `nsimplify(x)` replaces floating point numbers by nearby rationals,
#' square roots of rationals, or rational multiples of constants such as `pi`
#' and `E` (or of the constants listed in `constants`).
#'
#' @param e A number or expression.
#' @param constants A list of constants to try.
#' @param tolerance Relative tolerance.
#' @param rational If `TRUE`, only produce rationals.
#' @return An expression.
#' @examples
#' nsimplify(0.333333333333)
#' nsimplify(1.414213562373095)
#' nsimplify(6.283185307179586, constants = list(sym_pi))
#' @export
nsimplify <- function(e, constants = list(), tolerance = 1e-10, rational = FALSE) {
  e <- as_sym(e)
  .map_bottom_up(e, function(z) {
    if (z$type != "Float") return(z)
    if (.float_prec(z) > 53 || !is.finite(z$value) || (z$value == 0 && .float_mp(z)$m != 0))
      return(make_rational(.decimal_to_bigq(.fmt_float_sym(z, strip = TRUE))))
    # _real_to_rational: nsimplify(fl, rational=False) first, then a base 10
    # conversion when that is not a Rational (and rational = TRUE); a nonzero
    # float identified as zero is not identified (nsimplify_real raises)
    r <- .nsimplify_num(z$value, constants, tolerance, FALSE)
    if (r$type == "Float" || (is_zero_num(r) && z$value != 0) ||
        (isTRUE(rational) && !r$type %in% c("Integer", "Rational")))
      r <- make_rational(.decimal_to_bigq(.fmt_float_sym(z, strip = TRUE)))
    # _real_to_rational: a nonzero float never becomes zero
    if (is_zero_num(r) && .float_mp(z)$m != 0) r <- make_rational(.mp_to_q(.float_mp(z)))
    r
  })
}

.nsimplify_num <- function(v, constants, tol, rational) {
  if (v == 0) return(S_Zero)
  try_rat <- function(w, maxden = 10000) {
    r <- .best_rational(w, maxden)
    if (!is.null(r) && abs(as.double(r) - w) <= tol * max(1, abs(w))) return(make_rational(r))
    NULL
  }
  r <- try_rat(v)
  if (!is.null(r) || rational) return(if (is.null(r)) make_rational(gmp::as.bigq(v)) else r)
  for (cst in c(constants, list(S_Pi, S_E))) {
    cv <- .neval(as_sym(cst), list())
    r <- try_rat(v / cv, 1000)
    if (!is.null(r)) return(.mul(list(r, as_sym(cst))))
  }
  r2 <- try_rat(v^2, 1000)
  if (!is.null(r2)) return(.mul(list(if (v < 0) S_NegativeOne else S_One, .sqrt(r2))))
  for (k in 2:12) {
    r <- try_rat(v - sqrt(k), 1000)
    if (!is.null(r)) return(.add(list(r, .sqrt(.int(k)))))
    r <- try_rat(v + sqrt(k), 1000)
    if (!is.null(r)) return(.sub(r, .sqrt(.int(k))))
  }
  Float(v)
}

.best_rational <- function(w, maxden) {
  # continued fraction convergents
  a <- floor(w)
  h0 <- 1; h1 <- a
  k0 <- 0; k1 <- 1
  frac <- w - a
  best <- gmp::as.bigq(a)
  for (i in 1:40) {
    if (abs(frac) < 1e-15) break
    w <- 1 / frac
    a <- floor(w)
    frac <- w - a
    h2 <- a * h1 + h0
    k2 <- a * k1 + k0
    if (k2 > maxden) break
    best <- gmp::as.bigq(h2, k2)
    h0 <- h1; h1 <- h2
    k0 <- k1; k1 <- k2
  }
  best
}

#' Simplify signs
#'
#' `signsimp()` takes the sign out of sums inside products, as SymPy's
#' `signsimp()` (`-x*(-x + y) + y` becomes `x*(x - y) + y`), and `posify()`
#' replaces symbols with positive ones.
#'
#' @param e An expression.
#' @return An expression; `posify()` returns `list(expr, reps)`.
#' @export
signsimp <- function(e) {
  e <- as_sym(e)
  if (!length(e$args)) return(e)
  # sub_pre: a sum from which -1 can be extracted becomes -1*(-sum), and
  # re-evaluating the products folds that -1 into their coefficients
  adds <- Filter(could_extract_minus_sign, atoms(e, "Add"))
  reps <- list()
  for (a in adds) {
    na <- .neg(a)
    if (na$type == "Mul") next
    reps[[length(reps) + 1L]] <- list(a, new_sym("Mul", list(S_NegativeOne, na)))
  }
  if (length(reps)) e <- xreplace(e, reps)
  e <- .map_bottom_up(e, function(z) if (z$type == "Mul") .mul(z$args) else z)
  if (e$type == "Add") return(.add(lapply(e$args, signsimp)))
  e
}

#' @rdname signsimp
#' @export
posify <- function(e) {
  e <- as_sym(e)
  fs <- Filter(function(s) is.na(ask_fact(s, "positive")), free_symbols(e))
  reps <- lapply(fs, function(s) list(s, Dummy(s$name, positive = TRUE)))
  list(expr = xreplace(e, reps), reps = lapply(reps, rev))
}

