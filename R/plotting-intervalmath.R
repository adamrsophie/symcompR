# Interval arithmetic for plotting, after sympy/plotting/intervalmath:
# interval_arithmetic.py, interval_membership.py and lib_interval.py. An
# interval holds double bounds and a three valued validity flag (TRUE,
# FALSE, or NA for SymPy's None). The elementary functions of lib_interval
# are R's Math group generics on the class; comparisons return an
# intervalMembership pair of three valued logic values.

.pi_iv <- function(a, b = a, is_valid = TRUE) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (!(a < b)) { t <- a; a <- b; b <- t }
  structure(list(start = a, end = b, is_valid = is_valid), class = "sym_interval")
}
.pi_iv_all <- function(v = FALSE) .pi_iv(-Inf, Inf, is_valid = v)
.is_pi_iv <- function(x) inherits(x, "sym_interval")
.fuzzy_and2 <- function(a, b) if (isFALSE(a) || isFALSE(b)) FALSE else if (is.na(a) || is.na(b)) NA else TRUE
.fuzzy_or2 <- function(a, b) if (isTRUE(a) || isTRUE(b)) TRUE else if (is.na(a) || is.na(b)) NA else FALSE

#' Interval arithmetic for plotting
#'
#' SymPy's `sympy.plotting.intervalmath`. `interval(a, b)` is an interval of
#' doubles with a validity flag: `TRUE` when a function result is continuous
#' and in the domain, `FALSE` when the argument is outside the domain, and
#' `NA` (SymPy's `None`) when it is partly outside or the function is not
#' continuous there. Arithmetic, powers and the elementary functions
#' (`exp`, `log`, `log10`, `sqrt`, `abs`, `sin`, `cos`, `tan`, `asin`,
#' `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`,
#' `floor`, `ceiling`) work on intervals; `imin()`, `imax()` and `ceil()`
#' are SymPy's names. Comparisons return an `intervalMembership` pair
#' (truth, validity), combined with `&`, `|`, `xor()` and `!`.
#'
#' @param a,b The bounds; one value gives a point interval.
#' @param is_valid The validity flag.
#' @param x An interval or a number.
#' @param ... Intervals or numbers.
#' @return An interval, or an `intervalMembership`.
#' @examples
#' i <- interval(1, 2)
#' i + 3
#' i * interval(-1, 2)
#' sqrt(interval(-1, 4))
#' sin(interval(0, 3))
#' imax(interval(0, 1), interval(-1, 3))
#' i > 0
#' i$mid
#' @name intervalmath
#' @export
interval <- function(a, b = NULL, is_valid = TRUE) {
  if (.is_pi_iv(a)) return(.pi_iv(a$start, a$end, is_valid))
  if (is.null(b)) .pi_iv(a, a, is_valid) else .pi_iv(a, b, is_valid)
}

#' @export
`$.sym_interval` <- function(x, name) {
  x <- unclass(x)
  switch(name, mid = (x$start + x$end) / 2, width = x$end - x$start, x[[name]])
}
.pi_ivf <- function(v) if (is.infinite(v)) (if (v > 0) "inf" else "-inf") else sprintf("%f", v)
#' @export
format.sym_interval <- function(x, ...) sprintf("interval(%s, %s)", .pi_ivf(x$start), .pi_ivf(x$end))
#' @export
print.sym_interval <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
#' @export
as.character.sym_interval <- function(x, ...) sprintf("[%s, %s]", .pi_ivf(x$start), .pi_ivf(x$end))

.pi_iv_contains <- function(i, o) if (.is_pi_iv(o)) i$start <= o$start && o$end <= i$end else i$start <= o && i$end >= o

.pi_iv_lt <- function(s, o) {
  if (!.is_pi_iv(o)) {
    if (s$end < o) return(.pi_ivm(TRUE, s$is_valid))
    if (s$start > o) return(.pi_ivm(FALSE, s$is_valid))
    return(.pi_ivm(NA, s$is_valid))
  }
  v <- .fuzzy_and2(s$is_valid, o$is_valid)
  if (s$end < o$start) return(.pi_ivm(TRUE, v))
  if (s$start > o$end) return(.pi_ivm(FALSE, v))
  .pi_ivm(NA, v)
}
.pi_iv_le <- function(s, o) {
  if (!.is_pi_iv(o)) {
    if (s$end <= o) return(.pi_ivm(TRUE, s$is_valid))
    if (s$start > o) return(.pi_ivm(FALSE, s$is_valid))
    return(.pi_ivm(NA, s$is_valid))
  }
  v <- .fuzzy_and2(s$is_valid, o$is_valid)
  if (s$end <= o$start) return(.pi_ivm(TRUE, v))
  if (s$start > o$end) return(.pi_ivm(FALSE, v))
  .pi_ivm(NA, v)
}
.pi_iv_gt <- function(s, o) {
  if (.is_pi_iv(o)) return(.pi_iv_lt(o, s))
  if (s$start > o) return(.pi_ivm(TRUE, s$is_valid))
  if (s$end < o) return(.pi_ivm(FALSE, s$is_valid))
  .pi_ivm(NA, s$is_valid)
}
.pi_iv_ge <- function(s, o) {
  if (.is_pi_iv(o)) return(.pi_iv_le(o, s))
  if (s$start >= o) return(.pi_ivm(TRUE, s$is_valid))
  if (s$end < o) return(.pi_ivm(FALSE, s$is_valid))
  .pi_ivm(NA, s$is_valid)
}
.pi_iv_eq <- function(s, o, ne = FALSE) {
  flip <- function(m) if (ne) .pi_ivm(if (is.na(m$a)) NA else !m$a, m$b) else m
  if (!.is_pi_iv(o)) {
    if (s$start == o && s$end == o) return(flip(.pi_ivm(TRUE, s$is_valid)))
    return(flip(.pi_ivm(if (.pi_iv_contains(s, o)) NA else FALSE, s$is_valid)))
  }
  v <- .fuzzy_and2(s$is_valid, o$is_valid)
  if (s$start == o$start && s$end == o$end) return(flip(.pi_ivm(TRUE, v)))
  flip(.pi_ivm(if (!is.na(.pi_iv_lt(s, o)$a)) FALSE else NA, v))
}

.pi_iv_add <- function(s, o) {
  if (!.is_pi_iv(o)) return(.pi_iv(s$start + o, s$end + o, is_valid = s$is_valid))
  .pi_iv(s$start + o$start, s$end + o$end, is_valid = .fuzzy_and2(s$is_valid, o$is_valid))
}
.pi_iv_sub <- function(s, o) {
  if (!.is_pi_iv(o)) return(.pi_iv(s$start - o, s$end - o, is_valid = s$is_valid))
  .pi_iv(s$start - o$end, s$end - o$start, is_valid = .fuzzy_and2(s$is_valid, o$is_valid))
}
.pi_iv_neg <- function(s) .pi_iv(-s$end, -s$start, is_valid = s$is_valid)
.pi_iv_mul <- function(s, o) {
  if (!.is_pi_iv(o)) return(.pi_iv(s$start * o, s$end * o, is_valid = s$is_valid))
  if (isFALSE(s$is_valid) || isFALSE(o$is_valid)) return(.pi_iv_all(FALSE))
  if (is.na(s$is_valid) || is.na(o$is_valid)) return(.pi_iv_all(NA))
  p <- c(s$start * o$start, s$end * o$start, s$start * o$end, s$end * o$end)
  .pi_iv(min(p), max(p))
}
.pi_iv_div <- function(s, o) {
  if (!.is_pi_iv(s)) s <- .pi_iv(s)
  if (!isTRUE(s$is_valid)) return(.pi_iv_all(s$is_valid))
  if (!.is_pi_iv(o)) {
    if (o == 0) return(.pi_iv_all(FALSE))
    return(.pi_iv(s$start / o, s$end / o))
  }
  if (isFALSE(o$is_valid)) return(.pi_iv_all(FALSE))
  if (is.na(o$is_valid)) return(.pi_iv_all(NA))
  if (.pi_iv_contains(o, 0)) return(.pi_iv_all(NA))
  if (o$end < 0) { s <- .pi_iv_neg(s); o <- .pi_iv_neg(o) }
  p <- c(s$start / o$start, s$end / o$start, s$start / o$end, s$end / o$end)
  .pi_iv(max(p), min(p))
}
# the parity of the numerator and denominator of nsimplify(p)
.pi_iv_parity <- function(p) {
  r <- nsimplify(Float(p))
  nd <- fraction(r)
  c(num = as.numeric(N(nd[[1L]])) %% 2 == 0, den = as.numeric(N(nd[[2L]])) %% 2 == 0)
}
.pi_iv_pow <- function(s, o) {
  if (!isTRUE(s$is_valid)) return(s)
  if (.is_pi_iv(o)) return(exp(o * log(s)))
  if (o < 0) return(.pi_iv_div(1, .pi_iv_pow(s, abs(o))))
  if (o == round(o)) {
    if (o %% 2 == 1) return(.pi_iv(s$start^o, s$end^o))
    if (s$start < 0 && s$end > 0) return(.pi_iv(0, max(s$start^o, s$end^o)))
    return(.pi_iv(s$start^o, s$end^o))
  }
  par <- .pi_iv_parity(o)
  if (par[["num"]]) {
    a <- abs(s$start)^o; b <- abs(s$end)^o
    return(if (a < 0) .pi_iv(0, max(a, b)) else .pi_iv(a, b))
  }
  if (par[["den"]]) {
    if (s$end < 0) return(.pi_iv_all(FALSE))
    if (s$start < 0) return(.pi_iv(0, s$end^o, is_valid = NA))
    return(.pi_iv(s$start^o, s$end^o))
  }
  a <- if (s$start < 0) -abs(s$start)^o else s$start^o
  b <- if (s$end < 0) -abs(s$end)^o else s$end^o
  .pi_iv(a, b, is_valid = s$is_valid)
}
.pi_iv_rpow <- function(o, s) {
  if (!isTRUE(s$is_valid)) return(s)
  if (o < 0) {
    if (s$end - s$start > 0) return(.pi_iv_all(FALSE))
    if (.pi_iv_parity(s$start)[["den"]]) return(.pi_iv_all(FALSE))
    v <- -abs(o)^s$start
    return(.pi_iv(v, v))
  }
  .pi_iv(o^s$start, o^s$end)
}

#' @export
Ops.sym_interval <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.pi_iv_neg(e1))
    if (.Generic == "+") return(e1)
    stop("unary ", .Generic, " is not defined for intervals", call. = FALSE)
  }
  num <- function(v) is.numeric(v) && length(v) == 1L
  if (!((.is_pi_iv(e1) || num(e1)) && (.is_pi_iv(e2) || num(e2))))
    stop("NotImplementedError: interval arithmetic with ", class(if (.is_pi_iv(e1)) e2 else e1)[1L], call. = FALSE)
  l <- .is_pi_iv(e1)
  switch(.Generic,
    `+` = if (l) .pi_iv_add(e1, e2) else .pi_iv_add(e2, e1),
    `-` = if (l) .pi_iv_sub(e1, e2) else .pi_iv(e1 - e2$end, e1 - e2$start, is_valid = e2$is_valid),
    `*` = if (l) .pi_iv_mul(e1, e2) else .pi_iv_mul(e2, e1),
    `/` = .pi_iv_div(e1, e2),
    `^` = if (l) .pi_iv_pow(e1, e2) else .pi_iv_rpow(e1, e2),
    `<` = if (l) .pi_iv_lt(e1, e2) else .pi_iv_gt(e2, e1),
    `>` = if (l) .pi_iv_gt(e1, e2) else .pi_iv_lt(e2, e1),
    `<=` = if (l) .pi_iv_le(e1, e2) else .pi_iv_ge(e2, e1),
    `>=` = if (l) .pi_iv_ge(e1, e2) else .pi_iv_le(e2, e1),
    `==` = if (l) .pi_iv_eq(e1, e2) else .pi_iv_eq(e2, e1),
    `!=` = if (l) .pi_iv_eq(e1, e2, ne = TRUE) else .pi_iv_eq(e2, e1, ne = TRUE),
    stop("operator ", .Generic, " is not defined for intervals", call. = FALSE))
}

.pi_iv_mono <- function(f, x) .pi_iv(f(x$start), f(x$end), is_valid = x$is_valid)
.pi_iv_sin <- function(x, cosine = FALSE) {
  f <- if (cosine) base::cos else base::sin
  if (cosine) { if (!(is.finite(x$start) && is.finite(x$end))) return(.pi_iv(-1, 1, is_valid = x$is_valid)) }
  else if (!isTRUE(x$is_valid)) return(.pi_iv(-1, 1, is_valid = x$is_valid))
  na <- floor(x$start / (pi / 2)); nb <- floor(x$end / (pi / 2))
  a <- min(f(x$start), f(x$end)); b <- max(f(x$start), f(x$end))
  if (nb - na > 4) return(.pi_iv(-1, 1, is_valid = x$is_valid))
  if (na == nb) return(.pi_iv(a, b, is_valid = x$is_valid))
  sh <- if (cosine) c(0, 2) else c(1, 3)
  if ((na - sh[1L]) %/% 4 != (nb - sh[1L]) %/% 4) b <- 1
  if ((na - sh[2L]) %/% 4 != (nb - sh[2L]) %/% 4) a <- -1
  if (cosine) .pi_iv(a, b, is_valid = x$is_valid) else .pi_iv(a, b)
}
# a function defined on [lo, hi] (open at the ends when `open`)
.pi_iv_domain <- function(f, x, lo, hi, open = FALSE, check_valid = TRUE) {
  out <- if (open) x$start >= hi || x$end <= lo else x$start > hi || x$end < lo
  if ((check_valid && isFALSE(x$is_valid)) || out) return(.pi_iv_all(FALSE))
  part <- if (open) x$start <= lo || x$end >= hi else x$start < lo || x$end > hi
  if (part) return(.pi_iv_all(NA))
  .pi_iv_mono(f, x)
}
.pi_iv_step <- function(f, x) {
  if (isFALSE(x$is_valid)) return(.pi_iv_all(FALSE))
  a <- f(x$start); b <- f(x$end)
  .pi_iv(a, b, is_valid = if (a == b) x$is_valid else NA)
}
.pi_iv_logf <- function(f, x) {
  if (!isTRUE(x$is_valid)) return(.pi_iv_all(x$is_valid))
  if (x$end <= 0) return(.pi_iv_all(FALSE))
  if (x$start <= 0) return(.pi_iv_all(NA))
  .pi_iv(f(x$start), f(x$end))
}

#' @export
Math.sym_interval <- function(x, ...) {
  switch(.Generic,
    abs = if (x$start < 0 && x$end > 0) .pi_iv(0, max(abs(x$start), abs(x$end)), is_valid = x$is_valid) else .pi_iv(abs(x$start), abs(x$end)),
    exp = .pi_iv_mono(base::exp, x), sinh = .pi_iv_mono(base::sinh, x), tanh = .pi_iv_mono(base::tanh, x),
    atan = .pi_iv_mono(base::atan, x), asinh = .pi_iv_mono(base::asinh, x),
    cosh = if (x$start < 0 && x$end > 0) .pi_iv(1, max(base::cosh(x$start), base::cosh(x$end)), is_valid = x$is_valid) else .pi_iv_mono(base::cosh, x),
    log = .pi_iv_logf(base::log, x), log10 = .pi_iv_logf(base::log10, x),
    sqrt = if (x$end < 0) .pi_iv_all(FALSE) else if (x$start < 0) .pi_iv_all(NA) else .pi_iv_mono(base::sqrt, x),
    sin = .pi_iv_sin(x), cos = .pi_iv_sin(x, cosine = TRUE), tan = .pi_iv_div(.pi_iv_sin(x), .pi_iv_sin(x, cosine = TRUE)),
    asin = .pi_iv_domain(base::asin, x, -1, 1), acos = .pi_iv_domain(base::acos, x, -1, 1),
    atanh = .pi_iv_domain(base::atanh, x, -1, 1, open = TRUE),
    acosh = if (x$end < 1) .pi_iv_all(FALSE) else if (x$start < 1) .pi_iv_all(NA) else .pi_iv_mono(base::acosh, x),
    floor = .pi_iv_step(base::floor, x), ceiling = .pi_iv_step(base::ceiling, x),
    stop("NotImplementedError: ", .Generic, " for intervals", call. = FALSE))
}

.pi_iv_ext <- function(args, f) {
  ok <- vapply(args, function(a) (is.numeric(a) && length(a) == 1L) || .is_pi_iv(a), NA)
  if (!all(ok)) stop("NotImplementedError: interval minimum or maximum of non intervals", call. = FALSE)
  keep <- Filter(function(a) is.numeric(a) || isTRUE(a$is_valid), args)
  if (!length(keep)) return(.pi_iv_all(if (all(vapply(args, function(a) isFALSE(a$is_valid), NA))) FALSE else NA))
  .pi_iv(f(vapply(keep, function(a) if (is.numeric(a)) a else a$start, 0)),
      f(vapply(keep, function(a) if (is.numeric(a)) a else a$end, 0)))
}
#' @rdname intervalmath
#' @export
imin <- function(...) .pi_iv_ext(list(...), min)
#' @rdname intervalmath
#' @export
imax <- function(...) .pi_iv_ext(list(...), max)
#' @rdname intervalmath
#' @export
ceil <- function(x) if (.is_pi_iv(x)) ceiling(x) else .pi_iv(ceiling(x))

# intervalMembership: a pair of three valued logic values
.pi_ivm <- function(a, b) structure(list(a = a, b = b), class = "intervalMembership")
#' @rdname intervalmath
#' @export
intervalMembership <- function(a, b) .pi_ivm(a, b)
#' @export
format.intervalMembership <- function(x, ...) {
  f <- function(v) if (is.na(v)) "None" else if (v) "True" else "False"
  sprintf("intervalMembership(%s, %s)", f(x$a), f(x$b))
}
#' @export
print.intervalMembership <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
#' @export
`[[.intervalMembership` <- function(x, i) {
  if (!i %in% 1:2) stop("IndexError: ", i, " must be a valid indexing for the 2-tuple.", call. = FALSE)
  unclass(x)[[i]]
}
#' @export
Ops.intervalMembership <- function(e1, e2) {
  if (.Generic == "!") return(.pi_ivm(if (is.na(e1$a)) NA else !e1$a, e1$b))
  if (.Generic %in% c("==", "!=")) {
    eq <- inherits(e2, "intervalMembership") && identical(e1$a, e2$a) && identical(e1$b, e2$b)
    return(if (.Generic == "==") eq else !eq)
  }
  if (!inherits(e2, "intervalMembership")) stop("ValueError: The comparison is not supported for ", format(e2), ".", call. = FALSE)
  switch(.Generic,
    `&` = .pi_ivm(.fuzzy_and2(e1$a, e2$a), .fuzzy_and2(e1$b, e2$b)),
    `|` = .pi_ivm(.fuzzy_or2(e1$a, e2$a), .fuzzy_and2(e1$b, e2$b)),
    stop("operator ", .Generic, " is not defined for intervalMembership", call. = FALSE))
}
