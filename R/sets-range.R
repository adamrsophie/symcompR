# Range and DisjointUnion, after sympy/sets/fancysets.py and sets.py.

#' Ranges of integers and disjoint unions
#'
#' `Range(stop)` or `Range(start, stop, step)` is the set of integers
#' `start, start + step, ...` before `stop`, in SymPy's canonical form
#' (the stop moves to the first value past the last element, and an empty
#' range is `Range(0, 0, 1)`). It supports membership with [contains()],
#' intersections with intervals, ranges and the integers, and
#' `sequence_terms(r, n)` for its first elements. `DisjointUnion(...)`
#' tags the elements of each set with the set's index, like SymPy.
#'
#' @param start,stop,step Integers or infinities.
#' @param ... Sets.
#' @return A set.
#' @examples
#' Range(1, 10, 3)
#' contains(Range(1, 10, 3), 7)
#' Intersection(Range(1, 10, 3), Interval(3, 8))
#' @name Range
NULL

#' @rdname Range
#' @export
Range <- function(start, stop = NULL, step = 1) {
  if (is.null(stop)) { stop <- start; start <- 0 }
  start <- as_sym(start); stop <- as_sym(stop); step <- as_sym(step)
  empty <- new_sym("Range", list(S_Zero, S_Zero, S_One))
  if (is_zero_num(step)) stop("step cannot be 0", call. = FALSE)
  if (is_inf_const(start) && is_inf_const(stop)) {
    # equal ends or a step against the direction give the empty range
    if (identical(start$key, stop$key) || (is_exact_number(step) && is_const(start, "oo") == (num_sign(step) > 0))) return(empty)
    if (!(is_exact_number(step) && abs(step$value) == 1)) stop("step must be 1 or -1 for an unbounded Range", call. = FALSE)
    return(new_sym("Range", list(start, stop, step)))
  }
  if (is_inf_const(step)) {
    # one element when the first step reaches past stop in its direction
    up <- is_const(step, "oo")
    ok <- .is_true(if (up) Lt(start, stop) else Gt(start, stop))
    return(if (ok) new_sym("Range", list(start, .add(list(start, S_One)), S_One)) else empty)
  }
  if (is_exact_number(step) && (is_inf_const(stop) || is_inf_const(start))) {
    pos <- num_sign(step) > 0
    if (is_const(stop, "oo") && !pos || is_const(stop, "-oo") && pos ||
        is_const(start, "oo") && pos || is_const(start, "-oo") && !pos) return(empty)
  }
  span <- .sub(stop, start)
  if (is_exact_number(span) && is_exact_number(step)) {
    nb <- .qceil(gmp::as.bigq(span$value) / gmp::as.bigq(step$value))
    if (nb <= 0) return(empty)
    stop <- .add(list(start, .mul(list(.zs(nb), step))))
  }
  new_sym("Range", list(start, stop, step))
}

.range_parts <- function(r) list(start = r$args[[1L]], stop = r$args[[2L]], step = r$args[[3L]])
.range_size <- function(r) {
  p <- .range_parts(r)
  if (is_inf_const(p$start) || is_inf_const(p$stop)) return(S_Infinity)
  .div(.sub(p$stop, p$start), p$step)
}
.range_empty <- function(r) is_zero_num(.range_size(r))
.range_last <- function(r) { p <- .range_parts(r); .sub(p$stop, p$step) }
.range_inf <- function(r) {
  p <- .range_parts(r)
  if (num_sign(p$step) > 0) p$start else .range_last(r)
}
.range_sup <- function(r) {
  p <- .range_parts(r)
  if (num_sign(p$step) > 0) .range_last(r) else p$start
}
.range_item <- function(r, i) {
  p <- .range_parts(r)
  if (is_const(p$start, "-oo")) .sub(.range_last(r), .mul(list(.int(i), p$step)))
  else .add(list(p$start, .mul(list(.int(i), p$step))))
}

.range_contains <- function(r, x) {
  if (.range_empty(r)) return(S_false)
  if (identical(.facts(x)[["integer"]], FALSE)) return(S_false)
  if (!is_exact_number(x)) return(new_sym("Contains", list(x, r)))
  if (x$type != "Integer") return(S_false)
  p <- .range_parts(r)
  if (is_inf_const(p$start) && is_inf_const(p$stop)) return(S_true)
  ref <- if (is_inf_const(p$start)) .range_last(r) else p$start
  if (!is_exact_number(.sub(x, ref)) || !is_exact_number(p$step)) return(new_sym("Contains", list(x, r)))
  if (gmp::as.bigz(.Z(.sub(x, ref))) %% .zabs(.Z(p$step)) != 0) return(S_false)
  lo <- .range_inf(r); hi <- .range_sup(r)
  inb <- And(if (is_inf_const(lo)) S_true else Ge(x, lo), if (is_inf_const(hi)) S_true else Le(x, hi))
  if (.is_true(inb)) S_true else if (.is_false(inb)) S_false else new_sym("Contains", list(x, r))
}

# Range of the elements of r inside [lo, hi] (open ends allowed).
.range_clip <- function(r, lo, hi, lo_open = FALSE, hi_open = FALSE) {
  p <- .range_parts(r)
  st <- .Z(p$step)
  a <- .range_inf(r); b <- .range_sup(r)
  ast <- abs(st)
  ref <- if (!is_inf_const(a)) .Z(a) else if (!is_inf_const(b)) .Z(b) else gmp::as.bigz(0)
  first <- if (is_inf_const(lo)) a else {
    l <- if (lo_open) .exact_floor(lo) else .exact_ceil(lo)
    if (is.null(l)) return(NULL)
    if (lo_open) l <- l + 1
    l <- gmp::as.bigz(l)
    k <- l + .pmod(ref - l, ast)
    if (!is_inf_const(a) && k < .Z(a)) a else .zs(k)
  }
  last <- if (is_inf_const(hi)) b else {
    h <- if (hi_open) .exact_ceil(hi) else .exact_floor(hi)
    if (is.null(h)) return(NULL)
    if (hi_open) h <- h - 1
    h <- gmp::as.bigz(h)
    k <- h - .pmod(h - ref, ast)
    if (!is_inf_const(b) && k > .Z(b)) b else .zs(k)
  }
  if (!is_inf_const(first) && !is_inf_const(last) && .is_true(Gt(first, last))) return(S_EmptySet)
  if (is_inf_const(last)) return(Range(first, S_Infinity, .zs(ast)))
  if (is_inf_const(first)) return(Range(S_NegativeInfinity, .add(list(last, .zs(ast))), .zs(ast)))
  if (st > 0) Range(first, .add(list(last, .zs(ast))), .zs(ast)) else Range(last, .sub(first, .zs(ast)), .zs(-ast))
}
# Exact floor of a real constant: rationals, r*m^(1/k) radicals, or a
# value whose double bounds certify the floor; NULL otherwise.
.exact_floor <- function(v) {
  v <- as_sym(v)
  if (is_exact_number(v)) return(.bigq_floor(gmp::as.bigq(v$value)))
  if (length(free_symbols(v))) return(NULL)
  sp <- as_coeff_mul_split(v)
  r <- sp[[1L]]; p <- sp[[2L]]
  if (is_exact_number(r) && p$type == "Pow" && is_exact_number(p$args[[1L]]) && p$args[[2L]]$type == "Rational" &&
      gmp::numerator(gmp::as.bigq(p$args[[2L]]$value)) == 1 && p$args[[1L]]$value > 0) {
    kz <- gmp::denominator(gmp::as.bigq(p$args[[2L]]$value))
    q <- gmp::as.bigq(r$value)
    if (q == 0) return(gmp::as.bigz(0))
    m <- gmp::as.bigq(p$args[[1L]]$value)
    if (kz > 10000L) {
      # a huge root order: m^(1/k) lies in [1, 2) once m < 2^k
      if (abs(q) != 1) return(NULL)
      if (m == 1) return(if (q > 0) gmp::as.bigz(1) else gmp::as.bigz(-1))
      if (m < 1) return(if (q > 0) gmp::as.bigz(0) else gmp::as.bigz(-1))
      if (gmp::sizeinbase(.bigq_floor(m) + 1, 2) > kz) return(NULL)
      return(if (q > 0) gmp::as.bigz(1) else gmp::as.bigz(-2))
    }
    k <- as.integer(as.character(kz))
    X <- abs(q)^k * gmp::as.bigq(p$args[[1L]]$value)
    fl <- gmp::as.bigz(.iroot(gmp::as.bigz(.bigq_floor(X)), k)$root)
    if (q > 0) return(fl)
    exact <- gmp::as.bigq(fl)^k == X
    return(if (exact) -fl else -fl - 1)
  }
  tryCatch(.nt_certified_floor(v), error = function(e) NULL)
}
.exact_ceil <- function(v) { f <- .exact_floor(.neg(as_sym(v))); if (is.null(f)) NULL else -f }
.qceil <- function(q) -.bigq_floor(-gmp::as.bigq(q))

.range_intersect <- function(A, B) {
  if (A$type == "Range" && .range_empty(A) || B$type == "Range" && .range_empty(B)) return(S_EmptySet)
  if (B$type == "Range") {
    # the integers leave the other range as it is
    allz <- function(R) is_inf_const(R$args[[1L]]) && is_inf_const(R$args[[2L]])
    if (allz(A)) return(B)
    if (allz(B)) return(A)
    pa <- .range_parts(A); pb <- .range_parts(B)
    if (!all(vapply(list(pa$step, pb$step), is_exact_number, TRUE))) return(NULL)
    s1 <- abs(.Z(pa$step)); s2 <- abs(.Z(pb$step))
    rref <- function(R) if (!is_inf_const(.range_inf(R))) .Z(.range_inf(R)) else if (!is_inf_const(.range_sup(R))) .Z(.range_sup(R)) else gmp::as.bigz(0)
    r1 <- rref(A); r2 <- rref(B)
    g <- gmp::gcd.bigz(s1, s2)
    if (.pmod(r2 - r1, g) != 0) return(S_EmptySet)
    L <- s1 * s2 %/% g
    m2 <- s2 %/% g
    k <- if (m2 == 1) gmp::as.bigz(0) else .pmod((r2 - r1) %/% g * gmp::inv.bigz(.pmod(s1 %/% g, m2), m2), m2)
    x0 <- r1 + s1 * k
    lo <- .max_sym(.range_inf(A), .range_inf(B)); hi <- .min_sym(.range_sup(A), .range_sup(B))
    if (is_inf_const(lo)) {
      # the last common element at or below hi closes the range
      last <- .Z(hi) - .pmod(.Z(hi) - x0, L)
      return(Range(S_NegativeInfinity, .zs(last + L), .zs(L)))
    }
    first <- x0 + L * .qceil(gmp::as.bigq(.Z(lo) - x0, L))
    return(.range_clip(Range(.zs(first), S_Infinity, .zs(L)), lo, hi))
  }
  if (B$type == "Interval") {
    i <- .iv(B)
    if (length(free_symbols(i$a)) || length(free_symbols(i$b)) || length(free_symbols(A))) return(NULL)
    return(.range_clip(A, i$a, i$b, i$lo && !is_inf_const(i$a), i$ro && !is_inf_const(i$b)))
  }
  if (B$type == "SpecialSet") {
    if (B$name %in% c("Integers", "Reals", "Complexes", "Rationals")) return(A)
    if (B$name == "Naturals") return(.range_clip(A, S_One, S_Infinity))
    if (B$name == "Naturals0") return(.range_clip(A, S_Zero, S_Infinity))
  }
  NULL
}

register_str_printer("Range", function(e) paste0("Range(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
register_latex_printer("Range", function(e) {
  if (.range_empty(e)) return("\\emptyset")
  n <- .range_size(e)
  p <- .range_parts(e)
  els <- if (is_const(p$start, "-oo")) c("\\ldots", vapply(1:0, function(i) .latex(.range_item(e, i)), ""))
         else if (!is_exact_number(n)) c(vapply(0:1, function(i) .latex(.range_item(e, i)), ""), "\\ldots")
         else if (.int_val(n) > 4) c(vapply(0:1, function(i) .latex(.range_item(e, i)), ""), "\\ldots", .latex(.range_last(e)))
         else vapply(seq_len(.int_val(n)) - 1L, function(i) .latex(.range_item(e, i)), "")
  paste0("\\left\\{", paste(els, collapse = ", "), "\\right\\}")
})
register_rebuilder("Range", function(e, args) Range(args[[1L]], args[[2L]], args[[3L]]))

#' @rdname Range
#' @export
DisjointUnion <- function(...) new_sym("DisjointUnion", lapply(list(...), as_sym))
register_str_printer("DisjointUnion", function(e) paste0("DisjointUnion(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
.disjoint_contains <- function(s, x) {
  if (!(x$type == "Tuple" && length(x$args) == 2L)) return(S_false)
  i <- x$args[[2L]]
  if (i$type != "Integer" || .int_val(i) < 0 || .int_val(i) >= length(s$args)) return(S_false)
  .contains(s$args[[.int_val(i) + 1L]], x$args[[1L]])
}
