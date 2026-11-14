# Schur numbers and partitions (sympy/combinatorics/schur_number.py), the
# number of groups of a given order and the nilpotent, abelian and cyclic
# numbers (sympy/combinatorics/group_numbers.py), and common_assumptions()
# of sympy/core/assumptions.py.

def_function("SchurNumber", nargs = 1L,
  eval = function(args) {
    k <- args[[1L]]
    if (!is_number_type(k)) return(NULL)
    if (is_const(k, "oo")) return(S_Infinity)
    if (is_zero_num(k)) return(S_Zero)
    if (k$type != "Integer" || k$value < 0) stop("ValueError: k should be a positive integer", call. = FALSE)
    if (k$value <= 5) return(.int(c(1L, 4L, 13L, 44L, 160L)[as.integer(k$value)]))
    NULL
  })

#' Schur numbers and sum-free partitions
#'
#' `SchurNumber(k)` is the largest `n` for which `1, ..., n` splits into `k`
#' sum-free sets (known for `k <= 5`, otherwise unevaluated), and
#' `schur_lower_bound(k)` its lower bound. `schur_partition(n)` splits
#' `1, ..., n` into the fewest sum-free sets that SymPy's construction finds.
#'
#' @param k A nonnegative integer or symbol.
#' @param n A positive integer.
#' @return An expression, or a list of integer vectors.
#' @examples
#' SchurNumber(4)
#' schur_lower_bound(SchurNumber(8))
#' schur_partition(9)
#' @export
SchurNumber <- function(k) fcall("SchurNumber", list(as_sym(k)))

#' @rdname SchurNumber
#' @param s A `SchurNumber()` expression.
#' @export
schur_lower_bound <- function(s) {
  f <- if (is_sym(s) && s$type == "Function" && identical(s$name, "SchurNumber")) s$args[[1L]] else as_sym(s)
  if (f$type == "Integer") {
    v <- as.integer(f$value)
    if (v == 6L) return(.int(536L))
    if (v == 7L) return(.int(1680L))
    return(.sub(.mul(list(.int(3L), schur_lower_bound(.int(v - 1L)))), S_One))
  }
  .div(.sub(.pow(.int(3L), f), S_One), S_Two)
}

#' @rdname SchurNumber
#' @export
schur_partition <- function(n) {
  if (is_sym(n)) {
    if (!is_number_type(n)) stop("ValueError: Input value must be a number", call. = FALSE)
    if (is_const(n, "oo")) stop("ValueError: Input must be finite", call. = FALSE)
    n <- as.numeric(n$value)
  }
  if (n <= 0) stop("ValueError: n must be a non-zero positive integer.", call. = FALSE)
  k <- if (n <= 3) 1 else ceiling(log(2 * n + 1, 3))
  sets <- if (n == 1) list(1) else if (n == 2) list(c(1, 2)) else if (n == 3) list(c(1, 2, 3)) else list(c(1, 4), c(2, 3))
  while (length(sets) < k) {
    cur <- sets
    sets <- lapply(cur, function(it) c((it * 3)[it * 3 <= n], (it * 3 - 1)[it * 3 - 1 <= n]))
    last <- 3 * (0:length(cur)) + 1
    sets[[length(sets) + 1L]] <- last[last <= n]
    lo <- length(sets); hi <- (n - 1) %/% 3
    if (lo <= hi) sets[[length(sets)]] <- c(sets[[length(sets)]], 3 * (lo:hi) + 1)
  }
  lapply(sets, as.integer)
}

.group_factors <- function(n) {
  n <- as.numeric(n)
  if (n != round(n)) stop("ValueError: n must be an integer", call. = FALSE)
  if (n <= 0) stop("ValueError: n must be a positive integer, not ", n, call. = FALSE)
  f <- factorint(n)
  list(p = as.numeric(names(f)), e = as.numeric(f))
}
.is_nilpotent_factors <- function(fs) {
  for (p in fs$p) for (j in seq_along(fs$p)) {
    m <- 1
    for (r in seq_len(fs$e[j])) { m <- m * fs$p[j] %% p; if (m == 1) return(FALSE) }
  }
  TRUE
}

#' Numbers of groups of a given order
#'
#' `groups_count(n)` is the number of groups of order `n` up to isomorphism
#' where SymPy knows it; `is_nilpotent_number(n)`, `is_abelian_number(n)`
#' and `is_cyclic_number(n)` tell whether every group of order `n` is
#' nilpotent, abelian or cyclic.
#'
#' @param n A positive integer.
#' @return An integer or a logical.
#' @examples
#' groups_count(16)
#' groups_count(30)
#' is_cyclic_number(15)
#' @export
groups_count <- function(n) {
  fs <- .group_factors(n)
  n <- as.numeric(n)
  g <- function(a, b) as.numeric(gmp::gcd.bigz(a, b))
  if (length(fs$p) == 1L) {
    p <- fs$p; e <- fs$e
    if (p == 2) { A <- c(1, 1, 2, 5, 14, 51, 267, 2328, 56092, 10494213, 49487367289); if (e < length(A)) return(.int_or_big(A[e + 1])) }
    if (p == 3) { A <- c(1, 1, 2, 5, 15, 67, 504, 9310, 1396077, 5937876645); if (e < length(A)) return(.int_or_big(A[e + 1])) }
    if (e <= 2) return(e)
    if (e == 3) return(5)
    if (e == 4) return(15)
    if (e == 5) return(61 + 2 * p + 2 * g(p - 1, 3) + g(p - 1, 4))
    if (e == 6) return(3 * p^2 + 39 * p + 344 + 24 * g(p - 1, 3) + 11 * g(p - 1, 4) + 2 * g(p - 1, 5))
    if (e == 7) {
      if (p == 5) return(34297)
      return(3 * p^5 + 12 * p^4 + 44 * p^3 + 170 * p^2 + 707 * p + 2455 + (4 * p^2 + 44 * p + 291) * g(p - 1, 3) +
             (p^2 + 19 * p + 135) * g(p - 1, 4) + (3 * p + 31) * g(p - 1, 5) + 4 * g(p - 1, 7) + 5 * g(p - 1, 8) + g(p - 1, 9))
    }
  }
  if (any(fs$e > 1)) {
    small <- c(`12` = 5, `18` = 5, `20` = 5, `24` = 15, `28` = 4, `36` = 14, `40` = 14, `44` = 4, `45` = 2, `48` = 52,
               `50` = 5, `52` = 5, `54` = 15, `56` = 13, `60` = 13, `63` = 4, `68` = 5, `72` = 50, `75` = 3, `76` = 4,
               `80` = 52, `84` = 15, `88` = 12, `90` = 10, `92` = 4)
    k <- as.character(n)
    if (k %in% names(small)) return(unname(small[k]))
    stop("ValueError: Number of groups of order n is unknown or not implemented", call. = FALSE)
  }
  if (length(fs$p) == 2L) { p <- min(fs$p); q <- max(fs$p); return(if (q %% p == 1) 2 else 1) }
  .holder_formula(fs$p)
}
.int_or_big <- function(v) v
# _holder_formula: the number of groups of a squarefree order
.holder_formula <- function(P) {
  F <- P[vapply(P, function(p) all(P %% p != 1), NA)]
  M <- setdiff(P, F)
  s <- 0
  subsets <- unlist(lapply(0:length(M), function(r) if (r == 0) list(numeric()) else
    lapply(utils::combn(seq_along(M), r, simplify = FALSE), function(ix) M[ix])), recursive = FALSE)
  for (ps in subsets) {
    pr <- 1
    for (p in setdiff(M, ps)) {
      cnt <- sum(c(F, ps) %% p == 1)
      pr <- pr * (p^cnt - 1) %/% (p - 1)
      if (!pr) break
    }
    s <- s + pr
  }
  s
}

#' @rdname groups_count
#' @export
is_nilpotent_number <- function(n) .is_nilpotent_factors(.group_factors(n))
#' @rdname groups_count
#' @export
is_abelian_number <- function(n) { fs <- .group_factors(n); all(fs$e < 3) && .is_nilpotent_factors(fs) }
#' @rdname groups_count
#' @export
is_cyclic_number <- function(n) { fs <- .group_factors(n); all(fs$e == 1) && .is_nilpotent_factors(fs) }

#' Assumptions shared by expressions
#'
#' `common_assumptions(exprs, check)` returns the assumptions (among those
#' named in `check`, by default all of SymPy's) whose values are known and
#' equal for every expression, as SymPy's `common_assumptions()`.
#'
#' @param exprs A list of expressions.
#' @param check Optional names of assumptions to consider.
#' @return A named logical vector, sorted by name.
#' @examples
#' common_assumptions(list(-4, 0, sqrt(Integer(2)), 2, sym_pi, oo))
#' common_assumptions(list(0, 1, 2), c("positive", "integer"))
#' @export
common_assumptions <- function(exprs, check = NULL) {
  defined <- c("algebraic", "commutative", "complex", "composite", "even", "extended_negative",
               "extended_nonnegative", "extended_nonpositive", "extended_nonzero", "extended_positive",
               "extended_real", "finite", "imaginary", "infinite", "integer", "irrational", "negative",
               "noninteger", "nonnegative", "nonpositive", "nonzero", "odd", "polar", "positive", "prime",
               "rational", "real", "transcendental", "zero")
  check <- if (is.null(check)) defined else check
  if (!length(check) || !length(exprs)) return(logical())
  vals <- lapply(exprs, function(e) {
    e <- as_sym(e)
    v <- vapply(check, function(k) as.logical(ask_fact(e, k)), NA)
    v[!is.na(v)]
  })
  common <- Reduce(intersect, lapply(vals, names))
  common <- common[vapply(common, function(k) length(unique(vapply(vals, function(v) v[[k]], NA))) == 1L, NA)]
  out <- vals[[1L]][common]
  out[order(names(out))]
}
