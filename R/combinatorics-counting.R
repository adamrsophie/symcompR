# Counting permutations, combinations and derangements of multisets:
# nP(), nC() and nD() of sympy/functions/combinatorial/numbers.py, and the
# exported Gosper helpers of sympy/concrete/gosper.py.

# _multiset_histogram: the multiplicities of a multiset (a vector or list of
# items, or a named list or vector of counts)
.mh <- function(n) {
  if (!is.null(names(n)) && all(nzchar(names(n)))) {
    v <- unlist(n)
    if (!is.numeric(v) || any(v < 0) || any(v != round(v))) stop("ValueError", call. = FALSE)
    m <- as.integer(v[v > 0])
  } else {
    keys <- vapply(as.list(n), function(z) if (is_sym(z)) z$key else paste(class(z)[1L], format(z)), "")
    m <- as.integer(table(factor(keys, levels = unique(keys))))
  }
  list(m = m, items = length(m), N = sum(m))
}
.cnt_int <- function(n) is.numeric(n) && length(n) == 1L || inherits(n, "bigz") || (is_sym(n) && n$type == "Integer")
.cnt_z <- function(n) gmp::as.bigz(if (is_sym(n)) .Z(n) else n)
.cnt_out <- function(v) make_rational(gmp::as.bigq(v))

#' Permutations, combinations and derangements
#'
#' `nP(n, k)` counts the permutations of `k` items taken from `n` (a number
#' of distinct items, or a multiset given as a vector, list, or named list of
#' counts), `nC(n, k)` the combinations, with `replacement` allowing repeats,
#' and `nD(i)` the derangements of the multiset `i` (or `nD(n = 5)` of `n`
#' distinct items, or `nD(m = list(...))` of a multiset given by
#' multiplicities); all as SymPy's functions of the same names. With `k`
#' omitted, `nP()` and `nC()` sum over all sizes.
#'
#' @param n A nonnegative integer or a multiset.
#' @param k The number of items taken, or `NULL` for all sizes.
#' @param replacement Allow items to be taken again.
#' @param i A multiset (vector, list, or named list of counts).
#' @param m For `nD()`: a multiset given as counts of multiplicities.
#' @return An `Integer`.
#' @examples
#' nP(3, 2)
#' nP(c("a", "a", "b"), 2)
#' nC(c("a", "a", "b"), 2)
#' nC(5, 2, replacement = TRUE)
#' nD(c(1, 2, 3))
#' nD(n = 4)
#' @export
nP <- function(n, k = NULL, replacement = FALSE) {
  if (!is.null(k) && k < 0) stop("ValueError: k cannot be negative", call. = FALSE)
  if (.cnt_int(n)) return(.cnt_out(.np_int(.cnt_z(n), k, replacement)))
  .cnt_out(.np_mh(.mh(n), k, replacement))
}
.np_int <- function(n, k, replacement) {
  if (!is.null(k) && k == 0) return(gmp::as.bigz(1))
  if (is.null(k)) return(sum(do.call(c, lapply(0:as.integer(n), function(i) .np_int(n, i, replacement)))))
  k <- gmp::as.bigz(k)
  if (replacement) return(n^k)
  if (k > n) return(gmp::as.bigz(0))
  if (k == 1) return(n)
  prod(gmp::as.bigz(seq(as.numeric(n - k + 1), as.numeric(n))))
}
.np_mh <- function(h, k, replacement) {
  if (!is.null(k) && k == 0) return(gmp::as.bigz(1))
  if (is.null(k)) return(sum(do.call(c, lapply(0:h$N, function(i) .np_mh(h, i, replacement)))))
  if (replacement) return(gmp::as.bigz(h$items)^k)
  if (k == h$N) return(gmp::factorialZ(k) / prod(gmp::factorialZ(h$m[h$m > 1])))
  if (k > h$N) return(gmp::as.bigz(0))
  if (k == 1) return(gmp::as.bigz(h$items))
  memo <- new.env()
  rec <- function(m, k) {
    if (k == 0) return(gmp::as.bigz(1))
    m <- m[m > 0]
    key <- paste(c(k, sort(m)), collapse = ",")
    if (!is.null(memo[[key]])) return(memo[[key]])
    if (k == 1) return(gmp::as.bigz(length(m)))
    tot <- gmp::as.bigz(0)
    for (j in seq_along(m)) { m2 <- m; m2[j] <- m2[j] - 1L; tot <- tot + rec(m2, k - 1) }
    memo[[key]] <- tot
    tot
  }
  rec(h$m, k)
}

#' @rdname nP
#' @export
nC <- function(n, k = NULL, replacement = FALSE) {
  if (.cnt_int(n)) {
    n <- .cnt_z(n)
    if (is.null(k)) {
      if (!replacement) return(.cnt_out(gmp::as.bigz(2)^n))
      return(.cnt_out(sum(do.call(c, lapply(0:as.integer(n), function(i) .cnt_z(nC(n, i, TRUE)))))))
    }
    if (k < 0) stop("ValueError: k cannot be negative", call. = FALSE)
    if (replacement) return(.cnt_out(gmp::chooseZ(n + k - 1, k)))
    return(.cnt_out(gmp::chooseZ(n, k)))
  }
  h <- .mh(n)
  if (is.null(k)) {
    if (!replacement) return(.cnt_out(prod(gmp::as.bigz(h$m + 1L))))
    return(.cnt_out(sum(do.call(c, lapply(0:h$N, function(i) .cnt_z(nC(n, i, TRUE)))))))
  }
  if (replacement) return(nC(h$items, k, TRUE))
  if (k %in% c(1, h$N - 1)) return(.cnt_out(h$items))
  if (k %in% c(0, h$N)) return(S_One)
  if (k > h$N) return(S_Zero)
  .cnt_out(.aop_product(h$m)[[k + 1L]])
}

# _AOP_product: coefficients of prod((1 - x^(m+1))/(1 - x)), the number of
# k-subsets of the multiset for every k
.aop_product <- function(m) {
  ord <- sum(m)
  need <- (ord + 2L) %/% 2L
  m <- as.integer(m)
  rv <- c(rep(gmp::as.bigz(1), m[length(m)] + 1L))
  if (length(rv) < need) rv <- c(rv, rep(gmp::as.bigz(0), need - length(rv)))
  rv <- rv[seq_len(need)]
  m <- m[-length(m)]
  while (length(m)) {
    ni <- m[length(m)]; m <- m[-length(m)]
    N <- ni + 1L
    was <- rv
    if (min(N, length(rv)) > 1L) for (i in 2:min(N, length(rv))) rv[i] <- rv[i] + rv[i - 1L]
    if (N < need) for (i in (N + 1L):need) rv[i] <- rv[i] + rv[i - 1L] - was[i - N]
  }
  rev_ <- rev(rv)
  if (ord %% 2L) c(rv, rev_) else c(rv[-length(rv)], rev_)
}

#' @rdname nP
#' @export
nD <- function(i = NULL, n = NULL, m = NULL) {
  ok <- function(v) {
    if (!is.numeric(v) || any(v != round(v))) stop("TypeError: expecting integer values", call. = FALSE)
    if (any(v < 0)) stop("ValueError: value must not be negative", call. = FALSE)
  }
  if (sum(!vapply(list(i, n, m), is.null, NA)) != 1L) stop("ValueError: enter only 1 of i, n, or m", call. = FALSE)
  if (!is.null(n)) {
    ok(n)
    if (n == 0) return(S_Zero)
    return(subfactorial(as_sym(n)))
  }
  if (!is.null(i)) {
    if (.cnt_int(i) && is.null(names(i))) stop("TypeError: items must be a list or dictionary", call. = FALSE)
    if (!length(i)) return(S_Zero)
    if (!is.null(names(i)) && all(nzchar(names(i)))) ok(unlist(i))
    mult <- .mh(i)$m
  } else {
    if (!is.null(names(m)) && all(nzchar(names(m)))) {
      ks <- as.numeric(names(m)); vs <- unlist(m); ok(ks); ok(vs)
      mult <- rep(ks, vs)
    } else { ok(unlist(m)); mult <- unlist(m) }
    mult <- mult[mult > 0]
  }
  if (!length(mult)) return(S_Zero)
  N <- sum(mult); nkey <- length(mult); big <- max(mult)
  if (big == 1) return(subfactorial(.int(nkey)))
  if (big * 2 > N) return(S_Zero)
  counts <- table(mult)
  if (big * 2 == N) {
    if (nkey == 2L && length(counts) == 1L) return(S_One)
    if (nkey - 1L == big) return(.cnt_out(gmp::factorialZ(big)))
  }
  # integrate(exp(-x)*prod(laguerre(k, x)^count), (x, 0, oo)), exactly:
  # the integral of x^j*exp(-x) is j!
  x <- Dummy("x")
  pe <- expand(.mul(lapply(names(counts), function(k) .pow(fcall("laguerre", list(.int(as.integer(k)), x)), .int(counts[[k]])))))
  cl <- .coeff_list(pe, x)
  deg <- length(cl) - 1L
  tot <- .add(lapply(seq_along(cl), function(j) .mul(list(cl[[j]], .cnt_out(gmp::factorialZ(deg - j + 1L))))))
  fcall("Abs", list(tot))
}

#' Gosper's algorithm helpers
#'
#' `gosper_normal(f, g, n)` computes the Gosper normal form `(A, B, C)` of
#' the rational function `f/g` in `n`, and `gosper_term(f, n)` a
#' hypergeometric term `g` with `g(n + 1) - g(n) = f(n)`, or `NULL`, as
#' SymPy's functions of these names.
#'
#' @param f,g Polynomials in `n` (for `gosper_term()`, `f` a hypergeometric term).
#' @param n The variable.
#' @param polys If `TRUE`, return `Poly` objects.
#' @return A list of three polynomials, or an expression or `NULL`.
#' @examples
#' n <- Symbol("n")
#' gosper_normal(4*n + 5, 2*(4*n + 1)*(2*n + 3), n)
#' gosper_term(factorial(n), n)
#' @export
gosper_normal <- function(f, g, n, polys = TRUE) {
  n <- as_sym(n)
  r <- .gosper_normal(expand(as_sym(f)), expand(as_sym(g)), n)
  if (polys) lapply(r, function(p) Poly(p, n)) else r
}

#' @rdname gosper_normal
#' @export
gosper_term <- function(f, n) .gosper_term_f(as_sym(f), as_sym(n))
