# Linear recurrences, after sympy/discrete/recurrences.py.

#' Linear recurrences with constant coefficients
#'
#' `linrec(coeffs, init, n)` evaluates the `n`-th term of the sequence with
#' `a(n) = coeffs[1] a(n - 1) + ... + coeffs[k] a(n - k)` and initial terms
#' `init` (padded with zeros), by SymPy's square and reduce method in
#' `O(k^2 log n)` operations. Coefficients and initial terms may be symbolic.
#'
#' @param coeffs Coefficients of the recurrence.
#' @param init Initial terms.
#' @param n A nonnegative integer.
#' @return A `sym` expression.
#' @examples
#' linrec(c(1, 1), c(0, 1), 10)
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' linrec(list(x, y, z), c(1, 1, 1), 4)
#' @export
linrec <- function(coeffs, init, n) {
  if (!length(coeffs)) return(S_Zero)
  if (is_sym(coeffs))
    stop(.sympy_condition("TypeError", "Expected a sequence of coefficients for the recurrence"))
  if (is_sym(init))
    stop(.sympy_condition("TypeError", "Expected a sequence of values for the initialization of the recurrence"))
  n <- .as_int(n)
  if (n < 0) stop(.sympy_condition("ValueError", "Point of evaluation of recurrence must be a non-negative integer"))
  cf <- lapply(as.list(coeffs), as_sym)
  b <- lapply(as.list(init), as_sym)
  k <- length(cf)
  if (length(b) > k)
    stop(.sympy_condition("TypeError", "Count of initial values should not exceed the order of the recurrence"))
  b <- c(b, rep(list(S_Zero), k - length(b)))
  if (n < k) return(b[[as.integer(n) + 1L]])
  exact <- all(vapply(c(cf, b), function(e) e$type %in% c("Integer", "Rational"), TRUE))
  if (exact) {
    # rational coefficients: the same arithmetic on gmp rationals
    cq <- do.call(c, lapply(cf, function(e) gmp::as.bigq(e$value)))
    bq <- do.call(c, lapply(b, function(e) gmp::as.bigq(e$value)))
    co <- .linrec_coeffs_q(cq, n)
    return(make_rational(sum(co * bq)))
  }
  co <- .linrec_coeffs(cf, n)
  terms <- Map(function(u, v) u * v, co, b)
  Reduce(`+`, terms[-length(terms)], terms[[length(terms)]])
}

.linrec_coeffs <- function(c, n) {
  k <- length(c)
  square_reduce <- function(u, offset) {
    w <- rep(list(S_Zero), 2L * length(u) - 1L + offset)
    for (i in seq_along(u)) for (j in seq_along(u)) {
      idx <- offset + i + j - 1L
      w[[idx]] <- w[[idx]] + u[[i]] * u[[j]]
    }
    if (length(w) - 1L >= k) for (j in (length(w) - 1L):k) for (i in 0:(k - 1L))
      w[[j - i]] <- w[[j - i]] + w[[j + 1L]] * c[[i + 1L]]
    w[seq_len(k)]
  }
  final <- function(n) {
    if (n < k) {
      out <- rep(list(S_Zero), k)
      out[[as.integer(n) + 1L]] <- S_One
      return(out)
    }
    square_reduce(final(n %/% 2), as.integer(n %% 2))
  }
  final(n)
}

.linrec_coeffs_q <- function(c, n) {
  k <- length(c)
  square_reduce <- function(u, offset) {
    w <- gmp::as.bigq(rep(0, 2L * length(u) - 1L + offset))
    for (i in seq_along(u)) {
      idx <- offset + i + seq_along(u) - 1L
      w[idx] <- w[idx] + u[i] * u
    }
    if (length(w) - 1L >= k) for (j in (length(w) - 1L):k) {
      if (w[j + 1L] == 0) next
      idx <- j - (0:(k - 1L))
      w[idx] <- w[idx] + w[j + 1L] * c
    }
    w[seq_len(k)]
  }
  final <- function(n) {
    if (n < k) {
      out <- gmp::as.bigq(rep(0, k))
      out[as.integer(n) + 1L] <- gmp::as.bigq(1)
      return(out)
    }
    square_reduce(final(n %/% 2), as.integer(n %% 2))
  }
  final(n)
}
