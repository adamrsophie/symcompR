# Integer functions of sympy/core/intfunc.py and sympy/core/exprtools.py,
# number theoretic transforms of sympy/discrete/transforms.py, bool_map of
# sympy/logic/boolalg.py and UnevaluatedExpr of sympy/core/expr.py.

.int_args <- function(args) lapply(args, function(a) .nt_int(a))

#' Integer helpers
#'
#' `igcd()` and `ilcm()` are the gcd and lcm of integers. `mod_inverse(a, m)`
#' is the inverse of `a` modulo `m` (for a rational `a`, its inverse in the
#' field). `integer_nthroot(y, n)` returns `list(root, exact)` with `root` the
#' floor of the `n`-th root; `integer_log(y, x)` returns `list(e, exact)` with
#' `x^e <= y < x^(e + 1)`. `trailing(n)` counts trailing zero bits and
#' `num_digits(n, base)` counts digits. `arity(f)` is the number of arguments
#' of an R function.
#'
#' @param ... Integers.
#' @param a,m,y,x,n Integers.
#' @param base Base.
#' @param f An R function.
#' @return An integer, `sym` number or list.
#' @examples
#' igcd(12, 18, -8)
#' mod_inverse(3, 11)
#' integer_nthroot(26, 3)
#' integer_log(125, 5)
#' @name integer-helpers
NULL

#' @rdname integer-helpers
#' @export
igcd <- function(...) {
  a <- .int_args(list(...))
  if (!length(a)) stop("igcd() takes at least one argument", call. = FALSE)
  g <- abs(a[[1L]])
  for (v in a[-1L]) g <- gmp::gcd.bigz(g, v)
  .nt_out(abs(g))
}

#' @rdname integer-helpers
#' @export
ilcm <- function(...) {
  a <- .int_args(list(...))
  if (!length(a)) stop("ilcm() takes at least one argument", call. = FALSE)
  l <- abs(a[[1L]])
  for (v in a[-1L]) {
    if (l == 0 || v == 0) return(0)
    l <- abs(l * v) %/% gmp::gcd.bigz(l, v)
  }
  .nt_out(l)
}

#' @rdname integer-helpers
#' @export
mod_inverse <- function(a, m) {
  q <- .nt_rat(a)
  if (!is.null(q)) return(make_rational(1 / q))
  a <- .nt_int(a)
  m <- .nt_int(m)
  am <- abs(m)
  r <- if (am > 1) .zinv(a, am) else NA
  if (is.na(r)) stop("inverse of ", a, " (mod ", m, ") does not exist", call. = FALSE)
  # a negative modulus gives the residue in (m, 0], as SymPy
  if (m < 0 && r != 0) r <- r - am
  .nt_out(r)
}

#' @rdname integer-helpers
#' @export
integer_nthroot <- function(y, n) {
  y <- .nt_int(y)
  n <- as.numeric(.nt_int(n))
  if (y < 0) stop("y must be nonnegative", call. = FALSE)
  if (n < 1) stop("n must be positive", call. = FALSE)
  r <- .iroot(y, n)
  list(root = .nt_out(r$root), exact = r$exact)
}

#' @rdname integer-helpers
#' @export
integer_log <- function(y, x) {
  y <- .nt_int(y)
  x <- .nt_int(x)
  # signed cases follow sympy/core/intfunc.py
  if (x < 0) {
    r <- integer_log(abs(y), -x)
    return(list(e = r$e, exact = r$exact && (r$e %% 2 == 1) == (y < 0)))
  }
  if (x <= 1) stop(.sympy_condition("ValueError", "base must be 2 or more"))
  if (y < 0) {
    if (x != 2) return(list(e = 0, exact = FALSE))
    return(list(e = integer_log(-y, x)$e, exact = FALSE))
  }
  if (y == 0) stop(.sympy_condition("ValueError", "n cannot be 0"))
  e <- 0
  p <- gmp::as.bigz(1)
  while (p * x <= y) {
    p <- p * x
    e <- e + 1
  }
  list(e = e, exact = p == y)
}

#' @rdname integer-helpers
#' @export
trailing <- function(n) {
  n <- abs(.nt_int(n))
  if (n == 0) return(0)
  .v2(n)
}

#' @rdname integer-helpers
#' @export
num_digits <- function(n, base = 10) {
  n <- abs(.nt_int(n))
  base <- .nt_int(base)
  if (base < 2) stop("base must be at least 2", call. = FALSE)
  d <- 1
  while (n >= base) {
    n <- n %/% base
    d <- d + 1
  }
  d
}

#' @rdname integer-helpers
#' @export
arity <- function(f) {
  fm <- formals(args(f))
  if ("..." %in% names(fm)) return(NA_integer_)
  length(fm)
}

# ---------------------------------------------------------------------------
# Number theoretic transform

.ntt <- function(seq, prime, inverse) {
  p <- .nt_int(prime)
  if (!.is_prime(p)) stop("Expected prime modulus for number-theoretic transform", call. = FALSE)
  a <- .zc(lapply(as.list(seq), function(v) .nt_int(v) %% p))
  n <- length(a)
  if (n < 1L) return(numeric())
  m <- 2^ceiling(log2(n))
  if ((p - 1) %% m != 0) stop("Expected prime modulus of the form (m*2**k + 1)", call. = FALSE)
  if (m > n) a <- c(a, .zc(rep(list(gmp::as.bigz(0)), m - n)))
  rt <- gmp::powm(gmp::as.bigz(primitive_root(p)), (p - 1) %/% m, p)
  if (inverse) rt <- gmp::powm(rt, p - 2, p)
  w <- rep(gmp::as.bigz(1), m)
  for (i in seq_len(m - 1L)) w[i + 1L] <- (w[i] * rt) %% p
  # ponytail: direct O(n^2) transform; switch to the butterfly when n is large.
  out <- lapply(0:(m - 1), function(k) sum(a * w[((seq_len(m) - 1) * k) %% m + 1]) %% p)
  out <- .zc(out)
  if (inverse) out <- (out * gmp::powm(gmp::as.bigz(m), p - 2, p)) %% p
  .nt_out(out)
}

#' Number theoretic transform
#'
#' `ntt(seq, prime)` and `intt(seq, prime)` are the discrete Fourier transform
#' and its inverse over the integers modulo a prime of the form `m 2^k + 1`.
#' The sequence is padded with zeros to a power of two.
#'
#' @param seq A vector of integers.
#' @param prime A prime modulus.
#' @return A vector of integers.
#' @examples
#' ntt(c(1, 2, 3, 5), 3 * 2^8 + 1)
#' intt(ntt(c(1, 2, 3, 5), 769), 769)
#' @export
ntt <- function(seq, prime) .ntt(seq, prime, FALSE)

#' @rdname ntt
#' @export
intt <- function(seq, prime) .ntt(seq, prime, TRUE)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# bool_map

#' Map variables between equivalent boolean expressions
#'
#' `bool_map(b1, b2)` looks for a one to one renaming of the variables of
#' `b1` that makes it equivalent to `b2`. It returns `list(expr, mapping)`
#' with the simplified `b1` and a named list from `b1` variables to `b2`
#' variables, or `FALSE`.
#'
#' @param b1,b2 Boolean expressions.
#' @return A list or `FALSE`.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b"); x <- Symbol("x"); y <- Symbol("y")
#' bool_map(Or(And(a, b), Not(a)), Or(And(x, y), Not(x)))
#' @export
bool_map <- function(b1, b2) {
  b1 <- simplify_logic(as_sym(b1)); b2 <- simplify_logic(as_sym(b2))
  v1 <- .bool_atoms(b1); v2 <- .bool_atoms(b2)
  if (length(v1) != length(v2)) return(FALSE)
  n <- length(v1)
  if (n > 8L) stop("bool_map() is limited to 8 variables", call. = FALSE)
  tt2 <- truth_table(b2, v2)$value
  for (p in if (n) permutations(seq_len(n)) else list(integer())) {
    m <- stats::setNames(v2[p], vapply(v1, .k, ""))
    e <- xreplace(b1, lapply(seq_len(n), function(i) list(v1[[i]], v2[[p[i]]])))
    if (identical(truth_table(e, v2)$value, tt2)) {
      names(m) <- vapply(v1, .str, "")
      return(list(expr = b1, mapping = m))
    }
  }
  FALSE
}

# ---------------------------------------------------------------------------
# UnevaluatedExpr

#' Expression that does not evaluate
#'
#' `UnevaluatedExpr(e)` holds `e` as an atom in arithmetic, so
#' `x * UnevaluatedExpr(1 / x)` stays as it is. [doit()] releases it.
#'
#' @param e An expression.
#' @return A `sym` UnevaluatedExpr.
#' @examples
#' x <- Symbol("x")
#' x * UnevaluatedExpr(1 / x)
#' doit(x * UnevaluatedExpr(1 / x))
#' @export
UnevaluatedExpr <- function(e) new_sym("UnevaluatedExpr", list(as_sym(e)))

register_str_printer("UnevaluatedExpr", function(e) .str(e$args[[1L]]))
register_latex_printer("UnevaluatedExpr", function(e) .latex(e$args[[1L]]))
register_doit("UnevaluatedExpr", function(e) e$args[[1L]])
register_neval("UnevaluatedExpr", function(e, env) .neval(e$args[[1L]], env))
