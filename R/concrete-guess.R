# Guessing formulas and recurrences for sequences, after sympy/concrete/guess.py.

# a sequence of rationals as a bigq vector
.guess_q <- function(l) {
  l <- lapply(as.list(l), as_sym)
  if (!all(vapply(l, .is_rat, TRUE)))
    stop(.sympy_condition("TypeError", "Expected a sequence of rational numbers"))
  do.call(c, lapply(l, function(e) gmp::as.bigq(e$value)))
}

.guess_index_error <- function() stop(.sympy_condition("IndexError", "list index out of range"))

#' Guessing recurrences and generating functions
#'
#' `find_simple_recurrence_vector(l)` returns the coefficients of a linear
#' recurrence with constant coefficients satisfied by the rational sequence
#' `l`, or `0` when none is found; `find_simple_recurrence(v, A, N)` returns
#' the recurrence as an expression in `A(N), A(N + 1), ...`.
#' `guess_generating_function_rational(v, X)` guesses a rational ordinary
#' generating function, and `guess_generating_function(v, X, types,
#' maxsqrtn)` tries the generating function types `"ogf"`, `"egf"`, `"lgf"`,
#' `"hlgf"`, `"lgdogf"` and `"lgdegf"` and returns a named list.
#' `guess(l, all, evaluate, niter, variables)` guesses a closed form built
#' from rational functions and nested products, as SymPy's `guess()`.
#' `rationalize(x, maxcoeff)` recovers a rational number from a float by its
#' continued fraction, stopping at the first partial quotient of at least
#' `maxcoeff`; an expression that is not a number is passed to `radsimp()`.
#'
#' @param l,v A vector or list of rational numbers.
#' @param A A function made by [Function()].
#' @param N,X Symbols.
#' @param types Generating function types, or `"all"`.
#' @param maxsqrtn The largest root of a rational function tried.
#' @param all Return every formula found instead of the first ones.
#' @param evaluate Evaluate the products.
#' @param niter Number of iterations.
#' @param variables Symbols to use instead of `i0, i1, ...`.
#' @param x A number.
#' @param maxcoeff Bound on the partial quotients.
#' @return A list, an expression, a named list or `NULL`.
#' @examples
#' find_simple_recurrence_vector(c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34))
#' find_simple_recurrence(c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34))
#' guess_generating_function_rational(c(5, 8, 13, 21, 34, 55, 89, 144))
#' guess(c(1, 2, 6, 24, 120))
#' rationalize(pi, maxcoeff = 250)
#' @export
find_simple_recurrence_vector <- function(l) {
  l <- .guess_q(l)
  q1 <- gmp::as.bigq(0)
  q2 <- gmp::as.bigq(1)
  b <- 0L
  z <- length(l) %/% 2L
  while (length(q2) <= z) {
    while ({ if (b >= length(l)) .guess_index_error(); l[b + 1L] == 0 }) {
      b <- b + 1L
      if (b == length(l)) {
        cc <- gmp::as.bigz(1)
        for (i in seq_along(q2)) cc <- gmp::lcm.bigz(cc, gmp::denominator(q2[i]))
        if (q2[1L] * cc < 0) cc <- -cc
        return(lapply(seq_along(q2), function(i) make_rational(q2[i] * cc)))
      }
    }
    a <- 1 / l[b + 1L]
    m <- a
    if (b + 1L <= length(l) - 1L) for (k in (b + 1L):(length(l) - 1L)) {
      s <- gmp::as.bigq(0)
      for (j in b:(k - 1L)) s <- s + l[j + 2L] * m[k - j]
      m <- c(m, -s * a)
    }
    l <- m
    m <- gmp::as.bigq(rep(0L, max(length(q2), b + length(q1))))
    m[seq_along(q2)] <- a * q2
    idx <- b + seq_along(q1)
    m[idx] <- m[idx] + q1
    while ({ if (!length(m)) .guess_index_error(); m[length(m)] == 0 }) m <- m[-length(m)]
    q1 <- q2
    q2 <- m
    b <- 1L
  }
  list(S_Zero)
}

#' @rdname find_simple_recurrence_vector
#' @export
find_simple_recurrence <- function(v, A = Function("a"), N = Symbol("n")) {
  p <- find_simple_recurrence_vector(v)
  n <- length(p)
  if (n <= 1L) return(S_Zero)
  .add(lapply(seq_len(n), function(k) .mul(list(A(N + (n - k)), p[[k]]))))
}

#' @rdname find_simple_recurrence_vector
#' @export
rationalize <- function(x, maxcoeff = 10000) {
  # a SymPy Float never equals the Integer floor(x), so it runs until
  # x - floor(x) vanishes and 1/0 gives zoo
  isfloat <- FALSE
  if (is_sym(x)) {
    if (x$type == "Float") { x <- num_value_dbl(x); isfloat <- TRUE }
    else if (!.is_rat(x)) return(radsimp(x))
  }
  exact <- is_sym(x)
  if (exact) x <- gmp::as.bigq(x$value)
  fl <- function(y) if (exact) gmp::as.bigz(floor(y)) else floor(y)
  p0 <- gmp::as.bigz(0); p1 <- gmp::as.bigz(1)
  q0 <- gmp::as.bigz(1); q1 <- gmp::as.bigz(0)
  a <- fl(x)
  while (a < maxcoeff || q1 == 0) {
    ab <- gmp::as.bigz(a)
    p <- ab * p1 + p0
    q <- ab * q1 + q0
    p0 <- p1; p1 <- p
    q0 <- q1; q1 <- q
    if (!isfloat && x == a) break
    if (isfloat && x == a) stop(.sympy_condition("TypeError", "Invalid comparison of non-real zoo"))
    x <- 1 / (x - a)
    a <- fl(x)
  }
  make_rational(gmp::as.bigq(p, q))
}

#' @rdname find_simple_recurrence_vector
#' @export
guess_generating_function_rational <- function(v, X = Symbol("x")) {
  q <- find_simple_recurrence_vector(v)
  n <- length(q)
  if (n <= 1L) return(NULL)
  v <- lapply(as.list(v), as_sym)
  p <- lapply(seq_len(length(v) %/% 2L) - 1L, function(i)
    .add(c(list(S_Zero), lapply(seq_len(min(i + 1L, n)) - 1L, function(k) .mul(list(v[[i - k + 1L]], q[[k + 1L]]))))))
  num <- .add(c(list(S_Zero), lapply(seq_along(p), function(k) .mul(list(p[[k]], .pow(X, .int(k - 1L)))))))
  den <- .add(c(list(S_Zero), lapply(seq_len(n), function(k) .mul(list(q[[k]], .pow(X, .int(k - 1L)))))))
  .div(num, den)
}

# the convolution powers of w, as in each branch of guess_generating_function()
.ggf_try <- function(w, X, maxsqrtn) {
  t <- c(list(S_One), rep(list(S_Zero), length(w) - 1L))
  for (d in seq_len(max(1L, maxsqrtn)) - 1L) {
    t <- lapply(seq_along(w) - 1L, function(n)
      .add(c(list(S_Zero), lapply(0:n, function(i) .mul(list(t[[n - i + 1L]], w[[i + 1L]]))))))
    g <- guess_generating_function_rational(t, X = X)
    if (!is.null(g) && !is_zero_num(g)) return(.pow(g, make_rational(gmp::as.bigq(1L, d + 1L))))
  }
  NULL
}

# f'(x)/f(x) of the series with coefficients z
.ggf_logder <- function(z) {
  a <- z[[1L]]
  w <- list()
  for (n in seq_len(length(z) - 1L) - 1L) {
    s <- .add(c(list(S_Zero), lapply(seq_len(n) - 1L, function(i) .mul(list(w[[length(w) - i]], z[[i + 2L]])))))
    w[[n + 1L]] <- .div(.sub(.mul(list(z[[n + 2L]], .int(n + 1L))), s), a)
  }
  w
}

#' @rdname find_simple_recurrence_vector
#' @export
guess_generating_function <- function(v, X = Symbol("x"), types = "all", maxsqrtn = 2) {
  if ("all" %in% types) types <- c("ogf", "egf", "lgf", "hlgf", "lgdogf", "lgdegf")
  maxsqrtn <- as.integer(maxsqrtn)
  v <- lapply(as.list(v), as_sym)
  result <- list()
  facdiv <- function() {
    w <- list(); f <- S_One
    for (i in seq_along(v)) {
      if (i > 1L) f <- .mul(list(f, .int(i - 1L)))
      w[[i]] <- .div(v[[i]], f)
    }
    w
  }
  if ("ogf" %in% types) {
    g <- .ggf_try(v, X, maxsqrtn)
    if (!is.null(g)) result$ogf <- g
  }
  if ("egf" %in% types) {
    g <- .ggf_try(facdiv(), X, maxsqrtn)
    if (!is.null(g)) result$egf <- g
  }
  if ("lgf" %in% types) {
    w <- list(); f <- S_NegativeOne
    for (i in seq_along(v)) {
      f <- .neg(f)
      w[[i]] <- .div(.mul(list(f, v[[i]])), .int(i))
    }
    g <- .ggf_try(w, X, maxsqrtn)
    if (!is.null(g)) result$lgf <- g
  }
  if ("hlgf" %in% types) {
    g <- .ggf_try(lapply(seq_along(v), function(i) .div(v[[i]], .int(i))), X, maxsqrtn)
    if (!is.null(g)) result$hlgf <- g
  }
  nz <- !is_zero_num(v[[1L]])
  if (nz && ("lgdogf" %in% types || ("ogf" %in% types && is.null(result$ogf)))) {
    g <- .ggf_try(.ggf_logder(v), X, maxsqrtn)
    if (!is.null(g)) {
      result$lgdogf <- g
      if (is.null(result$ogf)) result$ogf <- exp(integrate(g, X))
    }
  }
  if (nz && ("lgdegf" %in% types || ("egf" %in% types && is.null(result$egf)))) {
    g <- .ggf_try(.ggf_logder(facdiv()), X, maxsqrtn)
    if (!is.null(g)) {
      result$lgdegf <- g
      if (is.null(result$egf)) result$egf <- exp(integrate(g, X))
    }
  }
  result
}

# the first nullspace vector of a rational matrix, normalized as SymPy's
# Matrix.nullspace(): from the reduced row echelon form, the first free
# variable is 1
.q_nullspace1 <- function(rows, ncol) {
  piv <- integer()
  r <- 1L
  for (j in seq_len(ncol)) {
    if (r > length(rows)) break
    k <- which(vapply(r:length(rows), function(i) rows[[i]][j] != 0, TRUE))
    if (!length(k)) next
    k <- k[1L] + r - 1L
    tmp <- rows[[r]]; rows[[r]] <- rows[[k]]; rows[[k]] <- tmp
    rows[[r]] <- rows[[r]] / rows[[r]][j]
    for (i in seq_along(rows)) if (i != r && rows[[i]][j] != 0) rows[[i]] <- rows[[i]] - rows[[i]][j] * rows[[r]]
    piv <- c(piv, j)
    r <- r + 1L
  }
  free <- setdiff(seq_len(ncol), piv)
  if (!length(free)) return(NULL)
  f <- free[1L]
  v <- gmp::as.bigq(rep(0L, ncol))
  v[f] <- gmp::as.bigq(1)
  for (i in seq_along(piv)) v[piv[i]] <- -rows[[i]][f]
  v
}

# SymPy's rational_interpolate() on the points (1, y1), (2, y2), ...,
# without the final cancellation
.guess_rinterp <- function(ydata, degnum, X) {
  yq <- do.call(c, lapply(ydata, function(e) gmp::as.bigq(e$value)))
  n <- length(yq)
  k <- n - degnum - 1L
  nr <- degnum + k + 1L
  nc <- degnum + k + 2L
  rows <- lapply(seq_len(nr), function(i) {
    xi <- gmp::as.bigq(i)
    v <- gmp::as.bigq(rep(1L, nc))
    if (max(degnum, k) >= 1L) for (j in seq_len(max(degnum, k))) v[j + 1L] <- v[j] * xi
    for (j in 0:k) v[degnum + k + 2L - j] <- -v[k - j + 1L] * yq[i]
    v
  })
  r <- .q_nullspace1(rows, nc)
  num <- .add(c(list(S_Zero), lapply(0:degnum, function(i) .mul(list(make_rational(r[i + 1L]), .pow(X, .int(i)))))))
  den <- .add(c(list(S_Zero), lapply(0:k, function(i) .mul(list(make_rational(r[i + degnum + 2L]), .pow(X, .int(i)))))))
  .div(num, den)
}

#' @rdname find_simple_recurrence_vector
#' @export
guess <- function(l, all = FALSE, evaluate = TRUE, niter = 2, variables = NULL) {
  l <- lapply(as.list(l), as_sym)
  if (length(l) > 1L && any(vapply(l[-length(l)], is_zero_num, TRUE))) return(list())
  N <- length(l)
  niter <- min(N - 1L, as.integer(niter))
  myprod <- if (evaluate) product else Product
  g <- list()
  res <- list()
  symb <- if (is.null(variables)) lapply(seq_len(max(niter, 0L)) - 1L, function(i) Symbol(paste0("i", i)))
          else as.list(variables)
  for (k in seq_along(symb) - 1L) {
    s <- symb[[k + 1L]]
    g[[k + 1L]] <- l
    n <- length(l)
    r <- list()
    if (n - 3L >= 0L) for (i in (n - 3L):0L) {
      ri <- .guess_rinterp(g[[k + 1L]][-n], i, s)
      ns <- .int(n)
      if (!is_zero_num(subs(denom(ri), s, ns)) &&
          is_zero_num(.sub(subs(ri, s, ns), g[[k + 1L]][[n]])) &&
          !any(vapply(r, function(e) identical(e$key, ri$key), TRUE)))
        r[[length(r) + 1L]] <- ri
    }
    if (length(r)) {
      if (k >= 1L) for (i in (k - 1L):0L)
        r <- lapply(r, function(v) .mul(list(g[[i + 1L]][[1L]],
          myprod(v, list(symb[[i + 2L]], S_One, .sub(symb[[i + 1L]], S_One))))))
      if (!all) return(r)
      res <- c(res, r)
    }
    l <- lapply(seq_len(N - k - 1L), function(i) .div(l[[i + 1L]], l[[i]]))
  }
  res
}
