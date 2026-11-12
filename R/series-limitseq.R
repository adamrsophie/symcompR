# limit_seq and its helpers, after sympy/series/limitseq.py.

#' Limits of sequences
#'
#' `limit_seq(expr, n)` is the limit of the sequence `expr` as the integer
#' `n` grows, as SymPy's `limit_seq`: it rewrites Fibonacci numbers and
#' factorials, tries the limit at infinity, and otherwise applies the
#' Stolz-Cesaro theorem to the dominant terms of numerator and denominator.
#' Oscillating sequences give `AccumBounds`; `NULL` means no answer was found.
#' `dominant(expr, n)` is the term of a sum that grows fastest in `n`.
#'
#' @param expr An expression in `n`.
#' @param n The index symbol.
#' @param trials How many times to apply Stolz-Cesaro.
#' @return An expression, `AccumBounds` or `NULL`.
#' @examples
#' n <- Symbol("n", integer = TRUE)
#' limit_seq(n / (n + 1), n)
#' limit_seq((-1)^n, n)
#' dominant(5 * n^3 + 4 * n^2 + n + 1, n)
#' @name limit_seq
NULL

# expr.rewrite(fibonacci, GoldenRatio) then expr.rewrite(factorial, subfactorial, gamma)
.limitseq_rewrite <- function(e) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  if (z$name == "fibonacci" && length(z$args) == 1L)
    return(.div(.sub(.pow(S_GoldenRatio, z$args[[1L]]), .inv(.pow(.neg(S_GoldenRatio), z$args[[1L]]))),
                .sub(.mul(list(S_Two, S_GoldenRatio)), S_One)))
  if (z$name %in% c("factorial", "subfactorial")) return(.rewrite_gamma(z))
  z
})

.limit_inf <- function(e, n) {
  tryCatch(.limit_doit(Limit(e, n, oo), deep = FALSE), NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
}

#' @rdname limit_seq
#' @export
dominant <- function(expr, n) {
  expr <- as_sym(expr); n <- as_sym(n)
  ex <- expand(expr, func = TRUE)
  terms <- if (ex$type == "Add") .sympy_args(ex) else list(ex)
  term0 <- terms[[length(terms)]]
  comp <- list(term0)
  for (t in terms[-length(terms)]) {
    r <- .div(term0, t)
    e <- gammasimp(r)
    if (identical(e$key, r$key)) e <- sym_factor(r)
    l <- limit_seq(e, n)
    if (is.null(l)) return(NULL)
    if (isTRUE(ask_fact(l, "zero"))) {
      term0 <- t; comp <- list(term0)
    } else if (!(is_const(l, "oo") || is_const(l, "-oo"))) {
      comp[[length(comp) + 1L]] <- t
    }
  }
  if (length(comp) > 1L) return(NULL)
  term0
}

.limit_seq_trials <- function(expr, n, trials) {
  for (i in seq_len(trials)) {
    if (!.has_type(expr, "Sum")) {
      r <- .limit_inf(expr, n)
      if (!is.null(r)) return(r)
    }
    nd <- .as_numer_denom(expr); num <- nd[[1L]]; den <- nd[[2L]]
    if (!.has1(den, n) || !.has1(num, n)) {
      r <- .limit_inf(doit(expr), n)
      return(r)
    }
    num <- .difference_delta(expand(num), n); den <- .difference_delta(expand(den), n)
    expr <- gammasimp(.div(num, den))
    if (!.has_type(expr, "Sum")) {
      r <- .limit_inf(expr, n)
      if (!is.null(r)) return(r)
    }
    nd <- .as_numer_denom(expr)
    num <- dominant(nd[[1L]], n)
    if (is.null(num)) return(NULL)
    den <- dominant(nd[[2L]], n)
    if (is.null(den)) return(NULL)
    expr <- gammasimp(.div(num, den))
  }
  NULL
}

#' @rdname limit_seq
#' @export
limit_seq <- function(expr, n = NULL, trials = 5L) {
  expr <- as_sym(expr)
  fs <- free_symbols(expr)
  if (is.null(n)) {
    if (length(fs) == 1L) n <- fs[[1L]]
    else if (!length(fs)) return(expr)
    else stop("Expression has more than one variable. Please specify a variable.", call. = FALSE)
  } else {
    n <- as_sym(n)
    if (!(n$key %in% vapply(fs, `[[`, "", "key"))) return(expr)
  }
  expr <- .limitseq_rewrite(expr)
  n_ <- Dummy("n", integer = TRUE, positive = TRUE)
  n1 <- Dummy("n", odd = TRUE, positive = TRUE)
  n2 <- Dummy("n", even = TRUE, positive = TRUE)
  neg_pow <- any(vapply(atoms(expr, "Pow"), function(p)
    isTRUE(ask_fact(p$args[[1L]], "negative")) && .has1(p$args[[2L]], n), NA))
  if (neg_pow || .has_function(expr, c("cos", "sin"))) {
    L1 <- .limit_seq_trials(xreplace(expr, list(list(n, n1))), n1, trials)
    if (!is.null(L1)) {
      L2 <- .limit_seq_trials(xreplace(expr, list(list(n, n2))), n2, trials)
      if (is.null(L2) || !identical(L1$key, L2$key)) {
        if (!is.null(L2) && .is_comparable(L1) && .is_comparable(L2))
          return(AccumBounds(fcall("Min", list(L1, L2)), fcall("Max", list(L1, L2))))
        return(NULL)
      }
    }
  } else {
    L1 <- .limit_seq_trials(xreplace(expr, list(list(n, n_))), n_, trials)
  }
  if (!is.null(L1)) return(L1)
  if (expr$type == "Add") {
    lims <- lapply(expr$args, limit_seq, n = n, trials = trials)
    if (any(vapply(lims, is.null, NA))) return(NULL)
    return(.add(lims))
  }
  if (!.has_type(expr, "Sum")) {
    lim <- .limit_seq_trials(fcall("Abs", list(xreplace(expr, list(list(n, n_))))), n_, trials)
    if (!is.null(lim) && isTRUE(ask_fact(lim, "zero"))) return(S_Zero)
  }
  NULL
}

# difference_delta of sympy/series/limitseq.py, with the _eval_difference_delta
# methods of Sum, Add, Mul (product rule over expr.args) and Pow (n only in the
# exponent)
.difference_delta <- function(expr, n, step = S_One) {
  np <- .add(list(n, step))
  r <- switch(expr$type,
    Sum = {
      lim <- expr$args[[length(expr$args)]]$args
      k <- lim[[1L]]; upper <- lim[[3L]]
      f <- if (length(expr$args) == 2L) expr$args[[1L]]
           else .sum_new(expr$args[[1L]], lapply(expr$args[2:(length(expr$args) - 1L)], function(t) t$args))
      doit(Sum(f, list(k, .add(list(upper, S_One)), subs(upper, n, np))))
    },
    Add = .add(lapply(expr$args, .difference_delta, n = n, step = step)),
    Mul = {
      a <- .sympy_args(expr); arg0 <- a[[1L]]; rest <- .mul(a[-1L])
      .add(list(.mul(list(subs(arg0, n, np), .difference_delta(rest, n, step))),
                .mul(list(.difference_delta(arg0, n, step), rest))))
    },
    Pow = {
      b <- expr$args[[1L]]; e <- expr$args[[2L]]
      if (.has1(e, n) && !.has1(b, n)) .mul(list(.sub(.pow(b, .sub(subs(e, n, np), e)), S_One), expr))
    },
    NULL)
  if (!is.null(r) && !is_zero_num(r)) return(r)
  .sub(subs(expr, n, np), expr)
}
