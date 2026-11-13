# rsolve() after sympy/solvers/recurr.py: the driver that collects the
# coefficients of y(n + k), clears denominators, shifts the lowest term to
# y(n), calls rsolve_hyper() and fits the initial values.

#' Solve linear recurrences
#'
#' `rsolve(f, y)` solves the linear recurrence `f = 0` (or an [Eq()]) with
#' polynomial or rational coefficients for the sequence term `y`, such as
#' `y(n)`, as SymPy's `rsolve()`. The inhomogeneous part must be a sum of
#' hypergeometric terms. Arbitrary constants are `C0`, `C1`, ... `init`
#' gives initial values as a list of `list(key, value)` pairs, where a key is
#' an index or a term such as `y(0)`, or as a list of values for the indices
#' `0, 1, ...`. The result is `NULL` when no solution is found.
#'
#' @param f The recurrence: an expression equal to zero or an equation.
#' @param y The unknown term, for example `y(n)`.
#' @param init Optional initial values.
#' @param ... For compatibility with the older calling form
#'   `rsolve(eq, y(n), n, init = ...)`: `init` may be given here by name.
#' @return An expression or `NULL`.
#' @examples
#' n <- Symbol("n", integer = TRUE); y <- Function("y")
#' rsolve(y(n + 2) - y(n + 1) - y(n), y(n), list(list(y(0), 0), list(y(1), 1)))
#' rsolve(y(n + 1) - 2 * y(n), y(n))
#' @export
rsolve <- function(f, y, init = NULL, ...) {
  dots <- list(...)
  if (is_sym(init) && is_symbol_like(init)) init <- dots$init
  f <- as_sym(f)
  y <- as_sym(y)
  if (f$type == "Equality") f <- .sub(f$args[[1L]], f$args[[2L]])
  n <- y$args[[1L]]
  yname <- y$name
  is_y <- function(z) z$type == "Function" && identical(z$name, yname)
  has_y <- function(z) {
    if (is_y(z)) return(TRUE)
    for (a in z$args) if (has_y(a)) return(TRUE)
    FALSE
  }
  hk <- integer()
  hv <- list()
  i_part <- list()
  for (g in .aargs(expand(f))) {
    fac <- .margs(g)
    dep <- Filter(has_y, fac)
    coeff <- .mul(Filter(Negate(has_y), fac))
    if (!length(dep)) {
      i_part[[length(i_part) + 1L]] <- coeff
      next
    }
    for (h in dep) {
      k <- NULL
      if (is_y(h)) {
        d <- expand(.sub(h$args[[1L]], n))
        if (!.has_free(d, list(n)) && length(free_symbols(d)) == 0L && is_exact_number(d))
          k <- as.integer(trunc(as.numeric(d$value)))
      }
      if (is.null(k))
        .sp_stop("ValueError", paste0("'", yname, "(", .sp_py_str(n), " + k)' expected, got '", .sp_py_str(h), "'"))
      j <- match(k, hk)
      if (is.na(j)) {
        hk <- c(hk, k)
        hv[[length(hv) + 1L]] <- list(coeff)
      } else hv[[j]] <- c(hv[[j]], list(coeff))
    }
  }
  hv <- lapply(hv, function(cs) simplify(.add(cs)))
  i_part <- .add(i_part)
  common <- S_One
  if (!is_zero_num(i_part) && !.is_hypergeometric(i_part, n) &&
      !(i_part$type == "Add" && all(vapply(.aargs(expand(i_part)), function(x) .is_hypergeometric(x, n), TRUE))))
    .sp_stop("ValueError", paste0("The independent term should be a sum of hypergeometric functions, got '",
                                  .sp_py_str(i_part), "'"))
  for (coeff in hv) {
    if (.sp_is_rational_function(coeff, n)) {
      if (!is_polynomial(coeff, n)) common <- lcm(common, .as_numer_denom(coeff)[[2L]], n)
    } else .sp_stop("ValueError", paste0("Polynomial or rational function expected, got '", .sp_py_str(coeff), "'"))
  }
  ind <- .as_numer_denom(i_part)
  if (is_polynomial(ind[[2L]], n)) common <- lcm(common, ind[[2L]], n)
  if (!is_one_num(common)) {
    hv <- lapply(hv, function(coeff) {
      nd <- .as_numer_denom(coeff)
      .mul(list(nd[[1L]], quo(common, nd[[2L]], n)))
    })
    i_part <- .mul(list(ind[[1L]], quo(common, ind[[2L]], n)))
  }
  K_min <- min(hk)
  if (K_min < 0L) {
    K <- abs(K_min)
    i_part <- expand(subs(i_part, n, .add(list(n, .int(K)))))
    hv <- lapply(hv, function(coeff) expand(subs(coeff, n, .add(list(n, .int(K))))))
    hk <- hk + K
  }
  K_max <- max(hk)
  coeffs <- lapply(0:K_max, function(i) { j <- match(i, hk); if (is.na(j)) S_Zero else hv[[j]] })
  res <- .rsolve_hyper(coeffs, .neg(i_part), n)
  if (is.null(res)) return(NULL)
  solution <- res$result
  symbols <- res$C
  if (!length(init)) init <- NULL
  if (length(symbols) && !is.null(init)) {
    pairs <- if (all(vapply(init, function(p) is.list(p) && !is_sym(p) && length(p) == 2L, TRUE))) init else
      lapply(seq_along(init), function(i) list(i - 1L, init[[i]]))
    equations <- list()
    for (p in pairs) {
      k <- p[[1L]]
      v <- as_sym(p[[2L]])
      if (is.numeric(k)) {
        i <- as.integer(k)
      } else {
        k <- as_sym(k)
        if (is_exact_number(k)) i <- as.integer(trunc(as.numeric(k$value)))
        else if (is_y(k)) i <- as.integer(trunc(as.numeric(k$args[[1L]]$value)))
        else .sp_stop("ValueError", paste0("Integer or term expected, got '", .sp_py_str(k), "'"))
      }
      eq <- .sub(subs(solution, n, .int(i)), v)
      if (.has_type(eq, "Constant") && .sp_has_nan(eq)) eq <- .sub(limit(solution, n, .int(i)), v)
      equations[[length(equations) + 1L]] <- eq
    }
    result <- .sp_solve(equations, symbols, list(dict = TRUE, "_raw" = TRUE))
    if (!length(result)) return(NULL)
    solution <- .sd_subs(solution, result[[1L]])
  }
  solution
}

.sp_has_nan <- function(e) {
  if (.sp_nan(e)) return(TRUE)
  for (a in e$args) if (.sp_has_nan(a)) return(TRUE)
  FALSE
}

# Expr.is_rational_function(n)
.sp_is_rational_function <- function(e, n) {
  if (!.has_free(e, list(n))) return(TRUE)
  nd <- .as_numer_denom(e)
  is_polynomial(nd[[1L]], n) && is_polynomial(nd[[2L]], n)
}
