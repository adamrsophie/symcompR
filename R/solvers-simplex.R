# Linear programming by the simplex method, after sympy/solvers/simplex.py.

.lp_err <- function(cls, msg) stop(.sympy_condition(cls, paste0("\n", msg)))

# the sign of a Rational or Float
.lp_sgn <- function(x) {
  if (x$type %in% c("Integer", "Rational")) return(sign(as.numeric(sign(x$value))))
  sign(num_value_dbl(x))
}
.lp_cmp <- function(a, b) .lp_sgn(.sub(a, b))

# Bland's rule key of a variable label list(is_y, index)
.lp_lab_cmp <- function(a, b) {
  if (a[[1L]] != b[[1L]]) return(if (a[[1L]]) 1L else -1L)
  sign(a[[2L]] - b[[2L]])
}

.lp_pivot <- function(M, i, j) {
  Mi <- M[i, ]; Mj <- M[, j]; Mij <- M[[i, j]]
  if (.lp_sgn(Mij) == 0) stop(.sympy_condition("ZeroDivisionError", "Tried to pivot about zero-valued entry."))
  q <- lapply(Mi, function(z) .div(z, Mij))
  A <- M
  for (r in seq_len(nrow(M))) for (cc in seq_len(ncol(M)))
    A[[r, cc]] <- .sub(M[[r, cc]], .mul(list(Mj[[r]], q[[cc]])))
  for (cc in seq_len(ncol(M))) A[[i, cc]] <- q[[cc]]
  for (r in seq_len(nrow(M))) A[[r, j]] <- .div(.neg(Mj[[r]]), Mij)
  A[[i, j]] <- .div(S_One, Mij)
  A
}

# a list matrix of numbers from a Matrix, a list of rows, a vector or NULL
.lp_mat <- function(x, col = FALSE) {
  if (inherits(x, "symmatrix")) return(matrix(x$data, x$nrow, x$ncol, byrow = TRUE))
  if (is.matrix(x)) return(matrix(lapply(x, as_sym), nrow(x), ncol(x)))
  if (is.null(x) || !length(x)) return(matrix(list(), 0L, 0L))
  if (is.list(x) && !is_sym(x) && any(vapply(x, function(r) (is.list(r) && !is_sym(r)) || (is.atomic(r) && length(r) > 1L), TRUE))) {
    rows <- lapply(x, function(r) lapply(as.list(r), as_sym))
    nc <- max(vapply(rows, length, 0L))
    if (nc == 0L) return(matrix(list(), length(rows), 0L))
    return(matrix(unlist(lapply(rows, identity), recursive = FALSE), length(rows), nc, byrow = TRUE))
  }
  v <- lapply(as.list(x), as_sym)
  matrix(v, length(v), 1L)
}

# Matrix(x) of an argument, which fails for None
.lp_matrix_arg <- function(x) {
  if (is.null(x))
    stop(.sympy_condition("TypeError", "\nData type not understood; expecting list of lists or lists of values."))
  .lp_mat(x)
}

.lp_t <- function(M) if (length(M)) t(M) else matrix(list(), ncol(M), nrow(M))
.lp_negm <- function(M) { if (length(M)) M[] <- lapply(M, .neg); M }

# sympy.solvers.simplex._simplex
.simplex <- function(A, B, C, D = NULL, dual = FALSE) {
  A <- .lp_mat(A); B <- .lp_mat(B); C <- .lp_mat(C)
  D <- if (is.null(D) || !length(D)) matrix(list(S_Zero), 1L, 1L) else .lp_mat(D)
  if (dual) {
    r <- .simplex(.lp_negm(.lp_t(A)), .lp_t(C), .lp_t(B), .lp_negm(D))
    return(list(.neg(r[[1L]]), r[[2L]], r[[3L]]))
  }
  if (length(A) && length(B)) {
    M <- rbind(cbind(A, B), cbind(C, D))
  } else {
    if (length(A) || length(B)) stop(.sympy_condition("ValueError", "must give A and B"))
    M <- cbind(C, D)
  }
  n <- ncol(M) - 1L
  m <- nrow(M) - 1L
  if (!all(vapply(M, function(z) z$type %in% c("Integer", "Rational", "Float"), TRUE)))
    .lp_err("TypeError", "Only rationals and floats are allowed.")
  X <- lapply(seq_len(n) - 1L, function(j) list(FALSE, j))
  Y <- lapply(seq_len(m) - 1L, function(i) list(TRUE, i))
  bland <- function(cols) {
    best <- cols[1L]
    for (cc in cols[-1L]) if (.lp_lab_cmp(X[[cc]], X[[best]]) < 0) best <- cc
    best
  }
  choose_row <- function(rows, cc) {
    best <- rows[1L]
    for (r in rows[-1L]) {
      d <- .lp_cmp(.div(M[[r, n + 1L]], M[[r, cc]]), .div(M[[best, n + 1L]], M[[best, cc]]))
      if (d < 0 || (d == 0 && .lp_lab_cmp(Y[[r]], Y[[best]]) < 0)) best <- r
    }
    best
  }
  repeat {
    neg <- which(vapply(seq_len(m), function(i) .lp_sgn(M[[i, n + 1L]]) < 0, TRUE))
    if (!length(neg)) break
    k <- neg[1L]
    piv_cols <- which(vapply(seq_len(n), function(j) .lp_sgn(M[[k, j]]) < 0, TRUE))
    if (!length(piv_cols)) .lp_err("InfeasibleLPError", "The constraint set is empty!")
    cc <- bland(piv_cols)
    piv_rows <- c(which(vapply(seq_len(m), function(i) .lp_sgn(M[[i, cc]]) > 0 && .lp_sgn(M[[i, n + 1L]]) >= 0, TRUE)), k)
    r <- choose_row(piv_rows, cc)
    M <- .lp_pivot(M, r, cc)
    tmp <- X[[cc]]; X[[cc]] <- Y[[r]]; Y[[r]] <- tmp
  }
  repeat {
    piv_cols <- which(vapply(seq_len(n), function(j) .lp_sgn(M[[m + 1L, j]]) < 0, TRUE))
    if (!length(piv_cols)) break
    cc <- bland(piv_cols)
    piv_rows <- which(vapply(seq_len(m), function(i) .lp_sgn(M[[i, cc]]) > 0, TRUE))
    if (!length(piv_rows)) .lp_err("UnboundedLPError", "Objective function can assume arbitrarily large values!")
    r <- choose_row(piv_rows, cc)
    M <- .lp_pivot(M, r, cc)
    tmp <- X[[cc]]; X[[cc]] <- Y[[r]]; Y[[r]] <- tmp
  }
  argmax <- vector("list", n)
  argmin_dual <- vector("list", m)
  for (i in seq_len(n)) {
    v <- X[[i]]
    if (!v[[1L]]) argmax[[v[[2L]] + 1L]] <- S_Zero else argmin_dual[[v[[2L]] + 1L]] <- M[[m + 1L, i]]
  }
  for (i in seq_len(m)) {
    v <- Y[[i]]
    if (v[[1L]]) argmin_dual[[v[[2L]] + 1L]] <- S_Zero else argmax[[v[[2L]] + 1L]] <- M[[i, n + 1L]]
  }
  list(.neg(M[[m + 1L, n + 1L]]), argmax, argmin_dual)
}

# lts - gts of a non-strict inequality
.lp_nonpos <- function(i) {
  if (i$type == "LessThan") return(.sub(i$args[[1L]], i$args[[2L]]))
  .sub(i$args[[2L]], i$args[[1L]])
}

# sympy.solvers.simplex._rel_as_nonpos; the symbols go in sorted order
.lp_rel_as_nonpos <- function(constr, syms) {
  r <- list(); rk <- list(); np <- list(); aux <- list(); unbound <- character()
  univariate <- list()
  zi <- 0L
  nextu <- function() { zi <<- zi + 1L; Dummy(paste0("z", zi)) }
  skeys <- vapply(syms, .k, "")
  for (i in constr) {
    if (.is_true(i)) next
    if (.is_false(i)) return(NULL)
    if (.has1(i, S_Infinity) || .has1(i, S_NegativeInfinity))
      stop(.sympy_condition("ValueError", "only finite bounds are permitted"))
    if (i$type %in% c("LessThan", "GreaterThan")) {
      e <- .lp_nonpos(i)
      freei <- free_symbols(e)
      if (length(freei) > 1L) np[[length(np) + 1L]] <- e
      else if (length(freei)) {
        x <- freei[[1L]]
        if (.k(x) %in% unbound) next
        ivl <- solve_univariate_inequality(Le(e, 0), x, relational = FALSE)
        univariate[[.k(x)]] <- if (is.null(univariate[[.k(x)]])) ivl else Intersection(univariate[[.k(x)]], ivl)
      } else if (!is_zero_num(e)) return(FALSE)
    } else {
      stop(.sympy_condition("TypeError", filldedent(sprintf(
        "only equalities like Eq(x, y) or non-strict\ninequalities like x >= y are allowed in lp, not %s", sstr(i)))))
    }
  }
  for (x in syms) {
    i <- univariate[[.k(x)]]
    if (!is.null(i) && i$type == "EmptySet") return(NULL)
    if (is.null(i)) { unbound <- c(unbound, .k(x)); next }
    a <- set_inf(i); b <- set_sup(i)
    if (is_inf_const(a)) {
      u <- nextu(); r[[length(r) + 1L]] <- list(x, .sub(b, u)); aux[[length(aux) + 1L]] <- u
    } else if (is_inf_const(b)) {
      if (!is_zero_num(a)) { u <- nextu(); r[[length(r) + 1L]] <- list(x, .add(list(a, u))); aux[[length(aux) + 1L]] <- u }
    } else {
      u <- nextu(); aux[[length(aux) + 1L]] <- u
      r[[length(r) + 1L]] <- list(x, .add(list(u, a)))
      np[[length(np) + 1L]] <- .sub(u, .sub(b, a))
    }
  }
  for (x in syms) if (.k(x) %in% unbound) {
    u <- nextu(); r[[length(r) + 1L]] <- list(x, .sub(u, x)); aux[[length(aux) + 1L]] <- u
  }
  list(np, r, aux)
}

# sympy.polys.matrices.linsolve._lin_eq2dict, for its errors only: whether
# a has a term in the symbols, signaling a nonlinear term
.lp_lin_terms <- function(a, keys) {
  if (a$key %in% keys) return(TRUE)
  if (a$type == "Add") return(any(vapply(a$args, .lp_lin_terms, TRUE, keys = keys)))
  if (a$type == "Mul") {
    got <- FALSE
    for (ai in a$args) if (.lp_lin_terms(ai, keys)) {
      if (got) stop(.sympy_condition("NonlinearError", filldedent(paste0("\n    nonlinear cross-term: ", sstr(a)))))
      got <- TRUE
    }
    return(got)
  }
  if (!any(vapply(free_symbols(a), .k, "") %in% keys)) return(FALSE)
  stop(.sympy_condition("NonlinearError", paste0("nonlinear term: ", sstr(a))))
}

.lp_coeffs <- function(exprs, xx) {
  keys <- vapply(xx, .k, "")
  for (e in exprs) .lp_lin_terms(e, keys)
  rows <- lapply(exprs, function(e) {
    e <- expand(e)
    cs <- lapply(xx, function(s) .linear_coeff(e, s))
    rest <- expand(.sub(e, .add(c(list(S_Zero), Map(function(c1, s) .mul(list(c1, s)), cs, xx)))))
    list(cs, .neg(rest))
  })
  A <- if (length(rows)) matrix(unlist(lapply(rows, `[[`, 1L), recursive = FALSE), length(rows), length(xx), byrow = TRUE)
       else matrix(list(), 0L, length(xx))
  B <- matrix(lapply(rows, `[[`, 2L), length(rows), 1L)
  list(A, B)
}

.lp <- function(how, f, constr) {
  F <- as_sym(f)
  np <- lapply(as.list(constr), as_sym)
  syms <- ordered(unique_syms(lapply(c(list(F), np), free_symbols)))
  n0 <- length(np)
  for (i in seq_len(n0)) if (np[[i]]$type == "Equality") {
    d <- .sub(np[[i]]$args[[1L]], np[[i]]$args[[2L]])
    np[[i]] <- Le(d, 0)
    np[[length(np) + 1L]] <- Le(.neg(d), 0)
  }
  res <- .lp_rel_as_nonpos(np, syms)
  if (is.null(res)) .lp_err("InfeasibleLPError", "Inconsistent/False constraint")
  if (isFALSE(res)) stop(.sympy_condition("TypeError", "cannot unpack non-iterable bool object"))
  r <- res[[2L]]; aux <- res[[3L]]
  F <- xreplace(F, r)
  np <- lapply(res[[1L]], function(e) xreplace(e, r))
  xx <- c(syms, aux)
  AB <- .lp_coeffs(np, xx)
  CD <- .lp_coeffs(list(F), xx)
  if (how == "max") {
    s <- .simplex(AB[[1L]], AB[[2L]], .lp_negm(CD[[1L]]), .lp_negm(CD[[2L]]))
    o <- .neg(s[[1L]])
  } else {
    s <- .simplex(AB[[1L]], AB[[2L]], CD[[1L]], CD[[2L]])
    o <- s[[1L]]
  }
  p <- s[[2L]]
  keys <- xx
  vals <- p
  if (length(r)) {
    rv <- lapply(r, function(kv) xreplace(kv[[2L]], Map(list, xx, p)))
    for (j in seq_along(r)) {
      hit <- which(vapply(keys, function(z) identical(z$key, r[[j]][[1L]]$key), TRUE))
      vals[[hit]] <- rv[[j]]
    }
    auxk <- vapply(aux, .k, "")
    keep <- which(!vapply(keys, .k, "") %in% auxk)
    keys <- keys[keep]; vals <- vals[keep]
    ord <- ordered(keys)
    vals <- vals[match(vapply(ord, .k, ""), vapply(keys, .k, ""))]
    keys <- ord
  }
  list(o, stats::setNames(vals, vapply(keys, sstr, "")))
}

#' Linear programming
#'
#' `lpmin(f, constr)` and `lpmax(f, constr)` minimize and maximize the
#' linear function `f` subject to the linear constraints `constr` (a list
#' of non-strict inequalities and equalities) by the simplex method, as
#' SymPy's functions of the same names; they return `list(optimum,
#' argopt)` with the values of the symbols in a named list. Infeasible
#' problems signal an `InfeasibleLPError` and unbounded ones an
#' `UnboundedLPError` condition. `linprog(c, A, b, A_eq, b_eq, bounds)`
#' minimizes `c x` subject to `A x <= b`, `A_eq x = b_eq` and the bounds
#' (by default `x >= 0`) in matrix form, and `show_linprog()` returns the
#' objective and constraints that `linprog()` would solve.
#'
#' @param f The objective, linear in its symbols.
#' @param constr A list of constraints.
#' @param c Objective coefficients.
#' @param A,b Inequality constraints `A x <= b` (lists of rows, vectors or
#'   matrices).
#' @param A_eq,b_eq Equality constraints.
#' @param bounds `NULL`, a pair `list(lo, hi)` for all variables, a list of
#'   pairs, or a named list of pairs keyed by 0 based index; `NULL` in a pair
#'   means no bound.
#' @return A list of the optimum and the solution.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' lpmax(x - y, list(Le(x, y + 2), Ge(x, 0), Ge(y, 0)))
#' lpmin(x, list(Ge(y, 1), Ge(x, y)))
#' linprog(c(1, -1), list(c(1, 1)), 2, bounds = list(`1` = list(NULL, NULL)))
#' @export
lpmin <- function(f, constr) .lp("min", f, constr)

#' @rdname lpmin
#' @export
lpmax <- function(f, constr) .lp("max", f, constr)

# sympy.solvers.simplex._handle_bounds
.lp_handle_bounds <- function(bounds) {
  make_list <- function(len, pairs) {
    li <- rep(list(S_Zero), len)
    for (pr in pairs) { idx <- if (pr[[1L]] < 0) len + pr[[1L]] + 1L else pr[[1L]] + 1L; li[[idx]] <- as_sym(pr[[2L]]) }
    li
  }
  unbound <- integer(); row <- list(); row2 <- list()
  b_len <- length(bounds)
  for (x in seq_along(bounds) - 1L) {
    a <- bounds[[x + 1L]][[1L]]; b <- bounds[[x + 1L]][[2L]]
    if (is.null(a) && is.null(b)) unbound <- c(unbound, x)
    else if (is.null(a)) {
      b_len <- b_len + 1L
      row <- c(row, list(make_list(b_len, list(list(x, 1), list(-1, 1))), make_list(b_len, list(list(x, -1), list(-1, -1)))))
      row2 <- c(row2, list(as_sym(b), .neg(as_sym(b))))
    } else if (is.null(b)) {
      if (!is_zero_num(as_sym(a))) {
        b_len <- b_len + 1L
        row <- c(row, list(make_list(b_len, list(list(x, 1), list(-1, -1))), make_list(b_len, list(list(x, -1), list(-1, 1)))))
        row2 <- c(row2, list(as_sym(a), .neg(as_sym(a))))
      }
    } else {
      b_len <- b_len + 1L
      row <- c(row, list(make_list(b_len, list(list(x, 1), list(-1, -1))), make_list(b_len, list(list(x, -1), list(-1, 1))),
                         make_list(b_len, list(list(-1, 1)))))
      row2 <- c(row2, list(as_sym(a), .neg(as_sym(a)), .sub(as_sym(b), as_sym(a))))
    }
  }
  for (x in unbound) {
    b_len <- b_len + 2L
    row <- c(row, list(make_list(b_len, list(list(x, 1), list(-1, 1), list(-2, -1))), make_list(b_len, list(list(x, -1), list(-1, -1), list(-2, 1)))))
    row2 <- c(row2, list(S_Zero, S_Zero))
  }
  rows <- lapply(row, function(rr) c(rr, rep(list(S_Zero), b_len - length(rr))))
  A <- if (length(rows)) matrix(unlist(rows, recursive = FALSE), length(rows), b_len, byrow = TRUE) else matrix(list(), 0L, 0L)
  list(A, matrix(row2, length(row2), 1L))
}

.lp_is_pair <- function(z) is.list(z) && length(z) == 2L && is.null(names(z))
.lp_default_bounds <- function(bounds) {
  is.null(bounds) || (is.list(bounds) && !length(bounds)) ||
    (.lp_is_pair(bounds) && !is.list(bounds[[1L]]) && !is.null(bounds[[1L]]) && is_zero_num(as_sym(bounds[[1L]])) && is.null(bounds[[2L]]))
}

.lp_bounds_list <- function(bounds, ncols) {
  if (.lp_is_pair(bounds) && !is.list(bounds[[1L]]) && !is.list(bounds[[2L]])) return(rep(list(bounds), ncols))
  if (!is.null(names(bounds)) && all(nzchar(names(bounds))) && all(vapply(bounds, .lp_is_pair, TRUE))) {
    out <- rep(list(list(0, NULL)), ncols)
    for (nm in rev(names(bounds))) out[[as.integer(nm) + 1L]] <- bounds[[nm]]
    return(out)
  }
  if (length(bounds) == ncols && all(vapply(bounds, .lp_is_pair, TRUE))) return(bounds)
  stop(.sympy_condition("ValueError", "unexpected bounds"))
}

.lp_row <- function(c) {
  C <- .lp_mat(c)
  if (nrow(C) != 1L && ncol(C) == 1L) C <- .lp_t(C)
  if (nrow(C) != 1L) stop(.sympy_condition("ValueError", "C must be a single row."))
  C
}

#' @rdname lpmin
#' @export
linprog <- function(c, A = NULL, b = NULL, A_eq = NULL, b_eq = NULL, bounds = NULL) {
  C <- .lp_row(c)
  if (!length(A)) {
    if (length(b)) stop(.sympy_condition("ValueError", "A and b must both be given"))
    A <- matrix(list(), 0L, ncol(C)); b <- matrix(list(), 0L, 1L)
  } else { A <- .lp_mat(A); b <- .lp_mat(b) }
  if (ncol(A) != ncol(C)) stop(.sympy_condition("ValueError", "number of columns in A and C must match"))
  if (is.null(A_eq)) {
    if (!is.null(b_eq)) stop(.sympy_condition("ValueError", "A_eq and b_eq must both be given"))
  } else {
    Ae <- .lp_matrix_arg(A_eq); be <- .lp_matrix_arg(b_eq)
    A <- rbind(A, Ae, .lp_negm(Ae)); b <- rbind(b, be, .lp_negm(be))
  }
  if (!.lp_default_bounds(bounds)) {
    bl <- .lp_bounds_list(bounds, ncol(A))
    hb <- .lp_handle_bounds(bl)
    aux <- ncol(hb[[1L]]) - ncol(A)
    if (length(A)) {
      A <- rbind(cbind(A, matrix(list(S_Zero), nrow(A), aux)), hb[[1L]])
      b <- rbind(b, hb[[2L]])
    } else { A <- hb[[1L]]; b <- hb[[2L]] }
    C <- cbind(C, matrix(list(S_Zero), 1L, aux))
  } else aux <- -ncol(A)
  s <- .simplex(A, b, C)
  p <- s[[2L]]
  keep <- if (aux > 0L) seq_len(length(p) - aux) else if (aux == 0L) integer() else seq_len(min(length(p), -aux))
  list(s[[1L]], p[keep])
}

#' @rdname lpmin
#' @export
show_linprog <- function(c, A = NULL, b = NULL, A_eq = NULL, b_eq = NULL, bounds = NULL) {
  C <- .lp_row(c)
  if (!length(A)) {
    if (length(b)) stop(.sympy_condition("ValueError", "A and b must both be given"))
    A <- matrix(list(), 0L, ncol(C)); b <- matrix(list(), 0L, 1L)
  } else {
    if (!length(b)) stop(.sympy_condition("ValueError", "A and b must both be given"))
    A <- .lp_mat(A); b <- .lp_mat(b)
  }
  if (ncol(A) != ncol(C)) stop(.sympy_condition("ValueError", "number of columns in A and C must match"))
  if (is.null(A_eq)) {
    if (!is.null(b_eq)) stop(.sympy_condition("ValueError", "A_eq and b_eq must both be given"))
    Ae <- matrix(list(), 0L, ncol(C)); be <- matrix(list(), 0L, 1L)
  } else { Ae <- .lp_matrix_arg(A_eq); be <- .lp_matrix_arg(b_eq) }
  bl <- if (.lp_default_bounds(bounds)) rep(list(list(0, NULL)), ncol(C)) else .lp_bounds_list(bounds, ncol(A))
  x <- lapply(seq_len(ncol(A)), function(i) Symbol(paste0("x", i)))
  dot <- function(row) .add(c(list(S_Zero), Map(function(a, s) .mul(list(a, s)), row, x)))
  f <- dot(C[1L, ])
  cons <- c(lapply(seq_len(nrow(A)), function(i) Le(dot(A[i, ]), b[[i, 1L]])),
            lapply(seq_len(nrow(Ae)), function(i) Eq(dot(Ae[i, ]), be[[i, 1L]])))
  for (i in seq_along(bl)) {
    if (!is.null(bl[[i]][[1L]])) cons[[length(cons) + 1L]] <- Ge(x[[i]], bl[[i]][[1L]])
    if (!is.null(bl[[i]][[2L]])) cons[[length(cons) + 1L]] <- Le(x[[i]], bl[[i]][[2L]])
  }
  list(f, cons)
}
