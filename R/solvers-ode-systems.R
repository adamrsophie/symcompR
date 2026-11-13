# Systems of linear ODEs, after sympy/solvers/ode/systems.py: the matrix
# form of a system, its canonical form, the solution of first order systems
# by matrix exponentials (linodesolve), the reduction of higher order
# systems to first order, the division of a system into its connected
# components, and dsolve_system().

.ode_order_err <- function(msg) .ode_stop("ODEOrderError", msg, "ValueError")
.ode_nonlin_err <- function(msg) .ode_stop("ODENonlinearError", msg, c("NonlinearError", "ValueError"))

# the number of nodes of an expression (sympy.core.sorting._nodes)
.ode_nodes <- function(e) {
  if (!is_sym(e)) return(1L)
  if (e$type == "Derivative") {
    n <- .ode_nodes(e$args[[1L]])
    for (vc in .ode_variable_count(e)) n <- n + vc[[2L]]
    return(n)
  }
  if (!length(e$args)) return(1L)
  1L + sum(vapply(e$args, .ode_nodes, 0L))
}

# Basic.subs with a dict: the larger keys first, then by default_sort_key
.ode_subs_unordered <- function(e, keys, values) {
  if (!length(keys)) return(e)
  nodes <- vapply(keys, .ode_nodes, 0L)
  o <- order(-nodes)
  ks <- keys[o]
  vs <- values[o]
  no <- nodes[o]
  idx <- integer()
  for (n in unique(no)) {
    grp <- which(no == n)
    if (length(grp) > 1L) {
      sk <- .key_sort(lapply(ks[grp], .sort_key), as.list(grp))
      grp <- unlist(sk)
    }
    idx <- c(idx, grp)
  }
  pairs <- lapply(idx, function(i) list(ks[[i]], vs[[i]]))
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) subs(z, pairs)))
  subs(e, pairs)
}

.ode_exp_atoms <- function(e) {
  a <- if (inherits(e, "symmatrix")) do.call(c, lapply(e$data, function(z) atoms(z, "Function"))) else atoms(e, "Function")
  a <- Filter(function(z) identical(z$name, "exp"), a)
  if (length(a)) a[!duplicated(vapply(a, .k, ""))] else a
}

# Matrix.is_zero_matrix: TRUE, FALSE, or NA when undecided
.ode_is_zero_matrix <- function(M) {
  if (is.null(M)) return(NA)
  z <- vapply(M$data, function(v) {
    r <- .sp_fact(v, "zero")
    if (is.null(r) || is.na(r)) NA else isTRUE(r)
  }, NA)
  if (any(!is.na(z) & !z)) return(FALSE)
  if (all(!is.na(z) & z)) return(TRUE)
  NA
}

# systems._get_func_order: the highest order of each function
.ode_get_func_order <- function(eqs, funcs) {
  o <- vapply(funcs, function(fn) max(vapply(eqs, .ode_order, 0L, func = fn)), 0L)
  stats::setNames(o, vapply(funcs, .k, ""))
}

# systems.simpsol
.ode_simpsol <- function(sol, wrt1, wrt2, doit = TRUE) {
  has_any <- function(e, ws) any(vapply(ws, function(w) .has1(e, w), TRUE))
  simpdep <- function(term, wrt1) {
    canonicalise <- function(a) {
      a <- factor_terms(a)
      nd <- .as_numer_denom(a)
      num <- collect(expand_mul(nd[[1L]]), wrt1)
      .div(num, nd[[2L]])
    }
    term <- powsimp(term)
    ea <- .ode_exp_atoms(term)
    .ode_subs_unordered(term, ea, lapply(ea, function(z) fcall("exp", list(canonicalise(z$args[[1L]])))))
  }
  simpcoeff <- function(coeff, wrt2) {
    coeff <- together(coeff)
    if (isTRUE(is_polynomial(coeff))) coeff <- ratsimp(radsimp(coeff))
    fs <- free_symbols(coeff)
    if (!is.null(wrt2)) {
      wk <- vapply(wrt2, .k, "")
      syms <- c(wrt2, .sp_sorted(Filter(function(s) !(s$key %in% wk), fs)))
    } else syms <- .sp_sorted(fs)
    if (length(syms)) coeff <- collect(coeff, syms)
    together(coeff)
  }
  simp_coeff_dep <- function(expr, wrt1, wrt2 = NULL) {
    add_dep <- function(e) e$type == "Add" && has_any(e, wrt1)
    expr <- .map_bottom_up(expr, function(z) {
      if (z$type == "Mul" && any(vapply(z$args, add_dep, TRUE))) expand_mul(z, deep = FALSE) else z
    })
    terms <- if (expr$type == "Add") .sympy_args(expr) else list(expr)
    dk <- character()
    dd <- list()
    dc <- list()
    for (term in terms) {
      cd <- .sp_indep(term, wrt1, as_add = FALSE)
      coeff <- cd[[1L]]
      dep <- simpdep(cd[[2L]], wrt1)
      if (!identical(dep$key, S_One$key)) {
        dep2 <- factor_terms(dep)
        if (!has_any(dep2, wrt1)) {
          coeff <- .mul(list(coeff, dep2))
          dep <- S_One
        }
      }
      j <- match(dep$key, dk)
      if (is.na(j)) {
        dk <- c(dk, dep$key)
        dd[[length(dd) + 1L]] <- dep
        dc[[length(dc) + 1L]] <- coeff
      } else dc[[j]] <- .add(list(dc[[j]], coeff))
    }
    out <- lapply(seq_along(dd), function(i) {
      c1 <- simpcoeff(dc[[i]], wrt2)
      if (!is.null(wrt2)) c1 <- simp_coeff_dep(c1, wrt2)
      .mul(list(c1, dd[[i]]))
    })
    .add(out)
  }
  simprhs <- function(rhs, rk, rv) {
    if (length(rk)) rhs <- .ode_subs_unordered(rhs, rk, rv)
    rhs <- factor_terms(rhs)
    rhs <- simp_coeff_dep(rhs, wrt1, wrt2)
    signsimp(rhs)
  }
  rk <- list()
  rv <- list()
  if (doit) {
    ints <- list()
    for (s in sol) ints <- c(ints, atoms(s, "Integral"))
    if (length(ints)) {
      ints <- ints[!duplicated(vapply(ints, .k, ""))]
      rk <- ints
      rv <- lapply(ints, function(i) doit(factor_terms(i)))
    }
  }
  lapply(sol, function(s) Eq(s$args[[1L]], simprhs(s$args[[2L]], rk, rv)))
}

# systems.linodesolve_type
.ode_linodesolve_type <- function(A, t, b = NULL) {
  is_non_constant <- !.ode_matrix_is_constant(A, t)
  is_non_homogeneous <- !(is.null(b) || isTRUE(.ode_is_zero_matrix(b)))
  type <- paste0("type", 2L * is_non_constant + is_non_homogeneous + 1L)
  match <- list(type_of_equation = type, antiderivative = NULL)
  if (is_non_constant) {
    ca <- .ode_is_commutative_anti_derivative(A, t)
    if (!isTRUE(ca[[2L]]))
      .ode_nie(.filldedent("The system does not have a commutative antiderivative, it cannot be solved by linodesolve."))
    match["antiderivative"] <- list(ca[[1L]])
    m2 <- .ode_first_order_type5_6_subs(A, t, b = b)
    for (nm in names(m2)) match[nm] <- list(m2[[nm]])
  }
  match
}

# systems._first_order_type5_6_subs
.ode_first_order_type5_6_subs <- function(A, t, b = NULL) {
  match <- list()
  ft <- .ode_factor_matrix(A, t)
  is_homogeneous <- is.null(b) || isTRUE(.ode_is_zero_matrix(b))
  if (!is.null(ft)) {
    t_ <- Symbol(paste0(t$name, "_"))
    F_t <- integrate(ft[[1L]], t)
    inverse <- solveset(Eq(t_, F_t), t)
    if (is_sym(inverse) && inverse$type == "FiniteSet" && !.has_type(inverse, "Piecewise") && length(inverse$args) == 1L) {
      A <- ft[[2L]]
      if (!is_homogeneous) {
        b <- .mat_map(b, function(z) .div(z, ft[[1L]]))
        b <- .mat_map(b, function(z) subs(z, t, inverse$args[[1L]]))
      }
      type <- paste0("type", 5L + (!is_homogeneous))
      match <- list(func_coeff = A, tau = F_t, t_ = t_, type_of_equation = type)
      match["rhs"] <- list(b)
    }
  }
  match
}

#' Matrix form of a system of linear ODEs
#'
#' `linear_ode_to_matrix(eqs, funcs, t, order)` writes the system of linear
#' ODEs `eqs` of the given order as
#' `A_n X^(n) = A_(n-1) X^(n-1) + ... + A_0 X + b`, as SymPy's
#' `linear_ode_to_matrix()`.
#'
#' @param eqs A list of linear ODEs.
#' @param funcs A list of the unknown functions.
#' @param t The independent variable.
#' @param order The order of the system.
#' @return A tuple of the list of matrices `A_n, ..., A_0` and the vector `b`.
#' @examples
#' t <- Symbol("t"); f <- Function("f"); g <- Function("g")
#' linear_ode_to_matrix(list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), f(t))),
#'                      list(f(t), g(t)), t, 1)
#' @export
linear_ode_to_matrix <- function(eqs, funcs, t, order) {
  eqs <- lapply(eqs, as_sym)
  funcs <- lapply(funcs, as_sym)
  t <- as_sym(t)
  order <- as.integer(.ode_num(order))
  for (eq in eqs) for (fn in funcs)
    if (.ode_order(eq, fn) > order) .ode_order_err(sprintf("Cannot represent system in %d-order form", order))
  As <- list()
  rhs <- NULL
  for (o in seq(order, 0L)) {
    syms <- lapply(funcs, function(fn) .ode_d(fn, t, o))
    Ab <- tryCatch(.pv_linear_eq_to_matrix(eqs, syms), NonlinearError = function(e) .ode_nonlin_err("The system of ODEs is nonlinear."))
    Ai <- .mat_map(Ab[[1L]], expand_mul)
    b <- Ab[[2L]]
    As[[length(As) + 1L]] <- if (o == order) Ai else .mat_map(Ai, .neg)
    if (o) eqs <- lapply(b$data, .neg)
    else rhs <- b
  }
  .pytuple(list(As, rhs))
}

# systems.matrix_exp
.ode_matrix_exp <- function(A, t) {
  pj <- matrix_exp_jordan_form(A, t)
  .matmul(.matmul(pj[[1L]], pj[[2L]]), inverse(pj[[1L]]))
}

#' Matrix exponential of a constant matrix times t
#'
#' `matrix_exp_jordan_form(A, t)` returns matrices `P` and `expJ` with
#' `exp(A*t) = P * expJ * P^-1`, as SymPy's `matrix_exp_jordan_form()`:
#' `expJ` is block diagonal, with real blocks for complex conjugate pairs of
#' eigenvalues of a real matrix.
#'
#' @param A A square matrix independent of `t`.
#' @param t The symbol.
#' @return A tuple `(P, expJ)`.
#' @examples
#' t <- Symbol("t")
#' matrix_exp_jordan_form(Matrix(list(list(1, 1), list(0, 1))), t)
#' @export
matrix_exp_jordan_form <- function(A, t) {
  A <- Matrix(A)
  t <- as_sym(t)
  N <- A$nrow
  M <- A$ncol
  if (N != M) .ode_stop("ValueError", sprintf("Needed square matrix but got shape (%d, %d)", N, M))
  if (any(vapply(A$data, function(z) .has1(z, t), TRUE))) .ode_stop("ValueError", "Matrix A should not depend on t")
  jf <- jordan_form(A)
  P0 <- jf$P
  J <- jf$J
  n <- J$nrow
  starts <- c(1L, which(vapply(seq_len(n - 1L), function(i) is_zero_num(.mget(J, i, i + 1L)), TRUE)) + 1L)
  if (n == 1L) starts <- 1L
  sizes <- diff(c(starts, n + 1L))
  basis <- lapply(seq_len(n), function(j) lapply(seq_len(n), function(i) .mget(P0, i, j)))
  ek <- character()
  ev <- list()
  chains <- list()
  for (bi in seq_along(starts)) {
    e <- .mget(J, starts[bi], starts[bi])
    j <- match(e$key, ek)
    if (is.na(j)) {
      ek <- c(ek, e$key)
      ev[[length(ev) + 1L]] <- e
      chains[[length(chains) + 1L]] <- list()
      j <- length(ev)
    }
    chains[[j]][[length(chains[[j]]) + 1L]] <- basis[starts[bi] + seq_len(sizes[bi]) - 1L]
  }
  ord <- unlist(.key_sort(lapply(ev, .sort_key), as.list(seq_along(ev))))
  isreal <- !any(vapply(A$data, function(z) .has1(z, S_I), TRUE))
  blocks <- list()
  vectors <- list()
  seen <- character()
  fact <- function(k) factorial(.int(k))
  for (ii in ord) {
    e <- ev[[ii]]
    for (chain in chains[[ii]]) {
      nn <- length(chain)
      ec <- fcall("conjugate", list(e))
      if (isreal && !identical(e$key, ec$key) && ec$key %in% ek) {
        if (e$key %in% seen) next
        seen <- c(seen, ec$key)
        exprt <- fcall("exp", list(.mul(list(fcall("re", list(e)), t))))
        imrt <- .mul(list(fcall("im", list(e)), t))
        imb <- list(list(fcall("cos", list(imrt)), fcall("sin", list(imrt))),
                    list(.neg(fcall("sin", list(imrt))), fcall("cos", list(imrt))))
        B <- zeros(2L * nn, 2L * nn)
        for (i in seq_len(2L * nn) - 1L) for (j in seq_len(2L * nn) - 1L) {
          bi_ <- i %/% 2L; bj_ <- j %/% 2L
          v <- if (bj_ >= bi_) .div(.mul(list(imb[[i %% 2L + 1L]][[j %% 2L + 1L]], .pow(t, .int(bj_ - bi_)))), fact(bj_ - bi_)) else S_Zero
          B$data[[i * 2L * nn + j + 1L]] <- v
        }
        blocks[[length(blocks) + 1L]] <- .mat_map(B, function(z) .mul(list(exprt, z)))
        for (i in seq_len(nn)) {
          vectors[[length(vectors) + 1L]] <- lapply(chain[[i]], function(z) fcall("re", list(z)))
          vectors[[length(vectors) + 1L]] <- lapply(chain[[i]], function(z) fcall("im", list(z)))
        }
      } else {
        vectors <- c(vectors, chain)
        B <- zeros(nn, nn)
        for (i in seq_len(nn) - 1L) for (j in seq_len(nn) - 1L)
          B$data[[i * nn + j + 1L]] <- if (j >= i) .div(.pow(t, .int(j - i)), fact(j - i)) else S_Zero
        ex <- fcall("exp", list(.mul(list(e, t))))
        blocks[[length(blocks) + 1L]] <- .mat_map(B, function(z) .mul(list(ex, z)))
      }
    }
  }
  expJ <- do.call(diag_matrix, blocks)
  P <- zeros(N, N)
  for (i in seq_len(N)) for (j in seq_len(N)) P$data[[(i - 1L) * N + j]] <- vectors[[j]][[i]]
  .pytuple(list(P, expJ))
}

#' Solve a system of first order linear ODEs
#'
#' `linodesolve(A, t, b)` solves `X' = A X + b` for the vector of unknown
#' functions `X`, as SymPy's `linodesolve()`: constant coefficient systems
#' by the Jordan form of `A`, and systems with a commutative antiderivative
#' `B` of `A` by `exp(B)`. The arbitrary constants are Dummy symbols.
#'
#' @param A The coefficient matrix.
#' @param t The independent variable.
#' @param b The vector of non-homogeneous terms, or `NULL`.
#' @param B A commutative antiderivative of `A`, or `NULL`.
#' @param type One of `"auto"` and `"type1"` to `"type6"`.
#' @param doit Whether to evaluate the integrals of the solution.
#' @param tau The substitution of the independent variable for the types
#'   `type5` and `type6`.
#' @return A list of the components of the solution.
#' @examples
#' t <- Symbol("t")
#' linodesolve(Matrix(list(list(1, 0), list(0, 2))), t)
#' @export
linodesolve <- function(A, t, b = NULL, B = NULL, type = "auto", doit = FALSE, tau = NULL) {
  if (!inherits(A, "symmatrix"))
    .ode_stop("ValueError", .filldedent("The coefficients of the system of ODEs should be of type Matrix"))
  if (A$nrow != A$ncol) .ode_stop("NonSquareMatrixError", .filldedent("The coefficient matrix must be a square"))
  if (!is.null(b)) {
    if (!inherits(b, "symmatrix"))
      .ode_stop("ValueError", .filldedent("The non-homogeneous terms of the system of ODEs should be of type Matrix"))
    if (A$nrow != b$nrow)
      .ode_stop("ValueError", .filldedent("The system of ODEs should have the same number of non-homogeneous terms and the number of equations"))
  }
  if (!is.null(B)) {
    if (!inherits(B, "symmatrix"))
      .ode_stop("ValueError", .filldedent("The antiderivative of coefficients of the system of ODEs should be of type Matrix"))
    if (B$nrow != B$ncol) .ode_stop("NonSquareMatrixError", .filldedent("The antiderivative of the coefficient matrix must be a square"))
    if (A$nrow != B$nrow)
      .ode_stop("ValueError", .filldedent("The coefficient matrix and its antiderivative should have same dimensions"))
  }
  if (!(type %in% c(paste0("type", 1:6), "auto"))) .ode_stop("ValueError", .filldedent("The input type should be a valid one"))
  t <- as_sym(t)
  n <- A$nrow
  Cvect <- Matrix(lapply(seq_len(n), function(i) Dummy()))
  if (is.null(b) && type %in% c("type2", "type4", "type6")) b <- zeros(n, 1L)
  is_transformed <- !is.null(tau)
  passed_type <- type
  if (type == "auto") {
    si <- .ode_linodesolve_type(A, t, b = b)
    type <- si$type_of_equation
    B <- si$antiderivative
  }
  if (type %in% c("type5", "type6")) {
    is_transformed <- TRUE
    if (passed_type != "auto" && is.null(tau)) {
      si <- .ode_first_order_type5_6_subs(A, t, b = b)
      if (!length(si)) .ode_stop("ValueError", .filldedent(sprintf("The system passed isn't %s.", type)))
      tau <- si$tau
      t <- si$t_
      A <- si$func_coeff
      b <- si$rhs
    }
  }
  intx <- function(z) if (.ode_true(z)) Integral(z, t) else S_Zero
  if (type %in% c("type1", "type2", "type5", "type6")) {
    pj <- matrix_exp_jordan_form(A, t)
    P <- .mat_map(pj[[1L]], simplify)
    J <- pj[[2L]]
    if (type %in% c("type1", "type5")) sol_vector <- .matmul(P, .matmul(J, Cvect))
    else {
      Jinv <- .mat_map(J, function(z) subs(z, t, .neg(t)))
      inner <- .mat_map(.matmul(.matmul(Jinv, inverse(P)), b), intx)
      sol_vector <- .matmul(.matmul(P, J), inner + Cvect)
    }
  } else {
    if (is.null(B)) B <- .ode_is_commutative_anti_derivative(A, t)[[1L]]
    if (type == "type3") sol_vector <- .matmul(matrix_exp(B), Cvect)
    else sol_vector <- .matmul(matrix_exp(B), .mat_map(.matmul(matrix_exp(.mat_map(B, .neg)), b), intx) + Cvect)
  }
  if (is_transformed) sol_vector <- .mat_map(sol_vector, function(z) subs(z, t, tau))
  gens <- .sp_sorted(.ode_exp_atoms(sol_vector))
  sv <- sol_vector$data
  if (type != "type1") sv <- lapply(sv, expand_mul)
  sv <- lapply(sv, function(s) if (length(gens)) collect(s, gens, exact = TRUE) else s)
  if (doit) sv <- lapply(sv, doit)
  sv
}

# systems._matrix_is_constant
.ode_matrix_is_constant <- function(M, t) {
  all(vapply(M$data, function(z) is_zero_num(.sp_indep(z, list(t), as_add = TRUE)[[2L]]), TRUE))
}

#' Canonical form of a system of ODEs
#'
#' `canonical_odes(eqs, funcs, t)` solves the system `eqs` for the highest
#' derivative of each function, as SymPy's `canonical_odes()`; there is one
#' system for each solution.
#'
#' @param eqs A list of ODEs.
#' @param funcs A list of the unknown functions.
#' @param t The independent variable.
#' @return A list of systems, each a list of equations.
#' @examples
#' x <- Symbol("x"); f <- Function("f"); g <- Function("g")
#' canonical_odes(list(Eq(Derivative(f(x), x), g(x)), Eq(Derivative(g(x), x) + f(x), 0)),
#'                list(f(x), g(x)), x)
#' @export
canonical_odes <- function(eqs, funcs, t) {
  eqs <- lapply(eqs, as_sym)
  funcs <- lapply(funcs, as_sym)
  t <- as_sym(t)
  order <- .ode_get_func_order(eqs, funcs)
  hd <- lapply(funcs, function(fn) .ode_d(fn, t, order[[fn$key]]))
  canon <- .sp_solve(eqs, hd, list(dict = TRUE, `_raw` = TRUE))
  lapply(canon, function(d) {
    dk <- vapply(d$k, .k, "")
    lapply(hd, function(h) {
      j <- match(h$key, dk)
      if (is.na(j)) .ode_stop("KeyError", .str(h))
      Eq(h, d$v[[j]])
    })
  })
}

# systems._is_commutative_anti_derivative
.ode_is_commutative_anti_derivative <- function(A, t) {
  B <- .mat_map(A, function(z) integrate(z, t))
  C <- .matmul(B, A) - .matmul(A, B)
  C <- .mat_map(.mat_map(C, expand), factor_terms)
  z <- .ode_is_zero_matrix(C)
  .pytuple(list(B, if (is.na(z)) FALSE else z))
}

# systems._factor_matrix
.ode_factor_matrix <- function(A, t) {
  term <- NULL
  for (el in A$data) {
    tt <- .sp_indep(el, list(t))[[2L]]
    if (.has1(tt, t)) {
      term <- tt
      break
    }
  }
  if (!is.null(term)) {
    Af <- .mat_map(A, function(z) ratsimp(.div(z, term)))
    term <- if (.ode_matrix_is_constant(Af, t)) list(term, Af) else NULL
  }
  term
}

# systems._get_poly_coeffs
.ode_get_poly_coeffs <- function(p, order) {
  cs <- rep(list(S_Zero), order + 1L)
  ms <- monoms(p)
  cf <- coeffs(p)
  for (i in seq_along(cf)) cs[[order - ms[[i]][1L] + 1L]] <- as_sym(cf[[i]])
  cs
}

# systems._is_second_order_type2
.ode_is_second_order_type2 <- function(A, t) {
  term <- .ode_factor_matrix(A, t)
  is_type2 <- FALSE
  if (!is.null(term)) {
    term <- .inv(term[[1L]])
    is_type2 <- isTRUE(is_polynomial(term))
  }
  if (is_type2) {
    p <- Poly(expand(term), t)
    ms <- monoms(p)
    if (ms[[1L]][1L] %in% c(2, 4)) {
      cs <- .ode_get_poly_coeffs(p, 4L)
      a1 <- powdenest(.sqrt(cs[[1L]]), force = TRUE)
      c1 <- powdenest(.sqrt(cs[[5L]]), force = TRUE)
      b1 <- powdenest(.sqrt(.sub(cs[[3L]], .mul(list(S_Two, a1, c1)))), force = TRUE)
      is_type2 <- identical(cs[[2L]]$key, .mul(list(S_Two, a1, b1))$key) && identical(cs[[4L]]$key, .mul(list(S_Two, b1, c1))$key)
      term <- .add(list(.mul(list(a1, .pow(t, S_Two))), .mul(list(b1, t)), c1))
    } else is_type2 <- FALSE
  }
  list(is_type2, term)
}

# systems._match_second_order_type
.ode_match_second_order_type <- function(A1, A0, t, b = NULL) {
  match <- list(type_of_equation = "type0")
  n <- A1$nrow
  if (.ode_matrix_is_constant(A1, t) && .ode_matrix_is_constant(A0, t)) return(match)
  if (isTRUE(.ode_is_zero_matrix(.mat_map(A1 + .mat_map(A0, function(z) .mul(list(z, t))), expand_mul)))) {
    match$type_of_equation <- "type1"
    match$A1 <- A1
  } else if (isTRUE(.ode_is_zero_matrix(A1)) && (is.null(b) || isTRUE(.ode_is_zero_matrix(b)))) {
    st <- .ode_is_second_order_type2(A0, t)
    if (isTRUE(st[[1L]])) {
      term <- st[[2L]]
      abc <- .ode_get_poly_coeffs(Poly(term, t), 2L)
      tsq <- expand(.pow(term, S_Two))
      A <- .mat_map(A0, function(z) ratsimp(.mul(list(z, tsq)))) +
        .mat_map(eye(n, n), function(z) .mul(list(.sub(.div(.pow(abc[[2L]], S_Two), .int(4)), .mul(list(abc[[1L]], abc[[3L]]))), z)))
      tau <- integrate(.inv(term), t)
      t_ <- Symbol(paste0(t$name, "_"))
      match$type_of_equation <- "type2"
      match$A0 <- A
      match[["g(t)"]] <- .sqrt(term)
      match$tau <- tau
      match$is_transformed <- TRUE
      match$t_ <- t_
    }
  }
  match
}

# systems._second_order_subs_type1
.ode_second_order_subs_type1 <- function(A, b, funcs, t) {
  U <- lapply(funcs, function(fn) .sub(.mul(list(t, .ode_d(fn, t))), fn))
  sol <- linodesolve(A, t, .mat_map(b, function(z) .mul(list(t, z))))
  reduced <- lapply(seq_along(U), function(i) Eq(U[[i]], sol[[i]]))
  canonical_odes(reduced, funcs, t)[[1L]]
}

# systems._second_order_subs_type2
.ode_second_order_subs_type2 <- function(A, funcs, t_) {
  new_funcs <- lapply(funcs, function(fn) fcall(paste0(fn$name, "_"), list(t_)))
  rhss <- .matmul(A, Matrix(new_funcs))$data
  new_eqs <- lapply(seq_along(new_funcs), function(i) Eq(.ode_d(new_funcs[[i]], t_, 2L), rhss[[i]]))
  list(new_eqs, new_funcs)
}

# systems._is_euler_system
.ode_is_euler_system <- function(As, t) {
  all(vapply(seq_along(As), function(i) .ode_matrix_is_constant(.mat_map(As[[i]], function(z) ratsimp(.mul(list(z, .pow(t, .int(i - 1L)))))), t), TRUE))
}

# a match dictionary as a pydict with string keys, its order entry keyed
# by the functions
.ode_match_dict <- function(m) {
  if (is.null(m)) return(NULL)
  vals <- lapply(names(m), function(nm) {
    v <- m[[nm]]
    if (nm == "order" && is.list(v) && !is.null(v$keys)) v
    else if (nm == "order" && !is.null(names(v))) v
    else v
  })
  .pydict(as.list(names(m)), vals)
}

# systems._classify_linear_system, as a named list
.ode_cls_linear <- function(eqs, funcs, t, is_canon = FALSE) {
  if (length(funcs) != length(eqs))
    .ode_stop("ValueError", paste0("Number of functions given is not equal to the number of equations ", .ode_list_str(funcs)))
  for (fn in funcs) if (length(fn$args) != 1L)
    .ode_stop("ValueError", paste0("dsolve() and classify_sysode() work with functions of one variable only, not ", .str(fn)))
  order <- .ode_get_func_order(eqs, funcs)
  system_order <- max(order)
  is_higher_order <- system_order > 1L
  is_second_order <- system_order == 2L && all(order == 2L)
  r <- tryCatch({
    canon_eqs <- if (!is_canon) canonical_odes(eqs, funcs, t) else list(eqs)
    if (length(canon_eqs) == 1L) linear_ode_to_matrix(canon_eqs[[1L]], funcs, t, system_order)
    else list(implicit = canon_eqs)
  }, ODENonlinearError = function(e) NULL)
  if (is.null(r)) return(NULL)
  if (!is.null(r$implicit)) return(list(is_implicit = TRUE, canon_eqs = r$implicit))
  As <- r[[1L]]
  b <- r[[2L]]
  is_homogeneous <- isTRUE(.ode_is_zero_matrix(b))
  ordd <- .pydict(funcs, lapply(unname(order), function(o) o))
  match <- list(no_of_equation = length(eqs), eq = eqs, func = funcs, order = ordd, is_linear = TRUE,
                is_homogeneous = is_homogeneous, is_general = TRUE)
  if (!is_homogeneous) match$rhs <- b
  is_constant <- all(vapply(As, .ode_matrix_is_constant, TRUE, t = t))
  if (!is_higher_order) {
    A <- As[[2L]]
    match$func_coeff <- A
    is_constant <- .ode_matrix_is_constant(A, t)
    match$is_constant <- is_constant
    si <- tryCatch(.ode_linodesolve_type(A, t, b = b), NotImplementedError = function(e) NULL)
    if (is.null(si)) return(NULL)
    for (nm in names(si)) match[nm] <- list(si[[nm]])
    antiderivative <- match$antiderivative
    match$antiderivative <- NULL
    if (!is_constant) match["commutative_antiderivative"] <- list(antiderivative)
    return(match)
  }
  match$type_of_equation <- "type0"
  if (is_second_order) {
    ms <- .ode_match_second_order_type(As[[2L]], As[[3L]], t)
    for (nm in names(ms)) match[nm] <- list(ms[[nm]])
    match$is_second_order <- TRUE
  }
  if (match$type_of_equation == "type0" && !is_constant) {
    if (.ode_is_euler_system(As, t)) {
      match$is_transformed <- TRUE
      match$type_of_equation <- "type1"
      match$t_ <- Symbol(paste0(t$name, "_"))
    } else {
      terms <- .ode_factor_matrix(As[[length(As)]], t)
      mids <- if (length(As) > 2L) As[2:(length(As) - 1L)] else list()
      is_jordan <- function(M) {
        JB <- zeros(M$nrow, M$nrow)
        for (i in seq_len(M$nrow)) {
          JB$data[[(i - 1L) * M$nrow + i]] <- .mget(M, 1L, 1L)
          if (i < M$nrow) JB$data[[(i - 1L) * M$nrow + i + 1L]] <- S_One
        }
        all(vapply(seq_along(M$data), function(k) identical(M$data[[k]]$key, JB$data[[k]]$key), TRUE))
      }
      if (all(vapply(mids, function(A) isTRUE(.ode_is_zero_matrix(A)), TRUE)) && !is.null(terms) && !is_jordan(terms[[2L]])) {
        jf <- jordan_form(terms[[2L]])
        match$type_of_equation <- "type2"
        match$J <- jf$J
        match[["f(t)"]] <- terms[[1L]]
        match$P <- jf$P
        match$is_transformed <- TRUE
      }
    }
    if (match$type_of_equation != "type0" && is_second_order) match$is_second_order <- NULL
  }
  match$is_higher_order <- is_higher_order
  match
}

.ode_list_str <- function(xs) paste0("[", paste(vapply(xs, .str, ""), collapse = ", "), "]")

.ode_classify_linear_system <- function(eqs, funcs, t, is_canon = FALSE) {
  m <- .ode_cls_linear(lapply(eqs, as_sym), lapply(funcs, as_sym), as_sym(t), is_canon)
  # the equations are kept as given, a tuple or a list
  if (!is.null(m$eq) && inherits(eqs, "pytuple")) m$eq <- .pytuple(m$eq)
  .ode_match_dict(m)
}

# systems._preprocess_eqs
.ode_preprocess_eqs <- function(eqs) lapply(eqs, function(eq) if (eq$type == "Equality") eq else Eq(eq, S_Zero))

# systems._eqs2dict and _dict2graph: the graph of the functions each
# equation depends on
.ode_eqs2graph <- function(eqs, funcs) {
  fk <- vapply(funcs, .k, "")
  nodes <- list()
  edges <- list()
  orig <- list()
  for (eq in eqs) {
    f1 <- .ode_undef_atoms(eq$args[[1L]])[[1L]]
    f2 <- Filter(function(g) g$key %in% fk && !identical(g$key, f1$key), .ode_undef_atoms(eq$args[[2L]]))
    f2 <- f2[order(match(vapply(f2, .k, ""), fk))]
    nodes[[length(nodes) + 1L]] <- f1
    for (g in f2) edges[[length(edges) + 1L]] <- list(f1, g)
    orig[[f1$key]] <- eq
  }
  list(nodes = nodes, edges = edges, orig = orig)
}

# systems._is_type1
.ode_is_type1 <- function(eqs, funcs, t) {
  r <- tryCatch(linear_ode_to_matrix(eqs, funcs, t, 1L), ODENonlinearError = function(e) NULL, ODEOrderError = function(e) NULL)
  if (is.null(r)) return(FALSE)
  A0 <- r[[1L]][[2L]]
  b <- r[[2L]]
  .ode_matrix_is_constant(A0, t) && isTRUE(.ode_is_zero_matrix(b))
}

# systems._combine_type1_subsystems
.ode_combine_type1_subsystems <- function(subsystem, funcs, t) {
  indices <- which(vapply(seq_along(subsystem), function(i) .ode_is_type1(subsystem[[i]], funcs[[i]], t), TRUE))
  remove <- integer()
  for (ip in seq_along(indices)) {
    i <- indices[ip]
    for (j in indices[-seq_len(ip)]) {
      if (any(vapply(subsystem[[j]], function(eq2) any(vapply(funcs[[i]], function(fn) .has1(eq2, fn), TRUE)), TRUE))) {
        subsystem[[j]] <- c(subsystem[[i]], subsystem[[j]])
        remove <- c(remove, i)
      }
    }
  }
  subsystem[setdiff(seq_along(subsystem), remove)]
}

# systems._component_division
.ode_component_division <- function(eqs, funcs, t) {
  g <- .ode_eqs2graph(eqs, funcs)
  subsystems <- list()
  for (cc in connected_components(list(g$nodes, g$edges))) {
    ck <- vapply(cc, .k, "")
    ed <- Filter(function(e) e[[1L]]$key %in% ck, g$edges)
    sccs <- strongly_connected_components(list(cc, ed))
    subsystem <- lapply(sccs, function(scc) lapply(scc, function(f) g$orig[[f$key]]))
    subsystem <- .ode_combine_type1_subsystems(subsystem, sccs, t)
    subsystems[[length(subsystems) + 1L]] <- subsystem
  }
  subsystems
}

# systems._linear_ode_solver
.ode_linear_ode_solver <- function(match) {
  t <- if (!is.null(match$t_)) match$t_ else match$t
  sv <- linodesolve(match$func_coeff, t, b = match$rhs, B = match$commutative_antiderivative,
                    type = match$type_of_equation, tau = match$tau)
  lapply(seq_along(match$func), function(i) Eq(match$func[[i]], sv[[i]]))
}

# systems._select_equations
.ode_select_equations <- function(eqs, funcs, key = function(f) f) {
  lk <- vapply(eqs, function(e) e$args[[1L]]$key, "")
  lapply(funcs, function(f) {
    j <- match(key(f)$key, lk)
    if (is.na(j)) .ode_stop("KeyError", .str(key(f)))
    Eq(f, eqs[[j]]$args[[2L]])
  })
}

# systems._higher_order_ode_solver
.ode_higher_order_ode_solver <- function(match) {
  eqs <- match$eq
  funcs <- match$func
  t <- match$t
  sysorder <- match$order
  if (inherits(sysorder, "pydict"))
    sysorder <- stats::setNames(vapply(sysorder$values, as.integer, 0L), vapply(sysorder$keys, .k, ""))
  type <- if (is.null(match$type_of_equation)) "type0" else match$type_of_equation
  is_second_order <- isTRUE(match$is_second_order)
  is_transformed <- isTRUE(match$is_transformed)
  is_euler <- is_transformed && type == "type1"
  is_higher_order_type2 <- is_transformed && type == "type2" && !is.null(match$P)
  if (is_second_order) {
    nf <- .ode_second_order_to_first_order(eqs, funcs, t, A1 = match$A1, A0 = match$A0, b = match$rhs, type = type, t_ = match$t_)
  } else {
    nf <- .ode_higher_order_to_first_order(eqs, sysorder, t, funcs = funcs, type = type, J = match$J,
                                           f_t = match[["f(t)"]], P = match$P, b = match$rhs)
  }
  new_eqs <- nf[[1L]]
  new_funcs <- nf[[2L]]
  if (is_transformed && !is.null(match$t_)) t <- match$t_
  if (!is_higher_order_type2) new_eqs <- .ode_select_equations(new_eqs, lapply(new_funcs, function(f) .ode_d(f, t)))
  sol <- NULL
  if (!is_higher_order_type2) sol <- tryCatch(.ode_strong_component_solver(new_eqs, new_funcs, t), NotImplementedError = function(e) NULL)
  if (is.null(sol)) sol <- tryCatch(.ode_component_solver(new_eqs, new_funcs, t), NotImplementedError = function(e) NULL)
  if (is.null(sol)) return(NULL)
  is_second_order_type2 <- is_second_order && type == "type2"
  us <- if (is_transformed) "__" else "_"
  sol <- .ode_select_equations(sol, funcs, key = function(f) fcall(paste0(f$name, us, "0"), list(t)))
  if (is_transformed) {
    if (is_second_order_type2) {
      g_t <- match[["g(t)"]]
      tau <- match$tau
      sol <- lapply(sol, function(s) Eq(s$args[[1L]], .mul(list(subs(s$args[[2L]], t, tau), g_t))))
    } else if (is_euler) {
      t0 <- match$t
      tau <- match$t_
      sol <- lapply(sol, function(s) subs(s, tau, fcall("log", list(t0))))
    } else if (is_higher_order_type2) {
      P <- match$P
      sv <- .matmul(P, Matrix(lapply(sol, function(s) s$args[[2L]])))$data
      sol <- lapply(seq_along(funcs), function(i) Eq(funcs[[i]], sv[[i]]))
    }
  }
  sol
}

# systems._strong_component_solver
.ode_strong_component_solver <- function(eqs, funcs, t) {
  match <- .ode_cls_linear(eqs, funcs, t, is_canon = TRUE)
  sol <- NULL
  if (length(match)) {
    match$t <- t
    if (isTRUE(match$is_higher_order)) sol <- .ode_higher_order_ode_solver(match)
    else if (isTRUE(match$is_linear)) sol <- .ode_linear_ode_solver(match)
    if (is.null(sol) && length(eqs) == 1L) {
      sol <- dsolve(eqs[[1L]], func = funcs[[1L]])
      variables <- free_symbols(eqs[[1L]])
      new_constants <- lapply(seq_len(.ode_order(eqs[[1L]], funcs[[1L]])), function(i) Dummy())
      sol <- constant_renumber(sol, variables = variables, newconstants = new_constants)
      sol <- list(sol)
    }
  }
  sol
}

.ode_get_funcs_from_canon <- function(eqs) lapply(eqs, function(eq) eq$args[[1L]]$args[[1L]])

# systems._weak_component_solver
.ode_weak_component_solver <- function(wcc, t) {
  eqs <- do.call(c, wcc)
  funcs <- .ode_get_funcs_from_canon(eqs)
  sol <- .ode_strong_component_solver(eqs, funcs, t)
  if (length(sol)) return(sol)
  sol <- list()
  for (scc in wcc) {
    funcs <- .ode_get_funcs_from_canon(scc)
    sk <- lapply(sol, function(s) s$args[[1L]])
    sv <- lapply(sol, function(s) s$args[[2L]])
    comp_eqs <- lapply(scc, function(eq) .ode_subs_unordered(eq, sk, sv))
    scc_sol <- .ode_strong_component_solver(comp_eqs, funcs, t)
    if (is.null(scc_sol)) .ode_nie(.filldedent("The system of ODEs passed cannot be solved by dsolve_system."))
    sol <- c(sol, scc_sol)
  }
  sol
}

# systems._component_solver
.ode_component_solver <- function(eqs, funcs, t) {
  sol <- list()
  for (wcc in .ode_component_division(eqs, funcs, t)) sol <- c(sol, .ode_weak_component_solver(wcc, t))
  sol
}

# systems._second_order_to_first_order
.ode_second_order_to_first_order <- function(eqs, funcs, t, type = "auto", A1 = NULL, A0 = NULL, b = NULL, t_ = NULL) {
  if ((type == "type1" && is.null(A1)) || (type == "type2" && is.null(A0)) || (type == "auto" && (is.null(A1) || is.null(A0)))) {
    r <- linear_ode_to_matrix(eqs, funcs, t, 2L)
    A2 <- r[[1L]][[1L]]
    A1 <- r[[1L]][[2L]]
    A0 <- r[[1L]][[3L]]
    b <- r[[2L]]
    if (!(A2 == eye(A2$nrow))) .ode_stop("ValueError", .filldedent("The system must be in its canonical form."))
  }
  if (type == "auto") {
    m <- .ode_match_second_order_type(A1, A0, t)
    type <- m$type_of_equation
    A1 <- m$A1
    A0 <- m$A0
  }
  sys_order <- stats::setNames(rep(2L, length(funcs)), vapply(funcs, .k, ""))
  if (type == "type1") {
    if (is.null(b)) b <- zeros(length(eqs), length(eqs))
    eqs <- .ode_second_order_subs_type1(A1, b, funcs, t)
    sys_order[] <- 1L
  }
  if (type == "type2") {
    if (is.null(t_)) t_ <- Symbol(paste0(t$name, "_"))
    t <- t_
    r <- .ode_second_order_subs_type2(A0, funcs, t_)
    eqs <- r[[1L]]
    funcs <- r[[2L]]
    sys_order <- stats::setNames(rep(2L, length(funcs)), vapply(funcs, .k, ""))
  }
  .ode_higher_order_to_first_order(eqs, sys_order, t, funcs = funcs)
}

# systems._higher_order_type2_to_sub_systems
.ode_higher_order_type2_to_sub_systems <- function(J, f_t, funcs, t, max_order, b = NULL, P = NULL) {
  if (is.null(J) || is.null(f_t) || !.ode_matrix_is_constant(J, t))
    .ode_stop("ValueError", .filldedent("Correctly input for args 'A' and 'f_t' for Linear, Higher Order, Type 2"))
  if (is.null(P) && !is.null(b) && !isTRUE(.ode_is_zero_matrix(b)))
    .ode_stop("ValueError", .filldedent("Provide the keyword 'P' for matrix P in A = P * J * P-1."))
  new_funcs <- Matrix(lapply(funcs, function(f) fcall(paste0(f$name, "__0"), list(t))))
  ne <- .mat_map(new_funcs, function(f) .ode_d(f, t, max_order)) - .mat_map(.matmul(J, new_funcs), function(z) .mul(list(f_t, z)))
  if (!is.null(b) && !isTRUE(.ode_is_zero_matrix(b))) ne <- ne - .matmul(inverse(P), b)
  new_eqs <- canonical_odes(ne$data, new_funcs$data, t)[[1L]]
  list(new_eqs, new_funcs$data)
}

# the coefficients of t^o d^o/dt^o f(log(t)) in the derivatives of f at
# log(t), lowest order first: those of the falling factorial u(u-1)...(u-o+1)
.ode_euler_coeffs <- function(o) {
  u <- Dummy("u")
  p <- expand(.mul(lapply(seq_len(o) - 1L, function(k) .sub(u, .int(k)))))
  lapply(0:o, function(i) coeff(p, u, i))
}

# systems._higher_order_to_first_order
.ode_higher_order_to_first_order <- function(eqs, sys_order, t, funcs = NULL, type = "type0", J = NULL, f_t = NULL, P = NULL, b = NULL) {
  if (type == "type1") {
    t_ <- Symbol(paste0(t$name, "_"))
    new_funcs <- lapply(funcs, function(f) fcall(paste0(f$name, "_"), list(t_)))
    max_order <- max(vapply(funcs, function(f) sys_order[[f$key]], 0L))
    sk <- c(funcs, list(t))
    sv <- c(new_funcs, list(fcall("exp", list(t_))))
    for (o in seq_len(max_order)) {
      cs <- .ode_euler_coeffs(o)
      for (i in seq_along(funcs)) {
        nf <- new_funcs[[i]]
        ex <- .div(.add(lapply(seq_along(cs), function(k) .mul(list(.ode_d(nf, t_, k - 1L), cs[[k]])))), .pow(t, .int(o)))
        sk[[length(sk) + 1L]] <- .ode_d(funcs[[i]], t, o)
        sv[[length(sv) + 1L]] <- ex
      }
    }
    new_eqs <- lapply(eqs, function(eq) .ode_subs_unordered(eq, sk, sv))
    new_sys_order <- stats::setNames(vapply(funcs, function(f) sys_order[[f$key]], 0L), vapply(new_funcs, .k, ""))
    new_eqs <- canonical_odes(new_eqs, new_funcs, t_)[[1L]]
    return(.ode_higher_order_to_first_order(new_eqs, new_sys_order, t_, funcs = new_funcs))
  }
  if (type == "type2") {
    max_order <- max(vapply(funcs, function(f) sys_order[[f$key]], 0L))
    return(.ode_higher_order_type2_to_sub_systems(J, f_t, funcs, t, max_order, P = P, b = b))
  }
  new_funcs <- list()
  for (prev in funcs) {
    nm <- prev$name
    fn <- fcall(paste0(nm, "_0"), list(t))
    new_funcs[[length(new_funcs) + 1L]] <- fn
    sk <- list(prev)
    sv <- list(fn)
    new_eqs <- list()
    so <- sys_order[[prev$key]]
    if (so > 1L) for (i in seq_len(so - 1L)) {
      nf <- fcall(paste0(nm, "_", i), list(t))
      sk[[length(sk) + 1L]] <- .ode_d(prev, t, i)
      sv[[length(sv) + 1L]] <- nf
      new_funcs[[length(new_funcs) + 1L]] <- nf
      prev_f <- sv[[i]]
      new_eqs[[length(new_eqs) + 1L]] <- Eq(.ode_d(prev_f, t), nf)
    }
    eqs <- c(lapply(eqs, function(eq) .ode_subs_unordered(eq, sk, sv)), new_eqs)
  }
  list(eqs, new_funcs)
}

# systems.dsolve_system
.ode_dsolve_system <- function(eqs, funcs = NULL, t = NULL, ics = NULL, doit = FALSE, simplify = TRUE) {
  if (!is.list(eqs) || is_sym(eqs)) .ode_stop("ValueError", .filldedent("List of equations should be passed. The input is not valid."))
  eqs <- .ode_preprocess_eqs(lapply(eqs, as_sym))
  if (!is.null(funcs) && (!is.list(funcs) || is_sym(funcs)))
    .ode_stop("ValueError", .filldedent("Input to the funcs should be a list of functions."))
  if (is.null(funcs)) funcs <- .ode_extract_funcs(eqs)
  funcs <- lapply(funcs, as_sym)
  if (any(vapply(funcs, function(f) length(f$args) != 1L, TRUE)))
    .ode_stop("ValueError", .filldedent("dsolve_system can solve a system of ODEs with only one independent variable."))
  if (length(eqs) != length(funcs)) .ode_stop("ValueError", .filldedent("Number of equations and number of functions do not match"))
  if (!is.null(t) && !(is_sym(t) && t$type == "Symbol")) .ode_stop("ValueError", .filldedent("The independent variable must be of type Symbol"))
  if (is.null(t)) {
    d <- NULL
    for (eq in eqs) if (length(dd <- atoms(eq, "Derivative"))) {
      d <- dd[[1L]]
      break
    }
    t <- Filter(function(s) s$type == "Symbol", atoms(d, "Symbol"))[[1L]]
  }
  sols <- list()
  for (ce in canonical_odes(eqs, funcs, t)) {
    sol <- tryCatch(.ode_strong_component_solver(ce, funcs, t), NotImplementedError = function(e) NULL)
    if (is.null(sol)) sol <- .ode_component_solver(ce, funcs, t)
    sols[[length(sols) + 1L]] <- sol
  }
  if (!length(sols)) return(sols)
  variables <- list()
  for (eq in eqs) variables <- c(variables, free_symbols(eq))
  vk <- unique(vapply(variables, .k, ""))
  variables <- variables[!duplicated(vapply(variables, .k, ""))]
  lapply(sols, function(sol) {
    sol <- .ode_select_equations(sol, funcs)
    sol <- constant_renumber(sol, variables = variables)
    consts <- function() {
      fs <- list()
      for (s in sol) fs <- c(fs, free_symbols(s))
      fs <- fs[!duplicated(vapply(fs, .k, ""))]
      Filter(function(s) !(s$key %in% vk), fs)
    }
    if (length(ics)) {
      sc <- .ode_solve_ics(sol, funcs, consts(), ics)
      sol <- lapply(sol, function(s) .ode_subs_unordered(s, sc$keys, sc$values))
    }
    if (simplify) sol <- .ode_simpsol(sol, list(t), .sp_sorted(consts()), doit = doit)
    sol
  })
}

# ExprWithLimits._eval_subs for Integral: a variable of an indefinite
# integral that is replaced by an expression becomes the upper limit, and
# the integrand is not changed where the variable is bound
.ode_integral_eval_subs <- function(e, old, new) {
  func <- e$args[[1L]]
  limits <- rev(e$args[-1L])
  fk <- vapply(free_symbols(e), .k, "")
  ofk <- vapply(free_symbols(old), .k, "")
  if (!(old$type %in% c("Symbol", "Dummy")) || length(intersect(ofk, fk))) {
    sub_into_func <- TRUE
    for (i in seq_along(limits)) {
      xab <- limits[[i]]$args
      if (length(xab) == 1L && identical(old$key, xab[[1L]]$key))
        xab <- if (.diff_wrt_ok(new)) list(new) else list(old, old)
      limits[[i]] <- do.call(Tuple, c(list(xab[[1L]]), lapply(xab[-1L], function(l) subs(l, old, new))))
      if (length(intersect(vapply(free_symbols(xab[[1L]]), .k, ""), ofk))) {
        sub_into_func <- FALSE
        break
      }
    }
    if (.is_undefined_application(old)) sub_into_func <- TRUE
    if (sub_into_func) func <- subs(func, old, new)
  } else {
    for (i in seq_along(limits)) {
      xab <- limits[[i]]$args
      if (length(xab) == 3L) {
        limits[[i]] <- do.call(Tuple, c(list(xab[[1L]]), lapply(xab[-1L], function(l) subs(l, old, new))))
        if (identical(old$key, xab[[1L]]$key)) break
      }
    }
  }
  new_sym("Integral", c(list(func), rev(limits)))
}
register_subs_hook("Integral", .ode_integral_eval_subs)

# A CRootOf has no free symbols and its index is not a subexpression: subs
# replaces it only as a whole (collect() substitutes a Dummy for the number 0
# when a real characteristic root gives the pattern sin(0))
register_subs_hook("CRootOf", function(e, old, new) if (identical(e$key, old$key)) new else e)
