# Methods of sums and products over their limits: ExprWithIntLimits
# (index, reorder, reorder_limit, change_index) of
# sympy/concrete/expr_with_intlimits.py, reverse_order, euler_maclaurin and
# the convergence tests of sympy/concrete/summations.py and products.py.

#' Limits of sums and products
#'
#' `sum_index(s, x)` is the position of the index `x` among the limits of a
#' `Sum` or `Product` (from 0, as SymPy's `index`). `reorder_limit(s, i, j)`
#' swaps two limits and `sym_reorder(s, ...)` (SymPy's `reorder`) applies swaps given as pairs of
#' positions or indices. `change_index(s, var, trafo, newvar)` substitutes a
#' linear transformation of an index. `reverse_order(s, ...)` reverses the
#' limits at the given positions or indices (Karr's convention).
#' `euler_maclaurin(s, m, n, eps, eval_integral)` returns the Euler-Maclaurin
#' approximation of a sum with `m` explicit terms and `n` correction terms,
#' and an error bound. `is_convergent()` and `is_absolutely_convergent()`
#' apply SymPy's convergence tests; they stop when no test decides.
#'
#' @param s A `Sum` or `Product`.
#' @param x,var An index symbol.
#' @param i,j Positions of limits, from 0.
#' @param ... Pairs (lists or vectors of two) of positions or indices, or
#'   positions or indices to reverse.
#' @param trafo A linear expression in `var`.
#' @param newvar The new index.
#' @param m,n Numbers of explicit terms and of correction terms.
#' @param eps Tolerance.
#' @param eval_integral Evaluate the integral.
#' @return An expression, a list or `TRUE`/`FALSE`.
#' @examples
#' n <- Symbol("n", integer = TRUE); a <- Symbol("a"); b <- Symbol("b")
#' x <- Symbol("x")
#' reverse_order(Sum(x * n, list(n, a, b)), 0)
#' change_index(Sum(x^2, list(x, a, b)), x, x + 1)
#' is_convergent(Sum(1 / n^2, list(n, 1, oo)))
#' @name sum-limits
NULL

.lims <- function(s) lapply(s$args[-1L], function(t) t$args)
.rebuild_limits <- function(s, f, lims) .sum_new(f, lims, s$type)

#' @rdname sum-limits
#' @export
sum_index <- function(s, x) {
  vars <- vapply(.lims(s), function(l) l[[1L]]$key, "")
  w <- which(vars == as_sym(x)$key)
  if (length(w) != 1L) stop("Number of instances of variable not equal to one", call. = FALSE)
  w - 1L
}

.limit_pos <- function(s, i) if (is.numeric(i)) as.integer(i) else sum_index(s, i)

#' @rdname sum-limits
#' @export
reorder_limit <- function(s, i, j) {
  s <- as_sym(s); lims <- .lims(s)
  vars <- vapply(lims, function(l) l[[1L]]$key, "")
  lx <- lims[[i + 1L]]; ly <- lims[[j + 1L]]
  bound <- function(e) any(vapply(free_symbols(e), `[[`, "", "key") %in% vars)
  if (any(vapply(list(lx[[2L]], lx[[3L]], ly[[2L]], ly[[3L]]), bound, NA)))
    stop(.str(s), " could not be reordered: could not interchange the two limits specified.", call. = FALSE)
  lims[[i + 1L]] <- ly; lims[[j + 1L]] <- lx
  .rebuild_limits(s, s$args[[1L]], lims)
}

#' @rdname sum-limits
#' @export
sym_reorder <- function(s, ...) {
  s <- as_sym(s); out <- s
  for (r in list(...)) {
    r <- as.list(r)
    if (length(r) != 2L) stop("Invalid number of arguments", call. = FALSE)
    out <- reorder_limit(out, .limit_pos(s, r[[1L]]), .limit_pos(s, r[[2L]]))
  }
  out
}

#' @rdname sum-limits
#' @export
change_index <- function(s, var, trafo, newvar = NULL) {
  s <- as_sym(s); var <- as_sym(var); trafo <- as_sym(trafo)
  newvar <- if (is.null(newvar)) var else as_sym(newvar)
  lims <- .lims(s)
  alpha <- beta <- NULL
  for (k in seq_along(lims)) {
    l <- lims[[k]]
    if (!identical(l[[1L]]$key, var$key)) next
    cl <- .coeff_list(trafo, var)
    if (length(cl) != 2L) stop("Index transformation is not linear", call. = FALSE)
    alpha <- cl[[1L]]; beta <- cl[[2L]]
    lo <- .add(list(.mul(list(alpha, l[[2L]])), beta)); hi <- .add(list(.mul(list(alpha, l[[3L]])), beta))
    lims[[k]] <- if (!length(free_symbols(alpha))) {
      if (is_one_num(alpha)) list(newvar, lo, hi)
      else if (is_negone_num(alpha)) list(newvar, hi, lo)
      else stop("Linear transformation results in non-linear summation stepsize", call. = FALSE)
    } else list(newvar, hi, lo)
  }
  f <- subs(s$args[[1L]], var, .div(.sub(var, beta), alpha))
  f <- subs(f, var, newvar)
  .rebuild_limits(s, f, lims)
}

#' @rdname sum-limits
#' @export
reverse_order <- function(s, ...) {
  s <- as_sym(s)
  idx <- vapply(list(...), function(i) .limit_pos(s, i), 0L)
  e <- 1L
  lims <- .lims(s)
  for (k in seq_along(lims)) if ((k - 1L) %in% idx) {
    e <- -e
    l <- lims[[k]]
    lims[[k]] <- list(l[[1L]], .add(list(l[[3L]], S_One)), .sub(l[[2L]], S_One))
  }
  f <- if (s$type == "Product") .pow(s$args[[1L]], .int(e)) else .mul(list(.int(e), s$args[[1L]]))
  .rebuild_limits(s, f, lims)
}

#' @rdname sum-limits
#' @export
euler_maclaurin <- function(s, m = 0, n = 0, eps = 0, eval_integral = TRUE) {
  s <- as_sym(s)
  m <- as.integer(m); n <- as.integer(n)
  f <- s$args[[1L]]
  lims <- .lims(s)
  if (length(lims) != 1L) stop("More than 1 limit", call. = FALSE)
  i <- lims[[1L]][[1L]]; a <- lims[[1L]][[2L]]; b <- lims[[1L]][[3L]]
  if (.is_true(Gt(a, b))) {
    if (is_one_num(.sub(a, b))) return(list(S_Zero, S_Zero))
    a0 <- a; a <- .add(list(b, S_One)); b <- .sub(a0, S_One); f <- .neg(f)
  }
  tot <- S_Zero
  ev3 <- function(t) N(fcall("Abs", list(t)), 3)
  if (m) {
    if (b$type == "Integer" && a$type == "Integer") m <- min(m, as.integer(as.character(b$value - a$value + 1)))
    if (!eps || isTRUE(is_polynomial(f, i))) {
      tot <- .add(lapply(seq_len(m) - 1L, function(k) subs(f, i, .add(list(a, .int(k))))))
    } else {
      term <- subs(f, i, a)
      if (!is_zero_num(term)) {
        test <- tryCatch(num_value_dbl(ev3(term)) < eps, error = function(e) NA)
        if (isTRUE(test)) return(list(tot, fcall("Abs", list(term))))
        if (is.na(test)) return(list(term, S_Zero))
      }
      tot <- term
      for (k in seq_len(m - 1L)) {
        term <- subs(f, i, .add(list(a, .int(k))))
        if (isTRUE(tryCatch(num_value_dbl(ev3(term)) < eps, error = function(e) NA)) && !is_zero_num(term))
          return(list(tot, fcall("Abs", list(term))))
        tot <- .add(list(tot, term))
      }
    }
    if (identical(.add(list(.sub(b, a), S_One))$key, .int(m)$key)) return(list(tot, S_Zero))
    a <- .add(list(a, .int(m)))
  }
  xd <- Dummy("x")
  I <- Integral(subs(f, i, xd), list(xd, a, b))
  if (eval_integral) I <- doit(I)
  tot <- .add(list(tot, I))
  fpoint <- function(e) if (is_const(b, "oo")) list(subs(e, i, a), S_Zero) else list(subs(e, i, a), subs(e, i, b))
  fab <- fpoint(f)
  iterm <- .div(.add(fab), S_Two)
  g <- diff(f, i)
  term <- S_Zero
  for (k in seq_len(n + 1L)) {
    gab <- fpoint(g)
    term <- .mul(list(.div(bernoulli(.int(2L * k)), fcall("factorial", list(.int(2L * k)))), .sub(gab[[2L]], gab[[1L]])))
    if (k > n) break
    if (eps && !is_zero_num(term)) {
      te <- N(term, 3)
      if (is_nan(te)) return(list(S_NaN, S_NaN))
      if (isTRUE(tryCatch(num_value_dbl(N(fcall("Abs", list(te)), 3)) < eps, error = function(e) NA))) break
    }
    tot <- .add(list(tot, term))
    g <- diff(g, i, 2)
  }
  list(.add(list(tot, iterm)), fcall("Abs", list(term)))
}

# ---------------------------------------------------------------------------
# Convergence

# Order(expr, (x, oo)).expr
.order_oo_expr <- function(expr, x) {
  o <- .so_new(expr, list(x), list(S_Infinity))
  if (o$type == "Order") o$args[[1L]] else o
}

# expr.as_independent(x, as_Add=False)[1]
.dep_factors <- function(e, x, list = FALSE) {
  fs <- if (e$type == "Mul") Filter(function(t) .has1(t, x), e$args) else if (.has1(e, x)) list(e) else list()
  if (list) fs else .mul(fs)
}

# cond.as_set() for a condition in one variable
.conv_cond_set <- function(cond, x) {
  if (identical(cond$key, S_true$key)) return(Reals)
  if (identical(cond$key, S_false$key)) return(S_EmptySet)
  solve_univariate_inequality(cond, x, relational = FALSE)
}

# Wild matching with the ported matcher: a named list or NULL
.conv_match <- function(e, pat) tryCatch(sym_match(e, pat), error = function(err) NULL)

.lim_seq_safe <- function(e, n) tryCatch(limit_seq(e, n), error = function(err) NULL)

#' @rdname sum-limits
#' @export
is_convergent <- function(s) {
  s <- as_sym(s)
  if (s$type == "Product") return(.product_is_convergent(s))
  if (s$type != "Sum") stop("is_convergent() needs a Sum or a Product", call. = FALSE)
  .sum_is_convergent(s)
}

#' @rdname sum-limits
#' @export
is_absolutely_convergent <- function(s) {
  s <- as_sym(s)
  is_convergent(.sum_new(fcall("Abs", list(s$args[[1L]])), .lims(s), "Sum"))
}

.product_is_convergent <- function(s) {
  term <- s$args[[1L]]; lims <- .lims(s)
  r <- tryCatch(list(v = is_convergent(.sum_new(fcall("log", list(term)), lims, "Sum"))), error = function(e) NULL)
  if (!is.null(r)) return(r$v)
  if (isTRUE(tryCatch(is_absolutely_convergent(.sum_new(.sub(term, S_One), lims, "Sum")), error = function(e) NA)))
    return(TRUE)
  stop("The algorithm to find the product convergence of ", .str(term), " is not yet implemented", call. = FALSE)
}

.sum_is_convergent <- function(s) {
  p <- Wild("p"); q <- Wild("q"); r <- Wild("r")
  lim0 <- s$args[[2L]]$args
  sym <- lim0[[1L]]; lower <- lim0[[2L]]; upper <- lim0[[3L]]
  term <- simplify(s$args[[1L]])
  if (length(free_symbols(term)) > 1L)
    stop("convergence checking for more than one symbol containing series is not handled", call. = FALSE)
  if (isTRUE(ask_fact(lower, "finite")) && isTRUE(ask_fact(upper, "finite"))) return(TRUE)
  if (is_const(lower, "-oo")) {
    if (is_const(upper, "oo"))
      return(.sum_is_convergent(Sum(term, list(sym, 0, oo))) && .sum_is_convergent(Sum(term, list(sym, -oo, 0))))
    term <- simplify(xreplace(term, list(list(sym, .neg(sym)))))
    lower <- .neg(upper); upper <- S_Infinity
  }
  sym_ <- Dummy(sym$name, integer = TRUE, positive = TRUE)
  term <- xreplace(term, list(list(sym, sym_)))
  sym <- sym_
  interval <- Interval(lower, upper)
  if (term$type == "Piecewise") {
    for (pc in term$args) {
      func <- pc$args[[1L]]; cond <- pc$args[[2L]]
      if (identical(cond$key, S_true$key) || is_const(tryCatch(set_sup(.conv_cond_set(cond, sym)), error = function(e) S_Zero), "oo"))
        return(.sum_is_convergent(Sum(func, list(sym, lower, upper))))
    }
    return(TRUE)
  }
  # divergence test
  lv <- .lim_seq_safe(term, sym)
  if (!is.null(lv) && identical(ask_fact(lv, "zero"), FALSE)) return(FALSE)
  lva <- .lim_seq_safe(fcall("Abs", list(term)), sym)
  if (!is.null(lva) && identical(ask_fact(lva, "zero"), FALSE)) return(FALSE)
  order <- .order_oo_expr(term, sym)
  # p-series test
  pt <- .conv_match(order, .pow(sym, p))
  if (!is.null(pt)) {
    pv <- pt[["p"]]
    if (.is_true(Lt(pv, -1))) return(TRUE)
    if (.is_true(Ge(pv, -1))) return(FALSE)
  }
  # comparison with 1/(n^p log(n)^q log(log(n))^r)
  li <- fcall("log", list(.inv(sym)))
  nl <- .conv_match(order, .inv(.mul(list(.pow(sym, p), .pow(li, q), .pow(fcall("log", list(.neg(li))), r)))))
  if (is.null(nl))
    nl <- .conv_match(order, .inv(.mul(list(.pow(sym, p), .pow(.neg(li), q), .pow(fcall("log", list(.neg(li))), r)))))
  if (!is.null(nl)) {
    P <- nl[["p"]]; Q <- nl[["q"]]; R <- nl[["r"]]
    gt1 <- function(v) .is_true(Gt(v, 1)); eq1 <- function(v) is_one_num(v)
    if (gt1(P) || (eq1(P) && gt1(Q)) || (eq1(P) && eq1(Q) && gt1(R))) return(TRUE)
    return(FALSE)
  }
  # limit comparison with 1/n
  lc <- .lim_seq_safe(.mul(list(sym, term)), sym)
  if (!is.null(lc) && !length(free_symbols(lc)) && .is_true(Gt(lc, 0))) return(FALSE)
  # ratio test
  nxt <- xreplace(term, list(list(sym, .add(list(sym, S_One)))))
  ratio <- combsimp(powsimp(.div(nxt, term)))
  lr <- .lim_seq_safe(ratio, sym)
  if (!is.null(lr) && !length(free_symbols(lr)) && !is_nan(lr) && lr$type != "AccumulationBounds") {
    if (.is_true(Gt(fcall("Abs", list(lr)), 1))) return(FALSE)
    if (.is_true(Lt(fcall("Abs", list(lr)), 1))) return(TRUE)
  }
  # Raabe's test
  if (!is.null(lr) && is_one_num(lr)) {
    tv <- gammasimp(.mul(list(sym, .sub(.div(term, subs(term, sym, .add(list(sym, S_One)))), S_One))))
    lv <- .lim_seq_safe(tv, sym)
    if (!is.null(lv) && !length(free_symbols(lv)) && lv$type != "AccumulationBounds") {
      if (.is_true(Gt(lv, 1))) return(TRUE)
      if (.is_true(Lt(lv, 1))) return(FALSE)
    }
  }
  # root test
  le <- .lim_seq_safe(.pow(fcall("Abs", list(term)), .inv(sym)), sym)
  if (!is.null(le) && !length(free_symbols(le)) && le$type != "AccumulationBounds") {
    if (.is_true(Lt(le, 1))) return(TRUE)
    if (.is_true(Gt(le, 1))) return(FALSE)
  }
  # alternating series test
  dv <- .conv_match(term, .mul(list(.pow(S_NegativeOne, .add(list(sym, p))), q)))
  if (!is.null(dv) && !.has1(dv[["p"]], sym)) {
    ok <- tryCatch({
      cd <- continuous_domain(dv[["q"]], sym, interval)
      isTRUE(is_subset(interval, cd)) && isTRUE(is_decreasing(dv[["q"]], interval, sym))
    }, error = function(e) FALSE)
    if (ok) return(TRUE)
  }
  # integral test
  check <- NULL
  maxima <- tryCatch(solveset(diff(term, sym), sym, interval), error = function(e) NULL)
  if (!is.null(maxima)) {
    if (maxima$type == "EmptySet") check <- interval
    else if (maxima$type == "FiniteSet") {
      sup <- set_sup(maxima)
      if (!length(free_symbols(sup))) check <- Interval(sup, set_sup(interval))
    }
  }
  if (!is.null(check)) {
    res <- tryCatch({
      cd <- continuous_domain(term, sym, check)
      if (isTRUE(is_subset(check, cd)) && (isTRUE(is_decreasing(term, check, sym)) || isTRUE(is_decreasing(.neg(term), check, sym)))) {
        iv <- doit(Integral(term, list(sym, lower, upper)))
        if (iv$type != "Integral" && !length(free_symbols(iv))) list(v = !isTRUE(is_infinite(iv)) && !is_nan(iv)) else NULL
      } else NULL
    }, error = function(e) NULL)
    if (!is.null(res)) return(res$v)
  }
  # Dirichlet and bounded times convergent tests
  if (order$type == "Mul") {
    args <- order$args
    m <- Dummy("m", integer = TRUE)
    dirichlet <- function(g) {
      iv <- .lim_seq_safe(doit(Sum(g, list(sym, set_inf(interval), m))), m)
      if (!is.null(iv) && isTRUE(ask_fact(iv, "finite"))) TRUE else NULL
    }
    bounded <- function(g1, g2) {
      lv <- .lim_seq_safe(g1, sym)
      if (!is.null(lv) && (isTRUE(ask_fact(lv, "finite")) ||
          (lv$type == "AccumulationBounds" && isTRUE(ask_fact(.sub(lv$args[[2L]], lv$args[[1L]]), "finite")))))
        if (isTRUE(tryCatch(is_absolutely_convergent(Sum(g2, list(sym, lower, upper))), error = function(e) NA))) return(TRUE)
      NULL
    }
    for (k in seq_len(length(args) - 1L)) {
      for (cmb in utils::combn(length(args), k, simplify = FALSE)) {
        a_n <- .mul(args[cmb]); b_n <- .mul(args[-cmb])
        dz <- tryCatch({
          cd <- continuous_domain(a_n, sym, interval)
          if (isTRUE(is_subset(interval, cd)) && isTRUE(is_decreasing(a_n, interval, sym))) dirichlet(b_n) else NULL
        }, error = function(e) NULL)
        if (!is.null(dz)) return(dz)
        bc <- tryCatch(bounded(a_n, b_n), error = function(e) NULL)
        if (!is.null(bc)) return(bc)
      }
    }
  }
  stop("The algorithm to find the Sum convergence of ", .str(xreplace(term, list(list(sym, lim0[[1L]])))),
       " is not yet implemented", call. = FALSE)
}

# continuous_domain of sympy/calculus/util.py
#' @rdname sum-limits
#' @param f An expression.
#' @param symbol The variable.
#' @param domain A real set.
#' @export
continuous_domain <- function(f, symbol, domain) {
  f <- as_sym(f); symbol <- as_sym(symbol); domain <- as_sym(domain)
  if (!isTRUE(is_subset(domain, Reals))) stop("Domain must be a subset of S.Reals.", call. = FALSE)
  implemented <- c("exp", "log", "Abs", "frac", "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan",
                   "acot", "asec", "acsc", "sinh", "cosh", "tanh", "coth", "sech", "csch", "asinh", "acosh",
                   "atanh", "acoth", "asech", "acsch")
  used <- unique(vapply(Filter(function(g) .has1(g, symbol), atoms(f, "Function")), `[[`, "", "name"))
  if (any(!(used %in% implemented))) stop("Unable to determine the domain of the given function.", call. = FALSE)
  sui <- function(rel) .conv_cond_set(rel, symbol)
  cont <- domain
  for (atom in atoms(f, "Pow")) {
    ex <- atom$args[[2L]]
    den <- .as_numer_denom(ex)[[2L]]
    if (!(isTRUE(ask_fact(ex, "rational")) && isTRUE(ask_fact(den, "odd"))))
      cont <- Intersection(sui(Ge(atom$args[[1L]], 0)), cont)
  }
  cons <- list(log = list(function(a) Gt(a, 0)), asin = list(function(a) Ge(a, -1), function(a) Le(a, 1)),
               acos = list(function(a) Ge(a, -1), function(a) Le(a, 1)), acosh = list(function(a) Ge(a, 1)),
               atanh = list(function(a) Gt(a, -1), function(a) Lt(a, 1)),
               asech = list(function(a) Gt(a, 0), function(a) Le(a, 1)))
  cons_u <- list(asec = list(function(a) Le(a, -1), function(a) Ge(a, 1)),
                 acsc = list(function(a) Le(a, -1), function(a) Ge(a, 1)),
                 acoth = list(function(a) Lt(a, -1), function(a) Gt(a, 1)))
  for (atom in atoms(f, "Function")) {
    a <- atom$args[[1L]]
    if (!is.null(cons[[atom$name]])) {
      for (cf in cons[[atom$name]]) cont <- Intersection(sui(cf(a)), cont)
    } else if (!is.null(cons_u[[atom$name]])) {
      cont <- Intersection(do.call(Union, lapply(cons_u[[atom$name]], function(cf) sui(cf(a)))), cont)
    } else if (atom$name == "acot") {
      cont <- Complement(cont, solveset(a, symbol, Reals))
    } else if (atom$name == "frac") {
      stop("continuous_domain of frac() is not implemented", call. = FALSE)
    }
  }
  Complement(cont, singularities(f, symbol, domain))
}
