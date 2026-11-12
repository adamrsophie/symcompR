# Expansion of hypergeometric and Meijer G functions into named special
# functions, a port of sympy/simplify/hyperexpand.py.
#
# A hypergeometric function is reduced in order, matched against a table of
# known closed forms (an "origin" reachable by integer shifts of the
# parameters), and the shifts are carried out with differential operators on
# a basis of functions. Meijer G functions use their own small table, then
# Slater's theorem (sums of hypergeometric functions).

.hx_x <- Dummy("x")
.hx_z0 <- Dummy("z0")

# _mod1: the fractional part of the numeric constant of x
.hx_mod1 <- function(x) {
  if (is_number_type(x)) return(fcall("Mod", list(x, S_One)))
  ca <- .coeff_add_split(x)
  .add(list(fcall("Mod", list(ca[[1L]], S_One)), ca[[2L]]))
}

# sift(params, _mod1): buckets in insertion order, named by the key of the
# fractional part
.hx_sift <- function(params, f = .hx_mod1) {
  out <- list()
  for (p in params) {
    m <- f(p); k <- paste0("k", m$key)
    if (is.null(out[[k]])) out[[k]] <- list(mod = m, vals = list())
    out[[k]]$vals[[length(out[[k]]$vals) + 1L]] <- p
  }
  out
}
.hx_bget <- function(b, mod) b[[paste0("k", mod$key)]]$vals

# list.sort() of parameters that differ by numbers
.hx_sort_num <- function(l, decreasing = FALSE) {
  if (length(l) < 2L) return(l)
  d <- vapply(l, function(x) {
    v <- .sub(x, l[[1L]])
    if (is_number_type(v)) as.double(.neval(v, list())) else NA_real_
  }, 0)
  if (anyNA(d)) return(.hx_sort_key(l))
  l[order(d, decreasing = decreasing, method = "radix")]
}
.hx_sort_key <- function(l) if (length(l) < 2L) l else .key_sort(lapply(l, .sort_key), l)

.hx_nonneg_false <- function(x) identical(ask_fact(x, "negative"), FALSE)

# ---------------------------------------------------------------------------
# Hyper_Function and G_Function: parameter lists of a hypergeometric and a
# Meijer G function

.HF <- function(ap, bq) structure(list(ap = lapply(ap, expand), bq = lapply(bq, expand)), class = "hx_hf")
.hf_call <- function(f, z) hyper(f$ap, f$bq, z)
.hf_gamma <- function(f) sum(vapply(f$ap, function(x) isTRUE(ask_fact(x, "integer")) && isTRUE(ask_fact(x, "negative")), TRUE))
.hf_sizes <- function(f) paste(length(f$ap), length(f$bq))
.hf_xreplace <- function(f, repl) .HF(lapply(f$ap, xreplace, rule = repl), lapply(f$bq, xreplace, rule = repl))
.hf_has <- function(f, x) any(vapply(c(f$ap, f$bq), function(p) .has1(p, x), TRUE))
.hf_str <- function(f) paste0("(", paste(vapply(f$ap, .str, ""), collapse = ", "), "; ", paste(vapply(f$bq, .str, ""), collapse = ", "), ")")

.hf_invariants <- function(f) {
  tr <- function(b) {
    mods <- lapply(b, `[[`, "mod")
    o <- .key_sort(lapply(mods, .sort_key), seq_along(mods))
    paste(vapply(o, function(i) paste0(b[[i]]$mod$key, "#", length(b[[i]]$vals)), ""), collapse = ";")
  }
  paste(.hf_gamma(f), tr(.hx_sift(f$ap)), tr(.hx_sift(f$bq)), sep = "|")
}

# Hyper_Function.difficulty: -1 if func is not reachable from f
.hf_difficulty <- function(f, func) {
  if (.hf_gamma(f) != .hf_gamma(func)) return(-1)
  oa <- .hx_sift(f$ap); ob <- .hx_sift(f$bq); a <- .hx_sift(func$ap); b <- .hx_sift(func$bq)
  diff <- S_Zero
  for (pr in list(list(a, oa), list(b, ob))) {
    bucket <- pr[[1L]]; obucket <- pr[[2L]]
    for (k in union(names(bucket), names(obucket))) {
      if (is.null(bucket[[k]]) || is.null(obucket[[k]]) || length(bucket[[k]]$vals) != length(obucket[[k]]$vals)) return(-1)
      l1 <- .hx_sort_num(bucket[[k]]$vals); l2 <- .hx_sort_num(obucket[[k]]$vals)
      for (i in seq_along(l1)) diff <- .add(list(diff, fcall("Abs", list(.sub(l1[[i]], l2[[i]])))))
    }
  }
  diff
}

.hf_suitable_origin <- function(f) {
  for (a in f$ap) for (b in f$bq) {
    d <- .sub(a, b)
    if (isTRUE(ask_fact(d, "integer")) && identical(ask_fact(d, "negative"), FALSE)) return(FALSE)
  }
  for (a in f$ap) if (is_zero_num(a)) return(FALSE)
  for (b in f$bq) if (isTRUE(ask_fact(b, "integer")) && isTRUE(ask_fact(b, "nonpositive"))) return(FALSE)
  TRUE
}

.GF <- function(an, ap, bm, bq) structure(list(an = lapply(an, expand), ap = lapply(ap, expand), bm = lapply(bm, expand), bq = lapply(bq, expand)), class = "hx_gf")
.gf_call <- function(f, z) meijerg(f$an, f$ap, f$bm, f$bq, z)
.gf_signature <- function(f) paste(length(f$an), length(f$ap), length(f$bm), length(f$bq))

# G_Function.compute_buckets: by fractional part, sorted by distance from the
# first element (descending for an and bq)
.gf_buckets <- function(f) {
  out <- list()
  lists <- list(f$an, f$ap, f$bm, f$bq); flips <- c(TRUE, FALSE, FALSE, TRUE)
  for (i in 1:4) {
    b <- .hx_sift(lists[[i]])
    for (k in names(b)) b[[k]]$vals <- .hx_sort_num(b[[k]]$vals, decreasing = flips[i])
    out[[i]] <- b
  }
  out
}

# ---------------------------------------------------------------------------
# Matrices of sym as R list matrices; B is a column (list), C a row (list)

.hx_zeros <- function(n, m = n) matrix(rep(list(S_Zero), n * m), n, m)
.hx_matrix <- function(rows) {
  n <- length(rows); m <- length(rows[[1L]])
  M <- .hx_zeros(n, m)
  for (i in seq_len(n)) for (j in seq_len(m)) M[[i, j]] <- as_sym(rows[[i]][[j]])
  M
}
.hx_mmap <- function(M, f) {
  if (is.matrix(M)) { for (i in seq_along(M)) M[[i]] <- f(M[[i]]); return(M) }
  lapply(M, f)
}
.hx_subs <- function(M, old, new) .hx_mmap(M, function(e) subs(e, old, new))
.hx_subs_d <- function(M, repl) .hx_mmap(M, function(e) subs(e, repl))
.hx_has_bad <- function(M) any(vapply(if (is.matrix(M)) as.list(M) else M, function(e)
  .has_any(e, list(S_NaN, S_Infinity, S_NegativeInfinity, S_ComplexInfinity)), TRUE))
.hx_rowmul <- function(C, M) lapply(seq_len(ncol(M)), function(j) .add(lapply(seq_along(C), function(i) .mul(list(C[[i]], M[[i, j]])))))
.hx_dot <- function(C, B) .add(lapply(seq_along(C), function(i) .mul(list(C[[i]], B[[i]]))))

.hx_eye <- function(n) { M <- .hx_zeros(n); for (i in seq_len(n)) M[[i, i]] <- S_One; M }

# ---------------------------------------------------------------------------
# Polynomials in .hx_x as coefficient lists, lowest degree first (Poly with
# the coefficients in their domain normal form)

.hx_coeffs <- function(expr, x = .hx_x) {
  e <- expand(as_sym(expr))
  cs <- list()
  for (t in .aargs(e)) {
    ind <- .as_indep(t, list(x), as_add = FALSE)
    dep <- ind[[2L]]
    k <- if (is_one_num(dep)) 0L else if (identical(dep$key, x$key)) 1L else
      if (dep$type == "Pow" && identical(dep$args[[1L]]$key, x$key) && dep$args[[2L]]$type == "Integer") as.integer(.int_val(dep$args[[2L]])) else
        stop("not a polynomial in ", .str(x), ": ", .str(expr), call. = FALSE)
    while (length(cs) < k + 1L) cs[[length(cs) + 1L]] <- list()
    cs[[k + 1L]][[length(cs[[k + 1L]]) + 1L]] <- ind[[1L]]
  }
  if (!length(cs)) return(list(S_Zero))
  out <- lapply(cs, function(l) { c <- .add(l); if (is_number_type(c)) c else .cancel(c) })
  # drop zero leading coefficients
  while (length(out) > 1L && is_zero_num(out[[length(out)]])) out <- out[-length(out)]
  out
}
.hx_poly_expr <- function(cs, x = .hx_x) .add(lapply(seq_along(cs), function(i) .mul(list(cs[[i]], .pow(x, .int(i - 1L))))))

# Operator.apply: sum_k c_k op^k(obj)
.hx_op <- function(cs, info = "") structure(list(coeffs = cs, info = info), class = "hx_op")
.hx_op_from <- function(expr, info = "") .hx_op(.hx_coeffs(expr), info)

.hx_vscale <- function(c, v) if (is.list(v) && !is_sym(v)) lapply(v, function(e) .mul(list(c, e))) else .mul(list(c, v))
.hx_vadd <- function(v, w) if (is.list(v) && !is_sym(v)) lapply(seq_along(v), function(i) .add(list(v[[i]], w[[i]]))) else .add(list(v, w))

.hx_apply <- function(o, obj, op) {
  cs <- o$coeffs
  diffs <- list(obj)
  for (c in cs[-1L]) diffs[[length(diffs) + 1L]] <- op(diffs[[length(diffs)]])
  r <- .hx_vscale(cs[[1L]], diffs[[1L]])
  for (i in seq_along(cs)[-1L]) r <- .hx_vadd(r, .hx_vscale(cs[[i]], diffs[[i]]))
  r
}
.hx_apply_operators <- function(obj, ops, op) {
  res <- obj
  for (o in rev(ops)) res <- .hx_apply(o, res, op)
  res
}

.hx_mult_op <- function(p) .hx_op_from(p, "mult")
.hx_shift_a <- function(ai) {
  if (is_zero_num(ai)) stop("Cannot increment zero upper index.", call. = FALSE)
  .hx_op_from(.add(list(.div(.hx_x, ai), S_One)), "ShiftA")
}
.hx_shift_b <- function(bi) {
  if (is_one_num(bi)) stop("Cannot decrement unit lower index.", call. = FALSE)
  .hx_op_from(.add(list(.div(.hx_x, .sub(bi, S_One)), S_One)), "ShiftB")
}

# Poly(n.all_coeffs()[:-1], A).as_expr(): (n - n(0))/A as a polynomial
.hx_drop_const <- function(cs, A) {
  if (length(cs) < 2L) return(S_Zero)
  .hx_poly_expr(cs[-1L], A)
}

.hx_unshift_a <- function(ap, bq, i, z) {
  ai <- .sub(ap[[i]], S_One); ap <- ap[-i]
  if (is_zero_num(ai)) stop("Cannot decrement unit upper index.", call. = FALSE)
  m <- .mul(c(list(z, ai), lapply(ap, function(a) .add(list(.hx_x, a)))))
  A <- Dummy("A")
  D <- .sub(.mul(list(ai, A)), ai)
  n <- .mul(c(list(D), lapply(bq, function(b) .add(list(D, b, S_NegativeOne)))))
  ncs <- .hx_coeffs(n, A)
  b0 <- .neg(ncs[[1L]])
  if (is_zero_num(b0)) stop("Cannot decrement upper index: cancels with lower", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, A), A, .add(list(.div(.hx_x, ai), S_One)))
  .hx_op_from(.div(.sub(n, m), b0), "UnShiftA")
}

.hx_unshift_b <- function(ap, bq, i, z) {
  bi <- .add(list(bq[[i]], S_One)); bq <- bq[-i]
  if (is_zero_num(bi)) stop("Cannot increment -1 lower index.", call. = FALSE)
  m <- .mul(c(list(.hx_x, .sub(bi, S_One)), lapply(bq, function(b) .add(list(.hx_x, b, S_NegativeOne)))))
  B <- Dummy("B")
  D <- .add(list(.mul(list(.sub(bi, S_One), B)), .neg(bi), S_One))
  n <- .mul(c(list(z), lapply(ap, function(a) .add(list(D, a)))))
  ncs <- .hx_coeffs(n, B)
  b0 <- ncs[[1L]]
  if (is_zero_num(b0)) stop("Cannot increment index: cancels with upper", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, B), B, .add(list(.div(.hx_x, .sub(bi, S_One)), S_One)))
  .hx_op_from(.div(.sub(m, n), b0), "UnShiftB")
}

.hx_mshift_a <- function(bi) .hx_op_from(.sub(bi, .hx_x), "MeijerShiftA")
.hx_mshift_b <- function(bi) .hx_op_from(.add(list(S_One, .neg(bi), .hx_x)), "MeijerShiftB")
.hx_mshift_c <- function(bi) .hx_op_from(.sub(.hx_x, bi), "MeijerShiftC")
.hx_mshift_d <- function(bi) .hx_op_from(.sub(.sub(bi, S_One), .hx_x), "MeijerShiftD")

.hx_munshift_a <- function(an, ap, bm, bq, i, z) {
  bi <- .sub(bm[[i]], S_One); bm <- bm[-i]
  m <- .mul(c(lapply(bm, function(b) .sub(b, .hx_x)), lapply(bq, function(b) .sub(.hx_x, b))))
  A <- Dummy("A")
  D <- .sub(bi, A)
  n <- .mul(c(list(z), lapply(an, function(a) .add(list(D, S_One, .neg(a)))), lapply(ap, function(a) .add(list(.neg(D), a, S_NegativeOne)))))
  ncs <- .hx_coeffs(n, A); b0 <- ncs[[1L]]
  if (is_zero_num(b0)) stop("Cannot decrement upper b index (cancels)", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, A), A, .sub(bi, .hx_x))
  .hx_op_from(.div(.sub(m, n), b0), "MeijerUnShiftA")
}

.hx_munshift_b <- function(an, ap, bm, bq, i, z) {
  ai <- .add(list(an[[i]], S_One)); an <- an[-i]
  m <- .mul(c(list(z), lapply(an, function(a) .add(list(S_One, .neg(a), .hx_x))), lapply(ap, function(a) .add(list(a, S_NegativeOne, .neg(.hx_x))))))
  B <- Dummy("B")
  D <- .add(list(B, ai, S_NegativeOne))
  n <- .mul(c(list(S_One), lapply(bm, function(b) .add(list(.neg(D), b))), lapply(bq, function(b) .sub(D, b))))
  ncs <- .hx_coeffs(n, B); b0 <- ncs[[1L]]
  if (is_zero_num(b0)) stop("Cannot increment upper a index (cancels)", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, B), B, .add(list(S_One, .neg(ai), .hx_x)))
  .hx_op_from(.div(.sub(m, n), b0), "MeijerUnShiftB")
}

.hx_munshift_c <- function(an, ap, bm, bq, i, z) {
  bi <- .sub(bq[[i]], S_One); bq <- bq[-i]
  m <- .mul(c(lapply(bm, function(b) .sub(b, .hx_x)), lapply(bq, function(b) .sub(.hx_x, b))))
  C <- Dummy("C")
  D <- .add(list(bi, C))
  n <- .mul(c(list(z), lapply(an, function(a) .add(list(D, S_One, .neg(a)))), lapply(ap, function(a) .add(list(.neg(D), a, S_NegativeOne)))))
  ncs <- .hx_coeffs(n, C); b0 <- ncs[[1L]]
  if (is_zero_num(b0)) stop("Cannot decrement lower b index (cancels)", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, C), C, .sub(.hx_x, bi))
  .hx_op_from(.div(.sub(m, n), b0), "MeijerUnShiftC")
}

.hx_munshift_d <- function(an, ap, bm, bq, i, z) {
  ai <- .add(list(ap[[i]], S_One)); ap <- ap[-i]
  m <- .mul(c(list(z), lapply(an, function(a) .add(list(S_One, .neg(a), .hx_x))), lapply(ap, function(a) .add(list(a, S_NegativeOne, .neg(.hx_x))))))
  B <- Dummy("B")
  D <- .add(list(ai, S_NegativeOne, .neg(B)))
  n <- .mul(c(list(S_One), lapply(bm, function(b) .add(list(.neg(D), b))), lapply(bq, function(b) .sub(D, b))))
  ncs <- .hx_coeffs(n, B); b0 <- ncs[[1L]]
  if (is_zero_num(b0)) stop("Cannot increment lower a index (cancels)", call. = FALSE)
  n <- subs(.hx_drop_const(ncs, B), B, .add(list(ai, S_NegativeOne, .neg(.hx_x))))
  .hx_op_from(.div(.sub(m, n), b0), "MeijerUnShiftD")
}

# ReduceOrder(ai, bj): cancel upper ai with lower bj when ai - bj is a
# nonnegative integer
.hx_reduce_order_op <- function(ai, bj) {
  n <- .sub(ai, bj)
  if (n$type != "Integer" || num_sign(n) < 0) return(NULL)
  if (isTRUE(ask_fact(bj, "integer")) && isTRUE(ask_fact(bj, "nonpositive"))) return(NULL)
  p <- S_One
  for (k in seq_len(as.integer(.int_val(n))) - 1L)
    p <- .mul(list(p, .div(.add(list(.hx_x, bj, .int(k))), .add(list(bj, .int(k))))))
  .hx_op_from(p, "ReduceOrder")
}

.hx_reduce_meijer <- function(b, a, sign) {
  n <- .sub(b, a)
  if (isTRUE(ask_fact(n, "negative")) || n$type != "Integer") return(NULL)
  p <- S_One
  for (k in seq_len(as.integer(.int_val(n))) - 1L) p <- .mul(list(p, .add(list(.mul(list(.int(sign), .hx_x)), a, .int(k)))))
  .hx_op_from(p, "ReduceOrderMeijer")
}
.hx_meijer_minus <- function(b, a) .hx_reduce_meijer(b, a, -1L)
.hx_meijer_plus <- function(a, b) .hx_reduce_meijer(.sub(S_One, a), .sub(S_One, b), 1L)

.hx_reduce_order_ <- function(ap, bq, gen, key) {
  o1 <- .key_sort(lapply(ap, key), seq_along(ap)); ap <- ap[o1]
  o2 <- .key_sort(lapply(bq, key), seq_along(bq)); bq <- bq[o2]
  nap <- list(); operators <- list()
  for (a in ap) {
    op <- NULL
    for (i in seq_along(bq)) {
      op <- gen(a, bq[[i]])
      if (!is.null(op)) { bq <- bq[-i]; break }
    }
    if (is.null(op)) nap[[length(nap) + 1L]] <- a else operators[[length(operators) + 1L]] <- op
  }
  list(nap, bq, operators)
}

.hx_reduce_order <- function(func) {
  r <- .hx_reduce_order_(func$ap, func$bq, .hx_reduce_order_op, .sort_key)
  list(.HF(r[[1L]], r[[2L]]), r[[3L]])
}

.hx_reduce_order_meijer <- function(func) {
  r1 <- .hx_reduce_order_(func$an, func$bq, .hx_meijer_plus, function(x) .sort_key(.neg(x)))
  r2 <- .hx_reduce_order_(func$bm, func$ap, .hx_meijer_minus, .sort_key)
  list(.GF(r1[[1L]], r2[[2L]], r2[[1L]], r1[[2L]]), c(r1[[3L]], r2[[3L]]))
}

# make_simp: cancel a rational function of z with z as the only generator
.hx_make_simp <- function(z) function(expr) {
  nd <- .as_numer_denom(expr)
  r <- .poly_cancel_in(expand(nd[[1L]]), nd[[2L]], z)
  .div(.mul(list(r[[1L]], r[[2L]])), r[[3L]])
}


.hx_make_derivative_operator <- function(M, z) function(C) {
  dC <- lapply(C, function(e) .mul(list(z, .diff1(e, z))))
  r <- .hx_vadd(dC, .hx_rowmul(C, M))
  lapply(r, .hx_make_simp(z))
}

# ---------------------------------------------------------------------------
# Formula: a hypergeometric function with a basis B, a row C with
# closed form C*B, and the matrix M with z*d/dz B = M*B

.hx_formula <- function(func, z, res, symbols, B = NULL, C = NULL, M = NULL) {
  symbols <- Filter(function(s) .hf_has(func, s), symbols)
  f <- list(func = func, z = z, symbols = symbols, B = B, C = C, M = M)
  if (!is.null(res)) f <- .hx_compute_basis(f, res)
  f
}

# Formula._compute_basis: B is closed_form and its successive z*d/dz
.hx_compute_basis <- function(f, closed_form) {
  afactors <- lapply(f$func$ap, function(a) .add(list(.hx_x, a)))
  bfactors <- lapply(f$func$bq, function(b) .add(list(.hx_x, b, S_NegativeOne)))
  expr <- .sub(.mul(c(list(.hx_x), bfactors)), .mul(c(list(f$z), afactors)))
  cs <- .hx_coeffs(expr)
  deg <- length(cs) - 1L
  n <- deg - 1L
  b <- list(closed_form)
  for (i in seq_len(n)) b[[length(b) + 1L]] <- .mul(list(f$z, .diff1(b[[length(b)]], f$z)))
  f$B <- b
  f$C <- c(list(S_One), rep(list(S_Zero), n))
  M <- .hx_zeros(n + 1L)
  for (i in seq_len(n)) M[[i, i + 1L]] <- S_One
  # last row: -l/LC with l the coefficients below the leading one, constant first
  lc <- cs[[deg + 1L]]
  for (j in seq_len(deg)) M[[n + 1L, j]] <- .neg(.div(cs[[j]], lc))
  f$M <- M
  f
}

.hx_closed_form <- function(f) .hx_dot(f$C, f$B)

# Formula.find_instantiations: the substitutions of the formula's symbols
# that make its parameters congruent to func's
.hx_find_instantiations <- function(f, func) {
  ap <- func$ap; bq <- func$bq
  if (length(ap) != length(f$func$ap) || length(bq) != length(f$func$bq))
    stop("Cannot instantiate other number of parameters", call. = FALSE)
  in_keys <- function(a, l) a$key %in% vapply(l, `[[`, "", "key")
  symbol_values <- lapply(f$symbols, function(a) {
    if (in_keys(a, f$func$ap)) ap
    else if (in_keys(a, f$func$bq)) bq
    else stop("At least one of the parameters of the formula must be equal to ", .str(a), call. = FALSE)
  })
  combos <- .hx_product(symbol_values)
  base_repl <- lapply(combos, function(vals) Map(list, f$symbols, vals))
  abuckets <- .hx_sift(ap); bbuckets <- .hx_sift(bq)
  critical <- rep(list(list(S_Zero)), length(f$symbols))
  result <- list()
  nn <- Dummy("n")
  for (repl in base_repl) {
    rmap <- function(x) xreplace(x, repl)
    symb_a <- .hx_sift(f$func$ap, function(x) .hx_mod1(rmap(x)))
    symb_b <- .hx_sift(f$func$bq, function(x) .hx_mod1(rmap(x)))
    for (pr in list(list(abuckets, symb_a), list(bbuckets, symb_b))) {
      bucket <- pr[[1L]]; obucket <- pr[[2L]]
      for (k in union(names(bucket), names(obucket))) {
        if (is.null(bucket[[k]]) || is.null(obucket[[k]]) || length(bucket[[k]]$vals) != length(obucket[[k]]$vals)) break
        for (si in seq_along(f$symbols)) {
          a <- f$symbols[[si]]
          a0 <- repl[[si]][[2L]]
          if (length(free_symbols(a0))) next
          exprs <- Filter(function(e) .has1(e, a), obucket[[k]]$vals)
          repl0 <- repl; repl0[[si]] <- list(a, .add(list(a0, nn)))
          for (e in exprs) for (target in bucket[[k]]$vals) {
            # solve(expr.xreplace(repl0) - target, _n), linear in _n
            g <- .sub(xreplace(e, repl0), target)
            g0 <- subs(g, nn, S_Zero); g1 <- .sub(subs(g, nn, S_One), g0)
            n0 <- simplify(.neg(.div(g0, g1)))
            if (length(free_symbols(n0))) stop("Value should not be true", call. = FALSE)
            critical[[si]][[length(critical[[si]]) + 1L]] <- n0
          }
        }
      }
    }
    values <- lapply(seq_along(f$symbols), function(si) {
      a0 <- repl[[si]][[2L]]
      vs <- vapply(critical[[si]], function(v) as.double(.neval(v, list())), 0)
      lo <- floor(min(vs)); hi <- ceiling(max(vs))
      lapply(lo:hi, function(n) .add(list(a0, .int(n))))
    })
    for (l in .hx_product(values)) result[[length(result) + 1L]] <- Map(list, f$symbols, l)
  }
  result
}

# itertools.product over a list of lists (last list varies fastest)
.hx_product <- function(lists) {
  out <- list(list())
  for (l in lists) {
    nxt <- list()
    for (o in out) for (v in l) nxt[[length(nxt) + 1L]] <- c(o, list(v))
    out <- nxt
  }
  out
}

.hx_collection <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      formulae <- .hx_add_formulae()
      symbolic <- list(); concrete <- list()
      for (f in formulae) {
        sz <- .hf_sizes(f$func)
        if (length(f$symbols)) symbolic[[sz]] <- c(symbolic[[sz]], list(f))
        else {
          if (is.null(concrete[[sz]])) concrete[[sz]] <- list()
          concrete[[sz]][[.hf_invariants(f$func)]] <- f
        }
      }
      cache <<- list(symbolic = symbolic, concrete = concrete)
    }
    cache
  }
})

# FormulaCollection.lookup_origin
.hx_lookup_origin <- function(func) {
  coll <- .hx_collection()
  inv <- .hf_invariants(func); sz <- .hf_sizes(func)
  if (!is.null(coll$concrete[[sz]]) && !is.null(coll$concrete[[sz]][[inv]])) return(coll$concrete[[sz]][[inv]])
  if (is.null(coll$symbolic[[sz]])) return(NULL)
  possible <- list()
  for (f in coll$symbolic[[sz]]) {
    for (repl in .hx_find_instantiations(f, func)) {
      func2 <- .hf_xreplace(f$func, repl)
      if (!.hf_suitable_origin(func2)) next
      d <- .hf_difficulty(func2, func)
      if (is_sym(d) || d != -1) possible[[length(possible) + 1L]] <- list(diff = d, repl = repl, f = f, func2 = func2)
    }
  }
  if (length(possible)) {
    dv <- vapply(possible, function(p) as.double(.neval(p$diff, list())), 0)
    possible <- possible[order(dv, method = "radix")]
  }
  for (p in possible) {
    f2 <- .hx_formula(p$func2, p$f$z, NULL, list(), .hx_subs_d(p$f$B, p$repl), .hx_subs_d(p$f$C, p$repl), .hx_subs_d(p$f$M, p$repl))
    if (!(.hx_has_bad(f2$B) || .hx_has_bad(f2$M) || .hx_has_bad(f2$C))) return(f2)
  }
  NULL
}

# MeijerFormula: an origin for Meijer G functions, with a matcher that
# instantiates it
.hx_meijer_formula <- function(an, ap, bm, bq, z, symbols, B, C, M, matcher)
  list(func = .GF(an, ap, bm, bq), z = z, symbols = symbols, B = B, C = C, M = M, matcher = matcher)

.hx_meijer_try_instantiate <- function(mf, func) {
  if (.gf_signature(func) != .gf_signature(mf$func)) return(NULL)
  res <- mf$matcher(func)
  if (is.null(res)) return(NULL)
  subs_ <- res[[1L]]; nf <- res[[2L]]
  .hx_meijer_formula(nf$an, nf$ap, nf$bm, nf$bq, mf$z, list(), .hx_subs_d(mf$B, subs_), .hx_subs_d(mf$C, subs_), .hx_subs_d(mf$M, subs_), NULL)
}

.hx_meijer_collection <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      out <- list()
      for (f in .hx_add_meijerg_formulae()) {
        s <- .gf_signature(f$func)
        out[[s]] <- c(out[[s]], list(f))
      }
      cache <<- out
    }
    cache
  }
})

.hx_meijer_lookup_origin <- function(func) {
  fs <- .hx_meijer_collection()[[.gf_signature(func)]]
  for (f in fs) {
    r <- .hx_meijer_try_instantiate(f, func)
    if (!is.null(r)) return(r)
  }
  NULL
}

# ---------------------------------------------------------------------------
# The formula tables (add_formulae and add_meijerg_formulae)

.hx_F <- function(nm, ...) fcall(nm, lapply(list(...), as_sym))
.hx_V <- function(...) lapply(list(...), as_sym)

.hx_add_formulae <- function() {
  a <- Dummy("a"); b <- Dummy("b"); c <- Dummy("c"); z <- Dummy("z")
  F <- .hx_F; V <- .hx_V; Mx <- function(...) .hx_matrix(list(...))
  h <- S_Half; R <- function(p, q) Rational(p, q); pi <- S_Pi; I <- S_I
  formulae <- list()
  add <- function(ap, bq, res)
    formulae[[length(formulae) + 1L]] <<- .hx_formula(.HF(lapply(ap, as_sym), lapply(bq, as_sym)), z, as_sym(res), list(a, b, c))
  addb <- function(ap, bq, B, C, M)
    formulae[[length(formulae) + 1L]] <<- .hx_formula(.HF(lapply(ap, as_sym), lapply(bq, as_sym)), z, NULL, list(a, b, c), lapply(B, as_sym), lapply(C, as_sym), M)
  scale <- function(s, v) lapply(v, function(e) .mul(list(as_sym(s), as_sym(e))))
  rt4 <- z^R(1, 4)

  add(list(), list(), exp(z))
  add(list(a), list(), F("HyperRep_power1", -a, z))
  addb(list(a, a - h), list(2 * a), V(F("HyperRep_power2", a, z), F("HyperRep_power2", a + h, z) / 2), V(1, 0),
       Mx(list((a - h) * z / (1 - z), (h - a) * z / (1 - z)), list(a / (1 - z), a * (z - 2) / (1 - z))))
  addb(list(1, 1), list(2), V(F("HyperRep_log1", z), 1), V(-1 / z, 0), Mx(list(0, z / (z - 1)), list(0, 0)))
  addb(list(h, 1), list(R(3, 2)), V(F("HyperRep_atanh", z), 1), V(1, 0), Mx(list(R(-1, 2), 1 / (1 - z) / 2), list(0, 0)))
  addb(list(h, h), list(R(3, 2)), V(F("HyperRep_asin1", z), F("HyperRep_power1", R(-1, 2), z)), V(1, 0),
       Mx(list(R(-1, 2), h), list(0, z / (1 - z) / 2)))
  addb(list(a, h + a), list(h), V(F("HyperRep_sqrts1", -a, z), -F("HyperRep_sqrts2", -a - h, z)), V(1, 0),
       Mx(list(0, -a), list(z * (-2 * a - 1) / 2 / (1 - z), h - z * (-2 * a - 1) / (1 - z))))
  addb(list(a, -a), list(h), V(F("HyperRep_cosasin", a, z), F("HyperRep_sinasin", a, z)), V(1, 0),
       Mx(list(0, -a), list(a * z / (1 - z), 1 / (1 - z) / 2)))
  addb(list(1, 1), list(3 * h), V(F("HyperRep_asin2", z), 1), V(1, 0), Mx(list((z - h) / (1 - z), 1 / (1 - z) / 2), list(0, 0)))
  addb(list(h, h), list(S_One), V(F("elliptic_k", z), F("elliptic_e", z)), V(2 / pi, 0),
       Mx(list(R(-1, 2), -1 / (2 * z - 2)), list(R(-1, 2), h)))
  addb(list(R(-1, 2), h), list(S_One), V(F("elliptic_k", z), F("elliptic_e", z)), V(0, 2 / pi),
       Mx(list(R(-1, 2), -1 / (2 * z - 2)), list(R(-1, 2), h)))
  addb(list(R(-1, 2), 1, 1), list(h, 2), V(z * F("HyperRep_atanh", z), F("HyperRep_log1", z), 1),
       V(R(-2, 3), -S_One / (3 * z), R(2, 3)),
       Mx(list(h, 0, z / (1 - z) / 2), list(0, 0, z / (z - 1)), list(0, 0, 0)))
  addb(list(R(-1, 2), 1, 1), list(2, 2), V(F("HyperRep_power1", h, z), F("HyperRep_log2", z), 1),
       V(R(4, 9) - 16 / (9 * z), 4 / (3 * z), 16 / (9 * z)),
       Mx(list(z / 2 / (z - 1), 0, 0), list(1 / (2 * (z - 1)), 0, h), list(0, 0, 0)))
  addb(list(1), list(b), V(z^(1 - b) * exp(z) * F("lowergamma", b - 1, z), 1), V(b - 1, 0), Mx(list(1 - b + z, 1), list(0, 0)))
  addb(list(a), list(2 * a),
       V(z^(h - a) * exp(z / 2) * F("besseli", a - h, z / 2) * gamma(a + h) / 4^(h - a),
         z^(h - a) * exp(z / 2) * F("besseli", a + h, z / 2) * gamma(a + h) / 4^(h - a)),
       V(1, 0), Mx(list(z / 2, z / 2), list(z / 2, z / 2 - 2 * a)))
  mz <- F("polar_lift", -1) * z
  addb(list(a), list(a + 1), V(mz^(-a) * a * F("lowergamma", a, mz), a * exp(z)), V(1, 0), Mx(list(-a, 1), list(0, z)))
  add(list(R(-1, 2)), list(h), exp(z) - sqrt(pi * z) * (-I) * F("erf", I * sqrt(z)))
  fre <- 2 * rt4 * exp(I * pi / 4) / sqrt(pi)
  addb(list(1), list(R(3, 4), R(5, 4)),
       V(sqrt(pi) * (I * sinh(2 * sqrt(z)) * F("fresnels", fre) + cosh(2 * sqrt(z)) * F("fresnelc", fre)) * exp(-I * pi / 4) / (2 * rt4),
         sqrt(pi) * rt4 * (sinh(2 * sqrt(z)) * F("fresnelc", fre) + I * cosh(2 * sqrt(z)) * F("fresnels", fre)) * exp(-I * pi / 4) / 2,
         1),
       V(1, 0, 0), Mx(list(R(-1, 4), 1, R(1, 4)), list(z, R(1, 4), 0), list(0, 0, 0)))
  addb(list(h, a), list(R(3, 2), a + 1),
       V(a / (2 * a - 1) * (-I) * sqrt(pi / z) * F("erf", I * sqrt(z)),
         a / (2 * a - 1) * (F("polar_lift", -1) * z)^(-a) * F("lowergamma", a, F("polar_lift", -1) * z),
         a / (2 * a - 1) * exp(z)),
       V(1, -1, 0), Mx(list(R(-1, 2), 0, 1), list(0, -a, 1), list(0, 0, z)))
  addb(list(1, 1), list(2, 2), V(F("Ei", z) - log(z), exp(z), 1, S_EulerGamma), V(1 / z, 0, 0, -1 / z),
       Mx(list(0, 1, -1, 0), list(0, z, 0, 0), list(0, 0, 0, 0), list(0, 0, 0, 0)))
  add(list(), list(h), cosh(2 * sqrt(z)))
  addb(list(), list(b),
       V(gamma(b) * z^((1 - b) / 2) * F("besseli", b - 1, 2 * sqrt(z)), gamma(b) * z^(1 - b / 2) * F("besseli", b, 2 * sqrt(z))),
       V(1, 0), Mx(list(0, 1), list(z, 1 - b)))
  x <- 4 * z^R(1, 4)
  fp <- function(a) F("besseli", a, x) + F("besselj", a, x)
  fm <- function(a) F("besseli", a, x) - F("besselj", a, x)
  B <- V(fp(2 * a - 1), fm(2 * a) * z^R(1, 4), fm(2 * a - 1) * sqrt(z), fp(2 * a) * z^R(3, 4))
  B <- scale(z^((1 - 2 * a) / 4), scale(gamma(2 * a), scale(2^(-2 * a), B)))
  addb(list(), list(h, a, a + h), B, V(1, 0, 0, 0),
       Mx(list(0, 1, 0, 0), list(0, h - a, 1, 0), list(0, 0, h, 1), list(z, 0, 0, 1 - a)))
  x <- 2 * (4 * z)^R(1, 4) * F("exp_polar", I * pi / 4)
  bj <- function(n) F("besselj", n, x); bi <- function(n) F("besseli", n, x)
  B <- V(bj(2 * a - 1) * bi(2 * a - 1), x * (bi(2 * a) * bj(2 * a - 1) - bi(2 * a - 1) * bj(2 * a)),
         x^2 * bi(2 * a) * bj(2 * a), x^3 * (bi(2 * a) * bj(2 * a - 1) + bi(2 * a - 1) * bj(2 * a)))
  B <- scale((2 * sqrt(F("polar_lift", -1) * z))^(1 - 2 * a) * gamma(2 * a)^2, B)
  addb(list(), list(a, a + h, 2 * a), B, V(1, 0, 0, 0),
       Mx(list(0, R(1, 4), 0, 0), list(0, (1 - 2 * a) / 2, R(-1, 2), 0), list(0, 0, 1 - 2 * a, R(1, 4)), list(-32 * z, 0, 0, 1 - a)))
  sz <- sqrt(z); bis <- function(n) F("besseli", n, sz)
  addb(list(a), list(a - h, 2 * a),
       V(z^(h - a) * bis(a - h)^2, z^(1 - a) * bis(a - h) * bis(a - R(3, 2)), z^(R(3, 2) - a) * bis(a - R(3, 2))^2),
       V(-gamma(a + h)^2 / 4^(h - a), 2 * gamma(a - h) * gamma(a + h) / 4^(1 - a), 0),
       Mx(list(1 - 2 * a, 1, 0), list(z / 2, h - a, h), list(0, z, 0)))
  B <- V(bis(1 - b) * bis(b - 1), sz * (bis(-b) * bis(b - 1) + bis(1 - b) * bis(b)), bis(-b) * bis(b))
  addb(list(h), list(b, 2 - b), scale(pi * (1 - b) / sin(pi * b), B), V(1, 0, 0),
       Mx(list(b - 1, h, 0), list(z, 0, z), list(0, h, -b)))
  addb(list(h), list(R(3, 2), R(3, 2)),
       V(F("Shi", 2 * sz) / 2 / sz, sinh(2 * sz) / 2 / sz, cosh(2 * sz)), V(1, 0, 0),
       Mx(list(R(-1, 2), h, 0), list(0, R(-1, 2), h), list(0, 2 * z, 0)))
  fr2 <- exp(pi * I / 4) * rt4 * 2 / sqrt(pi)
  addb(list(R(3, 4)), list(R(3, 2), R(7, 4)),
       V(F("fresnels", fr2) / (pi * fr2^3), sinh(2 * sz) / sz, cosh(2 * sz)), V(6, 0, 0),
       Mx(list(R(-3, 4), R(1, 16), 0), list(0, R(-1, 2), 1), list(0, z, 0)))
  addb(list(R(1, 4)), list(h, R(5, 4)),
       V(sqrt(pi) * exp(-I * pi / 4) * F("fresnelc", 2 * rt4 * exp(I * pi / 4) / sqrt(pi)) / (2 * rt4), cosh(2 * sz), sinh(2 * sz) * sz),
       V(1, 0, 0), Mx(list(R(-1, 4), R(1, 4), 0), list(0, 0, 1), list(0, z, h)))
  B <- V(bis(b - 1) * bis(2 * a - b), sz * bis(b) * bis(2 * a - b), sz * bis(b - 1) * bis(2 * a - b + 1), bis(b) * bis(2 * a - b + 1))
  addb(list(a, a + h), list(2 * a, b, 2 * a - b + 1), scale(gamma(b) * gamma(2 * a - b + 1) * (sz / 2)^(1 - 2 * a), B), V(1, 0, 0, 0),
       Mx(list(0, h, h, 0), list(z / 2, 1 - b, 0, z / 2), list(z / 2, 0, b - 2 * a, z / 2), list(0, h, h, -2 * a)))
  addb(list(1, 1), list(2, 2, R(3, 2)),
       V(F("Chi", 2 * sz) - log(2 * sz), cosh(2 * sz), sz * sinh(2 * sz), 1, S_EulerGamma), V(1 / z, 0, 0, 0, -1 / z),
       Mx(list(0, h, 0, R(-1, 2), 0), list(0, 0, 1, 0, 0), list(0, z, h, 0, 0), list(0, 0, 0, 0, 0), list(0, 0, 0, 0, 0)))
  addb(list(1, 1, a), list(2, 2, a + 1),
       V(a * (log(-z) + F("expint", 1, -z) + S_EulerGamma) / (z * (a^2 - 2 * a + 1)),
         a * (-z)^(-a) * (gamma(a) - F("uppergamma", a, -z)) / (a - 1)^2,
         a * exp(z) / (a^2 - 2 * a + 1), a / (z * (a^2 - 2 * a + 1))),
       V(1 - a, 1, -1 / z, 1),
       Mx(list(-1, 0, -1 / z, 1), list(0, -a, 1, 0), list(0, 0, z, 0), list(0, 0, 0, -1)))
  formulae
}

.hx_add_meijerg_formulae <- function() {
  a <- Dummy("a"); b <- Dummy("b"); c <- Dummy("c"); z <- Dummy("z"); rho <- Dummy("rho")
  V <- .hx_V; Mx <- function(...) .hx_matrix(list(...)); h <- S_Half; F <- .hx_F; pi <- S_Pi
  formulae <- list()
  add <- function(an, ap, bm, bq, B, C, M, matcher)
    formulae[[length(formulae) + 1L]] <<- .hx_meijer_formula(lapply(an, as_sym), lapply(ap, as_sym), lapply(bm, as_sym), lapply(bq, as_sym),
                                                             z, list(a, b, c, rho), B, C, M, matcher)
  m1 <- function(d) .hx_mod1(simplify(d))
  detect_uppergamma <- function(func) {
    x <- func$an[[1L]]; y <- func$bm[[1L]]; zz <- func$bm[[2L]]
    swapped <- FALSE
    if (is_zero_num(m1(.sub(x, y)))) { swapped <- TRUE; tmp <- y; y <- zz; zz <- tmp }
    if (!is_zero_num(m1(.sub(x, zz))) || .is_true(Gt(.sub(x, zz), 0))) return(NULL)
    l <- if (swapped) list(x, y) else list(y, x)
    list(list(list(rho, y), list(a, .sub(x, y))), .GF(list(x), list(), l, list()))
  }
  add(list(a + rho), list(), list(rho, a + rho), list(),
      V(gamma(1 - a) * z^rho * exp(z) * F("uppergamma", a, z), gamma(1 - a) * z^(a + rho)), V(1, 0),
      Mx(list(rho + z, -1), list(0, a + rho)), detect_uppergamma)
  detect_3113 <- function(func) {
    x <- func$an[[1L]]; u <- func$bm[[1L]]; v <- func$bm[[2L]]; w <- func$bm[[3L]]
    if (is_zero_num(m1(.sub(u, v)))) {
      if (is_zero_num(m1(.sub(v, w)))) return(NULL)
      sig <- list(h, h, S_Zero); x1 <- u; x2 <- v; y <- w
    } else if (is_zero_num(m1(.sub(x, u)))) {
      sig <- list(h, S_Zero, h); x1 <- u; y <- v; x2 <- w
    } else {
      sig <- list(S_Zero, h, h); y <- u; x1 <- v; x2 <- w
    }
    if (!is_zero_num(m1(.sub(x, x1))) || !is_zero_num(m1(.sub(x, x2))) || !identical(m1(.sub(x, y))$key, h$key) ||
        .is_true(Gt(.sub(x, x1), 0)) || .is_true(Gt(.sub(x, x2), 0))) return(NULL)
    list(list(list(a, x)), .GF(list(x), list(), lapply(sig, function(t) .add(list(x, -h, t))), list()))
  }
  s <- sin(2 * sqrt(z)); c_ <- cos(2 * sqrt(z)); S_ <- F("Si", 2 * sqrt(z)) - pi / 2; C <- F("Ci", 2 * sqrt(z))
  add(list(a), list(), list(a, a, a - h), list(),
      V(sqrt(pi) * z^(a - h) * (c_ * S_ - s * C), sqrt(pi) * z^a * (s * S_ + c_ * C), sqrt(pi) * z^a), V(-2, 0, 0),
      Mx(list(a - h, -1, 0), list(z, a, h), list(0, 0, a)), detect_3113)
  formulae
}

# ---------------------------------------------------------------------------
# Plans: the shift operators that turn an origin into the target

.hx_num_gt0 <- function(e) .is_true(Gt(e, 0))
.hx_num_le0 <- function(e) .is_true(Le(e, 0))
.hx_eq <- function(a, b) identical(a$key, b$key)

# devise_plan(target, origin, z)
.hx_devise_plan <- function(target, origin, z) {
  abuckets <- .hx_sift(target$ap); bbuckets <- .hx_sift(target$bq)
  nabuckets <- .hx_sift(origin$ap); nbbuckets <- .hx_sift(origin$bq)
  if (length(abuckets) != length(nabuckets) || length(bbuckets) != length(nbbuckets))
    stop(.hf_str(target), " not reachable from ", .hf_str(origin), call. = FALSE)
  ops <- list()
  # do_shifts mutates fro in place: return the operators and the final list
  do_shifts <- function(fro, to, inc, dec) {
    out <- list()
    for (i in seq_along(fro)) {
      if (.hx_num_gt0(.sub(to[[i]], fro[[i]]))) { sh <- inc; ch <- S_One } else { sh <- dec; ch <- S_NegativeOne }
      while (!.hx_eq(to[[i]], fro[[i]])) {
        out[[length(out) + 1L]] <- sh(fro, i)
        fro[[i]] <- .add(list(fro[[i]], ch))
      }
    }
    list(out, fro)
  }
  do_shifts_a <- function(nal, nbk, al, aother, bother)
    do_shifts(nal, al, function(p, i) .hx_shift_a(p[[i]]), function(p, i) .hx_unshift_a(c(p, aother), c(nbk, bother), i, z))
  do_shifts_b <- function(nal, nbk, bk, aother, bother)
    do_shifts(nbk, bk, function(p, i) .hx_unshift_b(c(nal, aother), c(p, bother), i, z), function(p, i) .hx_shift_b(p[[i]]))
  mods <- c(lapply(abuckets, `[[`, "mod"), lapply(bbuckets, `[[`, "mod"))
  mods <- .key_sort(lapply(mods, .sort_key), mods)
  others <- function(dic, key) unlist(lapply(names(dic)[names(dic) != key], function(k) dic[[k]]$vals), recursive = FALSE)
  for (r in mods) {
    k <- paste0("k", r$key)
    al <- list(); nal <- list(); bk <- list(); nbk <- list()
    if (!is.null(abuckets[[k]])) { al <- abuckets[[k]]$vals; nal <- nabuckets[[k]]$vals }
    if (!is.null(bbuckets[[k]])) { bk <- bbuckets[[k]]$vals; nbk <- nbbuckets[[k]]$vals }
    if (length(al) != length(nal) || length(bk) != length(nbk))
      stop(.hf_str(target), " not reachable from ", .hf_str(origin), call. = FALSE)
    al <- .hx_sort_key(al); nal <- .hx_sort_key(nal); bk <- .hx_sort_key(bk); nbk <- .hx_sort_key(nbk)
    aother <- others(nabuckets, k); bother <- others(nbbuckets, k)
    if (is.null(aother)) aother <- list()
    if (is.null(bother)) bother <- list()
    if (!length(al)) {
      ops <- c(ops, do_shifts_b(list(), nbk, bk, aother, bother)[[1L]])
    } else if (!length(bk)) {
      ops <- c(ops, do_shifts_a(nal, list(), al, aother, bother)[[1L]])
    } else {
      namax <- nal[[length(nal)]]; amax <- al[[length(al)]]
      if (.hx_num_le0(.sub(nbk[[1L]], namax)) || .hx_num_le0(.sub(bk[[1L]], amax)))
        stop("Non-suitable parameters.", call. = FALSE)
      if (.hx_num_gt0(.sub(namax, amax))) {
        r1 <- do_shifts_a(nal, nbk, al, aother, bother); ops <- c(ops, r1[[1L]]); nal <- r1[[2L]]
        r2 <- do_shifts_b(al, nbk, bk, aother, bother); ops <- c(ops, r2[[1L]]); nbk <- r2[[2L]]
      } else {
        r1 <- do_shifts_b(nal, nbk, bk, aother, bother); ops <- c(ops, r1[[1L]]); nbk <- r1[[2L]]
        r2 <- do_shifts_a(nal, bk, al, aother, bother); ops <- c(ops, r2[[1L]]); nal <- r2[[2L]]
      }
    }
    if (!is.null(nabuckets[[k]])) nabuckets[[k]]$vals <- al else if (length(al)) nabuckets[[k]] <- list(mod = r, vals = al)
    if (!is.null(nbbuckets[[k]])) nbbuckets[[k]]$vals <- bk else if (length(bk)) nbbuckets[[k]] <- list(mod = r, vals = bk)
  }
  rev(ops)
}

.hx_rf <- function(a, k) fcall("RisingFactorial", list(a, k))

# try_shifted_sum
.hx_try_shifted_sum <- function(func, z) {
  abuckets <- .hx_sift(func$ap); bbuckets <- .hx_sift(func$bq)
  a0 <- .hx_bget(abuckets, S_Zero)
  if (length(a0) != 1L) return(NULL)
  r <- a0[[1L]]
  if (.hx_num_le0(r)) return(NULL)
  b0 <- .hx_bget(bbuckets, S_Zero)
  if (is.null(b0)) return(NULL)
  l <- .hx_sort_num(b0)
  k <- l[[1L]]
  if (.hx_num_le0(k)) return(NULL)
  nap <- func$ap; nap <- nap[-which(vapply(nap, .hx_eq, TRUE, b = r))[1L]]
  nbq <- func$bq; nbq <- nbq[-which(vapply(nbq, .hx_eq, TRUE, b = k))[1L]]
  k <- .sub(k, S_One)
  nap <- lapply(nap, function(x) .sub(x, k)); nbq <- lapply(nbq, function(x) .sub(x, k))
  ops <- list()
  for (n in seq_len(as.integer(.int_val(r)) - 1L) - 1L) ops[[length(ops) + 1L]] <- .hx_shift_a(.int(n + 1L))
  ops <- rev(ops)
  fac <- .div(fcall("factorial", list(k)), .pow(z, k))
  fac <- .mul(list(fac, .mul(lapply(nbq, .hx_rf, k = k))))
  fac <- .div(fac, .mul(lapply(nap, .hx_rf, k = k)))
  ops[[length(ops) + 1L]] <- .hx_mult_op(fac)
  p <- S_Zero
  for (n in seq_len(as.integer(.int_val(k))) - 1L) {
    nn <- .int(n)
    m <- .div(.pow(z, nn), fcall("factorial", list(nn)))
    m <- .mul(list(m, .mul(lapply(nap, .hx_rf, k = nn))))
    m <- .div(m, .mul(lapply(nbq, .hx_rf, k = nn)))
    p <- .add(list(p, m))
  }
  list(.HF(nap, nbq), ops, .neg(p))
}

# try_polynomial: oo if the function is undefined, the polynomial if it
# terminates, NULL otherwise
.hx_try_polynomial <- function(func, z) {
  abuckets <- .hx_sift(func$ap); bbuckets <- .hx_sift(func$bq)
  a0 <- .hx_sort_num(.hx_bget(abuckets, S_Zero)); b0 <- .hx_sort_num(.hx_bget(bbuckets, S_Zero))
  al0 <- Filter(.hx_num_le0, a0); bl0 <- Filter(.hx_num_le0, b0)
  if (length(bl0) && all(vapply(al0, function(a) .is_true(Lt(a, bl0[[length(bl0)]])), TRUE))) return(S_Infinity)
  if (!length(al0)) return(NULL)
  a <- al0[[length(al0)]]
  fac <- S_One; res <- S_One
  for (n in seq_len(-as.integer(.int_val(a))) - 1L) {
    nn <- .int(n)
    fac <- .mul(list(fac, z))
    fac <- .div(fac, .add(list(nn, S_One)))
    fac <- .mul(list(fac, .mul(lapply(func$ap, function(x) .add(list(x, nn))))))
    fac <- .div(fac, .mul(lapply(func$bq, function(x) .add(list(x, nn)))))
    res <- .add(list(res, fac))
  }
  res
}

# try_lerchphi: functions whose parameters are rational-function multiples
# of a sum of Lerch transcendents
.hx_try_lerchphi <- function(func) {
  abuckets <- .hx_sift(func$ap); bbuckets <- .hx_sift(func$bq)
  paired <- list()
  for (k in names(abuckets)) {
    if (!is_zero_num(abuckets[[k]]$mod) && is.null(bbuckets[[k]])) return(NULL)
    bvalue <- if (is.null(bbuckets[[k]])) list() else bbuckets[[k]]$vals
    paired[[k]] <- list(abuckets[[k]]$vals, bvalue)
    bbuckets[[k]] <- NULL
  }
  if (length(bbuckets)) return(NULL)
  k0 <- paste0("k", S_Zero$key)
  if (is.null(abuckets[[k0]])) return(NULL)
  paired[[k0]][[2L]] <- c(paired[[k0]][[2L]], list(S_One))
  t <- Dummy("t")
  numer <- S_One; denom <- S_One
  for (k in names(paired)) {
    av <- paired[[k]][[1L]]; bv <- paired[[k]][[2L]]
    if (length(av) != length(bv)) return(NULL)
    for (i in seq_along(av)) {
      a <- av[[i]]; b <- bv[[i]]
      if (isTRUE(ask_fact(.sub(a, b), "positive"))) {
        kk <- .sub(a, b)
        numer <- .mul(list(numer, .hx_rf(.add(list(b, t)), kk))); denom <- .mul(list(denom, .hx_rf(b, kk)))
      } else {
        kk <- .sub(b, a)
        numer <- .mul(list(numer, .hx_rf(a, kk))); denom <- .mul(list(denom, .hx_rf(.add(list(a, t)), kk)))
      }
    }
  }
  part <- apart(.div(numer, denom), t)
  monomials <- list(); terms <- list(); term_keys <- list()
  for (arg in .aargs(part)) {
    nd <- .as_numer_denom(arg); nu <- nd[[1L]]; de <- nd[[2L]]
    if (!.has1(de, t)) {
      cs <- .hx_coeffs(nu, t)
      nz <- which(!vapply(cs, is_zero_num, TRUE))
      if (length(nz) != 1L) stop("p should be monomial", call. = FALSE)
      monomials[[length(monomials) + 1L]] <- list(.div(cs[[nz]], de), nz - 1L)
      next
    }
    if (.has1(nu, t)) stop("Need partial fraction decomposition with linear denominators", call. = FALSE)
    ind <- .as_indep(de, list(t), as_add = FALSE)
    indep <- ind[[1L]]; dep <- ind[[2L]]
    n <- S_One
    if (dep$type == "Pow") { n <- dep$args[[2L]]; dep <- dep$args[[1L]] }
    if (.hx_eq(dep, t)) a <- S_Zero
    else if (dep$type == "Add") {
      ia <- .as_indep(dep, list(t), as_add = TRUE); a <- ia[[1L]]; tmp <- ia[[2L]]
      b <- S_One
      if (!.hx_eq(tmp, t)) b <- .as_indep(tmp, list(t), as_add = FALSE)[[1L]]
      if (!.hx_eq(dep, .add(list(.mul(list(b, t)), a)))) stop("unrecognised form ", .str(dep), call. = FALSE)
      a <- .div(a, b)
      indep <- .mul(list(indep, .pow(b, n)))
    } else stop("unrecognised form of partial fraction", call. = FALSE)
    ka <- a$key
    if (is.null(terms[[ka]])) { terms[[ka]] <- list(); term_keys[[ka]] <- a }
    terms[[ka]][[length(terms[[ka]]) + 1L]] <- list(.div(nu, indep), n)
  }
  z <- Dummy("z")
  deriv <- list(); coeffs <- list()
  monomials <- monomials[order(vapply(monomials, `[[`, 0L, 2L), method = "radix")]
  mon <- list(.inv(.sub(S_One, z)))
  if (length(monomials)) for (k in seq_len(monomials[[length(monomials)]][[2L]]))
    mon[[k + 1L]] <- .mul(list(z, .diff1(mon[[k]], z)))
  put <- function(lst, key, item, val) {
    kk <- key$key
    if (is.null(lst[[kk]])) lst[[kk]] <- list(key = key, items = list())
    lst[[kk]]$items[[length(lst[[kk]]$items) + 1L]] <- val
    lst
  }
  for (mm in monomials) coeffs <- put(coeffs, S_One, NULL, .mul(list(mm[[1L]], mon[[mm[[2L]] + 1L]])))
  lerch <- function(k, a) fcall("lerchphi", list(z, k, a))
  for (ka in names(terms)) {
    a <- term_keys[[ka]]; l <- terms[[ka]]
    for (ck in l) coeffs <- put(coeffs, lerch(ck[[2L]], a), NULL, ck[[1L]])
    l <- l[order(vapply(l, function(ck) as.double(.neval(ck[[2L]], list())), 0), method = "radix")]
    kmax <- as.integer(.int_val(l[[length(l)]][[2L]]))
    for (k in seq_len(max(0L, kmax - 1L)) + 1L)
      deriv[[lerch(.int(k), a)$key]] <- list(key = lerch(.int(k), a), items = list(list(.neg(a), lerch(.int(k), a)), list(S_One, lerch(.int(k - 1L), a))))
    deriv[[lerch(S_One, a)$key]] <- list(key = lerch(S_One, a), items = list(list(.neg(a), lerch(S_One, a)), list(.inv(.sub(S_One, z)), S_One)))
  }
  basis_keys <- c(list(S_One), lapply(deriv, `[[`, "key"))
  trans <- stats::setNames(seq_along(basis_keys), vapply(basis_keys, `[[`, "", "key"))
  B <- lapply(basis_keys, expand_func)
  C <- rep(list(S_Zero), length(B))
  for (cc in coeffs) C[[trans[[cc$key$key]]]] <- .add(cc$items)
  M <- .hx_zeros(length(B))
  for (d in deriv) for (pr in d$items) M[[trans[[d$key$key]], trans[[pr[[2L]]$key]]]] <- pr[[1L]]
  .hx_formula(func, z, NULL, list(), B, C, M)
}

# build_hypergeometric_formula: a basis of contiguous functions when no
# origin is known
.hx_build_formula <- function(func) {
  z <- Dummy("z")
  if (length(func$ap)) {
    afactors <- lapply(func$ap, function(a) .add(list(.hx_x, a)))
    bfactors <- lapply(func$bq, function(b) .add(list(.hx_x, b, S_NegativeOne)))
    expr <- .sub(.mul(c(list(.hx_x), bfactors)), .mul(c(list(z), afactors)))
    cs <- .hx_coeffs(expr)
    n <- length(cs) - 1L
    basis <- list(); M <- .hx_zeros(n)
    for (k in seq_len(n) - 1L) {
      a <- .add(list(func$ap[[1L]], .int(k)))
      basis[[k + 1L]] <- hyper(c(list(a), func$ap[-1L]), func$bq, z)
      if (k < n - 1L) { M[[k + 1L, k + 1L]] <- .neg(a); M[[k + 1L, k + 2L]] <- a }
    }
    C <- c(list(S_One), rep(list(S_Zero), n - 1L))
    derivs <- list(.hx_eye(n))
    matmul <- function(A, B) {
      R <- .hx_zeros(nrow(A), ncol(B))
      for (i in seq_len(nrow(A))) for (j in seq_len(ncol(B)))
        R[[i, j]] <- .add(lapply(seq_len(ncol(A)), function(l) .mul(list(A[[i, l]], B[[l, j]]))))
      R
    }
    for (k in seq_len(n)) derivs[[k + 1L]] <- matmul(M, derivs[[k]])
    res <- rep(list(S_Zero), n)
    for (k in seq_along(cs)) {
      row <- .hx_rowmul(C, derivs[[k]])
      for (r in seq_len(n)) res[[r]] <- .add(list(res[[r]], .mul(list(cs[[k]], row[[r]]))))
    }
    lc <- cs[[length(cs)]]
    for (k in seq_len(n)) M[[n, k]] <- .neg(.div(.div(res[[k]], derivs[[n]][[1L, n]]), lc))
    return(.hx_formula(func, z, NULL, list(), basis, C, M))
  }
  basis <- list(); bq <- func$bq
  for (i in seq_along(bq)) {
    basis[[length(basis) + 1L]] <- hyper(list(), bq, z)
    bq[[i]] <- .add(list(bq[[i]], S_One))
  }
  basis[[length(basis) + 1L]] <- hyper(list(), bq, z)
  n <- length(basis)
  C <- c(list(S_One), rep(list(S_Zero), n - 1L))
  M <- .hx_zeros(n)
  M[[1L, n]] <- .div(z, .mul(func$bq))
  for (k in seq_len(n - 1L)) { M[[k + 1L, k]] <- func$bq[[k]]; M[[k + 1L, k + 1L]] <- .neg(func$bq[[k]]) }
  .hx_formula(func, z, NULL, list(), basis, C, M)
}

# hyperexpand_special: values at special points
.hx_hyperexpand_special <- function(ap, bq, z) {
  ap <- .tuple_args(ap); bq <- .tuple_args(bq)
  z_ <- z
  z <- unpolarify(z)
  if (is_zero_num(z)) return(S_One)
  if (length(ap) == 2L && length(bq) == 1L) {
    a <- ap[[1L]]; b <- ap[[2L]]; c <- bq[[1L]]
    g <- function(v) fcall("gamma", list(v))
    if (is_one_num(z)) return(.div(.div(.mul(list(g(.sub(.sub(c, a), b)), g(c))), g(.sub(c, a))), g(.sub(c, b))))
    mone <- identical(z$key, S_NegativeOne$key)
    if (mone && is_one_num(simplify(.add(list(b, .neg(a), c))))) { tmp <- b; b <- a; a <- tmp }
    if (mone && is_one_num(simplify(.add(list(a, .neg(b), c))))) {
      if (isTRUE(ask_fact(b, "integer")) && isTRUE(ask_fact(b, "negative")))
        return(.div(.div(.mul(list(S_Two, fcall("cos", list(.div(.mul(list(S_Pi, b)), S_Two))), g(.neg(b)), g(.add(list(b, .neg(a), S_One))))),
                         g(.neg(.div(b, S_Two)))), g(.add(list(.div(b, S_Two), .neg(a), S_One)))))
      return(.div(.div(.mul(list(g(.add(list(.div(b, S_Two), S_One))), g(.add(list(b, .neg(a), S_One))))), g(.add(list(b, S_One)))),
                  g(.add(list(.div(b, S_Two), .neg(a), S_One)))))
    }
  }
  hyper(ap, bq, z_)
}

# ---------------------------------------------------------------------------
# Expansion

.hx_zdiff <- function(z) function(f) .mul(list(z, .diff1(f, z)))

# _hyperexpand(func, z, ops0, z0, premult, prem, rewrite)
.hx_hyperexpand <- function(func, z, ops0 = list(), z0 = .hx_z0, premult = S_One, prem = S_Zero, rewrite = "default") {
  if (isTRUE(ask_fact(z, "zero"))) return(S_One)
  z <- polarify(z, subs = FALSE)
  if (identical(rewrite, "default")) rewrite <- "nonrepsmall"
  carryout_plan <- function(f, ops) {
    Mz <- .hx_subs(f$M, f$z, z0)
    C <- .hx_apply_operators(.hx_subs(f$C, f$z, z0), ops, .hx_make_derivative_operator(Mz, z0))
    Mp <- Mz
    if (!is_zero_num(prem)) for (i in seq_len(nrow(Mp))) Mp[[i, i]] <- .add(list(Mp[[i, i]], prem))
    C <- .hx_apply_operators(C, ops0, .hx_make_derivative_operator(Mp, z0))
    if (is_one_num(premult)) C <- lapply(C, .hx_make_simp(z0))
    r <- .mul(list(.hx_dot(C, .hx_subs(f$B, f$z, z0)), premult))
    res <- subs(r, z0, z)
    if (!is.null(rewrite)) res <- rewrite(res, rewrite)
    res
  }
  ro <- .hx_reduce_order(func)
  func <- ro[[1L]]; ops <- ro[[2L]]
  res <- .hx_try_polynomial(func, z0)
  if (!is.null(res)) {
    p <- .hx_apply_operators(res, ops, .hx_zdiff(z0))
    p <- .hx_apply_operators(.mul(list(p, premult)), ops0, .hx_zdiff(z0))
    return(unpolarify(subs(simplify(p), z0, z)))
  }
  p <- S_Zero
  res <- .hx_try_shifted_sum(func, z0)
  if (!is.null(res)) { func <- res[[1L]]; ops <- c(ops, res[[2L]]); p <- res[[3L]] }
  p <- .hx_apply_operators(p, ops, .hx_zdiff(z0))
  p <- .hx_apply_operators(.mul(list(p, premult)), ops0, .hx_zdiff(z0))
  p <- subs(simplify(p), z0, z)
  uz <- unpolarify(z)
  if ((is_one_num(uz) || identical(uz$key, S_NegativeOne$key)) && length(func$ap) == 2L && length(func$bq) == 1L) {
    f <- .hx_build_formula(func)
    r <- .replace_fun(carryout_plan(f, ops), "hyper", .hx_hyperexpand_special)
    if (!.has_function(r, "hyper")) return(.add(list(r, p)))
  }
  formula <- .hx_lookup_origin(func)
  if (is.null(formula)) formula <- .hx_try_lerchphi(func)
  if (is.null(formula)) formula <- .hx_build_formula(func)
  ops <- c(ops, .hx_devise_plan(func, formula$func, z0))
  r <- .add(list(carryout_plan(formula, ops), p))
  .replace_fun(powdenest(r, polar = TRUE), "hyper", .hx_hyperexpand_special)
}

# devise_plan_meijer(fro, to, z)
.hx_devise_plan_meijer <- function(fro, to, z) {
  st <- new.env()
  st$an <- fro$an; st$ap <- fro$ap; st$bm <- fro$bm; st$bq <- fro$bq
  try_shift <- function(slot, t, shifter, diff, counter) {
    f <- st[[slot]]
    for (idx in seq_len(min(length(f), length(t)))) {
      a <- f[[idx]]; b <- t[[idx]]
      d <- .sub(a, b)
      if (isTRUE(ask_fact(d, "integer")) && .hx_num_gt0(.div(.sub(b, a), diff)) &&
          all(vapply(counter, function(x) !.hx_eq(a, x), TRUE))) {
        sh <- shifter(idx)
        f[[idx]] <- .add(list(a, diff)); st[[slot]] <- f
        return(sh)
      }
    }
    NULL
  }
  ops <- list()
  repeat {
    op <- try_shift("an", to$an, function(i) .hx_munshift_b(st$an, st$ap, st$bm, st$bq, i, z), S_One, c(st$bm, st$bq))
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("ap", to$ap, function(i) .hx_munshift_d(st$an, st$ap, st$bm, st$bq, i, z), S_One, c(st$bm, st$bq))
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("bm", to$bm, function(i) .hx_munshift_a(st$an, st$ap, st$bm, st$bq, i, z), S_NegativeOne, c(st$an, st$ap))
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("bq", to$bq, function(i) .hx_munshift_c(st$an, st$ap, st$bm, st$bq, i, z), S_NegativeOne, c(st$an, st$ap))
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("an", to$an, function(i) .hx_mshift_b(st$an[[i]]), S_NegativeOne, list())
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("ap", to$ap, function(i) .hx_mshift_d(st$ap[[i]]), S_NegativeOne, list())
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("bm", to$bm, function(i) .hx_mshift_a(st$bm[[i]]), S_One, list())
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    op <- try_shift("bq", to$bq, function(i) .hx_mshift_c(st$bq[[i]]), S_One, list())
    if (!is.null(op)) { ops[[length(ops) + 1L]] <- op; next }
    break
  }
  same <- function(a, b) length(a) == length(b) && all(vapply(seq_along(a), function(i) .hx_eq(a[[i]], b[[i]]), TRUE))
  if (!same(st$an, to$an) || !same(st$ap, to$ap) || !same(st$bm, to$bm) || !same(st$bq, to$bq))
    stop("Could not devise plan.", call. = FALSE)
  rev(ops)
}

.hx_wcmp <- function(a, b) {
  for (i in seq_along(a)) if (a[i] != b[i]) return(if (a[i] < b[i]) -1L else 1L)
  0L
}

# _meijergexpand(func, z0, allow_hyper, rewrite, place)
.hx_meijergexpand <- function(func, z0, allow_hyper = FALSE, rewrite = "default", place = NULL) {
  if (identical(rewrite, "default")) rewrite <- NULL
  func0 <- func
  z <- Dummy("z")
  ro <- .hx_reduce_order_meijer(func)
  func <- ro[[1L]]; ops <- ro[[2L]]
  f <- .hx_meijer_lookup_origin(func)
  if (!is.null(f)) {
    ops <- c(ops, .hx_devise_plan_meijer(f$func, func, z))
    C <- .hx_apply_operators(.hx_subs(f$C, f$z, z), ops, .hx_make_derivative_operator(.hx_subs(f$M, f$z, z), z))
    C <- lapply(C, .hx_make_simp(z))
    r <- subs(.hx_dot(C, .hx_subs(f$B, f$z, z)), z, z0)
    return(powdenest(r, polar = TRUE))
  }
  t <- Dummy("t")
  can_do <- function(pbm, pap) {
    for (k in names(pbm)) if (length(pbm[[k]]$vals) > 1L) {
      l <- if (is.null(pap[[k]])) 0L else length(pap[[k]]$vals)
      if (l + 1L < length(pbm[[k]]$vals)) return(FALSE)
    }
    TRUE
  }
  g <- function(v) fcall("gamma", list(v))
  do_slater <- function(an, bm, ap, bq, z, zfinal, ops) {
    bk <- .gf_buckets(.GF(an, bm, ap, bq))
    pbm <- bk[[2L]]; pap <- bk[[3L]]
    if (!can_do(pbm, pap)) return(list(S_Zero, FALSE))
    cond <- length(an) + length(ap) < length(bm) + length(bq)
    if (length(an) + length(ap) == length(bm) + length(bq)) cond <- Lt(fcall("Abs", list(z)), S_One)
    if (identical(cond, FALSE)) return(list(S_Zero, FALSE))
    res <- S_Zero
    remove1 <- function(l, x) l[-which(vapply(l, .hx_eq, TRUE, b = x))[1L]]
    for (k in names(pbm)) {
      vals <- pbm[[k]]$vals
      if (length(vals) == 1L) {
        bh <- vals[[1L]]
        fac <- S_One
        bo <- remove1(bm, bh)
        for (bj in bo) fac <- .mul(list(fac, g(.sub(bj, bh))))
        for (aj in an) fac <- .mul(list(fac, g(.add(list(S_One, bh, .neg(aj))))))
        for (bj in bq) fac <- .div(fac, g(.add(list(S_One, bh, .neg(bj)))))
        for (aj in ap) fac <- .div(fac, g(.sub(aj, bh)))
        nap <- lapply(c(an, ap), function(a) .add(list(S_One, bh, .neg(a))))
        nbq <- lapply(c(bo, bq), function(b) .add(list(S_One, bh, .neg(b))))
        kk <- fcall("polar_lift", list(.pow(S_NegativeOne, .int(length(ap) - length(bm)))))
        harg <- .mul(list(kk, zfinal))
        premult <- .pow(.div(t, kk), bh)
        hyp <- .hx_hyperexpand(.HF(nap, nbq), harg, ops, t, premult, bh, rewrite = NULL)
        res <- .add(list(res, .mul(list(fac, hyp))))
      } else {
        b_ <- vals[[1L]]
        ki <- lapply(vals[-1L], function(bi) .sub(bi, b_))
        u <- length(ki)
        pv <- pap[[k]]$vals
        li <- lapply(pv[seq_len(min(u + 1L, length(pv)))], function(ai) .sub(ai, b_))
        bo <- bm; for (b in vals) bo <- remove1(bo, b)
        ao <- ap; for (a in pv[seq_len(min(u, length(pv)))]) ao <- remove1(ao, a)
        lu <- li[[length(li)]]
        di <- lapply(seq_len(u), function(i) .sub(li[[i]], ki[[i]]))
        s <- Dummy("s")
        integrand <- .pow(z, s)
        for (b in bm) {
          if (is_number_type(b) && is_zero_num(fcall("Mod", list(b, S_One)))) b <- .int(round(as.double(.neval(b, list()))))
          integrand <- .mul(list(integrand, g(.sub(b, s))))
        }
        for (a in an) integrand <- .mul(list(integrand, g(.add(list(S_One, .neg(a), s)))))
        for (b in bq) integrand <- .div(integrand, g(.add(list(S_One, .neg(b), s))))
        for (a in ap) integrand <- .div(integrand, g(.sub(a, s)))
        integrand <- expand_func(integrand)
        for (r in seq_len(round(as.double(.neval(lu, list())))) - 1L) {
          resid <- residue(integrand, s, .add(list(b_, .int(r))))
          resid <- .hx_apply_operators(resid, ops, .hx_zdiff(z))
          res <- .sub(res, resid)
        }
        au <- .add(list(b_, lu))
        kk <- fcall("polar_lift", list(.pow(S_NegativeOne, .int(length(ao) + length(bo) + 1L))))
        harg <- .mul(list(kk, zfinal))
        premult <- .pow(.div(t, kk), au)
        nap <- c(lapply(c(an, ap), function(a) .add(list(S_One, au, .neg(a)))), list(S_One))
        nbq <- lapply(c(bm, bq), function(b) .add(list(S_One, au, .neg(b))))
        hyp <- .hx_hyperexpand(.HF(nap, nbq), harg, ops, t, premult, au, rewrite = NULL)
        C <- .div(.pow(S_NegativeOne, lu), fcall("factorial", list(lu)))
        for (i in seq_len(u)) C <- .mul(list(C, .div(.pow(S_NegativeOne, di[[i]]), .hx_rf(.add(list(lu, .neg(li[[i]]), S_One)), di[[i]]))))
        for (a in an) C <- .mul(list(C, g(.add(list(S_One, .neg(a), au)))))
        for (b in bo) C <- .mul(list(C, g(.sub(b, au))))
        for (a in ao) C <- .div(C, g(.sub(a, au)))
        for (b in bq) C <- .div(C, g(.add(list(S_One, .neg(b), au))))
        res <- .add(list(res, .mul(list(C, hyp))))
      }
    }
    list(res, cond)
  }
  r1 <- do_slater(func$an, func$bm, func$ap, func$bq, z, z0, ops)
  slater1 <- r1[[1L]]; cond1 <- r1[[2L]]
  tr <- function(l) lapply(l, function(x) .sub(S_One, x))
  ops2 <- lapply(ops, function(o) .hx_op_from(subs(.hx_poly_expr(o$coeffs), list(list(z, .inv(t)), list(.hx_x, .neg(.hx_x))), simultaneous = TRUE), o$info))
  r2 <- do_slater(tr(func$bm), tr(func$an), tr(func$bq), tr(func$ap), t, .inv(z0), ops2)
  slater2 <- r2[[1L]]; cond2 <- r2[[2L]]
  slater1 <- powdenest(subs(slater1, z, z0), polar = TRUE)
  slater2 <- powdenest(subs(slater2, t, .inv(z0)), polar = TRUE)
  if (is_sym(cond2)) cond2 <- subs(cond2, t, .inv(z))
  m <- .gf_call(func, z)
  delta <- .mg_delta(m)
  if (.hx_num_gt0(delta) || (is_zero_num(delta) && length(.mg_ap(m)) == length(.mg_bq(m)) &&
                              !.is_false(Lt(fcall("re", list(.mg_nu(m))), -1)) &&
                              .hx_eq(fcall("polar_lift", list(z0)), fcall("polar_lift", list(S_One))))) {
    if (!identical(cond1, FALSE)) cond1 <- TRUE
    if (!identical(cond2, FALSE)) cond2 <- TRUE
  }
  slater1 <- rewrite(slater1, rewrite %||% (if (identical(cond1, TRUE)) "nonrep" else "nonrepsmall"))
  slater2 <- rewrite(slater2, rewrite %||% (if (identical(cond2, TRUE)) "nonrep" else "nonrepsmall"))
  if (!identical(cond1, FALSE) && !identical(cond2, FALSE)) {
    if (!is.null(place) && is_zero_num(as_sym(place))) cond2 <- FALSE
    if (!is.null(place) && .hx_eq(as_sym(place), S_ComplexInfinity)) cond1 <- FALSE
  }
  if (is_sym(cond1)) cond1 <- subs(cond1, z, z0)
  if (is_sym(cond2)) cond2 <- subs(cond2, z, z0)
  weight <- function(expr, cond) {
    c0 <- if (identical(cond, TRUE)) 0 else if (identical(cond, FALSE)) 1 else 2
    if (.has_any(expr, list(S_Infinity, S_ComplexInfinity, S_NegativeInfinity, S_NaN))) c0 <- 3
    co <- count_ops(expr)
    c(c0, .count_fun(expr, "hyper"), as.double(if (is_sym(co)) .neval(co, list()) else co))
  }
  w1 <- weight(slater1, cond1); w2 <- weight(slater2, cond2)
  wmin <- if (.hx_wcmp(w1, w2) <= 0) w1 else w2
  if (.hx_wcmp(wmin, c(0, 1, Inf)) <= 0) return(if (.hx_wcmp(w1, w2) < 0) slater1 else slater2)
  pw <- function() Piecewise(list(slater1, cond1), list(slater2, cond2), list(.gf_call(func0, z0), TRUE))
  if (max(w1[1L], w2[1L]) <= 1 && max(w1[2L], w2[2L]) <= 1) return(pw())
  r <- pw()
  if (!.has_function(r, "hyper") || allow_hyper) return(r)
  .gf_call(func0, z0)
}

#' Expand hypergeometric functions
#'
#' `hyperexpand(f)` rewrites the generalized hypergeometric functions
#' (`hyper()`) and Meijer G functions (`meijerg()`) in `f` in terms of named
#' special functions, as SymPy's `hyperexpand()`. Functions it cannot
#' express are left unchanged.
#'
#' @param f An expression.
#' @param allow_hyper Allow the result to contain simpler `hyper()` terms
#'   when a Meijer G function has no closed form.
#' @param rewrite The branch rewriting of the result: `"default"`,
#'   `"nonrepsmall"`, `"nonrep"`, or `NULL` for none.
#' @param place Where the expansion should be valid for a Meijer G function
#'   with two representations: `0` or `zoo`.
#' @return An expression.
#' @examples
#' z <- Symbol("z")
#' hyperexpand(hyper(list(1, 1), list(2), -z) * z)
#' hyperexpand(hyper(list(), list(Rational(1, 2)), -z^2 / 4))
#' @export
hyperexpand <- function(f, allow_hyper = FALSE, rewrite = "default", place = NULL) {
  if (inherits(f, "symmatrix")) return(.mat_map(f, function(e) hyperexpand(e, allow_hyper, rewrite, place)))
  f <- as_sym(f)
  do_replace <- function(ap, bq, z) {
    r <- .hx_hyperexpand(.HF(.tuple_args(ap), .tuple_args(bq)), z, rewrite = rewrite)
    if (is.null(r)) hyper(ap, bq, z) else r
  }
  do_meijer <- function(ap, bq, z) {
    r <- .hx_meijergexpand(.GF(.tuple_args(ap$args[[1L]]), .tuple_args(ap$args[[2L]]), .tuple_args(bq$args[[1L]]), .tuple_args(bq$args[[2L]])),
                           z, allow_hyper, rewrite = rewrite, place = place)
    if (.has_any(r, list(S_NaN, S_ComplexInfinity, S_Infinity, S_NegativeInfinity))) NULL else r
  }
  f <- .replace_fun(f, "hyper", do_replace)
  .map_bottom_up(f, function(z) {
    if (z$type == "Function" && z$name == "meijerg") {
      r <- do.call(do_meijer, z$args)
      if (!is.null(r)) return(r)
    }
    z
  })
}
