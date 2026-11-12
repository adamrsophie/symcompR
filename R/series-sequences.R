# Sequences, after sympy/series/sequences.py: SeqFormula, SeqPer, SeqAdd,
# SeqMul and EmptySequence as sym nodes.

.SEQ_TYPES <- c("SeqFormula", "SeqPer", "SeqAdd", "SeqMul", "EmptySequence")
.is_seq <- function(e) is_sym(e) && e$type %in% .SEQ_TYPES

S_EmptySequence <- new_sym("EmptySequence", key = "EmptySequence()")

# start and stop of a sequence as syms.
.seq_bounds <- function(s) {
  switch(s$type,
    SeqFormula = , SeqPer = list(s$args[[2L]]$args[[2L]], s$args[[2L]]$args[[3L]]),
    SeqAdd = , SeqMul = {
      b <- lapply(s$args, .seq_bounds)
      list(do.call(Max, lapply(b, `[[`, 1L)), do.call(Min, lapply(b, `[[`, 2L)))
    },
    EmptySequence = list(S_Infinity, S_NegativeInfinity))
}
.seq_var <- function(s) s$args[[2L]]$args[[1L]]

.seq_limits <- function(expr, limits, what) {
  find_x <- function() {
    fs <- free_symbols(expr)
    if (length(fs) == 1L) return(fs[[1L]])
    if (what == "SeqPer" || !length(fs)) return(Dummy("k"))
    stop(" specify dummy variables for ", .str(expr), ". If the formula contains more than one free symbol, ",
         "a dummy variable should be supplied explicitly e.g., SeqFormula(m*n**2, (n, 0, 5))", call. = FALSE)
  }
  if (is.null(limits)) return(list(find_x(), S_Zero, S_Infinity))
  limits <- lapply(if (is_sym(limits) && limits$type == "Tuple") limits$args else limits, as_sym)
  lim <- if (length(limits) == 3L) limits else if (length(limits) == 2L) c(list(find_x()), limits) else NULL
  if (is.null(lim) || !is_symbol_like(lim[[1L]])) stop("Invalid limits given: ", deparse(limits), call. = FALSE)
  if (is_const(lim[[2L]], "-oo") && is_const(lim[[3L]], "oo"))
    stop("Both the start and end value cannot be unbounded", call. = FALSE)
  lim
}

.seq_empty_interval <- function(a, b) identical(Interval(a, b)$type, "EmptySet")

#' Sequences
#'
#' `SeqFormula(formula, limits)` is the sequence of values of `formula` for
#' the variable running over `limits = list(n, start, stop)` (default
#' `list(n, 0, oo)`). `SeqPer(periodical, limits)` repeats a list of
#' values. `sym_sequence(seq, limits)` picks one of them, like SymPy's
#' `sequence()` (renamed so it does not mask [base::sequence()]). `+`, `-`
#' and `*` combine sequences termwise into `SeqAdd` and `SeqMul`, reducing
#' like terms; `-s` negates. `coeff(s, pt)` is the term at `pt` and
#' `sequence_terms(s, n)` the first `n` terms (SymPy's `s[:n]`).
#' `find_linear_recurrence(s, n, d, gfvar)` finds a linear recurrence among
#' the first `n` terms and, with `gfvar`, the generating function.
#'
#' @param formula,seq An expression, or a list of values for `SeqPer()`.
#' @param periodical A list of values repeated by the sequence.
#' @param limits `list(n, start, stop)` or `list(start, stop)`.
#' @param ... Sequences.
#' @param s A sequence.
#' @param n Number of terms.
#' @param d Largest recurrence order to try.
#' @param gfvar A symbol for the generating function.
#' @return A sequence, a list of terms, or for `find_linear_recurrence()`
#'   the coefficients (with `gfvar`, a list of the coefficients and the
#'   generating function).
#' @examples
#' n <- Symbol("n")
#' s <- SeqFormula(n^2, list(n, 0, Inf))
#' sequence_terms(s, 5)
#' sequence_terms(SeqPer(list(1, 2, 3)), 7)
#' SeqFormula(n^2) + SeqFormula(n^3)
#' find_linear_recurrence(sym_sequence(fibonacci(n)), 10)
#' @name sequences
NULL

#' @rdname sequences
#' @export
SeqFormula <- function(formula, limits = NULL) {
  formula <- as_sym(formula)
  lim <- .seq_limits(formula, limits, "SeqFormula")
  if (.seq_empty_interval(lim[[2L]], lim[[3L]])) return(S_EmptySequence)
  new_sym("SeqFormula", list(formula, do.call(Tuple, lim)))
}

#' @rdname sequences
#' @export
SeqPer <- function(periodical, limits = NULL) {
  if (is_sym(periodical) && periodical$type == "Tuple") periodical <- periodical$args
  if (!is.list(periodical) || is_sym(periodical))
    stop("invalid period ", if (is_sym(periodical)) .str(periodical) else deparse(periodical),
         " should be something like e.g (1, 2) ", call. = FALSE)
  flat <- function(v) if (is.list(v) && !is_sym(v)) do.call(c, lapply(v, flat)) else list(as_sym(v))
  per <- flat(periodical)
  lim <- .seq_limits(do.call(Tuple, per), limits, "SeqPer")
  if (.seq_empty_interval(lim[[2L]], lim[[3L]])) return(S_EmptySequence)
  new_sym("SeqPer", list(do.call(Tuple, per), do.call(Tuple, lim)))
}

#' @rdname sequences
#' @export
sym_sequence <- function(seq, limits = NULL) {
  if (is.list(seq) && !is_sym(seq) || (is_sym(seq) && seq$type == "Tuple")) SeqPer(seq, limits)
  else SeqFormula(seq, limits)
}

#' @rdname sequences
#' @export
EmptySequence <- S_EmptySequence

.seq_coeff <- function(s, pt) {
  pt <- as_sym(pt)
  b <- .seq_bounds(s)
  if (s$type == "EmptySequence" || .is_true(Lt(pt, b[[1L]])) || .is_true(Gt(pt, b[[2L]])))
    stop("Index ", .str(pt), " out of bounds ", .str(Interval(b[[1L]], b[[2L]])), call. = FALSE)
  switch(s$type,
    SeqFormula = subs(s$args[[1L]], .seq_var(s), pt),
    SeqPer = {
      per <- s$args[[1L]]$args
      k <- length(per)
      off <- if (is_const(b[[1L]], "-oo")) .sub(b[[2L]], pt) else .sub(pt, b[[1L]])
      idx <- as.integer(as.character(.Z(off) %% k))
      subs(per[[idx + 1L]], .seq_var(s), pt)
    },
    SeqAdd = .add(lapply(s$args, .seq_coeff, pt = pt)),
    SeqMul = .mul(lapply(s$args, .seq_coeff, pt = pt)))
}

.seq_ith_point <- function(s, i) {
  b <- .seq_bounds(s)
  if (is_const(b[[1L]], "-oo")) .sub(b[[2L]], .int(i)) else .add(list(b[[1L]], .int(i)))
}

.seq_length <- function(s) {
  if (s$type == "EmptySequence") return(S_Zero)
  b <- .seq_bounds(s)
  .add(list(.sub(b[[2L]], b[[1L]]), S_One))
}

#' @rdname sequences
#' @export
sequence_terms <- function(s, n = 10) {
  if (s$type == "Range") {
    L <- .range_size(s)
    if (is_exact_number(L)) n <- min(n, .int_val(L))
    return(lapply(seq_len(n) - 1L, function(i) .range_item(s, i)))
  }
  L <- .seq_length(s)
  if (is_exact_number(L)) n <- min(n, .int_val(L))
  lapply(seq_len(n) - 1L, function(i) .seq_coeff(s, .seq_ith_point(s, i)))
}

.seq_coeff_mul <- function(s, c) {
  c <- as_sym(c)
  switch(s$type,
    EmptySequence = s,
    SeqFormula = SeqFormula(.mul(list(s$args[[1L]], c)), s$args[[2L]]),
    SeqPer = SeqPer(lapply(s$args[[1L]]$args, function(v) .mul(list(v, c))), s$args[[2L]]),
    SeqAdd = new_sym("SeqAdd", lapply(s$args, .seq_coeff_mul, c = c)),
    SeqMul = new_sym("SeqMul", c(list(.seq_coeff_mul(s$args[[1L]], c)), s$args[-1L])))
}

# Common interval of two sequences.
.seq_intersect <- function(a, b) {
  I <- Intersection(Interval(.seq_bounds(a)[[1L]], .seq_bounds(a)[[2L]]), Interval(.seq_bounds(b)[[1L]], .seq_bounds(b)[[2L]]))
  if (I$type == "Interval") list(I$args[[1L]], I$args[[2L]]) else if (I$type == "FiniteSet") list(I$args[[1L]], I$args[[1L]]) else NULL
}

.seq_combine_periodic <- function(s, o, op) {
  per1 <- s$args[[1L]]$args; l1 <- length(per1)
  per2 <- o$args[[1L]]$args; l2 <- length(per2)
  L <- as.integer(as.character(gmp::lcm.bigz(l1, l2)))
  iv <- .seq_intersect(s, o)
  if (is.null(iv)) return(NULL)
  start <- iv[[1L]]; stop <- iv[[2L]]
  to_stop <- is_const(start, "-oo")
  as_int <- function(e) if (is_inf_const(e) || e$type != "Integer") NULL else .int_val(e)
  shift <- function(q) {
    b <- .seq_bounds(q)
    if (to_stop) {
      if (is_const(b[[1L]], "-oo")) list(as_int(.sub(b[[2L]], stop)), 1L) else list(as_int(.sub(stop, b[[1L]])), -1L)
    } else if (is_const(b[[1L]], "-oo")) list(as_int(.sub(b[[2L]], start)), -1L)
    else list(as_int(.sub(start, b[[1L]])), 1L)
  }
  s1 <- shift(s); s2 <- shift(o)
  if (is.null(s1[[1L]]) || is.null(s2[[1L]])) return(NULL)
  v1 <- .seq_var(s); v2 <- .seq_var(o)
  if (!identical(v1$key, v2$key)) per2 <- lapply(per2, function(p) subs(p, v2, v1))
  newp <- lapply(seq_len(L) - 1L, function(x) {
    i1 <- (s1[[1L]] + s1[[2L]] * x) %% l1
    i2 <- (s2[[1L]] + s2[[2L]] * x) %% l2
    op(per1[[i1 + 1L]], per2[[i2 + 1L]])
  })
  SeqPer(newp, list(v1, start, stop))
}

.seq_pair <- function(s, t, kind) {
  if (s$type == "SeqFormula" && t$type == "SeqFormula") {
    iv <- .seq_intersect(s, t)
    if (is.null(iv)) return(NULL)
    v1 <- .seq_var(s)
    f2 <- subs(t$args[[1L]], .seq_var(t), v1)
    f <- if (kind == "add") .add(list(s$args[[1L]], f2)) else .mul(list(s$args[[1L]], f2))
    return(SeqFormula(f, list(v1, iv[[1L]], iv[[2L]])))
  }
  if (s$type == "SeqPer" && t$type == "SeqPer")
    return(.seq_combine_periodic(s, t, if (kind == "add") function(a, b) .add(list(a, b)) else function(a, b) .mul(list(a, b))))
  NULL
}

.seq_start_key <- function(s) {
  b <- .seq_bounds(s)[[1L]]
  v <- if (is_const(b, "-oo")) -Inf else num_value_dbl(b)
  list(if (is.na(v)) Inf else v, .nodes(s), .sort_key(s))
}

# SeqAdd/SeqMul with SymPy's reduce(): combine the first pair that merges.
.seq_op <- function(type, args) {
  kind <- if (type == "SeqAdd") "add" else "mul"
  flat <- list()
  for (a in args) {
    if (!.is_seq(a)) stop("Input must be Sequences or  iterables of Sequences", call. = FALSE)
    flat <- c(flat, if (a$type == type) a$args else list(a))
  }
  if (kind == "add") flat <- Filter(function(a) a$type != "EmptySequence", flat)
  if (!length(flat)) return(S_EmptySequence)
  iv <- Reduce(function(I, a) Intersection(I, Interval(.seq_bounds(a)[[1L]], .seq_bounds(a)[[2L]])), flat,
               Interval(S_NegativeInfinity, S_Infinity))
  if (iv$type == "EmptySet") return(S_EmptySequence)
  repeat {
    merged <- FALSE
    for (i in seq_along(flat)) {
      for (j in seq_along(flat)) {
        if (i == j) next
        ns <- .seq_pair(flat[[i]], flat[[j]], kind)
        if (!is.null(ns)) {
          flat <- c(flat[-c(i, j)], list(ns))
          merged <- TRUE
          break
        }
      }
      if (merged) break
    }
    if (!merged) break
  }
  if (length(flat) == 1L) return(flat[[1L]])
  new_sym(type, .key_sort(lapply(flat, .seq_start_key), flat))
}

#' @rdname sequences
#' @export
SeqAdd <- function(...) .seq_op("SeqAdd", list(...))

#' @rdname sequences
#' @export
SeqMul <- function(...) .seq_op("SeqMul", list(...))

.seq_ops <- function(gen, e1, e2) {
  if (!(.is_seq(e1) && .is_seq(e2)))
    stop("cannot ", switch(gen, "+" = "add", "-" = "subtract", "*" = "multiply", "combine"), " sequence and ",
         if (.is_seq(e1)) class(e2)[1L] else class(e1)[1L], call. = FALSE)
  switch(gen,
    "+" = SeqAdd(e1, e2),
    "-" = SeqAdd(e1, .seq_coeff_mul(e2, S_NegativeOne)),
    "*" = SeqMul(e1, e2),
    "==" = identical(e1$key, e2$key),
    "!=" = !identical(e1$key, e2$key),
    stop("unsupported sequence operator ", gen, call. = FALSE))
}

#' @rdname sequences
#' @export
find_linear_recurrence <- function(s, n, d = NULL, gfvar = NULL) {
  x <- lapply(sequence_terms(s, n), function(t) simplify(expand(t)))
  lx <- length(x)
  r <- if (is.null(d)) lx %/% 2L else min(d, lx %/% 2L)
  coeffs <- list()
  for (l in seq_len(r)) {
    l2 <- 2L * l
    m <- Matrix(lapply(seq_len(l) - 1L, function(k) x[(k + 1L):(k + l)]))
    if (is_zero_num(simplify(sym_det(m)))) next
    y <- lapply(solve(m, Matrix(lapply(x[(l + 1L):l2], list)))$data, simplify)
    if (lx == l2) {
      coeffs <- rev(y)
      break
    }
    rows <- lapply(l:(lx - l - 1L), function(k) x[(k + 1L):(k + l)])
    ok <- all(vapply(seq_along(rows), function(i)
      is_zero_num(simplify(.sub(.add(Map(function(a, b) .mul(list(a, b)), rows[[i]], y)), x[[l2 + i]]))), TRUE))
    if (ok) {
      coeffs <- rev(y)
      break
    }
  }
  if (is.null(gfvar)) return(coeffs)
  gfvar <- as_sym(gfvar)
  l <- length(coeffs)
  if (!l) return(list(list(), NULL))
  num <- .mul(list(x[[l]], .pow(gfvar, .int(l - 1L))))
  den <- .sub(S_One, .mul(list(coeffs[[l]], .pow(gfvar, .int(l)))))
  if (l > 1L) for (i in seq_len(l - 1L)) {
    num <- .add(list(num, .mul(list(x[[i]], .pow(gfvar, .int(i - 1L))))))
    for (j in seq_len(l - i)) num <- .sub(num, .mul(list(coeffs[[i]], x[[j]], .pow(gfvar, .int(i + j - 1L)))))
    den <- .sub(den, .mul(list(coeffs[[i]], .pow(gfvar, .int(i)))))
  }
  list(coeffs, simplify(.div(factor_sym(num), factor_sym(den))))
}

for (.t in c("SeqFormula", "SeqPer")) {
  register_bound_vars(.t, function(e) e$args[[2L]]$args[[1L]]$key)
  register_str_printer(.t, local({
    tt <- .t
    function(e) paste0(tt, "(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")")
  }))
}
rm(.t)
for (.t in c("SeqAdd", "SeqMul")) register_str_printer(.t, local({
  tt <- .t
  function(e) paste0(tt, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
}))
rm(.t)
register_str_printer("EmptySequence", function(e) "EmptySequence")
register_latex_printer("EmptySequence", function(e) "\\emptyset")

.latex_seq <- function(e) {
  b <- .seq_bounds(e)
  if (length(free_symbols(b[[1L]])) || length(free_symbols(b[[2L]])))
    return(paste0("\\left\\{", .latex(e$args[[1L]]), "\\right\\}_{", .latex(.seq_var(e)), "=", .latex(b[[1L]]),
                  "}^{", .latex(b[[2L]]), "}"))
  L <- .seq_length(e)
  items <- if (is_const(b[[1L]], "-oo")) {
    c("\\ldots", vapply(3:0, function(k) .latex(.seq_coeff(e, .sub(b[[2L]], .int(k)))), ""))
  } else if (is_const(b[[2L]], "oo") || (is_exact_number(L) && .int_val(L) > 4)) {
    c(vapply(sequence_terms(e, 4L), .latex, ""), "\\ldots")
  } else vapply(sequence_terms(e, .int_val(L)), .latex, "")
  paste0("\\left[", paste(items, collapse = ", "), "\\right]")
}
for (.t in c("SeqFormula", "SeqPer", "SeqAdd", "SeqMul")) register_latex_printer(.t, .latex_seq)
rm(.t)

register_rebuilder("SeqFormula", function(e, args) SeqFormula(args[[1L]], args[[2L]]))
register_rebuilder("SeqPer", function(e, args) SeqPer(args[[1L]]$args, args[[2L]]))
register_rebuilder("SeqAdd", function(e, args) do.call(SeqAdd, args))
register_rebuilder("SeqMul", function(e, args) do.call(SeqMul, args))

#' Pade approximants of a series
#'
#' `approximants(l, X, simplify)` returns the successive Pade approximants
#' of the series whose coefficients are `l`, computed with the continued
#' fraction method of SymPy's `approximants()` generator.
#'
#' @param l A list of coefficients.
#' @param X The variable of the approximants.
#' @param simplify Simplify each approximant.
#' @return A list of rational functions.
#' @examples
#' approximants(list(1, 1, 2, 5, 14, 42, 132))
#' @export
approximants <- function(l, X = Symbol("x"), simplify = FALSE) {
  l <- lapply(l, as_sym)
  X <- as_sym(X)
  p1 <- list(S_One); q1 <- list(S_Zero)
  p2 <- list(S_Zero); q2 <- list(S_One)
  out <- list()
  while (length(l)) {
    b <- 0L
    while (is_zero_num(l[[b + 1L]])) {
      b <- b + 1L
      if (b == length(l)) return(out)
    }
    m <- list(.div(S_One, l[[b + 1L]]))
    if (b + 1L < length(l)) for (k in (b + 1L):(length(l) - 1L)) {
      s <- S_Zero
      # Python's m[b - j - 1] counts from the end
      for (j in b:(k - 1L)) s <- .sub(s, .mul(list(l[[j + 2L]], m[[length(m) + b - j]])))
      m[[length(m) + 1L]] <- .div(s, l[[b + 1L]])
    }
    l <- m
    a <- l[[1L]]
    l[[1L]] <- S_Zero
    p <- rep(list(S_Zero), max(length(p2), b + length(p1)))
    q <- rep(list(S_Zero), max(length(q2), b + length(q1)))
    for (k in seq_along(p2)) p[[k]] <- .mul(list(a, p2[[k]]))
    for (k in seq_along(p1)) p[[b + k]] <- .add(list(p[[b + k]], p1[[k]]))
    for (k in seq_along(q2)) q[[k]] <- .mul(list(a, q2[[k]]))
    for (k in seq_along(q1)) q[[b + k]] <- .add(list(q[[b + k]], q1[[k]]))
    while (length(p) && is_zero_num(p[[length(p)]])) p[[length(p)]] <- NULL
    while (length(q) && is_zero_num(q[[length(q)]])) q[[length(q)]] <- NULL
    p1 <- p2; p2 <- p
    q1 <- q2; q2 <- q
    c <- S_One
    for (v in c(p, q)) c <- lcm(c, fraction(v)[[2L]])
    num <- .add(lapply(seq_along(p), function(k) .mul(list(c, p[[k]], .pow(X, .int(k - 1L))))))
    den <- .add(lapply(seq_along(q), function(k) .mul(list(c, q[[k]], .pow(X, .int(k - 1L))))))
    r <- .div(num, den)
    out[[length(out) + 1L]] <- if (simplify) simplify(r) else r
  }
  out
}
