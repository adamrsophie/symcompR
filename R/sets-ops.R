# Set operations, after sympy/sets/sets.py (Interval, Union, Intersection,
# Complement, SymmetricDifference, simplify_union, simplify_intersection,
# _handle_finite_sets, _apply_operation, imageset) and the multiple
# dispatch handlers in sympy/sets/handlers/ (union, intersection, issubset,
# add, mul, power, functions), plus SetExpr from sympy/sets/setexpr.py.

.SET_TYPES <- c("Interval", "FiniteSet", "Union", "Intersection", "Complement",
                "ProductSet", "ImageSet", "ConditionSet", "SpecialSet", "EmptySet",
                "SymmetricDifference", "Range", "DisjointUnion",
                "CartesianComplexRegion", "PolarComplexRegion")

# ---------------------------------------------------------------------------
# Fuzzy logic and Python truth values

# fuzzy_bool of a SymPy value: TRUE, FALSE or NA
.fb <- function(e) {
  if (is.null(e)) return(NA)
  if (is.logical(e)) return(if (length(e) == 1L) e else NA)
  if (.is_true(e)) TRUE else if (.is_false(e)) FALSE else NA
}
# fuzzy_and / fuzzy_or over lazily computed values, as SymPy's generators
.f_and <- function(xs, f = identity) {
  rv <- TRUE
  for (v in xs) {
    a <- .fb(f(v))
    if (identical(a, FALSE)) return(FALSE)
    if (isTRUE(rv)) rv <- a
  }
  rv
}
.f_or <- function(xs, f = identity) {
  rv <- FALSE
  for (v in xs) {
    a <- .fb(f(v))
    if (isTRUE(a)) return(TRUE)
    if (identical(rv, FALSE)) rv <- a
  }
  rv
}
.f_not <- function(v) if (is.na(v)) NA else !v
.tfn <- function(v) if (is.na(v)) NULL else if (v) S_true else S_false
# bool() of a SymPy Boolean: relationals raise TypeError
.pb <- function(e) {
  if (is.logical(e)) return(isTRUE(e))
  if (.is_true(e)) return(TRUE)
  if (.is_false(e)) return(FALSE)
  .perr("TypeError", paste0("cannot determine truth value of Relational: ", .str(e)))
}
# Python ==: structural, but numbers compare by value
.py_eq <- function(a, b) {
  if (identical(a$key, b$key)) return(TRUE)
  if (is_number_type(a) && is_number_type(b) && !is_nan(a) && !is_nan(b)) return(isTRUE(.sort_num(a) == .sort_num(b)))
  FALSE
}
.lgl <- function(v) if (is_sym(v)) .is_true(v) else isTRUE(v)
.is_expr_s <- function(e) is_sym(e) && !is_set(e) && .is_expr_node(e) && !(e$type %in% c("Contains", "Lambda"))
.is_number_e <- function(e) .is_expr_s(e) && !length(free_symbols(e))

# Set equality: S.Reals == Interval(-oo, oo)
.reals_like <- function(s) (s$type == "SpecialSet" && s$name == "Reals") ||
  (s$type == "Interval" && is_const(s$args[[1L]], "-oo") && is_const(s$args[[2L]], "oo"))
.set_eq <- function(a, b) identical(a$key, b$key) || (is_set(a) && is_set(b) && .reals_like(a) && .reals_like(b))
.uniq_sets <- function(l) {
  out <- list()
  for (s in l) if (!any(vapply(out, .set_eq, NA, b = s))) out[[length(out) + 1L]] <- s
  out
}
.is_reals <- function(s) s$type == "SpecialSet" && s$name == "Reals"
.is_ivl <- function(s) s$type == "Interval" || .is_reals(s)

# Interval parts; Reals is the Interval(-oo, oo) singleton
.iv <- function(I) {
  if (I$type == "SpecialSet") return(list(a = S_NegativeInfinity, b = S_Infinity, lo = TRUE, ro = TRUE))
  list(a = I$args[[1L]], b = I$args[[2L]], lo = .is_true(I$args[[3L]]), ro = .is_true(I$args[[4L]]))
}

# Python `e in s`: a bool or TypeError
.py_in <- function(s, e) {
  c <- .contains(s, e)
  if (.is_true(c)) return(TRUE)
  if (.is_false(c)) return(FALSE)
  .perr("TypeError", paste0("did not evaluate to a bool: ", .str(c)))
}
# Python's < on comparable numbers: a numerically zero difference is equal.
# ponytail: R's Lt/Gt leave equal-valued radicals undecided where SymPy's
# _n2 evalf decides them; this double fallback (1e-12) is used only by the
# Interval/Interval intersection; fix in .rel_eval to cover every caller.
.cmpb <- function(op, a, b) {
  r <- switch(op, lt = Lt(a, b), gt = Gt(a, b), le = Le(a, b), ge = Ge(a, b))
  if (.is_true(r)) return(TRUE)
  if (.is_false(r)) return(FALSE)
  if (.is_comparable(a) && .is_comparable(b)) {
    d <- tryCatch(as.double(Re(.neval(.sub(a, b), list()))), error = function(e) NA_real_)
    if (!is.na(d)) {
      if (abs(d) < 1e-12) d <- 0
      return(switch(op, lt = d < 0, gt = d > 0, le = d <= 0, ge = d >= 0))
    }
  }
  .pb(r)
}
.py_max <- function(a, b) if (.pb(Gt(b, a))) b else a
.py_min <- function(a, b) if (.pb(Lt(b, a))) b else a

# ---------------------------------------------------------------------------
# Class chains and multiple dispatch (sympy.multipledispatch): the handler
# whose signature is most specific for both argument classes wins.

.SET_MRO <- list(Reals = c("Reals", "Interval"), Naturals0 = c("Naturals0", "Naturals"),
                 Complexes = c("Complexes", "CartesianComplexRegion", "ComplexRegion"),
                 CartesianComplexRegion = c("CartesianComplexRegion", "ComplexRegion"),
                 PolarComplexRegion = c("PolarComplexRegion", "ComplexRegion"))
.set_cls <- function(s) if (s$type == "SpecialSet") s$name else s$type
.cls_chain <- function(e) {
  if (is_set(e)) {
    t <- .set_cls(e)
    return(c(if (is.null(.SET_MRO[[t]])) t else .SET_MRO[[t]], "Set", "Basic"))
  }
  if (e$type == "Lambda") return(c("Lambda", "Expr", "Basic"))
  if (.is_expr_s(e) || e$type == "SetExpr") {
    if (e$type == "Integer") return(c(if (e$value == 0) "Zero", "Integer", "Rational", "Number", "Expr", "Basic"))
    if (is_const(e, "oo")) return(c("Infinity", "Number", "Expr", "Basic"))
    if (is_const(e, "-oo")) return(c("NegativeInfinity", "Number", "Expr", "Basic"))
    return(c("Expr", "Basic"))
  }
  "Basic"
}
.reg <- function(A, B, f) list(A, B, f)
.dispatch2 <- function(table, a, b) {
  ca <- .cls_chain(a); cb <- .cls_chain(b)
  best <- NULL; bi <- c(Inf, Inf)
  for (h in table) {
    ia <- match(h[[1L]], ca); ib <- match(h[[2L]], cb)
    if (is.na(ia) || is.na(ib)) next
    if (ia <= bi[1L] && ib <= bi[2L]) { best <- h[[3L]]; bi <- c(ia, ib) }
  }
  best
}

# ---------------------------------------------------------------------------
# Interval

#' @rdname sets
#' @export
Interval <- function(a, b, left_open = FALSE, right_open = FALSE) {
  start <- as_sym(a)
  end <- as_sym(b)
  lo <- .lgl(left_open)
  ro <- .lgl(right_open)
  d <- .sub(end, start)
  er <- function(e) if (is_nan(e)) NA else ask_fact(e, "extended_real")
  if (identical(.f_not(.f_and(list(start, end, d), er)), TRUE))
    .perr("ValueError", "Non-real intervals are not supported")
  if (isTRUE(.fb(tryCatch(Lt(end, start), error = function(e) NULL)))) return(S_EmptySet)
  if (!is_nan(d) && isTRUE(ask_fact(d, "negative"))) return(S_EmptySet)
  if (.py_eq(end, start)) {
    if (lo || ro || is_inf_const(start)) return(S_EmptySet)
    return(FiniteSet(end))
  }
  if (is_const(start, "-oo")) lo <- TRUE
  if (is_const(end, "oo")) ro <- TRUE
  if (is_const(start, "oo") || is_const(end, "-oo")) return(S_EmptySet)
  new_sym("Interval", list(start, end, if (lo) S_true else S_false, if (ro) S_true else S_false))
}

# Interval._is_comparable
.iv_comparable <- function(a, b) {
  A <- .iv(a); B <- .iv(b)
  all(vapply(list(A$a, A$b, B$a, B$b), .is_comparable, NA))
}

# ---------------------------------------------------------------------------
# Membership: Set.contains and the _contains methods

#' @rdname sets
#' @export
contains <- function(s, x) {
  if (inherits(s, "GeoEntity")) return(.geo_contains(s, x))
  if (.is_ab(s)) return(.ab_contains(s, x))
  s <- as_sym(s)
  x <- as_sym(x)
  .rel_to_r(.contains(s, x))
}

# Set.contains: the _contains result, or an unevaluated Contains
.contains <- function(s, x) {
  c <- .set_contains_raw(s, x)
  if (is.null(c)) new_sym("Contains", list(x, s)) else c
}

.set_contains_raw <- function(s, x) {
  switch(.set_cls(s),
         EmptySet = S_false,
         UniversalSet = S_true,
         Interval = , Reals = .interval_contains(s, x),
         FiniteSet = if (x$key %in% vapply(s$args, .k, "")) S_true else
           do.call(Or, lapply(s$args, function(e) Eq(e, x))),
         Union = do.call(Or, lapply(s$args, function(z) .contains(z, x))),
         Intersection = do.call(And, lapply(s$args, function(z) .contains(z, x))),
         Complement = And(.contains(s$args[[1L]], x), Not(.contains(s$args[[2L]], x))),
         Integers = if (!.is_expr_s(x)) S_false else .tfn(ask_fact(x, "integer")),
         Rationals = if (!.is_expr_s(x)) S_false else .tfn(ask_fact(x, "rational")),
         Naturals = , Naturals0 = {
           if (!.is_expr_s(x)) return(S_false)
           sg <- if (s$name == "Naturals") "positive" else "nonnegative"
           it <- ask_fact(x, "integer"); ps <- ask_fact(x, sg)
           if (isTRUE(ps) && isTRUE(it)) S_true else if (identical(it, FALSE) || identical(ps, FALSE)) S_false else NULL
         },
         Complexes = .cr_contains(.as_cregion(s), x),
         ProductSet = {
           if (x$type %in% c("Symbol", "Dummy")) return(NULL)
           if (x$type != "Tuple" || length(x$args) != length(s$args)) return(S_false)
           do.call(And, Map(.contains, s$args, x$args))
         },
         Range = .range_contains(s, x),
         CartesianComplexRegion = , PolarComplexRegion = .cr_contains(s, x),
         DisjointUnion = .disjoint_contains(s, x),
         ImageSet = .imageset_contains(s, x),
         ConditionSet = {
           v <- s$args[[1L]]
           And(subs(s$args[[2L]], v, x), .contains(s$args[[3L]], x))
         },
         NULL)
}

# ImageSet._contains for a one-variable Lambda: solve expr == other for the
# variable over the base set (_solveset_multi) and test the solutions
.imageset_contains <- function(s, other) {
  lam <- s$args[[1L]]; vs <- lam$args[[1L]]$args; expr <- lam$args[[2L]]
  if (length(vs) != 1L || vs[[1L]]$type == "Tuple" || expr$type == "Tuple") return(NULL)
  d <- Dummy(vs[[1L]]$name)
  e <- subs(expr, vs[[1L]], d)
  eq <- Eq(e, other)
  if (.is_false(eq)) return(S_false)
  if (!.has1(eq, d)) return(NULL)
  base <- s$args[[2L]]
  if (.is_reals(base)) {
    r <- Symbol(vs[[1L]]$name, real = TRUE)
    eq <- subs(eq, d, r); d <- r
  }
  sol <- tryCatch(solveset(eq, d, base), error = function(err) NULL)
  if (is.null(sol)) return(NULL)
  .tfn(.f_not(.set_is_empty(sol)))
}

# Interval._contains
.interval_contains <- function(s, x) {
  if (!.is_expr_s(x) || is_nan(x) || identical(ask_fact(x, "real"), FALSE) || .has1(x, S_ComplexInfinity))
    return(S_false)
  I <- .iv(s)
  if (is_const(I$a, "-oo") && is_const(I$b, "oo")) {
    r <- ask_fact(x, "real")
    if (!is.na(r)) return(.tfn(r))
  }
  And(if (I$lo) Lt(I$a, x) else Le(I$a, x), if (I$ro) Lt(x, I$b) else Le(x, I$b))
}

# Complexes is the CartesianComplexRegion of Reals x Reals
.as_cregion <- function(s) if (s$type == "SpecialSet" && s$name == "Complexes")
  new_sym("CartesianComplexRegion", list(new_sym("ProductSet", list(Reals, Reals)))) else s

# ---------------------------------------------------------------------------
# is_empty, is_finite_set

.set_is_empty <- function(s) {
  switch(.set_cls(s),
         Interval = , Reals = {
           I <- .iv(s)
           .fb(if (I$lo || I$ro) Ge(I$a, I$b) else Gt(I$a, I$b))
         },
         FiniteSet = FALSE, EmptySet = TRUE, UniversalSet = FALSE,
         Naturals = , Naturals0 = , Integers = , Rationals = , Complexes = FALSE,
         Union = .f_and(s$args, .set_is_empty),
         ProductSet = .f_or(s$args, .set_is_empty),
         DisjointUnion = .f_and(s$args, .set_is_empty),
         Range = tryCatch(ask_fact(.range_size(s), "zero"), error = function(e) NA),
         NA)
}

.set_is_finite <- function(s) {
  switch(.set_cls(s),
         Interval = , Reals = ask_fact(measure(s), "zero"),
         FiniteSet = TRUE, EmptySet = TRUE, UniversalSet = FALSE,
         Naturals = , Naturals0 = , Integers = , Rationals = , Complexes = FALSE,
         Union = .f_and(s$args, .set_is_finite),
         Intersection = if (isTRUE(.f_or(s$args, .set_is_finite))) TRUE else NA,
         Complement = {
           af <- .set_is_finite(s$args[[1L]])
           if (isTRUE(af)) TRUE else if (identical(af, FALSE) && isTRUE(.set_is_finite(s$args[[2L]]))) FALSE else NA
         },
         ProductSet = .f_or(list(.set_is_empty(s), .f_and(s$args, .set_is_finite))),
         Range = {
           p <- .range_parts(s)
           if (isTRUE(ask_fact(p$start, "integer")) && isTRUE(ask_fact(p$stop, "integer"))) TRUE
           else tryCatch(ask_fact(.range_size(s), "finite"), error = function(e) NA)
         },
         NA)
}

# ---------------------------------------------------------------------------
# is_subset: Set.is_subset, _eval_is_subset, _eval_is_superset and the
# is_subset_sets handlers

#' @rdname sets
#' @export
is_subset <- function(s, t) .is_subset(as_sym(s), as_sym(t))

.is_subset <- function(self, other) {
  if (!is_set(other)) .perr("ValueError", paste0("Unknown argument '", .str(other), "'"))
  if (.set_eq(self, other)) return(TRUE)
  ie <- .set_is_empty(self)
  if (isTRUE(ie)) return(TRUE)
  if (identical(ie, FALSE) && isTRUE(.set_is_empty(other))) return(FALSE)
  if (identical(.set_is_finite(self), FALSE) && isTRUE(.set_is_finite(other))) return(FALSE)
  r <- .eval_is_subset(self, other)
  if (!is.na(r)) return(r)
  r <- .eval_is_superset(other, self)
  if (!is.na(r)) return(r)
  r <- .dispatch2(.IS_SUBSET_SETS, self, other)(self, other)
  if (is.null(r)) NA else r
}

.eval_is_subset <- function(self, other) {
  switch(.set_cls(self),
         Union = .f_and(self$args, function(a) .is_subset(a, other)),
         FiniteSet = .f_and(self$args, function(e) .set_contains_raw(other, e)),
         Complement = if (isTRUE(.is_subset(self$args[[1L]], other))) TRUE else NA,
         Naturals = .is_subset(Range(1, S_Infinity), other),
         Naturals0 = .is_subset(Range(0, S_Infinity), other),
         Integers = .is_subset(Range(S_NegativeInfinity, S_Infinity), other),
         NA)
}

.eval_is_superset <- function(self, other) {
  switch(.set_cls(self),
         Union = if (isTRUE(.f_or(self$args, function(a) .is_subset(other, a)))) TRUE else NA,
         Naturals = .is_subset(other, Range(1, S_Infinity)),
         Naturals0 = .is_subset(other, Range(0, S_Infinity)),
         Integers = .is_subset(other, Range(S_NegativeInfinity, S_Infinity)),
         NA)
}

.ss_interval_interval <- function(a, b) {
  A <- .iv(a); B <- .iv(b)
  left_in <- if (B$lo && !A$lo) .fb(Gt(A$a, B$a)) else .fb(Ge(A$a, B$a))
  if (identical(left_in, FALSE)) return(FALSE)
  right_in <- if (B$ro && !A$ro) .fb(Lt(A$b, B$b)) else .fb(Le(A$b, B$b))
  .f_and(list(left_in, right_in))
}

.ss_interval_union <- function(a, b) {
  if (isTRUE(.f_or(b$args, function(s) .is_subset(a, s)))) return(TRUE)
  if (all(vapply(b$args, function(s) .is_ivl(s) || s$type == "FiniteSet", NA))) {
    A <- .iv(a)
    ivs <- lapply(Filter(.is_ivl, b$args), .iv)
    if (all(vapply(ivs, function(s) isTRUE(.fb(Lt(A$a, s$a))), NA))) return(FALSE)
    if (all(vapply(ivs, function(s) isTRUE(.fb(Gt(A$b, s$b))), NA))) return(FALSE)
    if (isTRUE(ask_fact(measure(a), "nonzero"))) {
      no_overlap <- function(s1, s2) .f_or(list(.fb(Le(s1$b, s2$a)), .fb(Ge(s1$a, s2$b))))
      if (all(vapply(ivs, function(s) isTRUE(no_overlap(s, A)), NA))) return(FALSE)
    }
  }
  NA
}

.ss_range_finiteset <- function(a, b) {
  n <- .range_size(a)
  if (!is_number_type(n) && !is_inf_const(n)) return(NA)
  if (.pb(Gt(n, .int(length(b$args))))) return(FALSE)
  if (any(vapply(a$args, function(z) length(free_symbols(z)) > 0L, NA)))
    return(.f_and(seq_len(as.integer(as.character(n$value))) - 1L, function(k) .contains(b, .range_item(a, k))))
  a_set <- vapply(seq_len(as.integer(as.character(n$value))) - 1L, function(k) .range_item(a, k)$key, "")
  b_remaining <- length(b$args)
  cnt <- 0L
  for (e in b$args) {
    if (e$type == "Integer") a_set <- setdiff(a_set, e$key)
    else if (!isTRUE(.f_not(ask_fact(e, "integer")))) cnt <- cnt + 1L
    b_remaining <- b_remaining - 1L
    if (length(a_set) > b_remaining + cnt) return(FALSE)
    if (!length(a_set)) return(TRUE)
  }
  NA
}

.ss_range_interval <- function(a, b) {
  p <- .range_parts(a)
  if (!isTRUE(ask_fact(p$step, "positive"))) return(NA)
  B <- .iv(b)
  inf <- .range_inf(a); sup <- .range_sup(a)
  cl <- if (B$lo && isTRUE(ask_fact(inf, "finite"))) Gt(inf, B$a) else Ge(inf, B$a)
  cr <- if (B$ro && isTRUE(ask_fact(sup, "finite"))) Lt(sup, B$b) else Le(sup, B$b)
  .f_and(list(cl, cr))
}

# f.expr.as_poly(*var, domain='QQ') is not None
.poly_qq <- function(e, vars) {
  vk <- vapply(vars, .k, "")
  if (!all(vapply(free_symbols(e), .k, "") %in% vk)) return(FALSE)
  if (!do.call(is_polynomial, c(list(e), vars))) return(FALSE)
  rat_ok <- function(t) {
    fs <- if (t$type == "Mul") t$args else list(t)
    all(vapply(fs, function(f) is_exact_number(f) || f$key %in% vk ||
                 (f$type == "Pow" && f$args[[1L]]$key %in% vk && f$args[[2L]]$type == "Integer"), NA))
  }
  ex <- expand(e)
  all(vapply(if (ex$type == "Add") ex$args else list(ex), rat_ok, NA))
}

.IS_SUBSET_SETS <- list(
  .reg("Set", "Set", function(a, b) NA),
  .reg("Interval", "Interval", function(a, b) .ss_interval_interval(a, b)),
  .reg("Interval", "FiniteSet", function(a, b) if (identical(.f_not(ask_fact(measure(a), "zero")), TRUE)) FALSE else NA),
  .reg("Interval", "Union", function(a, b) .ss_interval_union(a, b)),
  .reg("Range", "Range", function(a, b) {
    pa <- .range_parts(a); pb <- .range_parts(b)
    if (.py_eq(pa$step, S_One) && .py_eq(pb$step, S_One))
      return(.f_and(list(.fb(Ge(pa$start, pb$start)), .fb(Le(pa$stop, pb$stop)))))
    NA
  }),
  .reg("Range", "Interval", function(a, b) .ss_range_interval(a, b)),
  .reg("Range", "FiniteSet", function(a, b) .ss_range_finiteset(a, b)),
  .reg("ImageSet", "Reals", function(a, b) {
    lam <- a$args[[1L]]
    if (isTRUE(.is_subset(a$args[[2L]], b)) && .poly_qq(lam$args[[2L]], lam$args[[1L]]$args)) TRUE else NA
  }),
  .reg("Interval", "Range", function(a, b) if (isTRUE(ask_fact(measure(a), "extended_nonzero"))) FALSE else NA),
  .reg("Interval", "Rationals", function(a, b) if (isTRUE(ask_fact(measure(a), "extended_nonzero"))) FALSE else NA),
  .reg("Range", "Complexes", function(a, b) TRUE),
  .reg("Complexes", "Interval", function(a, b) FALSE),
  .reg("Complexes", "Range", function(a, b) FALSE),
  .reg("Complexes", "Rationals", function(a, b) FALSE),
  .reg("Rationals", "Complexes", function(a, b) TRUE),
  .reg("Reals", "Complexes", function(a, b) TRUE),
  .reg("Interval", "Reals", function(a, b) TRUE),
  .reg("Integers", "Rationals", function(a, b) TRUE),
  .reg("Rationals", "Reals", function(a, b) TRUE),
  .reg("Rationals", "Range", function(a, b) FALSE),
  .reg("Range", "Rationals", function(a, b) TRUE),
  .reg("ProductSet", "ProductSet", function(a, b) {
    if (length(a$args) != length(b$args)) return(FALSE)
    .f_and(seq_along(a$args), function(k) .is_subset(a$args[[k]], b$args[[k]]))
  }),
  .reg("ProductSet", "FiniteSet", function(a, b) NA))

# ---------------------------------------------------------------------------
# Equality of sets: Interval._eval_Eq and the _eval_is_eq handlers in
# sympy/sets/handlers/comparison.py, as is_eq runs them. NULL when undecided.

.set_is_eq <- function(lhs, rhs) {
  for (p in list(list(lhs, rhs), list(rhs, lhs))) {
    s1 <- p[[1L]]; s2 <- p[[2L]]
    if (is_set(s1) && .is_ivl(s1) && !(is_set(s2) && .is_ivl(s2))) {
      if (!is_set(s2) || s2$type == "FiniteSet") return(S_false)
    }
  }
  h1 <- .dispatch2(.EVAL_IS_EQ, lhs, rhs); h2 <- .dispatch2(.EVAL_IS_EQ, rhs, lhs)
  r <- if (is.null(h1)) NULL else h1(lhs, rhs)
  if (is.null(r) && !identical(h1, h2) && !is.null(h2)) r <- h2(rhs, lhs)
  if (!is.null(r)) return(r)
  if (is_set(lhs) && is_set(rhs) && .set_eq(lhs, rhs)) return(S_true)
  NULL
}

.EVAL_IS_EQ <- list(
  .reg("Interval", "FiniteSet", function(l, r) S_false),
  .reg("FiniteSet", "Interval", function(l, r) S_false),
  .reg("Interval", "Interval", function(l, r) {
    A <- .iv(l); B <- .iv(r)
    And(Eq(A$a, B$a), Eq(A$b, B$b), if (A$lo == B$lo) S_true else S_false, if (A$ro == B$ro) S_true else S_false)
  }),
  .reg("FiniteSet", "FiniteSet", function(l, r) {
    lk <- vapply(l$args, .k, ""); rk <- vapply(r$args, .k, "")
    .tfn(.f_and(list(function() .f_and(r$args[!(rk %in% lk)], function(e) .set_contains_raw(l, e)),
                     function() .f_and(l$args[!(lk %in% rk)], function(e) .set_contains_raw(r, e))),
                function(g) g()))
  }),
  .reg("ProductSet", "ProductSet", function(l, r) {
    if (length(l$args) != length(r$args)) return(S_false)
    .tfn(.f_and(seq_along(l$args), function(k) .fb(Eq(l$args[[k]], r$args[[k]]))))
  }),
  .reg("Set", "Basic", function(l, r) S_false),
  .reg("Set", "Set", function(l, r) .tfn(.f_and(list(list(l, r), list(r, l)), function(p) .is_subset(p[[1L]], p[[2L]])))))

# ---------------------------------------------------------------------------
# Union: Union.__new__, simplify_union and union_sets

#' @rdname sets
#' @export
Union <- function(...) {
  flat <- list()
  for (a in lapply(list(...), as_sym)) {
    if (is_set(a) && .set_eq(a, UniversalSet)) return(UniversalSet)
    if (a$type == "EmptySet") next
    if (a$type == "Union") flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  }
  .simplify_union(flat)
}

# Union(*args, evaluate=False): args ordered by Set._infimum_key
.union_raw <- function(args) new_sym("Union", .set_print_order(args))

.simplify_union <- function(args) {
  if (!length(args)) return(S_EmptySet)
  for (a in args) if (!is_set(a)) .perr("TypeError", "Input args to Union must be Sets")
  isfs <- vapply(args, function(a) a$type == "FiniteSet", NA)
  if (sum(isfs) > 1L)
    args <- c(list(FiniteSet(unlist(lapply(args[isfs], `[[`, "args"), recursive = FALSE))), args[!isfs])
  # ponytail: Python iterates a hash-ordered set here; argument order is used.
  args <- .uniq_sets(args)
  repeat {
    changed <- FALSE
    for (si in seq_along(args)) {
      for (ti in seq_along(args)) {
        if (ti == si) next
        new <- .dispatch2(.UNION_SETS, args[[si]], args[[ti]])(args[[si]], args[[ti]])
        if (!is.null(new)) {
          if (is_sym(new)) new <- list(new)
          args <- .uniq_sets(c(args[-c(si, ti)], new))
          changed <- TRUE
          break
        }
      }
      if (changed) break
    }
    if (!changed) break
  }
  if (length(args) == 1L) return(args[[1L]])
  .union_raw(args)
}

.us_interval_interval <- function(a, b) {
  if (!.iv_comparable(a, b)) return(NULL)
  A <- .iv(a); B <- .iv(b)
  end <- Min(A$b, B$b); start <- Max(A$a, B$a)
  if (.pb(Lt(end, start)) || (.py_eq(end, start) && !.py_in(a, end) && !.py_in(b, end))) return(NULL)
  start <- Min(A$a, B$a); end <- Max(A$b, B$b)
  lo <- (!.py_eq(A$a, start) || A$lo) && (!.py_eq(B$a, start) || B$lo)
  ro <- (!.py_eq(A$b, end) || A$ro) && (!.py_eq(B$b, end) || B$ro)
  Interval(start, end, lo, ro)
}

.us_interval_set <- function(a, b) {
  A <- .iv(a)
  olf <- A$lo && .is_true(.contains(b, A$a)) && isTRUE(ask_fact(A$a, "finite"))
  orf <- A$ro && .is_true(.contains(b, A$b)) && isTRUE(ask_fact(A$b, "finite"))
  if (olf || orf) {
    open_left <- A$lo && !.py_in(b, A$a)
    open_right <- A$ro && !.py_in(b, A$b)
    return(.uniq_sets(list(Interval(A$a, A$b, open_left, open_right), b)))
  }
  NULL
}

.UNION_SETS <- list(
  .reg("Naturals0", "Naturals", function(a, b) a),
  .reg("Rationals", "Naturals", function(a, b) a),
  .reg("Rationals", "Naturals0", function(a, b) a),
  .reg("Reals", "Naturals", function(a, b) a),
  .reg("Reals", "Naturals0", function(a, b) a),
  .reg("Reals", "Rationals", function(a, b) a),
  .reg("Integers", "Set", function(a, b) {
    i <- Intersection(a, b)
    if (.set_eq(i, a)) b else if (.set_eq(i, b)) a else NULL
  }),
  .reg("ComplexRegion", "Set", function(a, b) .cr_union(.as_cregion(a), .as_cregion(b))),
  .reg("EmptySet", "Set", function(a, b) b),
  .reg("UniversalSet", "Set", function(a, b) a),
  .reg("ProductSet", "ProductSet", function(a, b) {
    if (isTRUE(.is_subset(b, a))) return(a)
    if (length(b$args) != length(a$args)) return(NULL)
    if (length(a$args) == 2L) {
      if (.set_eq(a$args[[1L]], b$args[[1L]])) return(ProductSet(a$args[[1L]], Union(a$args[[2L]], b$args[[2L]])))
      if (.set_eq(a$args[[2L]], b$args[[2L]])) return(ProductSet(Union(a$args[[1L]], b$args[[1L]]), a$args[[2L]]))
    }
    NULL
  }),
  .reg("ProductSet", "Set", function(a, b) if (isTRUE(.is_subset(b, a))) a else NULL),
  .reg("Interval", "Interval", function(a, b) .us_interval_interval(a, b)),
  .reg("Interval", "UniversalSet", function(a, b) UniversalSet),
  .reg("Interval", "Set", function(a, b) .us_interval_set(a, b)),
  .reg("FiniteSet", "FiniteSet", function(a, b) FiniteSet(c(a$args, b$args))),
  .reg("FiniteSet", "Set", function(a, b) {
    inb <- vapply(a$args, function(x) .is_true(.contains(b, x)), NA)
    if (any(inb)) return(list(FiniteSet(a$args[!inb]), b))
    NULL
  }),
  .reg("Set", "Set", function(a, b) NULL))

# ---------------------------------------------------------------------------
# Intersection: Intersection.__new__, simplify_intersection,
# _handle_finite_sets and intersection_sets

#' @rdname sets
#' @export
Intersection <- function(...) {
  args <- .bool_print_order(.uniq_sets(lapply(list(...), as_sym)))
  flat <- list()
  for (a in args) {
    if (is_set(a) && a$type == "EmptySet") return(S_EmptySet)
    if (is_set(a) && .set_eq(a, UniversalSet)) next
    if (a$type == "Intersection") flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  }
  .simplify_intersection(flat)
}

.intersection_raw <- function(args) new_sym("Intersection", .set_print_order(args))

.simplify_intersection <- function(args) {
  if (!length(args)) return(UniversalSet)
  for (a in args) if (!is_set(a)) .perr("TypeError", "Input args to Union must be Sets")
  if (any(vapply(args, function(a) a$type == "EmptySet", NA))) return(S_EmptySet)
  rv <- .handle_finite_sets(args)
  if (!is.null(rv)) return(rv)
  for (k in seq_along(args)) if (args[[k]]$type == "Union") {
    s <- args[[k]]
    other <- .uniq_sets(args[-k])
    if (length(other)) {
      o <- do.call(Intersection, other)
      return(do.call(Union, lapply(s$args, function(a) Intersection(a, o))))
    }
    return(do.call(Union, s$args))
  }
  for (k in seq_along(args)) if (args[[k]]$type == "Complement") {
    s <- args[[k]]
    return(Complement(do.call(Intersection, c(args[-k], list(s$args[[1L]]))), s$args[[2L]]))
  }
  args <- .uniq_sets(args)
  repeat {
    changed <- FALSE
    for (si in seq_along(args)) {
      for (ti in seq_along(args)) {
        if (ti == si) next
        new <- .intersection_sets(args[[si]], args[[ti]])
        if (!is.null(new)) {
          args <- .uniq_sets(c(args[-c(si, ti)], list(new)))
          changed <- TRUE
          break
        }
      }
      if (changed) break
    }
    if (!changed) break
  }
  if (length(args) == 1L) return(args[[1L]])
  .intersection_raw(args)
}

# Intersection._handle_finite_sets
.fs_in <- function(S, e) .fb(.contains(S, e))
.uniq_keys <- function(xs) xs[!duplicated(vapply(xs, `[[`, "", "key"))]
.handle_finite_sets <- function(args) {
  isfs <- vapply(args, function(z) z$type == "FiniteSet", NA)
  if (!any(isfs)) return(NULL)
  fs_sets <- lapply(args[isfs], `[[`, "args")
  others <- args[!isfs]
  definite <- list()
  drop <- function(l, e) Filter(function(z) !identical(z$key, e$key), l)
  for (e in .uniq_keys(unlist(fs_sets, recursive = FALSE))) {
    inall <- .f_and(args, function(z) .contains(z, e))
    if (isTRUE(inall)) definite <- c(definite, list(e))
    if (!is.na(inall)) fs_sets <- lapply(fs_sets, drop, e = e)
  }
  fs_elements <- .uniq_keys(unlist(fs_sets, recursive = FALSE))
  fs_symsets <- lapply(fs_sets, FiniteSet)
  while (length(fs_elements)) {
    found <- FALSE
    for (e in fs_elements) {
      infs <- .f_and(fs_symsets, function(z) .contains(z, e))
      if (isTRUE(infs)) definite <- c(definite, list(e))
      if (!is.na(infs)) {
        for (n in seq_along(fs_sets)) if (e$key %in% vapply(fs_sets[[n]], `[[`, "", "key")) {
          fs_sets[[n]] <- drop(fs_sets[[n]], e); fs_symsets[[n]] <- FiniteSet(fs_sets[[n]])
        }
        fs_elements <- drop(fs_elements, e)
        found <- TRUE
        break
      }
    }
    if (!found) break
  }
  if (!all(lengths(fs_sets) > 0)) fs_sets <- list(list())
  if (length(definite)) fs_sets <- lapply(fs_sets, function(fs) .uniq_keys(c(fs, definite)))
  if (length(fs_sets) == 1L && !length(fs_sets[[1L]])) return(S_EmptySet)
  sets <- lapply(fs_sets, FiniteSet)
  allel <- .uniq_keys(unlist(fs_sets, recursive = FALSE))
  redundant <- function(o) all(vapply(allel, function(e) isTRUE(.fs_in(o, e)), NA))
  others <- Filter(Negate(redundant), others)
  if (length(others)) {
    rest <- do.call(Intersection, others)
    if (rest$type == "EmptySet") return(S_EmptySet)
    sets <- if (rest$type == "Intersection") c(sets, rest$args) else c(sets, list(rest))
  }
  if (length(sets) == 1L) return(sets[[1L]])
  .intersection_raw(.uniq_sets(sets))
}

.intersection_sets <- function(a, b) .dispatch2(.INTERSECTION_SETS, a, b)(a, b)

.is_interval_interval <- function(a, b) {
  if (.reals_like(a)) return(b)
  if (!.iv_comparable(a, b)) return(NULL)
  A <- .iv(a); B <- .iv(b)
  pick <- function(ea, eb, which) {
    if (identical(ea$key, eb$key)) return(ea)
    if (.has_float(eb) && !.has_float(ea)) return(eb)
    if (.has_float(ea) && !.has_float(eb)) return(ea)
    .iv(.bool_print_order(list(a, b))[[1L]])[[which]]
  }
  if (!(.cmpb("le", A$a, B$b) && .cmpb("le", B$a, A$b))) return(S_EmptySet)
  if (.cmpb("lt", A$a, B$a)) { start <- B$a; lo <- B$lo }
  else if (.cmpb("gt", A$a, B$a)) { start <- A$a; lo <- A$lo }
  else { start <- pick(A$a, B$a, "a"); lo <- A$lo || B$lo }
  if (.cmpb("lt", A$b, B$b)) { end <- A$b; ro <- A$ro }
  else if (.cmpb("gt", A$b, B$b)) { end <- B$b; ro <- B$ro }
  else { end <- pick(A$b, B$b, "b"); ro <- A$ro || B$ro }
  d <- .sub(end, start)
  if (is_number_type(d) && isTRUE(ask_fact(d, "zero")) && (lo || ro)) return(S_EmptySet)
  Interval(start, end, lo, ro)
}

# _intlike_interval
.intlike_interval <- function(a, b) {
  B <- .iv(b)
  if (is_const(B$a, "-oo") && is_const(B$b, "oo")) return(a)
  lo <- .py_max(set_inf(a), fcall("ceiling", list(B$a)))
  s <- tryCatch(Range(lo, .add(list(fcall("floor", list(B$b)), S_One))), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  .intersection_sets(s, b)
}

.is_range_interval <- function(a, b) {
  B <- .iv(b)
  if (!all(vapply(c(a$args, list(B$a, B$b)), function(v) .is_number_e(v) || is_inf_const(v), NA))) return(NULL)
  if (.range_empty(a)) return(S_EmptySet)
  start <- fcall("ceiling", list(.py_max(B$a, .range_inf(a))))
  if (!.py_in(b, start)) start <- .add(list(start, S_One))
  end <- fcall("floor", list(.py_min(B$b, .range_sup(a))))
  if (!.py_in(b, end)) end <- .sub(end, S_One)
  .intersection_sets(a, Range(start, .add(list(end, S_One))))
}

.is_range_range <- function(a, b) {
  if (!all(vapply(c(a$args, b$args), function(v) .is_number_e(v) || is_inf_const(v), NA))) return(NULL)
  if (.range_empty(a) || .range_empty(b)) return(S_EmptySet)
  if (.pb(Lt(.range_sup(b), .range_inf(a))) || .pb(Gt(.range_inf(b), .range_sup(a)))) return(S_EmptySet)
  .range_intersect(a, b)
}

.is_imageset_set <- function(self, other) {
  lam <- self$args[[1L]]; vs <- lam$args[[1L]]$args
  if (length(vs) > 1L || vs[[1L]]$type == "Tuple") return(NULL)
  base <- self$args[[2L]]
  if (!(base$type == "SpecialSet" && base$name == "Integers")) return(NULL)
  gm <- NULL
  if (other$type == "ImageSet" && other$args[[2L]]$type == "SpecialSet" && other$args[[2L]]$name == "Integers") {
    m <- Dummy("m")
    gm <- subs(other$args[[1L]]$args[[2L]], other$args[[1L]]$args[[1L]]$args[[1L]], m)
  } else if (other$type == "SpecialSet" && other$name == "Integers") {
    m <- Dummy("m")
    gm <- m
  }
  if (is.null(gm)) return(NULL)
  fn <- lam$args[[2L]]; n <- vs[[1L]]
  solns <- tryCatch(diophantine(.sub(fn, gm), syms = list(n, m), permute = TRUE), error = function(e) NULL)
  if (is.null(solns)) return(NULL)
  sl <- if (solns$type == "FiniteSet") solns$args else list()
  if (!length(sl)) return(S_EmptySet)
  if (any(vapply(sl, function(t) length(free_symbols(t)) > 0L, NA))) {
    if (length(sl) != 1L) return(NULL)
    soln <- sl[[1L]]$args[[1L]]
    t <- free_symbols(soln)[[1L]]
    return(imageset(Lambda(n, expand(subs(fn, n, subs(soln, t, n)))), Integers))
  }
  FiniteSet(lapply(sl, function(s) subs(fn, n, s$args[[1L]])))
}

.solution_union <- function(exprs, sym) {
  do.call(Union, lapply(exprs, function(i) {
    r <- tryCatch(solve_linear(i, 0, sym), error = function(e) NULL)
    if (!is.null(r) && identical(r[[1L]]$key, sym$key)) FiniteSet(r[[2L]]) else ConditionSet(sym, Eq(i, 0))
  }))
}

.is_imageset_reals <- function(self, other) {
  lam <- self$args[[1L]]; vs <- lam$args[[1L]]$args
  if (length(vs) > 1L || vs[[1L]]$type == "Tuple") return(NULL)
  base <- self$args[[2L]]
  f <- lam$args[[2L]]; n <- vs[[1L]]
  n_ <- Dummy(n$name, real = TRUE)
  ri <- .re_im(subs(f, n, n_))
  re <- subs(ri[[1L]], n_, n)
  im <- subs(expand_complex(ri[[2L]]), n_, n)
  z <- ask_fact(im, "zero")
  if (identical(z, FALSE)) return(S_EmptySet)
  if (!isTRUE(z)) {
    fsy <- free_symbols(im)
    if (!(length(fsy) == 1L && identical(fsy[[1L]]$key, n$key))) return(NULL)
    nu <- numer(im)
    base <- Intersection(base, .solution_union(if (nu$type == "Mul") nu$args else list(nu), n))
  }
  base <- Complement(base, .solution_union(denoms(f), n))
  imageset(Lambda(n, re), base)
}

.is_imageset_interval <- function(self, other) {
  lam <- self$args[[1L]]; vs <- lam$args[[1L]]$args
  if (length(vs) > 1L || vs[[1L]]$type == "Tuple") return(NULL)
  base <- self$args[[2L]]
  f <- lam$args[[2L]]; n <- vs[[1L]]
  if (!is_polynomial(f, n)) return(NULL)
  cl <- .coeff_list(f, n)
  if (length(cl) != 2L) return(NULL)
  a <- cl[[1L]]; b <- cl[[2L]]
  I <- .iv(other)
  if (!all(vapply(list(a, b, I$a, I$b), .is_comparable, NA))) return(NULL)
  new_int <- if (.is_true(Gt(a, 0))) Interval(.div(.sub(I$a, b), a), .div(.sub(I$b, b), a), I$lo, I$ro)
             else if (.is_true(Lt(a, 0))) Interval(.div(.sub(I$b, b), a), .div(.sub(I$a, b), a), I$ro, I$lo)
  if (is.null(new_int)) return(NULL)
  nv <- Intersection(new_int, base)
  if (nv$type == "Range" && !is_inf_const(.range_size(nv)))
    return(FiniteSet(lapply(seq_len(as.integer(as.character(.range_size(nv)$value))) - 1L,
                            function(j) .add(list(.mul(list(a, .range_item(nv, j))), b)))))
  if (nv$type != "Intersection") return(imageset(Lambda(n, .add(list(.mul(list(a, n)), b))), nv))
  NULL
}

.INTERSECTION_SETS <- list(
  .reg("ConditionSet", "ConditionSet", function(a, b) NULL),
  .reg("ConditionSet", "Set", function(a, b) ConditionSet(a$args[[1L]], a$args[[2L]], Intersection(a$args[[3L]], b))),
  .reg("Naturals", "Integers", function(a, b) a),
  .reg("Naturals", "Naturals", function(a, b) if (a$name == "Naturals") a else b),
  .reg("Interval", "Naturals", function(a, b) .intersection_sets(b, a)),
  .reg("ComplexRegion", "Set", function(a, b) .cr_intersect(.as_cregion(a), .as_cregion(b))),
  .reg("Integers", "Reals", function(a, b) a),
  .reg("Range", "Interval", function(a, b) .is_range_interval(a, b)),
  .reg("Range", "Naturals", function(a, b) .intersection_sets(a, Interval(set_inf(b), S_Infinity))),
  .reg("Range", "Range", function(a, b) .is_range_range(a, b)),
  .reg("Range", "Integers", function(a, b) a),
  .reg("Range", "Rationals", function(a, b) a),
  .reg("ImageSet", "Set", function(a, b) .is_imageset_set(a, b)),
  .reg("ImageSet", "Reals", function(a, b) .is_imageset_reals(a, b)),
  .reg("ImageSet", "Interval", function(a, b) .is_imageset_interval(a, b)),
  .reg("ProductSet", "ProductSet", function(a, b) {
    if (length(b$args) != length(a$args)) return(S_EmptySet)
    do.call(ProductSet, Map(Intersection, a$args, b$args))
  }),
  .reg("Interval", "Interval", function(a, b) .is_interval_interval(a, b)),
  .reg("EmptySet", "Set", function(a, b) S_EmptySet),
  .reg("UniversalSet", "Set", function(a, b) b),
  .reg("FiniteSet", "FiniteSet", function(a, b) {
    bk <- vapply(b$args, .k, "")
    FiniteSet(Filter(function(e) e$key %in% bk, a$args))
  }),
  .reg("FiniteSet", "Set", function(a, b) tryCatch(FiniteSet(Filter(function(el) .py_in(b, el), a$args)),
                                                   TypeError = function(e) NULL)),
  .reg("Set", "Set", function(a, b) NULL),
  .reg("Integers", "Rationals", function(a, b) a),
  .reg("Naturals", "Rationals", function(a, b) a),
  .reg("Rationals", "Reals", function(a, b) a),
  .reg("Integers", "Interval", function(a, b) .intlike_interval(a, b)),
  .reg("Naturals", "Interval", function(a, b) .intlike_interval(a, b)))

# ---------------------------------------------------------------------------
# Complement and SymmetricDifference

#' @rdname sets
#' @export
Complement <- function(s, t) .complement_reduce(as_sym(s), as_sym(t))

.complement_raw <- function(A, B) new_sym("Complement", list(A, B))

# Complement.reduce
.complement_reduce <- function(A, B) {
  if ((is_set(B) && .set_eq(B, UniversalSet)) || isTRUE(.is_subset(A, B))) return(S_EmptySet)
  if (B$type == "Union") return(do.call(Intersection, lapply(B$args, function(s) Complement(A, s))))
  r <- .set_complement(B, A)
  if (is.null(r)) .complement_raw(A, B) else r
}

# self._complement(other): other - self
.set_complement <- function(self, other) {
  switch(.set_cls(self),
         Interval = , Reals = {
           I <- .iv(self)
           if (.set_eq(other, Reals))
             return(Union(Interval(S_NegativeInfinity, I$a, TRUE, !I$lo), Interval(I$b, S_Infinity, !I$ro, TRUE)))
           if (other$type == "FiniteSet" && !any(vapply(other$args, .is_number_e, NA))) return(NULL)
           .set_complement_base(self, other)
         },
         Union = do.call(Intersection, lapply(self$args, function(s) Complement(other, s))),
         EmptySet = other,
         UniversalSet = S_EmptySet,
         FiniteSet = .fs_complement(self, other),
         .set_complement_base(self, other))
}

# Set._complement
.set_complement_base <- function(self, other) {
  if (self$type == "ProductSet" && other$type == "ProductSet") {
    if (length(self$args) != length(other$args)) return(other)
    k <- length(self$args)
    return(do.call(Union, lapply(seq_len(k), function(n)
      do.call(ProductSet, lapply(seq_len(k), function(i)
        if (i != n) other$args[[i]] else Complement(other$args[[i]], self$args[[i]]))))))
  }
  if (.is_ivl(other)) {
    if (.is_ivl(self) || self$type == "FiniteSet") return(Intersection(other, Complement(Reals, self)))
    return(NULL)
  }
  if (other$type == "Union") return(do.call(Union, lapply(other$args, function(o) Complement(o, self))))
  if (other$type == "Complement") return(.complement_raw(other$args[[1L]], Union(other$args[[2L]], self)))
  if (other$type == "EmptySet") return(S_EmptySet)
  if (other$type == "FiniteSet") {
    v <- vapply(other$args, function(x) .fb(.contains(self, x)), NA)
    return(Union(FiniteSet(other$args[v %in% FALSE]),
                 if (any(is.na(v))) .complement_raw(FiniteSet(other$args[is.na(v)]), self) else S_EmptySet))
  }
  NULL
}

# FiniteSet._complement
.fs_complement <- function(self, other) {
  if (.is_ivl(other)) {
    nums <- list(); syms <- list()
    for (m in self$args) {
      r <- if (.is_expr_s(m)) ask_fact(m, "real") else NA
      if (.is_number_e(m) && isTRUE(r)) nums[[length(nums) + 1L]] <- m
      else if (!identical(r, FALSE)) syms[[length(syms) + 1L]] <- m
    }
    if (.set_eq(other, Reals) && length(nums)) {
      nums <- .sort_set_elems(nums)
      ivs <- list(Interval(S_NegativeInfinity, nums[[1L]], TRUE, TRUE))
      for (k in seq_len(length(nums) - 1L)) ivs[[length(ivs) + 1L]] <- Interval(nums[[k]], nums[[k + 1L]], TRUE, TRUE)
      ivs[[length(ivs) + 1L]] <- Interval(nums[[length(nums)]], S_Infinity, TRUE, TRUE)
      u <- if (length(ivs) == 1L) ivs[[1L]] else .union_raw(ivs)
      if (length(syms)) return(.complement_raw(u, FiniteSet(syms)))
      return(u)
    }
    if (!length(nums)) return(if (length(syms)) .complement_raw(other, FiniteSet(syms)) else other)
  } else if (other$type == "FiniteSet") {
    unk <- Filter(function(i) { c <- .contains(other, i); !.is_true(c) && !.is_false(c) }, self$args)
    unk <- FiniteSet(unk)
    if (.set_eq(unk, self)) return(NULL)
    not_true <- Filter(function(i) !.is_true(.contains(self, i)), other$args)
    return(Complement(FiniteSet(not_true), unk))
  }
  .set_complement_base(self, other)
}

#' @rdname sets
#' @export
SymmetricDifference <- function(s, t) {
  A <- as_sym(s); B <- as_sym(t)
  switch(.set_cls(B),
         EmptySet = , UniversalSet = A,
         Union(Complement(B, A), Complement(A, B)))
}

# ---------------------------------------------------------------------------
# measure, boundary, inf and sup

#' @rdname sets
#' @export
measure <- function(s) {
  s <- as_sym(s)
  switch(.set_cls(s),
         Interval = , Reals = { I <- .iv(s); .sub(I$b, I$a) },
         FiniteSet = S_Zero, EmptySet = S_Zero,
         UniversalSet = S_Infinity,
         Union = .union_measure(s$args),
         ProductSet = .mul(lapply(s$args, measure)),
         CartesianComplexRegion = , PolarComplexRegion = measure(s$args[[1L]]),
         Integers = , Naturals = , Naturals0 = , Rationals = , Complexes = S_Zero,
         .perr("NotImplementedError", paste0("(", .str(s), ")._measure")))
}

#' @rdname sets
#' @export
boundary <- function(s) {
  s <- as_sym(s)
  switch(.set_cls(s),
         Interval = , Reals = { I <- .iv(s); FiniteSet(Filter(Negate(is_inf_const), list(I$a, I$b))) },
         FiniteSet = s, EmptySet = s, UniversalSet = S_EmptySet,
         Union = {
           a <- s$args
           do.call(Union, lapply(seq_along(a), function(i) {
             b <- boundary(a[[i]])
             for (j in seq_along(a)) if (j != i) b <- Complement(b, Complement(a[[j]], boundary(a[[j]])))
             b
           }))
         },
         Integers = , Naturals = , Naturals0 = s, Rationals = Reals, Complexes = S_EmptySet,
         .perr("NotImplementedError", paste0("(", .str(s), ")._boundary")))
}

#' @rdname sets
#' @export
set_inf <- function(s) {
  s <- as_sym(s)
  switch(.set_cls(s), Interval = , Reals = .iv(s)$a, FiniteSet = fcall("Min", s$args),
         Union = fcall("Min", lapply(s$args, set_inf)),
         Naturals = S_One, Naturals0 = S_Zero, Integers = , Rationals = S_NegativeInfinity,
         Range = .range_inf(s),
         .perr("NotImplementedError", paste0("(", .str(s), ")._inf")))
}

#' @rdname sets
#' @export
set_sup <- function(s) {
  s <- as_sym(s)
  switch(.set_cls(s), Interval = , Reals = .iv(s)$b, FiniteSet = fcall("Max", s$args),
         Union = fcall("Max", lapply(s$args, set_sup)),
         Naturals = , Naturals0 = , Integers = , Rationals = S_Infinity,
         Range = .range_sup(s),
         .perr("NotImplementedError", paste0("(", .str(s), ")._sup")))
}

# ---------------------------------------------------------------------------
# ImageSet and imageset

#' @rdname sets
#' @param lam A [Lambda()].
#' @param base A set.
#' @export
ImageSet <- function(lam, base) {
  base <- as_sym(base)
  if (!is_sym(lam) || lam$type != "Lambda") .perr("ValueError", "First argument must be a Lambda")
  if (!is_set(base)) .perr("TypeError", "Set arguments to ImageSet should of type Set")
  vs <- lam$args[[1L]]$args; expr <- lam$args[[2L]]
  if (length(vs) == 1L && identical(vs[[1L]]$key, expr$key)) return(base)
  if (!any(vapply(vs, function(v) .has1(expr, v), NA))) {
    ie <- .set_is_empty(base)
    if (isTRUE(ie)) return(S_EmptySet)
    if (identical(ie, FALSE)) return(FiniteSet(expr))
  }
  new_sym("ImageSet", list(lam, base))
}

# sympy.sets.sets.imageset for a Lambda and one set
.imageset <- function(lam, base) {
  set <- as_sym(base)
  if (!is_set(set)) .perr("ValueError", paste0("arguments after mapping should be sets, not ", set$type))
  f <- lam
  r <- tryCatch(.set_function(f, set), TypeError = function(e) NULL)
  if (is.null(r)) r <- ImageSet(f, set)
  else if (r$type == "EmptySet") return(r)
  if (r$type == "ImageSet") { f <- r$args[[1L]]; set <- r$args[[2L]] }
  vs <- f$args[[1L]]$args
  if (identical(vs[[1L]]$key, f$args[[2L]]$key)) return(set)
  if (set$type == "ImageSet") {
    inner <- set$args[[1L]]
    if (length(inner$args[[1L]]$args) == 1L && length(vs) == 1L)
      return(.imageset(Lambda(inner$args[[1L]]$args[[1L]], subs(f$args[[2L]], vs[[1L]], inner$args[[2L]])), set$args[[2L]]))
  }
  r
}

.lam_call <- function(f, v) apply_lambda(f, v)

.set_function <- function(f, x) {
  h <- .dispatch2(.SET_FUNCTION, f, x)
  if (is.null(h)) NULL else h(f, x)
}

# is_function_invertible_in_set for a Lambda
.fn_invertible <- function(f) {
  u <- Dummy("u")
  fd <- .diff1(.lam_call(f, u), u)
  if (.is_true(Gt(fd, 0)) || .is_true(Lt(fd, 0))) TRUE else NA
}

.sf_lambda_interval <- function(f, x) {
  expr <- f$args[[2L]]; vs <- f$args[[1L]]$args
  if (length(free_symbols(expr)) > 1L || length(vs) != 1L) return(NULL)
  var <- vs[[1L]]
  if (!isTRUE(ask_fact(var, "real")) && identical(ask_fact(subs(expr, var, Dummy("x", real = TRUE)), "real"), FALSE))
    return(NULL)
  if (expr$type == "Piecewise") {
    result <- S_EmptySet
    domain_set <- x
    for (pc in expr$args) {
      p_expr <- pc$args[[1L]]; p_cond <- pc$args[[2L]]
      intrvl <- if (.is_true(p_cond)) domain_set else Intersection(domain_set, .conv_cond_set(p_cond, var))
      image <- if (is_number_type(p_expr)) FiniteSet(p_expr) else .imageset(Lambda(var, p_expr), intrvl)
      result <- Union(result, image)
      domain_set <- Complement(domain_set, intrvl)
      if (domain_set$type == "EmptySet") break
    }
    return(result)
  }
  I <- .iv(x)
  if (!.is_comparable(I$a) || !.is_comparable(I$b)) return(NULL)
  sing <- tryCatch({
    s <- singularities(expr, var, x)
    if (s$type == "FiniteSet") s$args else if (s$type == "EmptySet") list() else stop("unsupported")
  }, error = function(e) NULL)
  if (is.null(sing)) return(NULL)
  if (length(sing) > 1L) sing <- .sort_set_elems(sing)
  in_list <- function(p, l) any(vapply(l, .py_eq, NA, b = p))
  .start <- NULL; .end <- NULL
  if (I$lo) .start <- limit(expr, var, I$a, "+") else if (!in_list(I$a, sing)) .start <- .lam_call(f, I$a)
  if (I$ro) .end <- limit(expr, var, I$b, "-") else if (!in_list(I$b, sing)) .end <- .lam_call(f, I$b)
  if (!length(sing)) {
    sol <- solveset(.diff1(expr, var), var)
    if (!(sol$type %in% c("FiniteSet", "EmptySet"))) return(NULL)
    solns <- if (sol$type == "FiniteSet") sol$args else list()
    inside <- Filter(function(i) isTRUE(ask_fact(i, "real")) && .py_in(x, i), solns)
    extr <- c(list(.start, .end), lapply(inside, function(i) .lam_call(f, i)))
    start <- fcall("Min", extr); end <- fcall("Max", extr)
    lo <- FALSE; ro <- FALSE
    if (.pb(Le(.start, .end))) {
      if (.py_eq(start, .start) && !in_list(start, solns)) lo <- I$lo
      if (.py_eq(end, .end) && !in_list(end, solns)) ro <- I$ro
    } else {
      if (.py_eq(start, .end) && !in_list(start, solns)) lo <- I$ro
      if (.py_eq(end, .start) && !in_list(end, solns)) ro <- I$lo
    }
    if (.py_eq(start, end) && is_inf_const(start)) return(FiniteSet(start))
    if (is_const(start, "-oo") && is_const(end, "oo") && (.has1(expr, S_Infinity) || .has1(expr, S_NegativeInfinity)))
      return(NULL)
    return(Interval(start, end, lo, ro))
  }
  k <- length(sing)
  mid <- if (k > 1L) do.call(Union, lapply(seq_len(k - 1L), function(i)
    .imageset(f, Interval(sing[[i]], sing[[i + 1L]], TRUE, TRUE)))) else S_EmptySet
  Union(Union(.imageset(f, Interval(I$a, sing[[1L]], I$lo, TRUE)), mid),
        .imageset(f, Interval(sing[[k]], I$b, TRUE, I$ro)))
}

.sf_range <- function(f, self) {
  if (.range_empty(self)) return(S_EmptySet)
  expr <- f$args[[2L]]
  if (!.is_expr_s(expr)) return(NULL)
  n <- .range_size(self)
  if (.py_eq(n, S_One)) return(FiniteSet(.lam_call(f, .range_item(self, 0L))))
  x <- f$args[[1L]]$args[[1L]]
  if (!.has1(expr, x) || .has1(.diff1(expr, x), x)) return(NULL)
  p <- .range_parts(self)
  F <- if (isTRUE(ask_fact(p$start, "finite"))) .lam_call(f, .add(list(.mul(list(p$step, x)), p$start)))
       else .lam_call(f, .add(list(.mul(list(.neg(p$step), x)), .range_last(self))))
  F <- expand_mul(F)
  if (!identical(F$key, expr$key)) return(.imageset(Lambda(x, F), Range(n)))
  NULL
}

# Wild-free match of a*n + b with a, b free of n: list(a, b) or NULL
.match_linear <- function(e, n) {
  if (!.has1(e, n)) return(list(S_Zero, e))
  term_a <- function(t) {
    if (identical(t$key, n$key)) return(S_One)
    if (t$type == "Mul") {
      k <- which(vapply(t$args, function(z) identical(z$key, n$key), NA))
      rest <- t$args[-k]
      if (length(k) == 1L && !any(vapply(rest, .has1, NA, tg = n))) return(.mul(rest))
    }
    NULL
  }
  if (e$type == "Add") {
    a <- list(); b <- list()
    for (t in e$args) {
      if (!.has1(t, n)) { b[[length(b) + 1L]] <- t; next }
      c <- term_a(t)
      if (is.null(c)) return(NULL)
      a[[length(a) + 1L]] <- c
    }
    return(list(.add(a), .add(b)))
  }
  c <- term_a(e)
  if (is.null(c)) NULL else list(c, S_Zero)
}

# Mod for comparable numbers, as Mod.eval's ratio branch
.num_mod <- function(p, q) {
  r <- .div(p, q)
  if (isTRUE(ask_fact(r, "integer"))) return(S_Zero)
  d <- trunc(num_value_dbl(N(r)))
  rv <- .sub(p, .mul(list(.int(d), q)))
  if (.is_true(Lt(.mul(list(rv, q)), 0))) rv <- .add(list(rv, q))
  rv
}

.sf_integers <- function(f, self) {
  expr <- f$args[[2L]]
  if (!.is_expr_s(expr)) return(NULL)
  n <- f$args[[1L]]$args[[1L]]
  if (identical(expr$key, fcall("Abs", list(n))$key)) return(Naturals0)
  c <- .lam_call(f, S_Zero)
  fx <- .sub(.lam_call(f, n), c)
  f_x <- .sub(.lam_call(f, .neg(n)), c)
  neg_count <- function(e) sum(vapply(if (e$type == "Add") e$args else list(e), could_extract_minus_sign, NA))
  if (neg_count(f_x) < neg_count(fx)) expr <- .add(list(f_x, c))
  mt <- .match_linear(expr, n)
  if (!is.null(mt) && !is_zero_num(mt[[1L]]) && !.has_float(mt[[1L]]) && !.has_float(mt[[2L]])) {
    a <- mt[[1L]]; b <- mt[[2L]]
    if (.py_eq(a, S_One) || .py_eq(a, S_NegativeOne))
      b <- .add(Filter(function(bi) !isTRUE(ask_fact(bi, "integer")), if (b$type == "Add") b$args else list(b)))
    if (.is_number_e(b) && isTRUE(ask_fact(a, "real"))) {
      rim <- .match_real_imag(b)
      if (!is.null(rim) && !is_zero_num(rim[[1L]]) && .is_comparable(rim[[1L]]) && .is_comparable(a))
        b <- .add(list(.num_mod(rim[[1L]], a), .mul(list(S_I, rim[[2L]]))))
    } else if (.is_number_e(b) && isTRUE(ask_fact(a, "imaginary"))) {
      rim <- .match_real_imag(b)
      ai <- .div(a, S_I)
      if (!is.null(rim) && !is_zero_num(rim[[2L]]) && .is_comparable(rim[[2L]]) && .is_comparable(ai))
        b <- .add(list(rim[[1L]], .mul(list(S_I, .num_mod(rim[[2L]], ai)))))
    }
    expr <- .add(list(.mul(list(a, n)), b))
  }
  if (!identical(expr$key, f$args[[2L]]$key)) return(ImageSet(Lambda(n, expr), Integers))
  NULL
}

# match_real_imag: list(re, im) or NULL
.match_real_imag <- function(e) {
  terms <- if (e$type == "Add") e$args else list(e)
  hasI <- vapply(terms, function(t) .has1(t, S_I), NA)
  r_ <- .add(terms[!hasI]); i_ <- .add(terms[hasI])
  if (is_zero_num(i_) && isTRUE(ask_fact(r_, "real"))) return(list(r_, i_))
  ic <- .div(i_, S_I)
  if (!is_zero_num(i_) && !.has1(ic, S_I) && isTRUE(ask_fact(ic, "real")) && isTRUE(ask_fact(r_, "real"))) return(list(r_, ic))
  NULL
}

.sf_naturals <- function(f, self) {
  expr <- f$args[[2L]]
  if (!.is_expr_s(expr)) return(NULL)
  x <- f$args[[1L]]$args[[1L]]
  if (length(setdiff(vapply(free_symbols(expr), .k, ""), x$key))) return(NULL)
  if (identical(expr$key, fcall("Abs", list(x))$key)) return(if (self$name == "Naturals") self else Naturals0)
  step <- coeff(expr, x)
  c <- subs(expr, x, S_Zero)
  if (c$type == "Integer" && step$type == "Integer" && identical(expr$key, .add(list(.mul(list(step, x)), c))$key)) {
    if (self$name == "Naturals") c <- .add(list(c, step))
    if (.pb(Gt(step, 0))) {
      if (.py_eq(step, S_One)) {
        if (is_zero_num(c)) return(Naturals0)
        if (.py_eq(c, S_One)) return(Naturals)
      }
      return(Range(c, S_Infinity, step))
    }
    return(Range(c, S_NegativeInfinity, step))
  }
  NULL
}

.sf_generic <- function(f, x) {
  X <- Symbol("x")
  ImageSet(Lambda(X, .lam_call(f, X)), x)
}

.SET_FUNCTION <- list(
  .reg("Lambda", "FiniteSet", function(f, x) FiniteSet(lapply(x$args, function(e) .lam_call(f, e)))),
  .reg("Lambda", "Interval", function(f, x) .sf_lambda_interval(f, x)),
  .reg("Lambda", "Union", function(f, x) do.call(Union, lapply(x$args, function(a) .imageset(f, a)))),
  .reg("Lambda", "Intersection", function(f, x)
    if (isTRUE(.fn_invertible(f))) do.call(Intersection, lapply(x$args, function(a) .imageset(f, a))) else .sf_generic(f, x)),
  .reg("Lambda", "EmptySet", function(f, x) x),
  .reg("Lambda", "Set", function(f, x) .sf_generic(f, x)),
  .reg("Lambda", "Range", function(f, x) .sf_range(f, x)),
  .reg("Lambda", "Integers", function(f, x) .sf_integers(f, x)),
  .reg("Lambda", "Naturals", function(f, x) .sf_naturals(f, x)),
  .reg("Lambda", "Reals", function(f, x)
    if (!.is_expr_s(f$args[[2L]])) NULL else .set_function(f, Interval(S_NegativeInfinity, S_Infinity))))

# ---------------------------------------------------------------------------
# Set arithmetic: _apply_operation, set_add, ... and the handlers in
# sympy/sets/handlers/add.py, mul.py and power.py

.handle_finite_sets_op <- function(op, x, y, commutative) {
  fx <- is_set(x) && x$type == "FiniteSet"; fy <- is_set(y) && y$type == "FiniteSet"
  if (fx && fy) return(FiniteSet(unlist(lapply(x$args, function(i) lapply(y$args, function(j) op(i, j))), recursive = FALSE)))
  if (fx || fy) {
    fs <- if (fx) x else y; other <- if (fx) y else x
    sets <- lapply(fs$args, function(i) .apply_operation(op, other, i, commutative))
    return(do.call(Union, lapply(sets, function(i) if (is_set(i)) i else FiniteSet(i))))
  }
  NULL
}

# ImageSet(lam, base).doit()
.imageset_doit <- function(lam, base) {
  r <- ImageSet(lam, base)
  if (r$type != "ImageSet") return(r)
  res <- .set_function(lam, base)
  if (is.null(res)) r else res
}

.apply_operation <- function(op, x, y, commutative) {
  d <- Dummy("d")
  out <- .handle_finite_sets_op(op, x, y, commutative)
  if (is.null(out)) out <- op(x, y)
  if (is.null(out) && commutative) out <- op(y, x)
  if (is.null(out)) {
    if (is_set(x) && !is_set(y)) out <- .imageset_doit(Lambda(d, op(d, y)), x)
    else if (!is_set(x) && is_set(y)) out <- .imageset_doit(Lambda(d, op(x, d)), y)
    else .perr("NotImplementedError", "ImageSet of two sets is not supported")
  }
  out
}

.op_dispatch <- function(table) function(x, y) {
  h <- .dispatch2(table, x, y)
  if (is.null(h)) NULL else h(x, y)
}

.SET_ADD <- list(
  .reg("Basic", "Basic", function(x, y) NULL),
  .reg("Expr", "Expr", function(x, y) .add(list(x, y))),
  .reg("Interval", "Interval", function(x, y) {
    A <- .iv(x); B <- .iv(y)
    Interval(.add(list(A$a, B$a)), .add(list(A$b, B$b)), A$lo || B$lo, A$ro || B$ro)
  }),
  .reg("Interval", "Infinity", function(x, y)
    if (is_const(.iv(x)$a, "-oo")) Interval(S_NegativeInfinity, S_Infinity) else FiniteSet(S_Infinity)),
  .reg("Interval", "NegativeInfinity", function(x, y)
    if (is_const(.iv(x)$b, "oo")) Interval(S_NegativeInfinity, S_Infinity) else FiniteSet(S_NegativeInfinity)))

.SET_SUB <- list(
  .reg("Basic", "Basic", function(x, y) NULL),
  .reg("Expr", "Expr", function(x, y) .sub(x, y)),
  .reg("Interval", "Interval", function(x, y) {
    A <- .iv(x); B <- .iv(y)
    Interval(.sub(A$a, B$b), .sub(A$b, B$a), A$lo || B$ro, A$ro || B$lo)
  }),
  .reg("Interval", "Infinity", function(x, y)
    if (is_const(.iv(x)$a, "-oo")) Interval(S_NegativeInfinity, S_Infinity) else FiniteSet(S_NegativeInfinity)),
  .reg("Interval", "NegativeInfinity", function(x, y)
    if (is_const(.iv(x)$a, "-oo")) Interval(S_NegativeInfinity, S_Infinity) else FiniteSet(S_NegativeInfinity)))

# Python min/max over (value, open) tuples
.tuple_lt <- function(p, q) {
  if (!.py_eq(p[[1L]], q[[1L]])) {
    if (is_nan(p[[1L]]) || is_nan(q[[1L]])) .perr("TypeError", "Invalid NaN comparison")
    return(.pb(Lt(p[[1L]], q[[1L]])))
  }
  p[[2L]] < q[[2L]]
}
.tuple_minmax <- function(vals, max = FALSE) {
  best <- vals[[1L]]
  for (v in vals[-1L]) if (if (max) .tuple_lt(best, v) else .tuple_lt(v, best)) best <- v
  best
}

.mul_interval_inf <- function(x, y) {
  I <- .iv(x)
  negative <- ask_fact(I$a, "extended_negative")
  positive <- ask_fact(I$b, "extended_positive")
  zero <- .contains(x, S_Zero)
  if (is.na(negative) || is.na(positive) || !(.is_true(zero) || .is_false(zero))) return(NULL)
  vals <- list()
  if (negative) vals[[length(vals) + 1L]] <- .neg(y)
  if (.is_true(zero)) vals[[length(vals) + 1L]] <- S_NaN
  if (positive) vals[[length(vals) + 1L]] <- y
  FiniteSet(vals)
}

.SET_MUL <- list(
  .reg("Basic", "Basic", function(x, y) NULL),
  .reg("Set", "Set", function(x, y) NULL),
  .reg("Expr", "Expr", function(x, y) .mul(list(x, y))),
  .reg("Interval", "Interval", function(x, y) {
    A <- .iv(x); B <- .iv(y)
    cv <- list(list(.mul(list(A$a, B$a)), A$lo || B$lo), list(.mul(list(A$a, B$b)), A$lo || B$ro),
               list(.mul(list(A$b, B$a)), A$ro || B$lo), list(.mul(list(A$b, B$b)), A$ro || B$ro))
    mn <- .tuple_minmax(cv); mx <- .tuple_minmax(cv, max = TRUE)
    Interval(mn[[1L]], mx[[1L]], mn[[2L]], mx[[2L]])
  }),
  .reg("Interval", "Infinity", function(x, y) .mul_interval_inf(x, y)),
  .reg("Interval", "NegativeInfinity", function(x, y) .mul_interval_inf(x, y)))

.SET_DIV <- list(
  .reg("Basic", "Basic", function(x, y) NULL),
  .reg("Expr", "Expr", function(x, y) .div(x, y)),
  .reg("Set", "Set", function(x, y) NULL),
  .reg("Interval", "Interval", function(x, y) {
    B <- .iv(y)
    if (isTRUE(ask_fact(.mul(list(B$a, B$b)), "negative"))) return(Interval(S_NegativeInfinity, S_Infinity))
    s2 <- if (.py_eq(B$a, S_Zero)) S_Infinity else .div(S_One, B$a)
    s1 <- if (.py_eq(B$b, S_Zero)) S_NegativeInfinity else .div(S_One, B$b)
    .set_mul(x, Interval(s1, s2, B$ro, B$lo))
  }))

.pow_interval_integer <- function(x, exponent) {
  I <- .iv(x)
  s1 <- .pow(I$a, exponent); s2 <- .pow(I$b, exponent)
  cond <- if (.pb(Gt(exponent, 0))) Gt(s2, s1) else Gt(I$b, .neg(I$a))
  if (.is_true(cond)) { lo <- I$lo; ro <- I$ro; sleft <- s2 }
  else { lo <- I$ro; ro <- I$lo; sleft <- s1 }
  if (isTRUE(ask_fact(I$a, "positive")) || isTRUE(ask_fact(I$b, "negative")))
    return(Interval(Min(s1, s2), Max(s1, s2), lo, ro))
  if (isTRUE(ask_fact(exponent, "odd"))) {
    if (isTRUE(ask_fact(exponent, "negative"))) {
      if (isTRUE(ask_fact(I$a, "zero"))) return(Interval(s2, S_Infinity, I$ro))
      if (isTRUE(ask_fact(I$b, "zero"))) return(Interval(S_NegativeInfinity, s1, TRUE, I$lo))
      return(Union(Interval(S_NegativeInfinity, s1, TRUE, I$lo), Interval(s2, S_Infinity, I$ro)))
    }
    return(Interval(s1, s2, I$lo, I$ro))
  }
  if (isTRUE(ask_fact(exponent, "even"))) {
    if (isTRUE(ask_fact(exponent, "negative"))) {
      if (isTRUE(ask_fact(I$a, "zero"))) return(Interval(s2, S_Infinity, I$ro))
      if (isTRUE(ask_fact(I$b, "zero"))) return(Interval(s1, S_Infinity, I$lo))
      return(Interval(0, S_Infinity))
    }
    return(Interval(S_Zero, sleft, !.py_in(x, S_Zero), lo))
  }
  NULL
}

.pow_interval_inf <- function(b, e) {
  I <- .iv(b)
  if (isTRUE(ask_fact(I$a, "nonnegative"))) {
    if (.pb(Lt(I$b, 1))) return(FiniteSet(S_Zero))
    if (.pb(Gt(I$a, 1))) return(FiniteSet(S_Infinity))
    return(Interval(0, S_Infinity))
  }
  if (isTRUE(ask_fact(I$b, "negative"))) {
    if (.pb(Gt(I$a, -1))) return(FiniteSet(S_Zero))
    if (.pb(Lt(I$b, -1))) return(FiniteSet(S_NegativeInfinity, S_Infinity))
    return(Interval(S_NegativeInfinity, S_Infinity))
  }
  if (.pb(Gt(I$a, -1))) {
    if (.pb(Lt(I$b, 1))) return(FiniteSet(S_Zero))
    return(Interval(0, S_Infinity))
  }
  Interval(S_NegativeInfinity, S_Infinity)
}

.SET_POW <- list(
  .reg("Basic", "Basic", function(x, y) NULL),
  .reg("Expr", "Expr", function(x, y) .pow(x, y)),
  .reg("Interval", "Zero", function(x, z) FiniteSet(S_One)),
  .reg("Interval", "Integer", function(x, e) .pow_interval_integer(x, e)),
  .reg("Interval", "Infinity", function(x, e) .pow_interval_inf(x, e)),
  .reg("Interval", "NegativeInfinity", function(x, e) .op_dispatch(.SET_POW)(.set_div(S_One, x), S_Infinity)))

.set_add <- function(x, y) .apply_operation(.op_dispatch(.SET_ADD), as_sym(x), as_sym(y), TRUE)
.set_sub <- function(x, y) .apply_operation(.op_dispatch(.SET_SUB), as_sym(x), as_sym(y), FALSE)
.set_mul <- function(x, y) .apply_operation(.op_dispatch(.SET_MUL), as_sym(x), as_sym(y), TRUE)
.set_div <- function(x, y) .apply_operation(.op_dispatch(.SET_DIV), as_sym(x), as_sym(y), FALSE)
.set_pow <- function(x, y) .apply_operation(.op_dispatch(.SET_POW), as_sym(x), as_sym(y), FALSE)

#' Set expressions
#'
#' `SetExpr(set)` is an expression whose value ranges over `set`, like
#' SymPy's `SetExpr`. The arithmetic operators `+`, `-`, `*`, `/` and `^`
#' apply interval arithmetic and act elementwise on finite sets; the
#' resulting set is in `$args[[1]]`.
#'
#' @param set A set.
#' @return A `sym` of type `SetExpr`.
#' @examples
#' SetExpr(Interval(0, 1)) + 1
#' SetExpr(Interval(1, 2)) * SetExpr(Interval(10, 20))
#' SetExpr(FiniteSet(1, 2)) + SetExpr(Interval(0, 10))
#' @export
SetExpr <- function(set) {
  s <- as_sym(set)
  if (!is_set(s)) .perr("TypeError", "SetExpr requires a Set")
  new_sym("SetExpr", list(s))
}
register_str_printer("SetExpr", function(e) paste0("SetExpr(", .str(e$args[[1L]]), ")"))
.is_setexpr <- function(e) is_sym(e) && e$type == "SetExpr"
.setexpr_ops <- function(gen, e1, e2) {
  op <- switch(gen, "+" = .set_add, "-" = .set_sub, "*" = .set_mul, "/" = .set_div, "^" = .set_pow, NULL)
  if (is.null(op)) return(NULL)
  x <- as_sym(e1); y <- as_sym(e2)
  if (.is_setexpr(x)) x <- x$args[[1L]]
  if (.is_setexpr(y)) y <- y$args[[1L]]
  SetExpr(op(x, y))
}

# ---------------------------------------------------------------------------
# Printers: a finite set that contains finite sets prints as
# FiniteSet(...); an interval prints its endpoints with str(); ConditionSet
# omits a UniversalSet base.

register_str_printer("FiniteSet", function(e) {
  els <- .key_sort(lapply(e$args, .sort_key), e$args)
  body <- paste(vapply(els, .str, ""), collapse = ", ")
  nested <- any(vapply(els, function(z) .has_type(z, "FiniteSet"), NA))
  if (nested) paste0("FiniteSet(", body, ")") else paste0("{", body, "}")
})

register_str_printer("Interval", function(e) {
  i <- .iv(e)
  base <- paste0(sstr(i$a), ", ", sstr(i$b))
  lo <- i$lo && !is_inf_const(i$a)
  ro <- i$ro && !is_inf_const(i$b)
  if (i$lo && i$ro && (lo || ro)) return(paste0("Interval.open(", base, ")"))
  if (lo) return(paste0("Interval.Lopen(", base, ")"))
  if (ro) return(paste0("Interval.Ropen(", base, ")"))
  paste0("Interval(", base, ")")
})

register_str_printer("ConditionSet", function(e) {
  b <- e$args[[3L]]
  if (b$type == "SpecialSet" && b$name == "UniversalSet")
    return(paste0("ConditionSet(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))
  paste0("ConditionSet(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", ", .str(b), ")")
})
