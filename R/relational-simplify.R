# Relational._eval_simplify, Relational.canonical and _canonical_coeff of
# sympy/core/relational.py, and BooleanFunction._eval_simplify of
# sympy/logic/boolalg.py.

.as_coeff_dict <- function(e) {
  out_k <- list(); out_v <- list()
  add1 <- function(m, c) {
    k <- m$key
    if (is.null(out_k[[k]])) { out_k[[k]] <<- m; out_v[[k]] <<- list(c) } else out_v[[k]][[length(out_v[[k]]) + 1L]] <<- c
  }
  if (e$type == "Add") {
    for (a in e$args) { sp <- as_coeff_mul_split(a); add1(sp[[2L]], sp[[1L]]) }
  } else {
    sp <- as_coeff_mul_split(e)
    if (!.is_rat(sp[[1L]])) add1(e, S_One) else add1(sp[[2L]], sp[[1L]])
  }
  list(k = out_k, v = lapply(out_v, function(v) if (length(v) == 1L) v[[1L]] else .add(v)))
}

.rel_reversedsign <- function(r) {
  a <- r$args[[1L]]; b <- r$args[[2L]]
  if (a$type == "BooleanAtom" || b$type == "BooleanAtom") return(r)
  op <- .rel_op(r)
  .rel_make(switch(op, "==" = "==", ">" = "<", ">=" = "<=", "<" = ">", "<=" = ">=", "!=" = "!="), .neg(a), .neg(b))
}

.is_number_expr <- function(z) .is_expr_node(z) && !length(free_symbols(z)) && !.has_type(z, "Symbol")

# Relational.canonical
.rel_canonical <- function(r) {
  args <- lapply(r$args, function(i) if (is_relational(i)) .rel_canonical(i) else i)
  if (!identical(vapply(args, .k, ""), vapply(r$args, .k, ""))) {
    r <- .rel_make(.rel_op(r), args[[1L]], args[[2L]])
    if (!is_relational(r)) return(r)
  }
  lhs <- r$args[[1L]]; rhs <- r$args[[2L]]
  if (.is_number_expr(rhs)) {
    if (is_number_type(rhs) && is_number_type(lhs) && .is_true(Gt(lhs, rhs))) r <- .rel_reversed(r)
  } else {
    ord <- ordered(list(lhs, rhs))
    if (.is_number_expr(lhs) || !identical(ord[[1L]]$key, lhs$key)) r <- .rel_reversed(r)
  }
  if (!is_relational(r)) return(r)
  lhs <- r$args[[1L]]; rhs <- r$args[[2L]]
  if (lhs$type == "BooleanAtom" || rhs$type == "BooleanAtom") return(r)
  if (.is_expr_node(lhs) && could_extract_minus_sign(lhs)) return(.rel_reversedsign(r))
  if (!.is_number_expr(rhs) && .is_expr_node(rhs) && could_extract_minus_sign(rhs)) {
    e1 <- ordered(list(lhs, .neg(rhs)))[[1L]]
    if (!identical(e1$key, lhs$key)) return(.rel_reversedsign(.rel_reversed(r)))
  }
  r
}

# _canonical_coeff
.canonical_coeff <- function(rel) {
  rel <- .rel_canonical(rel)
  if (!is_relational(rel) || rel$args[[2L]]$type %in% .BOOL_TYPES) return(rel)
  if (!.is_expr_node(rel$args[[1L]])) return(.rel_reversed(rel))
  ca <- .coeff_add_split(rel$args[[1L]]); b <- ca[[1L]]; l <- ca[[2L]]
  if (!.is_rat(b)) { b <- S_Zero; l <- rel$args[[1L]] }
  if (is_zero_num(l)) return(.rel_make(.rel_op(rel), rel$args[[1L]], rel$args[[2L]]))
  cm <- .coeff_mul_rat(l); m <- cm[[1L]]; lhs <- cm[[2L]]
  rhs <- .div(.sub(rel$args[[2L]], b), m)
  if (num_sign(m) < 0) return(.rel_make(.rel_rev_op(.rel_op(rel)), lhs, rhs))
  .rel_make(.rel_op(rel), lhs, rhs)
}

.may_be_indeterminate <- function(terms) {
  terms <- Filter(function(t) !isTRUE(ask_fact(t, "finite")), terms)
  if (length(terms) < 2L) return(FALSE)
  direction <- NULL
  for (t in terms) {
    sg <- if (isTRUE(ask_fact(t, "extended_nonnegative"))) 1L else if (isTRUE(ask_fact(t, "extended_nonpositive"))) -1L else return(TRUE)
    if (is.null(direction)) direction <- sg else if (sg != direction) return(TRUE)
  }
  FALSE
}

.relational_simplify <- function(self, kw) {
  measure <- function(z) .simp_measure(z, kw)
  op0 <- .rel_op(self)
  finish <- function(...) {
    ch <- list(...)
    res <- ch[[which.min(vapply(ch, measure, 0))]]
    if (measure(res) < kw$ratio * measure(self)) res else self
  }
  is_eqne <- function(op) op %in% c("==", "!=")
  simplify_mul <- function(rel) {
    if (is_zero_num(rel$args[[1L]]) && !is_zero_num(rel$args[[2L]])) rel <- .rel_reversed(rel)
    if (!is_relational(rel)) return(rel)
    L <- rel$args[[1L]]; R <- rel$args[[2L]]
    d <- Dummy("d")
    factored <- factor_terms(.sub(L, .mul(list(d, R))))
    if (.has1(factored, d)) {
      ci <- .as_indep(factored, list(d), FALSE); common <- ci[[1L]]; core <- ci[[2L]]
    } else {
      if (factored$type != "Mul") return(rel)
      common <- factored; core <- S_One
    }
    if (is_one_num(common) && common$type == "Integer") return(rel)
    keep <- list(); reverse <- FALSE; changed <- FALSE
    op <- .rel_op(rel)
    for (f in .margs(common)) {
      if (!.is_comm(f)) { keep[[length(keep) + 1L]] <- f; next }
      fin <- isTRUE(ask_fact(f, "finite"))
      safe <- if (is_eqne(op)) fin && isFALSE(ask_fact(f, "zero"))
              else fin && (isTRUE(ask_fact(f, "positive")) || isTRUE(ask_fact(f, "negative")))
      if (!safe) { keep[[length(keep) + 1L]] <- f; next }
      changed <- TRUE
      if (!is_eqne(op) && isTRUE(ask_fact(f, "negative"))) reverse <- !reverse
    }
    if (!changed) return(rel)
    keep <- .mul(keep)
    L <- .mul(list(keep, subs(core, d, S_Zero)))
    R <- .neg(.mul(list(keep, coeff(core, d))))
    newop <- if (reverse) .rel_rev_op(op) else op
    switch(newop, "<" = Lt(L, R, evaluate = FALSE), "<=" = Le(L, R, evaluate = FALSE), ">" = Gt(L, R, evaluate = FALSE),
           ">=" = Ge(L, R, evaluate = FALSE), "==" = Eq(L, R, evaluate = FALSE), "!=" = Ne(L, R, evaluate = FALSE))
  }
  unhollow <- function(z) {
    if (z$type == "Mul" && length(z$args) == 2L && .is_rat(z$args[[1L]]) && z$args[[2L]]$type == "Add")
      return(.add(lapply(z$args[[2L]]$args, function(t) .mul(list(z$args[[1L]], t)))))
    z
  }
  dummy_like <- function(z) {
    as <- list()
    for (nm in c("commutative", "finite", "real", "positive", "negative", "zero", "integer", "rational")) {
      v <- ask_fact(z, nm)
      if (!is.na(v)) as[[nm]] <- v
    }
    do.call(Dummy, c(list("d"), as))
  }
  additive_move_is_safe <- function(term, op) isTRUE(ask_fact(term, "finite")) && (is_eqne(op) || isTRUE(ask_fact(term, "real")))
  same_sign <- function(a, b) (isTRUE(ask_fact(b, "positive")) && isTRUE(ask_fact(a, "positive"))) ||
    (isTRUE(ask_fact(b, "negative")) && isTRUE(ask_fact(a, "negative")))
  r <- .rel_make(op0, .simplify(self$args[[1L]], kw), .simplify(self$args[[2L]], kw))
  if (!is_relational(r)) return(r)
  if (any(vapply(r$args, function(a) inherits(a, "symmatrix") || a$type %in% c("ImmutableMatrix", "MatrixSymbol"), TRUE))) return(r)
  if (!all(vapply(r$args, .is_expr_node, TRUE))) return(r)
  r <- simplify_mul(r)
  if (!is_relational(r)) return(r)
  base <- .rel_canonical(r)
  op <- .rel_op(r)
  pmask <- list(); prestore <- list()
  mask_power <- function(p) {
    if (!is.null(pmask[[p$key]])) return(pmask[[p$key]])
    d <- dummy_like(p); pmask[[p$key]] <<- d; prestore[[length(prestore) + 1L]] <<- list(d, p)
    d
  }
  hide <- function(e) .map_bottom_up(e, function(p) if (p$type == "Pow" && p$args[[1L]]$type == "Add") mask_power(p) else p)
  restore <- function(e) if (length(prestore)) xreplace(e, prestore) else e
  L <- expand_mul(hide(r$args[[1L]])); R <- expand_mul(hide(r$args[[2L]]))
  terms <- list()
  for (pr in list(list(L, FALSE), list(R, TRUE))) for (t in .aargs(pr[[1L]])) {
    t <- restore(t)
    if (isTRUE(ask_fact(t, "finite"))) next
    terms[[length(terms) + 1L]] <- if (pr[[2L]]) .neg(t) else t
  }
  track <- if (length(terms) < 2L) TRUE else .may_be_indeterminate(terms)
  ld <- .as_coeff_dict(L); rd <- .as_coeff_dict(R)
  shared <- list(); pinned <- list(); reduced_shared <- FALSE
  put <- function(dd, k, v) { dd$k[[k$key]] <- k; dd$v[[k$key]] <- v; dd }
  pop <- function(dd, key) { dd$k[[key]] <- NULL; dd$v[[key]] <- NULL; dd }
  for (key in unique(c(names(ld$k), names(rd$k)))) {
    kobj <- if (!is.null(ld$k[[key]])) ld$k[[key]] else rd$k[[key]]
    if (is_one_num(kobj) && kobj$type == "Integer") next
    original <- restore(kobj)
    if (additive_move_is_safe(original, op) || (!isTRUE(ask_fact(original, "finite")) && !track)) next
    in_l <- !is.null(ld$k[[key]]); in_r <- !is.null(rd$k[[key]])
    if (in_l && in_r) {
      lc <- ld$v[[key]]; rc <- rd$v[[key]]; ld <- pop(ld, key); rd <- pop(rd, key)
      if (!(.is_rat(lc) && .is_rat(rc))) return(finish(base))
      sl <- if (num_sign(lc) > 0) S_One else S_NegativeOne
      sr <- if (num_sign(rc) > 0) S_One else S_NegativeOne
      dl <- dummy_like(original); dr <- dummy_like(original); d <- dummy_like(original)
      ld <- put(ld, dl, sl); rd <- put(rd, dr, sr)
      if (!identical(lc$key, sl$key)) { ld <- put(ld, d, .sub(lc, sl)); reduced_shared <- TRUE }
      if (!identical(rc$key, sr$key)) { rd <- put(rd, d, .sub(rc, sr)); reduced_shared <- TRUE }
      shared[[length(shared) + 1L]] <- list(original, dl, dr, d)
    } else if (in_l) {
      c <- ld$v[[key]]; ld <- pop(ld, key)
      d <- dummy_like(original); ld <- put(ld, d, c)
      pinned[[length(pinned) + 1L]] <- list(d, original, TRUE)
    } else {
      c <- rd$v[[key]]; rd <- pop(rd, key)
      d <- dummy_like(original); rd <- put(rd, d, c)
      pinned[[length(pinned) + 1L]] <- list(d, original, FALSE)
    }
  }
  L <- restore(.add(lapply(names(ld$k), function(k) .mul(list(ld$v[[k]], ld$k[[k]])))))
  R <- restore(.add(lapply(names(rd$k), function(k) .mul(list(rd$v[[k]], rd$k[[k]])))))
  dif <- .simplify(.sub(L, R), kw)
  zrel <- simplify_mul(.rel_make(op, dif, S_Zero))
  if (!is_relational(zrel)) return(finish(base, zrel))
  op <- .rel_op(zrel)
  dif <- expand_mul(.sub(zrel$args[[1L]], zrel$args[[2L]]))
  trackers <- base::c(lapply(pinned, `[[`, 1L), do.call(base::c, lapply(shared, function(s) s[2:4])))
  has_tr <- function(z) any(vapply(trackers, function(t) .has1(z, t), TRUE))
  left <- list(); right <- list(); core <- dif
  for (p in pinned) {
    d <- p[[1L]]; original <- p[[2L]]
    c <- coeff(core, d)
    if (has_tr(c)) return(finish(base))
    core <- .sub(core, .mul(list(c, d)))
    if (.has1(core, d)) return(finish(base))
    restored <- .mul(list(c, original))
    if (p[[3L]]) left[[length(left) + 1L]] <- restored else right[[length(right) + 1L]] <- .neg(restored)
  }
  for (sh in shared) {
    original <- sh[[1L]]; dl <- sh[[2L]]; dr <- sh[[3L]]; d <- sh[[4L]]
    cl <- coeff(core, dl); cr <- coeff(core, dr); cd <- coeff(core, d)
    if (any(vapply(list(cl, cr, cd), has_tr, TRUE))) return(finish(base))
    core <- .sub(core, .add(list(.mul(list(cl, dl)), .mul(list(cr, dr)), .mul(list(cd, d)))))
    if (.has1(core, dl) || .has1(core, dr) || .has1(core, d)) return(finish(base))
    al <- cl; ar <- .neg(cr)
    choices <- list()
    if (same_sign(.add(list(al, cd)), al)) choices[[length(choices) + 1L]] <- list(.add(list(al, cd)), ar)
    if (same_sign(.sub(ar, cd), ar)) choices[[length(choices) + 1L]] <- list(al, .sub(ar, cd))
    if (!length(choices)) return(finish(base))
    ms <- vapply(choices, function(c) measure(.mul(list(c[[1L]], original))) + measure(.mul(list(c[[2L]], original))), 0)
    ch <- choices[[which.min(ms)]]
    left[[length(left) + 1L]] <- .mul(list(ch[[1L]], original))
    right[[length(right) + 1L]] <- .mul(list(ch[[2L]], original))
  }
  for (term in .aargs(expand_mul(core))) {
    c <- as_coeff_mul_split(term)[[1L]]
    if (isTRUE(ask_fact(c, "negative"))) right[[length(right) + 1L]] <- .neg(term) else left[[length(left) + 1L]] <- term
  }
  lhs <- gcd_terms(.add(left), clear = FALSE); rhs <- gcd_terms(.add(right), clear = FALSE)
  derived <- simplify_mul(switch(op, "<" = Lt(lhs, rhs, evaluate = FALSE), "<=" = Le(lhs, rhs, evaluate = FALSE),
                                 ">" = Gt(lhs, rhs, evaluate = FALSE), ">=" = Ge(lhs, rhs, evaluate = FALSE),
                                 "==" = Eq(lhs, rhs, evaluate = FALSE), "!=" = Ne(lhs, rhs, evaluate = FALSE)))
  if (.contains_key(derived, S_NaN$key)) return(finish(base))
  if (is_relational(derived)) derived <- .canonical_coeff(derived)
  if (is_relational(derived)) {
    a2 <- lapply(derived$args, unhollow)
    derived <- .rel_canonical(switch(.rel_op(derived), "<" = Lt(a2[[1L]], a2[[2L]], evaluate = FALSE),
      "<=" = Le(a2[[1L]], a2[[2L]], evaluate = FALSE), ">" = Gt(a2[[1L]], a2[[2L]], evaluate = FALSE),
      ">=" = Ge(a2[[1L]], a2[[2L]], evaluate = FALSE), "==" = Eq(a2[[1L]], a2[[2L]], evaluate = FALSE),
      "!=" = Ne(a2[[1L]], a2[[2L]], evaluate = FALSE)))
  }
  if (!length(trackers)) {
    v <- NULL
    if (!length(free_symbols(dif)) && .is_expr_node(dif)) {
      sg <- tryCatch(.mp_num_sign(dif), error = function(err) NA)
      if (!is.na(sg)) v <- .int(sg)
    } else if (isTRUE(tryCatch(.expr_equals(dif, S_Zero), error = function(err) FALSE))) v <- S_Zero
    if (!is.null(v)) derived <- .rel_make(op, v, S_Zero)
  }
  if (reduced_shared || !track) return(finish(derived, base))
  finish(base, derived)
}

.boolean_simplify <- function(e, kw) {
  rv <- simplify_univariate(e)
  if (!(rv$type %in% setdiff(.BOOL_TYPES, "BooleanAtom"))) return(.simplify(rv, kw))
  rv <- sym_rebuild(rv, lapply(rv$args, function(a) .simplify(a, kw)))
  simplify_logic(rv)
}
