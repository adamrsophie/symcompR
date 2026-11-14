# And._eval_simplify of sympy/logic/boolalg.py: equalities in one symbol
# are solved and substituted into the other relationals, then the two and
# three term patterns of _simplify_patterns_and and _simplify_patterns_and3
# replace pairs and triples of relationals (_apply_patternbased_simplification).

.fs_rel_tuple <- function(...) new_sym("Tuple", list(...))

# Relational.reversed, built without evaluation
.fs_rel_rev <- function(r) {
  rt <- switch(r$type, StrictLessThan = "StrictGreaterThan", LessThan = "GreaterThan",
               StrictGreaterThan = "StrictLessThan", GreaterThan = "LessThan", r$type)
  new_sym(rt, list(r$args[[2L]], r$args[[1L]]))
}

.fs_and_patterns2 <- function() {
  a <- Wild("a"); b <- Wild("b"); c <- Wild("c")
  T2 <- .fs_rel_tuple
  list(
    list(T2(Eq(a, b), Lt(a, b)), S_false),
    list(T2(Lt(b, a), Lt(a, b)), S_false),
    list(T2(Eq(a, b), Le(b, a)), Eq(a, b)),
    list(T2(Le(b, a), Le(a, b)), Eq(a, b)),
    list(T2(Le(a, b), Lt(a, b)), Lt(a, b)),
    list(T2(Le(a, b), Ne(a, b)), Lt(a, b)),
    list(T2(Lt(a, b), Ne(a, b)), Lt(a, b)),
    list(T2(Eq(a, b), Eq(a, .neg(b))), And(Eq(a, S_Zero), Eq(b, S_Zero))),
    list(T2(Le(b, a), Le(c, a)), Ge(a, Max(b, c))),
    list(T2(Le(b, a), Lt(c, a)), ITE(Gt(b, c), Ge(a, b), Gt(a, c))),
    list(T2(Lt(b, a), Lt(c, a)), Gt(a, Max(b, c))),
    list(T2(Le(a, b), Le(a, c)), Le(a, Min(b, c))),
    list(T2(Le(a, b), Lt(a, c)), ITE(Lt(b, c), Le(a, b), Lt(a, c))),
    list(T2(Lt(a, b), Lt(a, c)), Lt(a, Min(b, c))),
    list(T2(Le(a, b), Le(c, a)), ITE(Eq(b, c), Eq(a, b), ITE(Lt(b, c), S_false, And(Le(a, b), Ge(a, c))))),
    list(T2(Le(c, a), Le(a, b)), ITE(Eq(b, c), Eq(a, b), ITE(Lt(b, c), S_false, And(Le(a, b), Ge(a, c))))),
    list(T2(Lt(a, b), Lt(c, a)), ITE(Lt(b, c), S_false, And(Lt(a, b), Gt(a, c)))),
    list(T2(Lt(c, a), Lt(a, b)), ITE(Lt(b, c), S_false, And(Lt(a, b), Gt(a, c)))),
    list(T2(Le(a, b), Lt(c, a)), ITE(Le(b, c), S_false, And(Le(a, b), Gt(a, c)))),
    list(T2(Le(c, a), Lt(a, b)), ITE(Le(b, c), S_false, And(Lt(a, b), Ge(a, c)))),
    list(T2(Eq(a, b), Eq(a, c)), ITE(Eq(b, c), Eq(a, b), S_false)),
    list(T2(Lt(a, b), Lt(.neg(b), a)), ITE(Gt(b, 0), Lt(fcall("Abs", list(a)), b), S_false)),
    list(T2(Le(a, b), Le(.neg(b), a)), ITE(Ge(b, 0), Le(fcall("Abs", list(a)), b), S_false))
  )
}

.fs_and_patterns3 <- function() {
  a <- Wild("a"); b <- Wild("b"); c <- Wild("c")
  T3 <- .fs_rel_tuple
  list(
    list(T3(Ge(a, b), Ge(b, c), Gt(c, a)), S_false),
    list(T3(Ge(a, b), Gt(b, c), Gt(c, a)), S_false),
    list(T3(Gt(a, b), Gt(b, c), Gt(c, a)), S_false),
    list(T3(Ge(a, b), Ge(a, c), Ge(b, c)), And(Ge(a, b), Ge(b, c))),
    list(T3(Ge(a, b), Ge(a, c), Gt(b, c)), And(Ge(a, b), Gt(b, c))),
    list(T3(Ge(a, b), Gt(a, c), Gt(b, c)), And(Ge(a, b), Gt(b, c))),
    list(T3(Ge(a, c), Gt(a, b), Gt(b, c)), And(Gt(a, b), Gt(b, c))),
    list(T3(Ge(b, c), Gt(a, b), Gt(a, c)), And(Gt(a, b), Ge(b, c))),
    list(T3(Gt(a, b), Gt(a, c), Gt(b, c)), And(Gt(a, b), Gt(b, c))),
    list(T3(Ge(b, a), Ge(c, a), Ge(b, c)), And(Ge(c, a), Ge(b, c))),
    list(T3(Ge(b, a), Ge(c, a), Gt(b, c)), And(Ge(c, a), Gt(b, c))),
    list(T3(Ge(b, a), Gt(c, a), Gt(b, c)), And(Gt(c, a), Gt(b, c))),
    list(T3(Ge(c, a), Gt(b, a), Gt(b, c)), And(Ge(c, a), Gt(b, c))),
    list(T3(Ge(b, c), Gt(b, a), Gt(c, a)), And(Gt(c, a), Ge(b, c))),
    list(T3(Gt(b, a), Gt(c, a), Gt(b, c)), And(Gt(c, a), Gt(b, c))),
    list(T3(Ge(a, b), Ge(b, c), Ge(c, a)), And(Eq(a, b), Eq(b, c)))
  )
}

.fs_is_ineq <- function(r) r$type %in% c("StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan")

# itertools.product of the candidate tuples
.fs_product <- function(lists) {
  out <- list(list())
  for (l in lists) out <- unlist(lapply(out, function(p) lapply(l, function(z) c(p, list(z)))), recursive = FALSE)
  out
}

# itertools.permutations(seq_len(n), 3) and combinations(seq_len(n), 2)
.fs_perm3 <- function(n) {
  out <- list()
  for (i in seq_len(n)) for (j in seq_len(n)) for (k in seq_len(n))
    if (i != j && j != k && i != k) out[[length(out) + 1L]] <- c(i, j, k)
  out
}
.fs_comb2 <- function(n) {
  out <- list()
  for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) out[[length(out) + 1L]] <- c(i, j)
  out
}

# One pass of _apply_patternbased_twoterm_simplification (nterm 2) or
# _apply_patternbased_threeterm_simplification (nterm 3)
.fs_pattern_terms <- function(Rel, patterns, func, dominating, replacement, measure, nterm) {
  flip <- if (nterm == 2L) c("GreaterThan", "StrictGreaterThan") else c("LessThan", "StrictLessThan")
  changed <- TRUE
  while (changed && length(Rel) >= nterm) {
    changed <- FALSE
    Rel <- lapply(Rel, function(r) if (r$type %in% flip) .fs_rel_rev(r) else r)
    Rel <- ordered(Rel)
    rtmp <- lapply(Rel, function(r) if (.fs_is_ineq(r)) list(r) else list(r, .fs_rel_rev(r)))
    results <- list()
    idxs <- if (nterm == 2L) .fs_comb2(length(Rel)) else .fs_perm3(length(Rel))
    for (ix in idxs) {
      for (pt in patterns) {
        res <- list()
        for (ps in .fs_product(rtmp[ix])) {
          old <- do.call(.fs_rel_tuple, ps)
          m <- .matches(pt[[1L]], old, list())
          if (length(m)) res[[length(res) + 1L]] <- list(m, old)
        }
        for (r in res) {
          rule <- unname(lapply(r[[1L]], function(q) list(q$w, q$v)))
          np <- xreplace(pt[[2L]], rule)
          if (identical(np$key, dominating$key)) return(list(replacement))
          if (!(np$type == "ITE") && !.has_function(np, c("Min", "Max"))) {
            saving <- measure(do.call(func, r[[2L]]$args)) - measure(np)
            if (saving > 0) results[[length(results) + 1L]] <- list(saving, ix, np)
          }
        }
      }
    }
    if (length(results)) {
      o <- order(-vapply(results, `[[`, 0, 1L), method = "radix")
      best <- results[[o[1L]]]
      newrel <- best[[3L]]
      Rel <- Rel[-best[[2L]]]
      if (!identical(newrel$key, Not(dominating)$key)) {
        if (newrel$type == "And") Rel <- c(Rel, newrel$args) else Rel[[length(Rel) + 1L]] <- newrel
      }
      changed <- TRUE
    }
  }
  Rel
}

# sympy.core.relational._canonical: every relational replaced by its canonical form
.fs_canonical <- function(cond) {
  if (is_relational(cond)) return(.rel_canonical(cond))
  if (!length(cond$args) || !(cond$type %in% .BOOL_TYPES)) return(cond)
  sym_rebuild(cond, lapply(cond$args, .fs_canonical))
}

.fs_apply_patternbased <- function(rv, patterns, measure, dominating, threeterm = NULL) {
  replacement <- dominating
  isrel <- vapply(rv$args, is_relational, NA)
  Rel <- rv$args[isrel]
  nonRel <- rv$args[!isrel]
  if (length(Rel) <= 1L) return(rv)
  realok <- vapply(Rel, function(i) !any(vapply(free_symbols(i), function(s) identical(ask_fact(s, "real"), FALSE), NA)), NA)
  nonRealRel <- Rel[!realok]
  Rel <- lapply(Rel[realok], .rel_canonical)
  func <- if (rv$type == "And") And else Or
  if (!is.null(threeterm) && length(Rel) >= 3L)
    Rel <- .fs_pattern_terms(Rel, threeterm, func, dominating, replacement, measure, 3L)
  Rel <- .fs_pattern_terms(Rel, patterns, func, dominating, replacement, measure, 2L)
  do.call(func, c(lapply(ordered(Rel), .fs_canonical), nonRel, nonRealRel))
}

# the equality substitution part of And._eval_simplify
.fs_and_eqs <- function(rv, kw, measure) {
  isrel <- vapply(rv$args, is_relational, NA)
  if (!any(isrel)) return(rv)
  Rel <- rv$args[isrel]
  nonRel <- rv$args[!isrel]
  iseq <- vapply(Rel, function(i) i$type == "Equality", NA)
  eqs0 <- Rel[iseq]
  other <- Rel[!iseq]
  if (!length(eqs0)) return(rv)
  reps <- list()
  rep_of <- function(x) reps[[x$key]]
  eqs <- list()
  nonlineqs <- list()
  sifted <- list()
  for (e in ordered(eqs0)) {
    k <- as.character(length(free_symbols(e)))
    sifted[[k]] <- c(sifted[[k]], list(e))
  }
  while (!is.null(sifted[["1"]])) {
    one <- sifted[["1"]]
    sifted[["1"]] <- NULL
    for (e in one) {
      x <- free_symbols(e)[[1L]]
      lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
      skip <- FALSE
      if ((!identical(lhs$key, x$key) || .has1(rhs, x)) && is.null(rep_of(x))) {
        mb <- tryCatch(linear_coeffs(new_sym("Add", list(lhs, .neg(rhs))), x), error = function(err) NULL)
        if (!is.null(mb)) {
          enew <- Eq(x, .div(.neg(mb[[2L]]), mb[[1L]]))
          if (measure(enew) <= kw$ratio * measure(e)) {
            e <- enew
            lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
          } else {
            eqs[[length(eqs) + 1L]] <- e
            skip <- TRUE
          }
        }
      }
      if (skip) next
      if (!is.null(rep_of(x))) {
        eqs[[length(eqs) + 1L]] <- subs(e, x, rep_of(x))
      } else if (is_relational(e) && identical(lhs$key, x$key) && !.has1(rhs, x)) {
        reps[[x$key]] <- list(x, rhs)
        eqs[[length(eqs) + 1L]] <- e
      } else {
        nonlineqs[[length(nonlineqs) + 1L]] <- e
      }
    }
    resifted <- list()
    for (k in names(sifted)) for (e in sifted[[k]]) {
      if (length(reps)) e <- xreplace(e, unname(reps))
      kk <- as.character(length(free_symbols(e)))
      resifted[[kk]] <- c(resifted[[kk]], list(e))
    }
    sifted <- resifted
  }
  for (k in names(sifted)) eqs <- c(eqs, sifted[[k]])
  sb <- function(z) if (length(reps)) subs(z, unname(reps)) else z
  nonlineqs <- lapply(nonlineqs, sb)
  other <- lapply(other, sb)
  canon <- function(i) if (is_relational(i)) .rel_canonical(i) else i
  do.call(And, c(lapply(c(eqs, nonlineqs, other), canon), nonRel))
}

.fs_and_simplify <- function(e, kw) {
  rv <- .boolean_simplify(e, kw)
  if (!is_sym(rv) || rv$type != "And") return(rv)
  measure <- function(z) .simp_measure(z, kw)
  rv <- .fs_and_eqs(rv, kw, measure)
  if (!is_sym(rv) || rv$type != "And") return(rv)
  if (is.null(.fs_memo$and2)) {
    .fs_memo$and2 <- .fs_and_patterns2()
    .fs_memo$and3 <- .fs_and_patterns3()
  }
  .fs_apply_patternbased(rv, .fs_memo$and2, measure, S_false, .fs_memo$and3)
}
