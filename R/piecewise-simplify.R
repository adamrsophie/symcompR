# piecewise_simplify and Piecewise._intervals after
# sympy/functions/elementary/piecewise.py, simplify_univariate of
# sympy/logic/boolalg.py, and _solve_inequality of
# sympy/solvers/inequalities.py.

.rel_op <- function(r) switch(r$type, StrictLessThan = "<", LessThan = "<=", StrictGreaterThan = ">",
                              GreaterThan = ">=", Equality = "==", Unequality = "!=")
.rel_make <- function(op, a, b) switch(op, "<" = Lt(a, b), "<=" = Le(a, b), ">" = Gt(a, b), ">=" = Ge(a, b),
                                       "==" = Eq(a, b), "!=" = Ne(a, b))
.rel_rev_op <- function(op) switch(op, "<" = ">", "<=" = ">=", ">" = "<", ">=" = "<=", "==" = "==", "!=" = "!=")
# Relational.reversed: the sides swapped, the operator mirrored
.rel_reversed <- function(r) .rel_make(.rel_rev_op(.rel_op(r)), r$args[[2L]], r$args[[1L]])
.rel_lts <- function(r) if (r$type %in% c("GreaterThan", "StrictGreaterThan")) r$args[[2L]] else r$args[[1L]]
.rel_gts <- function(r) if (r$type %in% c("GreaterThan", "StrictGreaterThan")) r$args[[1L]] else r$args[[2L]]
.is_true <- function(e) is_sym(e) && identical(e$key, S_true$key)
.is_false <- function(e) is_sym(e) && identical(e$key, S_false$key)

# as_independent(*syms, as_Add = TRUE / FALSE)
.as_indep <- function(e, syms, as_add = TRUE) {
  parts <- if (as_add) .aargs(e) else .margs(e)
  dep <- vapply(parts, function(t) any(vapply(syms, function(s) .has1(t, s), TRUE)), TRUE)
  op <- if (as_add) .add else .mul
  list(op(parts[!dep]), op(parts[dep]))
}

#' Denominators of an expression
#'
#' `denoms(e)` returns the factors of the denominators of all
#' subexpressions of `e`, as SymPy's `denoms()`; with `symbols` only those
#' that contain one of them.
#'
#' @param e An expression.
#' @param ... Symbols.
#' @return A list of expressions.
#' @examples
#' x <- Symbol("x")
#' denoms(1/x + 1/(x + 1))
#' @export
denoms <- function(e, ...) {
  e <- as_sym(e)
  syms <- lapply(list(...), as_sym)
  out <- list()
  walk <- function(p) {
    if (.is_expr_node(p)) {
      den <- denom(p)
      if (!(is_one_num(den) && den$type == "Integer")) for (f in .margs(den)) out[[f$key]] <<- f
    }
    for (a in p$args) walk(a)
  }
  walk(e)
  out <- unname(out)
  if (!length(syms)) return(out)
  Filter(function(d) any(vapply(syms, function(s) .has1(d, s), TRUE)), out)
}

# _solve_inequality: a relational in s solved for s
.solve_inequality <- function(ie, s, linear = FALSE) {
  if (!.has1(ie, s)) return(ie)
  if (identical(ie$args[[2L]]$key, s$key)) ie <- .rel_reversed(ie)
  if (identical(ie$args[[1L]]$key, s$key) && !.has1(ie$args[[2L]], s)) return(ie)
  bden <- c(denoms(ie$args[[1L]]), denoms(ie$args[[2L]]))
  bden <- bden[!duplicated(vapply(bden, .k, ""))]
  classify <- function(r, i) {
    v <- tryCatch(subs(r, s, i), error = function(err) S_NaN)
    if (is_nan(v)) return(v)
    if (!(.is_true(v) || .is_false(v))) return(NULL)
    v
  }
  istrue <- function(v) !is.null(v) && .is_true(v)
  isfalse <- function(v) !is.null(v) && .is_false(v)
  rv <- NULL
  oo <- S_Infinity
  op <- .rel_op(ie)
  expr <- .sub(ie$args[[1L]], ie$args[[2L]])
  polyok <- isTRUE(is_polynomial(expr, s))
  deg <- if (polyok) degree(Poly(expr, s)) else NA
  if (polyok && deg == 0) rv <- .rel_make(op, expr, S_Zero)
  else if (!polyok || (!linear && deg > 1)) {
    if (linear) stop("NotImplementedError")
    rv <- tryCatch(reduce_rational_inequalities(list(list(list(expr, op))), s),
                   error = function(err) solve_univariate_inequality(ie, s))
    okoo <- classify(ie, oo)
    if (istrue(okoo) && isfalse(classify(rv, oo))) rv <- subs(rv, Lt(s, oo), S_true)
    oknoo <- classify(ie, .neg(oo))
    if (istrue(oknoo) && isfalse(classify(rv, .neg(oo)))) {
      rv <- subs(rv, Lt(.neg(oo), s), S_true)
      rv <- subs(rv, Gt(s, .neg(oo)), S_true)
    }
    if (.is_true(rv)) {
      rv <- if (istrue(okoo)) Le(s, oo) else Lt(s, oo)
      if (!istrue(oknoo)) rv <- And(Lt(.neg(oo), s), rv)
    }
  }
  conds <- list()
  if (is.null(rv)) {
    e <- expand(expr); rhs <- S_Zero
    bi <- .as_indep(e, list(s), TRUE); b <- bi[[1L]]
    e <- .sub(e, b); rhs <- .sub(rhs, b)
    ef <- factor_terms(e)
    ai <- .as_indep(ef, list(s), FALSE); a <- ai[[1L]]; e <- ai[[2L]]
    az <- ask_fact(a, "zero"); an <- ask_fact(a, "negative"); ap <- ask_fact(a, "positive")
    if (!isFALSE(az) || (is.na(an) && is.na(ap) && !(op %in% c("!=", "==")))) { e <- ef; a <- S_One }
    rhs <- .div(rhs, a)
    rv <- if (isTRUE(ask_fact(a, "positive"))) .rel_make(op, e, rhs) else .rel_make(.rel_rev_op(op), e, rhs)
    for (i in list(.neg(oo), oo)) {
      if (istrue(classify(rv, i)) && !istrue(classify(ie, i)))
        conds[[length(conds) + 1L]] <- if (is_const(i, "oo")) Lt(s, i) else Lt(i, s)
    }
  }
  cden <- vapply(denoms(rv), .k, "")
  for (d in bden) {
    if (d$key %in% cden) next
    c <- .solve_inequality(Eq(d, S_Zero), s, linear = linear)
    if (c$type == "Equality" && identical(c$args[[1L]]$key, s$key)) {
      if (istrue(classify(rv, c$args[[2L]]))) conds[[length(conds) + 1L]] <- Not(c)
    } else conds[[length(conds) + 1L]] <- Ne(d, S_Zero)
  }
  conds[[length(conds) + 1L]] <- rv
  do.call(And, conds)
}

#' @rdname piecewise_simplify
#' @export
distribute_or_over_and <- function(expr) .distribute(as_sym(expr), "Or", "And")

# Piecewise._intervals: list(ok, intervals) with intervals a list of
# list(lower, upper, expr, index), or list(FALSE, message)
.pw_intervals <- function(pw, sym, err_on_Eq = FALSE) {
  nonsymfail <- function(cond) list(FALSE, paste("A condition not involving", sstr(sym), "appeared:", sstr(cond)))
  solve_rel <- function(r) {
    if (!.has1(r, sym)) return(nonsymfail(r))
    rv <- tryCatch(.solve_inequality(r, sym), error = function(err) NULL)
    if (is.null(rv)) return(list(FALSE, paste("Unable to solve relational", sstr(r))))
    if (is_relational(rv)) {
      if (!identical(rv$args[[1L]]$key, sym$key) || .has1(rv$args[[2L]], sym))
        return(list(FALSE, paste("Unable to solve relational", sstr(r))))
      if (rv$type == "Equality") rv <- S_false
      else if (rv$type == "Unequality") rv <- Or(Lt(sym, rv$args[[2L]]), Gt(sym, rv$args[[2L]]))
    } else if (identical(rv$key, And(Lt(.neg(S_Infinity), sym), Lt(sym, S_Infinity))$key)) rv <- S_true
    list(TRUE, rv)
  }
  rels <- Filter(is_relational, unique(atoms(pw, .REL_TYPES)))
  reps <- list()
  for (r in rels) {
    o <- solve_rel(r)
    if (!isTRUE(o[[1L]])) return(list(FALSE, o[[2L]]))
    reps[[length(reps) + 1L]] <- list(r, o[[2L]])
  }
  args <- lapply(pw$args, function(i) if (length(reps)) xreplace(i, reps) else i)
  expr_cond <- list(); default <- NULL; idefault <- NULL
  for (i in seq_along(args)) {
    expr <- args[[i]]$args[[1L]]; cond <- args[[i]]$args[[2L]]
    if (.is_false(cond)) next
    if (.is_true(cond)) { default <- expr; idefault <- i; break }
    if (cond$type == "Equality") {
      if (err_on_Eq) return(list(FALSE, paste("encountered Eq condition:", sstr(cond))))
      next
    }
    cond <- to_cnf(cond)
    if (cond$type == "And") cond <- distribute_or_over_and(cond)
    if (cond$type == "Or") {
      for (o in cond$args) if (o$type != "Equality") expr_cond[[length(expr_cond) + 1L]] <- list(i, expr, o)
    } else if (!.is_false(cond)) expr_cond[[length(expr_cond) + 1L]] <- list(i, expr, cond)
  }
  int_expr <- list()
  k <- 1L
  while (k <= length(expr_cond)) {
    iarg <- expr_cond[[k]][[1L]]; expr <- expr_cond[[k]][[2L]]; cond <- expr_cond[[k]][[3L]]
    k <- k + 1L
    if (cond$type == "And") {
      lower <- .neg(S_Infinity); upper <- S_Infinity; exclude <- list()
      for (cond2 in cond$args) {
        if (!is_relational(cond2)) return(list(FALSE, "expecting only Relationals"))
        if (cond2$type == "Equality") {
          lower <- upper
          if (err_on_Eq) return(list(FALSE, "encountered secondary Eq condition"))
          break
        } else if (cond2$type == "Unequality") {
          l <- cond2$args[[1L]]; r <- cond2$args[[2L]]
          if (identical(l$key, sym$key)) exclude[[length(exclude) + 1L]] <- r
          else if (identical(r$key, sym$key)) exclude[[length(exclude) + 1L]] <- l
          else return(nonsymfail(cond2))
          next
        } else if (identical(.rel_lts(cond2)$key, sym$key)) upper <- Min(.rel_gts(cond2), upper)
        else if (identical(.rel_gts(cond2)$key, sym$key)) lower <- Max(.rel_lts(cond2), lower)
        else return(nonsymfail(cond2))
      }
      if (length(exclude)) {
        exclude <- ordered(exclude)
        newcond <- list()
        for (e in exclude) {
          if (.is_true(Lt(e, lower)) || .is_true(Gt(e, upper))) next
          if (!length(newcond)) newcond[[1L]] <- list(NULL, lower)
          newcond[[length(newcond) + 1L]] <- list(newcond[[length(newcond)]][[2L]], e)
        }
        if (!length(newcond)) newcond[[1L]] <- list(NULL, lower)
        newcond[[length(newcond) + 1L]] <- list(newcond[[length(newcond)]][[2L]], upper)
        newcond <- newcond[-1L]
        for (nc in newcond) expr_cond[[length(expr_cond) + 1L]] <- list(iarg, expr, And(Lt(nc[[1L]], sym), Lt(sym, nc[[2L]])))
        next
      }
    } else if (is_relational(cond) && cond$type != "Unequality") {
      lower <- .rel_lts(cond); upper <- .rel_gts(cond)
      if (identical(.rel_lts(cond)$key, sym$key)) lower <- .neg(S_Infinity)
      else if (identical(.rel_gts(cond)$key, sym$key)) upper <- S_Infinity
      else return(nonsymfail(cond))
    } else return(list(FALSE, paste("unrecognized condition:", sstr(cond))))
    upper <- Max(lower, upper)
    if (err_on_Eq && identical(lower$key, upper$key)) return(list(FALSE, "encountered Eq condition"))
    if (!.is_true(Ge(lower, upper))) int_expr[[length(int_expr) + 1L]] <- list(lower, upper, expr, iarg)
  }
  if (!is.null(default)) int_expr[[length(int_expr) + 1L]] <- list(.neg(S_Infinity), S_Infinity, default, idefault)
  ks <- vapply(int_expr, function(t) paste(t[[1L]]$key, t[[2L]]$key, t[[3L]]$key, t[[4L]], sep = "|"), "")
  list(TRUE, int_expr[!duplicated(ks)])
}

# _blessed: an Eq that can be substituted
.blessed <- function(e) {
  l <- e$args[[1L]]; r <- e$args[[2L]]
  dw <- function(z) z$type %in% c("Symbol", "Dummy") || (z$type == "Function" && !(z$name %in% names(.FUNS)))
  dw(l) && (dw(r) || r$type %in% c("Integer", "Rational") || (r$type == "Constant" && r$name %in% c("pi", "E", "EulerGamma", "Catalan", "GoldenRatio", "TribonacciConstant")))
}

.pw_pairs <- function(pw) lapply(pw$args, function(p) list(p$args[[1L]], p$args[[2L]]))
.pw_build <- function(pairs) do.call(Piecewise, lapply(pairs, function(p) list(p[[1L]], p[[2L]])))

.pw_simplify_eq_and <- function(args) {
  for (i in seq_along(args)) {
    expr <- args[[i]][[1L]]; cond <- args[[i]][[2L]]
    if (cond$type == "And") {
      isq <- vapply(cond$args, function(a) a$type == "Equality", TRUE)
      eqs <- cond$args[isq]; other <- cond$args[!isq]
    } else if (cond$type == "Equality") { eqs <- list(cond); other <- list() }
    else { eqs <- list(); other <- list() }
    if (length(eqs)) {
      eqs <- ordered(eqs)
      for (j in seq_along(eqs)) {
        e <- eqs[[j]]
        if (e$type == "Equality" && .blessed(e)) {
          expr <- subs(expr, e$args[[1L]], e$args[[2L]])
          if (j < length(eqs)) eqs[(j + 1L):length(eqs)] <- lapply(eqs[(j + 1L):length(eqs)], function(ei) subs(ei, e$args[[1L]], e$args[[2L]]))
          other <- lapply(other, function(ei) subs(ei, e$args[[1L]], e$args[[2L]]))
        }
      }
      cond <- do.call(And, c(eqs, other))
      args[[i]] <- list(expr, cond)
    }
  }
  args
}

.pw_simplify_equal_to_next <- function(args) {
  prevexpr <- NULL
  for (i in rev(seq_along(args))) {
    expr <- args[[i]][[1L]]; cond <- args[[i]][[2L]]
    if (!is.null(prevexpr)) {
      if (cond$type == "And") {
        isq <- vapply(cond$args, function(a) a$type == "Equality", TRUE)
        eqs <- cond$args[isq]; other <- cond$args[!isq]
      } else if (cond$type == "Equality") { eqs <- list(cond); other <- list() }
      else { eqs <- list(); other <- list() }
      pe <- prevexpr; ee <- expr
      if (length(eqs) && !length(other)) {
        for (e in ordered(eqs)) {
          if (length(args) == 2L || .blessed(e)) {
            pe <- subs(pe, e$args[[1L]], e$args[[2L]]); ee <- subs(ee, e$args[[1L]], e$args[[2L]])
          }
        }
      }
      if (identical(pe$key, ee$key)) args[[i]] <- list(args[[i + 1L]][[1L]], cond)
      else prevexpr <- expr
    } else prevexpr <- expr
  }
  args
}

#' Simplify a Piecewise expression
#'
#' `piecewise_simplify(expr)` simplifies the expressions and conditions of a
#' Piecewise, rewrites univariate conditions as disjoint intervals, uses
#' equality conditions to simplify their expressions, and merges segments
#' that agree, as SymPy's `piecewise_simplify()`. `simplify_univariate()`
#' writes a boolean in one variable as a union of intervals.
#' `distribute_or_over_and()` puts a boolean into disjunctive form.
#'
#' @param expr A Piecewise or boolean expression.
#' @param deep Unused; kept for SymPy's signature.
#' @param kw Options of `simplify()` (internal).
#' @param ... Ignored.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' piecewise_simplify(Piecewise(list(x^2, x < 0), list(x^2, x > 0), list(0, TRUE)))
#' simplify_univariate((x > 0) & (x < 2) | (x > 1))
#' @export
piecewise_simplify <- function(expr, deep = TRUE, kw = .simp_kw(), ...) {
  expr <- .pw_simplify_arguments(as_sym(expr), kw)
  if (expr$type != "Piecewise") return(expr)
  args <- .pw_pairs(expr)
  args <- .pw_simplify_eq_and(args)
  args <- .pw_simplify_equal_to_next(args)
  .pw_build(args)
}

.pw_simplify_arguments <- function(expr, kw) {
  f1 <- free_symbols(expr$args[[1L]]$args[[2L]])
  args <- NULL
  if (length(f1) == 1L && !length(atoms(expr, "Equality"))) {
    x <- f1[[1L]]
    iv <- .pw_intervals(expr, x, err_on_Eq = TRUE)
    include <- function(c, a) tryCatch(.is_true(subs(c, x, a)), error = function(err) FALSE)
    if (isTRUE(iv[[1L]])) {
      args <- list()
      covered <- S_EmptySet
      for (t in iv[[2L]]) {
        a <- t[[1L]]; b <- t[[2L]]; e <- t[[3L]]; i <- t[[4L]]
        c <- expr$args[[i]]$args[[2L]]
        incl_a <- include(c, a); incl_b <- include(c, b)
        ivl <- Interval(a, b, !incl_a, !incl_b)
        cset <- Complement(ivl, covered)
        if (identical(cset$key, S_EmptySet$key)) next
        a2 <- tryCatch(set_inf(cset), error = function(err) NULL)
        if (!is.null(a2)) { a <- a2; incl_a <- include(c, a) }
        inc <- function(v) isTRUE(tryCatch(contains(covered, v), error = function(err) FALSE)) || .is_true(tryCatch(contains(covered, v), error = function(err) S_false))
        ainf <- isTRUE(ask_fact(a, "infinite")); binf <- isTRUE(ask_fact(b, "infinite"))
        if (incl_a && incl_b) {
          c <- if (ainf && binf) S_true else if (binf) (if (inc(a)) Gt(x, a) else Ge(x, a)) else if (ainf) Le(x, b)
               else if (inc(a)) And(Lt(a, x), Le(x, b)) else And(Le(a, x), Le(x, b))
        } else if (incl_a) {
          c <- if (ainf) Lt(x, b) else if (inc(a)) And(Lt(a, x), Lt(x, b)) else And(Le(a, x), Lt(x, b))
        } else if (incl_b) {
          c <- if (binf) Gt(x, a) else And(Lt(a, x), Le(x, b))
        } else if (inc(a)) {
          c <- Lt(x, b)
        } else c <- And(Lt(a, x), Lt(x, b))
        covered <- Union(covered, ivl)
        if (is_const(a, "-oo") && incl_a) covered <- Union(covered, FiniteSet(.neg(S_Infinity)))
        if (is_const(b, "oo") && incl_b) covered <- Union(covered, FiniteSet(S_Infinity))
        args[[length(args) + 1L]] <- list(e, c)
      }
      if (!isTRUE(is_subset(Reals, covered))) args[[length(args) + 1L]] <- list(S_NaN, S_true)
    }
  }
  if (is.null(args)) {
    args <- .pw_pairs(expr)
    for (i in seq_along(args)) args[[i]][[2L]] <- .simplify(args[[i]][[2L]], kw)
  }
  kw2 <- kw; kw2$doit <- FALSE
  for (i in seq_along(args)) {
    e <- args[[i]][[1L]]
    newe <- .simplify(e, kw2)
    if (!identical(newe$key, e$key)) args[[i]][[1L]] <- newe
  }
  .pw_build(args)
}

#' @rdname piecewise_simplify
#' @export
simplify_univariate <- function(expr) {
  expr <- as_sym(expr)
  if (!(expr$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE"))) return(expr)
  if (length(atoms(expr, c("Equality", "Unequality")))) return(expr)
  free <- free_symbols(expr)
  if (length(free) != 1L) return(expr)
  x <- free[[1L]]
  iv <- .pw_intervals(fnode_pw(list(S_Zero, expr)), x, err_on_Eq = TRUE)
  if (!isTRUE(iv[[1L]])) return(expr)
  if (!length(iv[[2L]])) return(S_false)
  args <- list()
  for (t in iv[[2L]]) {
    a <- t[[1L]]; b <- t[[2L]]
    if (is_const(a, "-oo")) {
      c <- if (is_const(b, "oo")) S_true else if (.is_true(subs(expr, x, b))) Le(x, b) else Lt(x, b)
    } else {
      incl_a <- .is_true(subs(expr, x, a)); incl_b <- .is_true(subs(expr, x, b))
      binf <- isTRUE(ask_fact(b, "infinite"))
      c <- if (incl_a && incl_b) (if (binf) Ge(x, a) else And(Le(a, x), Le(x, b)))
           else if (incl_a) And(Le(a, x), Lt(x, b))
           else if (incl_b) (if (binf) Gt(x, a) else And(Lt(a, x), Le(x, b)))
           else And(Lt(a, x), Lt(x, b))
    }
    args[[length(args) + 1L]] <- c
  }
  do.call(Or, args)
}

# An unevaluated Piecewise of the given (expr, cond) pairs
fnode_pw <- function(...) {
  ps <- list(...)
  new_sym("Piecewise", lapply(ps, function(p) new_sym("ExprCondPair", list(as_sym(p[[1L]]), as_sym(p[[2L]])))))
}
