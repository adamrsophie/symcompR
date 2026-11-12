# ask() and its helpers, after sympy/assumptions/ask.py, with the handler
# dispatch of sympy/multipledispatch: each predicate's handler tries its
# registered signatures in SymPy's order (.ASK_HANDLER_ORDER), takes the
# first whose classes match the arguments' SymPy classes, and moves on to
# the next one when a handler raises MDNotImplementedError (here, returns
# .MDNI). No match leaves the answer unknown.

.MDNI <- structure("MDNotImplementedError", class = "ask_mdni")

# ---------------------------------------------------------------------------
# SymPy classes of expressions

.CONST_PYCLASS <- c(oo = "Infinity", `-oo` = "NegativeInfinity", zoo = "ComplexInfinity", nan = "NaN",
                    I = "ImaginaryUnit", pi = "Pi", E = "Exp1")

.ask_py_class <- function(e) {
  if (!is_sym(e)) {
    if (is.logical(e)) return("bool")
    if (is.numeric(e)) return("int")
    if (inherits(e, "symmatrix")) return("MutableDenseMatrix")
    return("object")
  }
  t <- e$type
  switch(t,
    Integer = {
      v <- e$value
      if (v == 0) "Zero" else if (v == 1) "One" else if (v == -1) "NegativeOne" else "Integer"
    },
    Rational = if (e$value == gmp::as.bigq(1, 2)) "Half" else "Rational",
    Float = "Float",
    Constant = if (!is.na(.CONST_PYCLASS[e$name])) .CONST_PYCLASS[[e$name]] else e$name,
    BooleanAtom = if (e$name == "true") "BooleanTrue" else "BooleanFalse",
    Function = if (!is.null(.ASK_MRO[[e$name]])) e$name else "AppliedUndef",
    ImmutableMatrix = "ImmutableDenseMatrix",
    Predicate = if (!is.null(e$extra$sympy_class)) e$extra$sympy_class else "Predicate",
    t)
}

.mro_cache <- new.env(parent = emptyenv())
.py_mro <- function(cls) {
  m <- .mro_cache[[cls]]
  if (!is.null(m)) return(m)
  m <- .ASK_MRO[[cls]]
  if (is.null(m)) m <- c(cls, "Expr", "Basic", "object")
  assign(cls, m, envir = .mro_cache)
  m
}
.issubclass <- function(cls, sig) sig %in% .py_mro(cls)
.isinstance <- function(e, cls) .issubclass(.ask_py_class(e), cls)

# ---------------------------------------------------------------------------
# Handler registry

.ask_handlers <- new.env(parent = emptyenv())

# register(), and register_many() when sigs has several entries; a
# signature of several classes is written "Integer,Integer"
.ask_register <- function(pred, sigs, fn) {
  h <- .ask_handlers[[pred]]
  if (is.null(h)) {
    h <- new.env(parent = emptyenv())
    assign(pred, h, envir = .ask_handlers)
  }
  for (s in sigs) assign(s, fn, envir = h)
  invisible(NULL)
}

.sig_cache <- new.env(parent = emptyenv())
.split_sig <- function(s) {
  v <- .sig_cache[[s]]
  if (is.null(v)) {
    v <- strsplit(s, ",", fixed = TRUE)[[1L]]
    assign(s, v, envir = .sig_cache)
  }
  v
}

# A handler result as TRUE, FALSE or NA (Python's True, False and None).
.ask_fz <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA)
  if (is_sym(x)) return(if (.is_true(x)) TRUE else if (.is_false(x)) FALSE else NA)
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  NA
}

# Dispatcher.__call__ through Predicate.eval: NA when nothing matches.
.dispatch <- function(order, funs, args, assumptions) {
  classes <- vapply(args, .ask_py_class, "")
  n <- length(args)
  for (sig in order) {
    sc <- .split_sig(sig)
    if (length(sc) != n) next
    ok <- TRUE
    for (i in seq_len(n)) if (!.issubclass(classes[[i]], sc[[i]])) { ok <- FALSE; break }
    if (!ok) next
    fn <- funs[[sig]]
    if (is.null(fn)) next
    r <- do.call(fn, c(args, list(assumptions)))
    if (inherits(r, "ask_mdni")) next
    return(.ask_fz(r))
  }
  NA
}

# Predicate.eval(args, assumptions)
.pred_eval <- function(pred, args, assumptions) {
  if (pred$type == "UndefinedPredicate") return(NA)
  custom <- .ask_custom[[pred$key]]
  if (isTRUE(pred$extra$custom)) {
    if (is.null(custom)) return(NA)
    return(.dispatch(custom$order, custom$funs, args, assumptions))
  }
  nm <- pred$name
  if (nm %in% .BINREL_NAMES) return(.binrel_eval(nm, args, assumptions))
  order <- .ASK_HANDLER_ORDER[[nm]]
  funs <- .ask_handlers[[nm]]
  if (!is.null(custom)) {
    order <- c(custom$order, order)
    funs <- c(custom$funs, if (!is.null(funs)) as.list(funs))
  }
  if (is.null(order) || is.null(funs)) return(NA)
  .dispatch(order, funs, args, assumptions)
}

# AppliedPredicate._eval_ask
.eval_ask <- function(e, assumptions) {
  if (e$type == "AppliedBinaryRelation") return(.binrel_eval_ask(e, assumptions))
  .pred_eval(.ap_fun(e), .ap_args(e), assumptions)
}

# ---------------------------------------------------------------------------
# Small helpers used by the handlers

.QA <- function(name, ...) new_sym("AppliedPredicate", c(list(.qp(name)), list(...)))
.QN <- function(name, ...) Not(.QA(name, ...))
.ask_recursive_q <- function(name, e, assumptions) .ask_recursive(.QA(name, e), assumptions)
.askr <- .ask_recursive_q

# old assumptions handler: expr.is_<fact>, or MDNotImplementedError
.old_fact_handler <- function(fact) {
  force(fact)
  function(expr, assumptions) {
    ret <- .fact(expr, fact)
    if (is.na(ret)) .MDNI else ret
  }
}

# old assumptions: expr.is_<fact>
.fact <- function(e, f) {
  if (!is_sym(e)) return(NA)
  r <- tryCatch(ask_fact(e, f), error = function(err) NA)
  if (length(r) != 1L) NA else r
}

.is_number <- function(e) {
  if (!is_sym(e)) return(is.numeric(e))
  if (e$type %in% c(.NUMBER_TYPES, "Constant")) return(TRUE)
  if (e$type %in% c("Symbol", "Dummy", "Wild") || .is_matexpr(e) || e$type %in% c(.PRED_TYPES, .APPLIED_TYPES, .BOOL_TYPES, .REL_TYPES))
    return(FALSE)
  .expr_is_number(e)
}

# Python truthiness of a fuzzy value
.tt <- function(x) isTRUE(x)

.fuzzy_not <- function(x) if (is.na(x)) NA else !x
.fuzzy_and_vals <- function(vals) {
  rv <- TRUE
  for (v in vals) {
    v <- .ask_fz(v)
    if (identical(v, FALSE)) return(FALSE)
    if (isTRUE(rv)) rv <- v
  }
  rv
}
.fuzzy_or_vals <- function(vals) .fuzzy_not(.fuzzy_and_vals(lapply(vals, function(v) .fuzzy_not(.ask_fz(v)))))

# fuzzy_and and fuzzy_or over lazily evaluated items: f(i) for i in 1..n
.fuzzy_and_lazy <- function(n, f) {
  rv <- TRUE
  for (i in seq_len(n)) {
    v <- .ask_fz(f(i))
    if (identical(v, FALSE)) return(FALSE)
    if (isTRUE(rv)) rv <- v
  }
  rv
}
.fuzzy_or_lazy <- function(n, f) .fuzzy_not(.fuzzy_and_lazy(n, function(i) .fuzzy_not(.ask_fz(f(i)))))

# _fuzzy_group(..., quick_exit = TRUE) over f(i)
.fuzzy_group_lazy <- function(n, f, quick_exit = FALSE) {
  saw_other <- FALSE
  for (i in seq_len(n)) {
    a <- .ask_fz(f(i))
    if (isTRUE(a)) next
    if (is.na(a)) return(NA)
    if (quick_exit && saw_other) return(NA)
    saw_other <- TRUE
  }
  !saw_other
}

# test_closed_group, ask_all and ask_any of handlers/common.py
.test_closed_group <- function(expr, assumptions, key) {
  a <- if (expr$type %in% c("Add", "Mul")) .sympy_args(expr) else expr$args
  .fuzzy_group_lazy(length(a), function(i) .askr(key, a[[i]], assumptions), quick_exit = TRUE)
}
.ask_all <- function(queries, assumptions) .fuzzy_and_lazy(length(queries), function(i) .ask_recursive(queries[[i]], assumptions))
.ask_any <- function(queries, assumptions) .fuzzy_or_lazy(length(queries), function(i) .ask_recursive(queries[[i]], assumptions))

# Python all() and any() over lazily evaluated truthiness
.py_all <- function(n, f) { for (i in seq_len(n)) if (!isTRUE(.ask_fz(f(i)))) return(FALSE); TRUE }
.py_any <- function(n, f) { for (i in seq_len(n)) if (isTRUE(.ask_fz(f(i)))) return(TRUE); FALSE }

# conjuncts(assumptions) as keys, for "Q.x(expr) in conjuncts(assumptions)"
.ask_conjuncts <- function(a) {
  a <- .ask_bool(a)
  if (a$type == "And") a$args else list(a)
}
.in_conjuncts <- function(p, assumptions) p$key %in% vapply(.ask_conjuncts(assumptions), .k, "")

.ask_bool <- function(a) {
  if (inherits(a, "sym_predicate_fn")) return(attr(a, "predicate"))
  if (is.logical(a)) return(if (isTRUE(a)) S_true else S_false)
  .ask_sympify(a)
}

# ---------------------------------------------------------------------------
# Known facts

.known_cache <- new.env(parent = emptyenv())

.facts_to_cnf <- function(strs) {
  cls <- lapply(strsplit(strs, " ", fixed = TRUE), function(lits)
    .clause(lapply(lits, function(l) {
      neg <- startsWith(l, "!")
      .Literal(.qp(sub("^!", "", l)), neg)
    })))
  .cnf(.named_clauses(cls))
}

.get_all_known_facts <- function() {
  v <- .known_cache$all
  if (is.null(v)) { v <- .facts_to_cnf(.ASK_ALL_FACTS); .known_cache$all <- v }
  v
}
.get_all_known_matrix_facts <- function() {
  v <- .known_cache$matrix
  if (is.null(v)) { v <- .facts_to_cnf(.ASK_MATRIX_FACTS); .known_cache$matrix <- v }
  v
}
.get_all_known_number_facts <- function() {
  v <- .known_cache$number
  if (is.null(v)) { v <- .facts_to_cnf(.ASK_NUMBER_FACTS); .known_cache$number <- v }
  v
}

# get_known_facts_dict: implied and rejected predicate names
.known_facts_entry <- function(pred) {
  if (pred$type != "Predicate" || isTRUE(pred$extra$custom)) return(NULL)
  .ASK_FACTS_DICT[[pred$name]]
}

.get_known_facts_dict <- function() {
  out <- lapply(.ASK_FACTS_DICT, function(v) list(structure(lapply(v[[1L]], .qp), pyset = TRUE),
                                                   structure(lapply(v[[2L]], .qp), pyset = TRUE)))
  out
}

# The integer encoding of the known facts, made once.
.known_enc <- function() {
  v <- .known_cache$enc
  if (is.null(v)) {
    v <- .enc_cnf()
    .enc_from_cnf(v, .get_all_known_facts())
    .known_cache$enc <- v
  }
  v
}

# ---------------------------------------------------------------------------
# ask

.cnf_memo <- new.env(parent = emptyenv())
.cnf_memo_n <- 0L
# CNF.from_prop for assumptions, remembered by key (handlers ask with the
# same assumptions many times)
.assump_cnf <- function(assumptions) {
  k <- assumptions$key
  v <- .cnf_memo[[k]]
  if (!is.null(v)) return(v)
  v <- .CNF_from_prop(assumptions)
  if (length(ls(.cnf_memo)) > 500L) rm(list = ls(.cnf_memo), envir = .cnf_memo)
  assign(k, v, envir = .cnf_memo)
  v
}

# _extract_all_facts
.extract_all_facts <- function(assump, exprs) {
  ekeys <- vapply(exprs, .k, "")
  facts <- list()
  for (clause in .subset2(assump, "clauses")) {
    args <- list()
    ok <- TRUE
    for (literal in clause) {
      l <- literal$lit
      if (.is_applied(l) && length(l$args) == 2L && l$args[[2L]]$key %in% ekeys) {
        args[[length(args) + 1L]] <- .Literal(.ap_fun(l), literal$neg)
      } else {
        ok <- FALSE
        break
      }
    }
    if (ok && length(args)) facts[[length(facts) + 1L]] <- .clause(args)
  }
  .cnf(.named_clauses(facts))
}

# _ask_single_fact
.ask_single_fact <- function(key, local_facts) {
  key <- .as_pred(key)
  if (!length(.subset2(local_facts, "clauses"))) return(NA)
  get_facts <- function(k) {
    v <- .known_facts_entry(k)
    if (is.null(v)) list(character(), character()) else v
  }
  kname <- if (key$type == "Predicate" && !isTRUE(key$extra$custom)) key$name else ""
  for (clause in .subset2(local_facts, "clauses")) {
    if (length(clause) != 1L) next
    f <- clause[[1L]]
    pred <- f$lit
    negated <- f$neg
    key_req <- get_facts(key)[[1L]]
    pname <- if (pred$type == "Predicate" && !isTRUE(pred$extra$custom)) pred$name else ""
    if (negated && pname %in% key_req) return(FALSE)
    if (!negated) {
      rr <- get_facts(pred)
      if (kname %in% rr[[1L]]) return(TRUE)
      if (kname %in% rr[[2L]]) return(FALSE)
    }
  }
  NA
}

.is_boolean_kind <- function(e) {
  is_sym(e) && e$type %in% c(.APPLIED_TYPES, .BOOL_TYPES, .REL_TYPES, "Nand", "Nor", "Xnor",
                             "Exclusive", "Contains")
}

# _normalize_expr: relations become Q.lt, Q.le, Q.eq and Q.ne
.normalize_expr <- function(e) {
  .map_bottom_up(e, function(z) {
    t <- z$type
    if (t %in% .REL_TYPES && length(z$args) == 2L) {
      a <- z$args[[1L]]; b <- z$args[[2L]]
      z <- switch(t,
        StrictGreaterThan = .apply_pred(.qp("lt"), list(b, a)),
        GreaterThan = .apply_pred(.qp("le"), list(b, a)),
        StrictLessThan = .apply_pred(.qp("lt"), list(a, b)),
        LessThan = .apply_pred(.qp("le"), list(a, b)),
        Equality = .apply_pred(.qp("eq"), list(a, b)),
        Unequality = .apply_pred(.qp("ne"), list(a, b)),
        z)
    }
    if (z$type == "AppliedBinaryRelation") {
      nm <- .ap_name(z)
      if (nm == "gt") z <- .apply_pred(.qp("lt"), list(z$args[[3L]], z$args[[2L]]))
      else if (nm == "ge") z <- .apply_pred(.qp("le"), list(z$args[[3L]], z$args[[2L]]))
    }
    z
  })
}

#' Ask whether a proposition holds under assumptions
#'
#' `ask(proposition, assumptions)` evaluates a proposition built from the
#' predicates of [Q] to `TRUE`, `FALSE` or `NA` (unknown), using the
#' assumptions given, the extra assumptions in `...` (joined with `&`) and
#' those of `context`, which is [global_assumptions] unless another
#' [AssumptionsContext()] is given. This is SymPy's `ask()`: the predicates'
#' handlers answer first, then a SAT based search over the facts known
#' about the expressions involved, then the linear real arithmetic solver for
#' inequalities. Relations such as `x > 0` may appear in the proposition and
#' the assumptions. A proposition or assumptions that are not boolean raise
#' a `TypeError`, and inconsistent assumptions a `ValueError`.
#'
#' @param proposition A boolean expression such as `Q$positive(x)`,
#'   `Q$even(x) | Q$odd(x)` or `x > 1`.
#' @param assumptions Assumptions, such as `Q$positive(x) & Q$integer(y)`, or
#'   `TRUE` for none.
#' @param ... Further assumptions, joined to `assumptions` with `&`.
#' @param context An [AssumptionsContext()] whose assumptions are also used.
#' @return `TRUE`, `FALSE` or `NA`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' ask(Q$rational(sym_pi))
#' ask(Q$even(x * y), Q$even(x) & Q$integer(y))
#' ask(Q$prime(4 * x), Q$integer(x))
#' ask(Q$odd(3 * x))
#' ask(x > 0, Q$positive(x))
#' @export
ask <- function(proposition, assumptions = TRUE, ..., context = global_assumptions) {
  if (inherits(proposition, "sym_predicate_fn"))
    .ask_type_error("proposition must be a valid logical expression")
  proposition <- .ask_bool(proposition)
  if (.is_predicate(proposition) || !.is_boolean_kind(proposition))
    .ask_type_error("proposition must be a valid logical expression")
  if (inherits(assumptions, "sym_predicate_fn"))
    .ask_type_error("assumptions must be a valid logical expression")
  parts <- c(list(.ask_bool(assumptions)), lapply(list(...), .ask_bool), .ctx_preds(context))
  assumptions <- if (length(parts) == 1L) parts[[1L]] else do.call(And, parts)
  if (.is_predicate(assumptions) || !.is_boolean_kind(assumptions))
    .ask_type_error("assumptions must be a valid logical expression")
  proposition <- .normalize_expr(proposition)
  assumptions <- .normalize_expr(assumptions)
  if (.is_applied(proposition)) {
    key <- .ap_fun(proposition)
    args <- .ap_args(proposition)
  } else {
    key <- .qp("is_true")
    args <- list(proposition)
  }
  assump_cnf <- .assump_cnf(assumptions)
  local_facts <- .extract_all_facts(assump_cnf, args)
  if (length(.subset2(local_facts, "clauses"))) {
    kf <- .known_enc()
    enc <- .enc_cnf()
    enc$symbols <- kf$symbols
    enc$encoding <- list2env(as.list(kf$encoding), parent = emptyenv())
    enc$data <- kf$data
    .enc_add_from_cnf(enc, local_facts)
    if (!.enc_satisfiable(enc)) .ask_value_error(paste0("inconsistent assumptions ", .str(assumptions)))
  }
  res <- .ask_single_fact(key, local_facts)
  if (!is.na(res)) return(res)
  res <- .eval_ask(.apply_pred(key, args), assumptions)
  if (!is.na(res)) return(res)
  res <- satask(proposition, assumptions = assumptions)
  if (!is.na(res)) return(res)
  tryCatch(lra_satask(proposition, assumptions = assumptions), UnhandledInput = function(e) NA)
}

# _ask_recursive: the handlers alone, without the SAT search
.ask_recursive <- function(proposition, assumptions = TRUE) {
  proposition <- .ask_bool(proposition)
  assumptions <- .ask_bool(assumptions)
  if (.is_applied(proposition)) {
    key <- .ap_fun(proposition)
    args <- .ap_args(proposition)
  } else {
    key <- .qp("is_true")
    args <- list(proposition)
  }
  local_facts <- .extract_all_facts(.assump_cnf(assumptions), args)
  res <- .ask_single_fact(key, local_facts)
  if (!is.na(res)) return(res)
  .eval_ask(.apply_pred(key, args), assumptions)
}
