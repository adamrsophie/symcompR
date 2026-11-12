# The SAT based part of ask(), after sympy/assumptions/satask.py,
# sathandlers.py, reasoning_engine.py and facts.py, and the wrapper
# functions of sympy/assumptions/wrapper.py.

.pyset <- function(lst) {
  lst <- unname(lst)
  if (length(lst)) lst <- lst[!duplicated(vapply(lst, .k, ""))]
  structure(lst, pyset = TRUE)
}

# ---------------------------------------------------------------------------
# sathandlers.py

.fact_args <- function(expr) if (expr$type %in% c("Add", "Mul")) .sympy_args(expr) else expr$args

# allargs, anyarg and exactlyonearg
.allargs <- function(symbol, fact, expr) do.call(And, lapply(.fact_args(expr), function(a) subs(fact, symbol, a)))
.anyarg <- function(symbol, fact, expr) do.call(Or, lapply(.fact_args(expr), function(a) subs(fact, symbol, a)))
.exactlyonearg <- function(symbol, fact, expr) {
  pa <- lapply(.fact_args(expr), function(a) subs(fact, symbol, a))
  do.call(Or, lapply(seq_along(pa), function(i) do.call(And, c(pa[i], lapply(pa[-i], Not)))))
}

# The placeholder symbol of the registered facts
.sat_x <- function() Symbol("x")

# _old_assump_getters, in their order
.OLD_ASSUMP_GETTERS <- c("positive", "zero", "negative", "rational", "irrational", "even", "odd",
                         "imaginary", "prime", "composite")

# The facts of class_fact_registry, by class; each function returns a list
# of facts (register() and multiregister() handlers alike).
.CLASS_FACTS <- list(
  Abs = function(expr) {
    arg <- expr$args[[1L]]
    list(.QA("nonnegative", expr),
         Equivalent(.QN("zero", arg), .QN("zero", expr)),
         Implies(.QA("even", arg), .QA("even", expr)),
         Implies(.QA("odd", arg), .QA("odd", expr)),
         Implies(.QA("integer", arg), .QA("integer", expr)))
  },
  Add = function(expr) {
    x <- .sat_x()
    allq <- function(p) .allargs(x, .QA(p, x), expr)
    list(Implies(allq("positive"), .QA("positive", expr)),
         Implies(allq("negative"), .QA("negative", expr)),
         Implies(allq("real"), .QA("real", expr)),
         Implies(allq("rational"), .QA("rational", expr)),
         Implies(allq("integer"), .QA("integer", expr)),
         Implies(.exactlyonearg(x, .QN("integer", x), expr), .QN("integer", expr)),
         Implies(allq("real"), Implies(.exactlyonearg(x, .QA("irrational", x), expr), .QA("irrational", expr))))
  },
  Mul = function(expr) {
    x <- .sat_x()
    allq <- function(p) .allargs(x, .QA(p, x), expr)
    allnz <- .allargs(x, .QN("zero", x), expr)
    list(Equivalent(.QA("zero", expr), .anyarg(x, .QA("zero", x), expr)),
         Implies(allq("positive"), .QA("positive", expr)),
         Implies(allq("real"), .QA("real", expr)),
         Implies(allq("rational"), .QA("rational", expr)),
         Implies(allq("integer"), .QA("integer", expr)),
         Implies(And(allnz, .exactlyonearg(x, .QN("rational", x), expr)), .QN("integer", expr)),
         Implies(allq("commutative"), .QA("commutative", expr)),
         Implies(allq("prime"), .QN("prime", expr)),
         Implies(And(.allargs(x, Or(.QA("imaginary", x), .QA("real", x)), expr), allnz),
                 Implies(.exactlyonearg(x, .QA("imaginary", x), expr), .QA("imaginary", expr))),
         Implies(And(allq("real"), allnz),
                 Implies(.exactlyonearg(x, .QA("irrational", x), expr), .QA("irrational", expr))),
         Implies(allq("integer"), Equivalent(.anyarg(x, .QA("even", x), expr), .QA("even", expr))))
  },
  MatMul = function(expr) {
    x <- .sat_x()
    list(Implies(.allargs(x, .QA("square", x), expr),
                 Equivalent(.QA("invertible", expr), .allargs(x, .QA("invertible", x), expr))))
  },
  Pow = function(expr) {
    base <- expr$args[[1L]]; e <- expr$args[[2L]]
    list(Implies(And(.QA("real", base), .QA("even", e), .QA("nonnegative", e)), .QA("nonnegative", expr)),
         Implies(And(.QA("nonnegative", base), .QA("odd", e), .QA("nonnegative", e)), .QA("nonnegative", expr)),
         Implies(And(.QA("nonpositive", base), .QA("odd", e), .QA("nonnegative", e)), .QA("nonpositive", expr)),
         Equivalent(.QA("zero", expr), And(.QA("zero", base), .QA("positive", e))))
  },
  Number = function(expr) .old_assump_facts(expr),
  NumberSymbol = function(expr) .old_assump_facts(expr),
  ImaginaryUnit = function(expr) .old_assump_facts(expr)
)

.old_assump_facts <- function(expr) {
  ret <- list()
  for (p in .OLD_ASSUMP_GETTERS) {
    prop <- .fact(expr, p)
    if (!is.na(prop)) ret[[length(ret) + 1L]] <- if (prop) .QA(p, expr) else .QN(p, expr)
  }
  ret
}

# class_fact_registry(expr): the set of facts registered for the classes of
# expr
.class_fact_registry <- function(expr) {
  cls <- .ask_py_class(expr)
  out <- list()
  for (k in names(.CLASS_FACTS)) if (.issubclass(cls, k)) out <- c(out, .CLASS_FACTS[[k]](expr))
  .pyset(out)
}

# ---------------------------------------------------------------------------
# satask.py

# find_symbols: the symbols of a CNF's predicates, or of an expression
.find_symbols <- function(pred) {
  if (inherits(pred, "ask_cnf")) {
    out <- list()
    for (a in .cnf_all_predicates(pred)) out <- c(out, .find_symbols(a))
    return(.pyset(out))
  }
  .pyset(atoms(pred, c("Symbol", "Dummy", "Wild")))
}

.keys_of <- function(lst) vapply(lst, .k, "")

# extract_predargs: the arguments of the predicates of proposition, and of
# the assumptions that share symbols with them
.extract_predargs <- function(proposition, assumptions = NULL) {
  req_keys <- .find_symbols(proposition)
  keys <- .cnf_all_predicates(proposition)
  lkeys <- list()
  if (!is.null(assumptions)) lkeys <- .cnf_all_predicates(assumptions)
  lkeys <- Filter(function(l) !(.is_true(l) || .is_false(l)), lkeys)
  lsyms <- lapply(lkeys, .find_symbols)
  repeat {
    rk <- .keys_of(req_keys)
    tmp <- list()
    for (s in lsyms) if (any(.keys_of(s) %in% rk)) tmp <- c(tmp, s)
    tmp <- .pyset(tmp)
    tmp_keys <- tmp[!(.keys_of(tmp) %in% rk)]
    if (!length(tmp_keys)) break
    req_keys <- .pyset(c(req_keys, tmp_keys))
  }
  rk <- .keys_of(req_keys)
  for (i in seq_along(lkeys)) if (any(.keys_of(lsyms[[i]]) %in% rk)) keys[[length(keys) + 1L]] <- lkeys[[i]]
  exprs <- list()
  for (key in keys) {
    if (.is_applied(key)) exprs <- c(exprs, .ap_args(key)) else exprs[[length(exprs) + 1L]] <- key
  }
  .pyset(exprs)
}

# get_relevant_clsfacts: the facts of the classes of exprs, and the
# expressions those facts are about
.get_relevant_clsfacts <- function(exprs, relevant_facts = NULL) {
  if (is.null(relevant_facts)) relevant_facts <- .CNF()
  newexprs <- list()
  for (expr in exprs) {
    for (fact in .class_fact_registry(expr)) {
      newfact <- .CNF_to_CNF(fact)
      relevant_facts <- .cnf_and(relevant_facts, newfact)
      for (key in .cnf_all_predicates(newfact))
        if (.is_applied(key)) newexprs <- c(newexprs, .ap_args(key))
    }
  }
  newexprs <- .pyset(newexprs)
  list(.pyset(newexprs[!(.keys_of(newexprs) %in% .keys_of(exprs))]), relevant_facts)
}

.is_matrix_kind <- function(e) .is_matexpr(e) || .is_matrixbase(e)

# The encoded known facts of the given kinds, made once per kind
.known_facts_template <- function(matrix, number) {
  k <- paste0("tmpl", matrix, number)
  v <- .known_cache[[k]]
  if (is.null(v)) {
    cnf <- .CNF()
    if (matrix) cnf <- .cnf_and(cnf, .get_all_known_matrix_facts())
    if (number) cnf <- .cnf_and(cnf, .get_all_known_number_facts())
    v <- .enc_cnf()
    .enc_from_cnf(v, cnf)
    assign(k, v, envir = .known_cache)
  }
  v
}

# get_all_relevant_facts: the facts of the classes of every expression
# involved, with the known facts repeated for each expression
.get_all_relevant_facts <- function(proposition, assumptions, use_known_facts = TRUE, iterations = Inf) {
  i <- 0L
  relevant_facts <- .CNF()
  all_exprs <- list()
  repeat {
    if (i == 0L) exprs <- .extract_predargs(proposition, assumptions)
    all_exprs <- c(all_exprs, exprs)
    r <- .get_relevant_clsfacts(exprs, relevant_facts)
    exprs <- r[[1L]]
    relevant_facts <- r[[2L]]
    i <- i + 1L
    if (i >= iterations) break
    if (!length(exprs)) break
  }
  all_exprs <- .pyset(all_exprs)
  ctx <- .enc_cnf()
  if (use_known_facts) {
    mat <- any(vapply(all_exprs, .is_matrix_kind, TRUE))
    num <- any(vapply(all_exprs, function(e) !.is_matrix_kind(e) && !.is_boolean_kind(e), TRUE))
    kf <- .known_facts_template(mat, num)
    n_lit <- length(kf$symbols)
    data <- list()
    symbols <- list()
    for (j in seq_along(all_exprs)) {
      expr <- all_exprs[[j]]
      delta <- (j - 1L) * n_lit
      symbols <- c(symbols, lapply(kf$symbols, function(p) new_sym("AppliedPredicate", list(p, expr))))
      data <- c(data, lapply(kf$data, function(cl) cl + sign(cl) * delta))
    }
    ctx$symbols <- symbols
    if (length(symbols)) {
      ks <- .keys_of(symbols)
      for (j in seq_along(ks)) assign(ks[[j]], j, envir = ctx$encoding)
    }
    ctx$data <- data
  }
  .enc_add_from_cnf(ctx, relevant_facts)
  ctx
}

.enc_copy <- function(enc) {
  e <- .enc_cnf()
  e$symbols <- enc$symbols
  e$encoding <- list2env(as.list(enc$encoding), parent = emptyenv())
  e$data <- enc$data
  e
}

# The clauses of cnf encoded with an extra literal: activation_literal
# implies each clause.
.encode_with_selector <- function(cnf, enc, activation_literal) {
  new <- lapply(unname(.subset2(cnf, "clauses")), function(cl) {
    v <- .enc_clause(enc, cl)
    unique(c(v[v != 0L], -activation_literal))
  })
  enc$data <- c(enc$data, new)
  invisible(enc)
}

#' Decide a proposition by satisfiability
#'
#' `satask()` is the search over facts that [ask()] uses once the handlers
#' of the predicates give no answer: it collects the facts known about every
#' expression in the proposition and the assumptions, such as
#' `Q$nonnegative(Abs(x))`, and asks whether the proposition, or its
#' negation, can hold together with them.
#'
#' @param proposition A boolean expression of predicates.
#' @param assumptions Assumptions, or `TRUE` for none.
#' @param use_known_facts Whether the facts between the predicates are used.
#' @param iterations How many times facts are collected from the
#'   expressions found; `Inf` until no new expression appears.
#' @param early_return Whether to answer from unit propagation alone when it
#'   decides the proposition, trusting the assumptions to be consistent.
#' @return `TRUE`, `FALSE` or `NA`. Inconsistent assumptions raise a
#'   `ValueError`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' satask(Q$zero(Abs(x)), Q$zero(x))
#' satask(Q$zero(x * y), Q$zero(x))
#' @export
satask <- function(proposition, assumptions = TRUE, use_known_facts = TRUE, iterations = Inf,
                   early_return = FALSE) {
  proposition <- .ask_bool(proposition)
  assumptions <- .ask_bool(assumptions)
  props <- .CNF_from_prop(proposition)
  nprops <- .CNF_from_prop(Not(proposition))
  acnf <- .CNF_from_prop(assumptions)
  sat <- .get_all_relevant_facts(props, acnf, use_known_facts = use_known_facts, iterations = iterations)
  .enc_add_from_cnf(sat, acnf)
  # ReasoningEngine
  if (any(vapply(sat$data, function(cl) length(cl) == 1L && cl == 0L, TRUE)))
    .ask_value_error("Inconsistent assumptions")
  guarded <- .enc_copy(sat)
  selector <- length(guarded$symbols) + 1L
  guarded$symbols[[selector]] <- S_true
  .encode_with_selector(props, guarded, selector)
  .encode_with_selector(nprops, guarded, -selector)
  n <- length(guarded$symbols)
  fixed <- .sat_fixed(guarded$data, n)
  if (is.null(fixed)) .ask_value_error("Inconsistent assumptions")
  if (early_return && fixed[selector] != 0L) return(fixed[selector] > 0L)
  can_true <- !is.null(.sat_solve(guarded$data, n, selector))
  can_false <- !is.null(.sat_solve(guarded$data, n, -selector))
  if (!can_true && !can_false) .ask_value_error("Inconsistent assumptions")
  if (can_true && can_false) return(NA)
  can_true
}

# ---------------------------------------------------------------------------
# facts.py

# ask_full_inference
.ask_full_inference <- function(proposition, assumptions, known_facts) {
  sat <- function(e) {
    enc <- .enc_cnf()
    .enc_from_cnf(enc, .CNF_from_prop(e))
    .enc_satisfiable(enc)
  }
  if (!sat(And(known_facts, assumptions, proposition))) return(FALSE)
  if (!sat(And(known_facts, assumptions, Not(proposition)))) return(TRUE)
  NA
}

# single_fact_lookup: for each key, the keys and negated keys it implies
.single_fact_lookup <- function(known_facts_keys, known_facts) {
  out <- list()
  for (key in known_facts_keys) {
    m <- list(key)
    for (other in known_facts_keys) {
      if (.keq(other, key)) next
      if (isTRUE(.ask_full_inference(other, key, known_facts))) m[[length(m) + 1L]] <- other
      if (isTRUE(.ask_full_inference(Not(other), key, known_facts))) m[[length(m) + 1L]] <- Not(other)
    }
    out[[length(out) + 1L]] <- list(key, .pyset(m))
  }
  structure(out, pydict = TRUE)
}

# generate_known_facts_dict: for each predicate of the keys, the predicates
# it implies and those it rejects (the fact itself stands for its CNF)
.generate_known_facts_dict <- function(keys, fact) {
  mapping <- .single_fact_lookup(keys, fact)
  out <- list()
  for (kv in mapping) {
    key <- kv[[1L]]
    implied <- list(); rejected <- list()
    for (v in kv[[2L]]) {
      if (v$type == "Not") rejected[[length(rejected) + 1L]] <- .ap_fun(v$args[[1L]])
      else implied[[length(implied) + 1L]] <- .ap_fun(v)
    }
    out[[length(out) + 1L]] <- list(.ap_fun(key), structure(list(.pyset(implied), .pyset(rejected)), pytuple = TRUE))
  }
  structure(out, pydict = TRUE)
}

# ---------------------------------------------------------------------------
# wrapper.py

.WRAPPER_FACTS <- c("algebraic", "antihermitian", "commutative", "complex", "composite", "even",
                    "extended_negative", "extended_nonnegative", "extended_nonpositive",
                    "extended_nonzero", "extended_positive", "extended_real", "finite",
                    "hermitian", "imaginary", "infinite", "integer", "irrational", "negative",
                    "noninteger", "nonnegative", "nonpositive", "nonzero", "odd", "polar",
                    "positive", "prime", "rational", "real", "transcendental", "zero")

#' Query predicates of an expression under assumptions
#'
#' `AssumptionsWrapper(expr, assumptions)` answers `$is_positive`,
#' `$is_integer` and the other `is_<fact>` queries of `expr`: from the
#' expression's own assumptions when `assumptions` is `NULL`, and through
#' [ask()] with the assumptions otherwise. `is_extended_real(obj,
#' assumptions)` and `is_extended_nonnegative(obj, assumptions)` are the
#' single queries.
#'
#' @param expr,obj An expression.
#' @param assumptions Assumptions such as `Q$positive(x)`, or `NULL`.
#' @return `AssumptionsWrapper()` returns an object whose `$is_<fact>`
#'   fields are `TRUE`, `FALSE` or `NA`; the others return one of those.
#' @examples
#' y <- Symbol("y")
#' AssumptionsWrapper(y)$is_positive
#' AssumptionsWrapper(y, Q$positive(y))$is_positive
#' is_extended_real(y, Q$extended_real(y))
#' @export
AssumptionsWrapper <- function(expr, assumptions = NULL) {
  structure(list(expr = .ask_sympify(expr), assumptions = assumptions), class = "AssumptionsWrapper")
}

#' @export
`$.AssumptionsWrapper` <- function(x, name) {
  if (name %in% c("expr", "assumptions")) return(.subset2(x, name))
  if (!startsWith(name, "is_") || !(substring(name, 4L) %in% .WRAPPER_FACTS)) return(NULL)
  .awf(.subset2(x, "expr"), substring(name, 4L), .subset2(x, "assumptions"))
}

#' @export
print.AssumptionsWrapper <- function(x, ...) {
  a <- .subset2(x, "assumptions")
  cat("AssumptionsWrapper(", .str(.subset2(x, "expr")),
      if (!is.null(a)) paste0(", ", .str(.ask_bool(a))), ")\n", sep = "")
  invisible(x)
}

#' @rdname AssumptionsWrapper
#' @export
is_extended_real <- function(obj, assumptions = NULL) .awf(.ask_sympify(obj), "extended_real", assumptions)

#' @rdname AssumptionsWrapper
#' @export
is_extended_nonnegative <- function(obj, assumptions = NULL)
  .awf(.ask_sympify(obj), "extended_nonnegative", assumptions)
