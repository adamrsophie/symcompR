# Propositional inference, after sympy/logic/inference.py and the DPLL
# algorithm of sympy/logic/algorithms/dpll.py.

# the iteration order of a CPython set of small integers built by adding
# vals in turn (open addressing with linear probes and perturbation)
.pyset_iter <- function(vals) {
  vals <- as.integer(vals)
  vals <- vals[!duplicated(vals)]
  ins <- function(tab, k) {
    mask <- length(tab) - 1L
    h <- if (k == -1L) -2L else k
    i <- bitwAnd(h, mask)
    perturb <- if (h >= 0L) gmp::as.bigz(h) else gmp::as.bigz(2)^64 + h
    repeat {
      if (is.na(tab[i + 1L])) { tab[i + 1L] <- k; return(tab) }
      if (i + 9L <= mask) for (j in 1:9) if (is.na(tab[i + j + 1L])) { tab[i + j + 1L] <- k; return(tab) }
      perturb <- perturb %/% 32L
      i <- as.integer((i * 5L + 1L + perturb) %% (mask + 1L))
    }
  }
  tab <- rep(NA_integer_, 8L)
  used <- 0L
  for (k in vals) {
    tab <- ins(tab, k)
    used <- used + 1L
    mask <- length(tab) - 1L
    if (used * 5L >= mask * 3L) {
      size <- 8L
      while (size <= used * 4L) size <- size * 2L
      old <- tab[!is.na(tab)]
      tab <- rep(NA_integer_, size)
      for (o in old) tab <- ins(tab, o)
    }
  }
  tab[!is.na(tab)]
}

# a model is a named list; its names are the strings of the symbols
.model_has <- function(model, s) !is.null(model[[.str(s)]])
.model_set <- function(model, s, v) { model[[.str(s)]] <- v; model }

#' The symbol of a literal
#'
#' `literal_symbol(literal)` is the symbol of a literal (a symbol or its
#' negation), as SymPy's `literal_symbol()`; `TRUE` and `FALSE` are returned
#' unchanged and anything else is an error.
#'
#' @param literal A literal.
#' @return A symbol or a logical value.
#' @examples
#' A <- Symbol("A")
#' literal_symbol(Not(A))
#' @export
literal_symbol <- function(literal) {
  if (isTRUE(literal) || isFALSE(literal)) return(literal)
  literal <- as_sym(literal)
  if (is_symbol_like(literal)) return(literal)
  if (literal$type == "Not") return(literal_symbol(literal$args[[1L]]))
  stop(.sympy_condition("ValueError", "Argument must be a boolean literal."))
}

# Or.make_args and And.make_args
.in_args <- function(x, args) any(vapply(args, function(a) identical(a$key, x$key), TRUE))

#' Evaluate a proposition in a model
#'
#' `pl_true(expr, model, deep)` returns `TRUE` or `FALSE` when the model (a
#' named list of logical values keyed by symbol names) decides `expr`, and
#' `NA` (SymPy's `None`) otherwise. With `deep = TRUE` an undecided result
#' is checked for validity and satisfiability, as SymPy's `pl_true()`.
#'
#' @param expr A boolean expression.
#' @param model A named list of logical values.
#' @param deep Look further when the model does not decide.
#' @return `TRUE`, `FALSE` or `NA`.
#' @examples
#' A <- Symbol("A"); B <- Symbol("B")
#' pl_true(Or(A, B), list(A = TRUE))
#' pl_true(And(A, B), list(A = TRUE))
#' @export
pl_true <- function(expr, model = NULL, deep = FALSE) {
  if (isTRUE(expr) || isFALSE(expr)) return(expr)
  expr <- as_sym(expr)
  if (.is_true(expr)) return(TRUE)
  if (.is_false(expr)) return(FALSE)
  valid_bool <- function(e) {
    if (is_symbol_like(e) || e$type == "BooleanAtom") return(TRUE)
    if (!(e$type %in% .BOOL_TYPES)) return(FALSE)
    all(vapply(e$args, valid_bool, TRUE))
  }
  if (!valid_bool(expr))
    stop(.sympy_condition("ValueError", sprintf("%s is not a valid boolean expression", sstr(expr))))
  model <- Filter(function(v) isTRUE(v) || isFALSE(v), as.list(model))
  pairs <- lapply(names(model), function(nm) list(nm, if (model[[nm]]) S_true else S_false))
  result <- if (length(pairs)) xreplace(expr, pairs) else expr
  if (.is_true(result)) return(TRUE)
  if (.is_false(result)) return(FALSE)
  if (deep) {
    at <- .bool_atoms(result)
    m <- stats::setNames(rep(list(TRUE), length(at)), vapply(at, .str, ""))
    if (isTRUE(pl_true(result, m))) {
      if (valid(result)) return(TRUE)
    } else if (isFALSE(satisfiable(result))) return(FALSE)
  }
  NA
}

#' Entailment and knowledge bases
#'
#' `entails(expr, formula_set)` tells whether the formulas imply `expr`.
#' `PropKB()` is a propositional knowledge base: `kb$tell(sentence)` adds
#' the clauses of a sentence, `kb$ask(query)` asks whether they entail the
#' query, `kb$retract(sentence)` removes clauses and `kb$clauses()` lists
#' them in SymPy's order. `KB()` is the abstract base, whose methods signal
#' `NotImplementedError`. These follow `sympy.logic.inference`.
#'
#' @param expr A boolean expression.
#' @param formula_set A list of boolean expressions.
#' @param sentence A sentence to tell first.
#' @return `TRUE` or `FALSE`; a knowledge base.
#' @examples
#' A <- Symbol("A"); B <- Symbol("B")
#' entails(B, list(Implies(A, B), A))
#' kb <- PropKB()
#' kb$tell(Implies(A, B)); kb$tell(A)
#' kb$ask(B)
#' @export
entails <- function(expr, formula_set = NULL) {
  if (inherits(formula_set, "sym")) formula_set <- list(formula_set)
  fs <- if (length(formula_set)) lapply(as.list(formula_set), .as_bool) else list()
  fs[[length(fs) + 1L]] <- Not(expr)
  !.sat_truthy(satisfiable(do.call(And, fs)))
}

.sat_truthy <- function(m) !(isFALSE(m) || (is.list(m) && !length(m)))

.kb_new <- function(cls, sentence, tell, ask, retract) {
  self <- new.env(parent = emptyenv())
  self$clauses_ <- list()
  self$tell <- function(sentence) tell(self, sentence)
  self$ask <- function(query) ask(self, query)
  self$retract <- function(sentence) retract(self, sentence)
  self$clauses <- function() if (length(self$clauses_)) ordered(unname(self$clauses_)) else list()
  class(self) <- c(cls, "KB")
  if (!is.null(sentence) && !isFALSE(sentence)) self$tell(sentence)
  self
}

#' @rdname entails
#' @export
KB <- function(sentence = NULL) {
  nie <- function(...) stop(.sympy_condition("NotImplementedError", ""))
  .kb_new("KB", sentence, nie, nie, nie)
}

#' @rdname entails
#' @export
PropKB <- function(sentence = NULL) {
  tell <- function(self, sentence) {
    for (c in .conjuncts(to_cnf(.as_bool(sentence)))) self$clauses_[[c$key]] <- c
    invisible(NULL)
  }
  ask <- function(self, query) entails(query, unname(self$clauses_))
  retract <- function(self, sentence) {
    for (c in .conjuncts(to_cnf(.as_bool(sentence)))) self$clauses_[[c$key]] <- NULL
    invisible(NULL)
  }
  .kb_new("PropKB", sentence, tell, ask, retract)
}

#' @export
print.KB <- function(x, ...) {
  cat("<", class(x)[1L], " with ", length(x$clauses_), " clauses>\n", sep = "")
  invisible(x)
}

# --- the DPLL algorithm on expressions ---------------------------------

#' The DPLL algorithm
#'
#' The building blocks of the Davis Putnam Logemann Loveland algorithm in
#' `sympy.logic.algorithms.dpll`: `dpll_satisfiable(expr)` returns a model
#' (a named list) or `FALSE`; `dpll(clauses, symbols, model)` searches from
#' a partial model; `find_pure_symbol()`, `find_unit_clause()` and
#' `unit_propagate()` are its steps. The `_int_repr` variants work on
#' clauses given as integer vectors (positive and negative literals), and
#' `to_int_repr(clauses, symbols)` encodes clauses so; `pl_true_int_repr()`
#' evaluates one such clause.
#'
#' @param expr A boolean expression.
#' @param clauses A list of clauses.
#' @param symbols A list of symbols, or an integer vector.
#' @param model A named list of logical values (keyed by symbol names, or by
#'   integers for the `_int_repr` functions).
#' @param unknown_clauses A list of clauses.
#' @param symbol A literal.
#' @param s An integer literal.
#' @param clause An integer clause.
#' @return A model, `FALSE`, a pair `list(symbol, value)` or a list of
#'   clauses.
#' @examples
#' A <- Symbol("A"); B <- Symbol("B")
#' dpll_satisfiable(And(A, Not(B)))
#' find_pure_symbol(list(A, B), list(Or(A, Not(B)), Or(Not(B), A)))
#' unit_propagate(list(Or(A, B), Or(Not(A), B)), A)
#' @export
dpll_satisfiable <- function(expr) {
  expr <- .as_bool(expr)
  clauses <- .conjuncts(to_cnf(expr))
  if (any(vapply(clauses, .is_false, TRUE))) return(FALSE)
  syms <- .bool_atoms(expr)
  syms <- .key_sort(lapply(syms, .sort_key), syms)
  cl <- to_int_repr(clauses, syms)
  result <- .dpll_int(cl, seq_along(syms), list())
  if (isFALSE(result) || !length(result)) return(result)
  out <- list()
  for (k in names(result)) out[[.str(syms[[as.integer(k)]])]] <- result[[k]]
  out
}

#' @rdname dpll_satisfiable
#' @export
dpll <- function(clauses, symbols, model) {
  clauses <- lapply(clauses, .as_bool)
  symbols <- lapply(symbols, as_sym)
  model <- as.list(model)
  rm_sym <- function(symbols, P) {
    i <- which(vapply(symbols, function(s) identical(s$key, P$key), TRUE))[1L]
    if (is.na(i)) stop(.sympy_condition("ValueError", "list.remove(x): x not in list"))
    symbols[-i]
  }
  for (finder in list(function() find_unit_clause(clauses, model), function() find_pure_symbol(symbols, clauses))) {
    pv <- finder()
    while (!is.null(pv[[1L]])) {
      P <- pv[[1L]]; value <- pv[[2L]]
      model <- .model_set(model, P, value)
      symbols <- rm_sym(symbols, P)
      if (!value) P <- Not(P)
      clauses <- unit_propagate(clauses, P)
      pv <- finder()
    }
  }
  unknown <- list()
  for (c in clauses) {
    val <- pl_true(c, model)
    if (isFALSE(val)) return(FALSE)
    if (!isTRUE(val)) unknown[[length(unknown) + 1L]] <- c
  }
  if (!length(unknown) || !length(clauses)) return(model)
  P <- symbols[[length(symbols)]]
  symbols <- symbols[-length(symbols)]
  model_copy <- .model_set(model, P, FALSE)
  model <- .model_set(model, P, TRUE)
  r <- dpll(unit_propagate(unknown, P), symbols, model)
  if (.sat_truthy(r)) return(r)
  dpll(unit_propagate(unknown, Not(P)), symbols, model_copy)
}

#' @rdname dpll_satisfiable
#' @export
find_pure_symbol <- function(symbols, unknown_clauses) {
  for (s in lapply(symbols, as_sym)) {
    pos <- FALSE; neg <- FALSE
    ns <- Not(s)
    for (c in unknown_clauses) {
      d <- .disjuncts(.as_bool(c))
      if (!pos && .in_args(s, d)) pos <- TRUE
      if (!neg && .in_args(ns, d)) neg <- TRUE
    }
    if (pos != neg) return(list(s, pos))
  }
  list(NULL, NULL)
}

#' @rdname dpll_satisfiable
#' @export
find_unit_clause <- function(clauses, model) {
  for (clause in clauses) {
    n <- 0L
    for (literal in .disjuncts(.as_bool(clause))) {
      s <- literal_symbol(literal)
      if (!.model_has(model, s)) {
        n <- n + 1L
        P <- s; value <- literal$type != "Not"
      }
    }
    if (n == 1L) return(list(P, value))
  }
  list(NULL, NULL)
}

#' @rdname dpll_satisfiable
#' @export
unit_propagate <- function(clauses, symbol) {
  symbol <- .as_bool(symbol)
  nsym <- Not(symbol)
  out <- list()
  for (c in clauses) {
    c <- .as_bool(c)
    if (c$type != "Or") { out[[length(out) + 1L]] <- c; next }
    args <- .bool_args(c)
    keep <- TRUE
    for (a in args) {
      if (identical(a$key, nsym$key)) {
        out[[length(out) + 1L]] <- do.call(Or, Filter(function(x) !identical(x$key, nsym$key), args))
        keep <- FALSE
        break
      }
      if (identical(a$key, symbol$key)) { keep <- FALSE; break }
    }
    if (keep) out[[length(out) + 1L]] <- c
  }
  out
}

#' @rdname dpll_satisfiable
#' @export
to_int_repr <- function(clauses, symbols) {
  keys <- vapply(symbols, function(s) as_sym(s)$key, "")
  lapply(clauses, function(c) {
    v <- vapply(.disjuncts(.as_bool(c)), function(a) {
      b <- if (a$type == "Not") a$args[[1L]] else a
      i <- match(b$key, keys)
      if (is.na(i)) stop(.sympy_condition("KeyError", sstr(b)))
      if (a$type == "Not") -i else i
    }, 0L)
    .pyset_iter(v)
  })
}

# int models: named lists keyed by the integer as a string
.im_get <- function(model, k) model[[as.character(k)]]

#' @rdname dpll_satisfiable
#' @export
pl_true_int_repr <- function(clause, model = list()) {
  result <- FALSE
  for (lit in clause) {
    p <- .im_get(model, abs(lit))
    if (!is.null(p) && lit < 0) p <- !p
    if (isTRUE(p)) return(TRUE)
    if (is.null(p)) result <- NA
  }
  result
}

#' @rdname dpll_satisfiable
#' @export
unit_propagate_int_repr <- function(clauses, s) {
  out <- list()
  for (cl in clauses) if (!(s %in% cl)) out[[length(out) + 1L]] <- .pyset_iter(cl[cl != -s])
  out
}

# the pure literal search, with symbols iterated in the given order
.find_pure_int <- function(symbols, unknown_clauses) {
  all_symbols <- .pyset_iter(unlist(unknown_clauses))
  found_pos <- .pyset_iter(symbols[symbols %in% all_symbols])
  negs <- -symbols
  found_neg <- .pyset_iter(negs[negs %in% all_symbols])
  for (p in found_pos) if (!(-p %in% found_neg)) return(list(p, TRUE))
  for (p in found_neg) if (!(-p %in% found_pos)) return(list(-p, FALSE))
  list(NULL, NULL)
}

#' @rdname dpll_satisfiable
#' @export
find_pure_symbol_int_repr <- function(symbols, unknown_clauses) {
  r <- .find_pure_int(as.integer(symbols), lapply(unknown_clauses, as.integer))
  if (!is.null(r[[1L]])) r[[1L]] <- as.integer(r[[1L]])
  r
}

#' @rdname dpll_satisfiable
#' @export
find_unit_clause_int_repr <- function(clauses, model) {
  bound <- as.integer(names(model))
  bound <- c(bound, -bound)
  for (cl in clauses) {
    unbound <- setdiff(as.integer(cl), bound)
    if (length(unbound) == 1L) {
      p <- unbound
      return(if (p < 0L) list(-p, FALSE) else list(p, TRUE))
    }
  }
  list(NULL, NULL)
}

#' @rdname dpll_satisfiable
#' @export
dpll_int_repr <- function(clauses, symbols, model) {
  .dpll_int(lapply(clauses, as.integer), sort(as.integer(symbols)), as.list(model))
}

# symbols is a set of positive integers, iterated (and popped) in increasing
# order as CPython does for small integers
.dpll_int <- function(clauses, symbols, model) {
  for (pure in c(FALSE, TRUE)) {
    pv <- if (pure) .find_pure_int(symbols, clauses) else find_unit_clause_int_repr(clauses, model)
    while (!is.null(pv[[1L]])) {
      P <- pv[[1L]]; value <- pv[[2L]]
      model[[as.character(P)]] <- value
      symbols <- symbols[symbols != P]
      if (!value) P <- -P
      clauses <- unit_propagate_int_repr(clauses, P)
      pv <- if (pure) .find_pure_int(symbols, clauses) else find_unit_clause_int_repr(clauses, model)
    }
  }
  unknown <- list()
  for (c in clauses) {
    val <- pl_true_int_repr(c, model)
    if (isFALSE(val)) return(FALSE)
    if (!isTRUE(val)) unknown[[length(unknown) + 1L]] <- c
  }
  if (!length(unknown)) return(model)
  P <- symbols[1L]
  symbols <- symbols[-1L]
  model_copy <- model
  model[[as.character(P)]] <- TRUE
  model_copy[[as.character(P)]] <- FALSE
  r <- .dpll_int(unit_propagate_int_repr(unknown, P), symbols, model)
  if (.sat_truthy(r)) return(r)
  .dpll_int(unit_propagate_int_repr(unknown, -P), symbols, model_copy)
}

#' Read DIMACS CNF
#'
#' `dimacs_load(s)` reads a boolean expression in the DIMACS CNF format
#' (variable `n` becomes the symbol `cnf_n`), and `dimacs_load_file(location)`
#' reads it from a file, as `load()` and `load_file()` in
#' `sympy.logic.utilities.dimacs`.
#'
#' @param s A string.
#' @param location A file path.
#' @return A boolean expression.
#' @examples
#' dimacs_load("1 2 \n 3")
#' @export
dimacs_load <- function(s) {
  clauses <- list()
  for (line in strsplit(s, "\n", fixed = TRUE)[[1L]]) {
    if (startsWith(line, "c")) next
    if (grepl("^p\\s*cnf\\s*(\\d*)\\s*(\\d*)", line, perl = TRUE)) next
    lits <- list()
    for (lit in strsplit(line, " ", fixed = TRUE)[[1L]]) {
      if (!nzchar(lit)) next
      v <- suppressWarnings(as.integer(trimws(lit)))
      if (is.na(v))
        stop(.sympy_condition("ValueError", sprintf("invalid literal for int() with base 10: '%s'", lit)))
      if (v == 0L) next
      s1 <- Symbol(paste0("cnf_", abs(v)))
      lits[[length(lits) + 1L]] <- if (v < 0L) Not(s1) else s1
    }
    if (length(lits)) clauses[[length(clauses) + 1L]] <- do.call(Or, lits)
  }
  do.call(And, clauses)
}

#' @rdname dimacs_load
#' @export
dimacs_load_file <- function(location) dimacs_load(paste(readLines(location, warn = FALSE), collapse = "\n"))
