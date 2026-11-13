# Boolean algebra algorithms, after sympy/logic/boolalg.py and
# sympy/logic/algorithms/dpll.py.

.bool_atoms <- function(e) {
  out <- list()
  walk <- function(x) {
    if (x$type %in% .BOOL_TYPES) {
      for (a in x$args) walk(a)
    } else if (x$type != "BooleanAtom") {
      out[[x$key]] <<- x
    }
  }
  walk(e)
  .sort_syms_any(unname(out))
}

.sort_syms_any <- function(lst) {
  if (length(lst) < 2L) return(lst)
  nm <- vapply(lst, function(s) if (is_symbol_like(s)) s$name else .str(s), "")
  lst[order(.natural_key(nm), method = "radix")]
}

#' Eliminate implications and convert to normal forms
#'
#' `eliminate_implications()` rewrites `Implies`, `Equivalent`, `Xor` and
#' `ITE` with `And`, `Or` and `Not`. `to_nnf()` gives negation normal form,
#' `to_cnf()` conjunctive normal form and `to_dnf()` disjunctive normal
#' form; with `simplify = TRUE` the result is minimized with
#' [simplify_logic()]. `is_cnf()`, `is_dnf()` and `is_nnf()` test the forms.
#'
#' @param e A boolean expression.
#' @param simplify For `to_cnf()` and `to_dnf()`, minimize the result; for
#'   `to_nnf()`, drop complementary arguments (the default).
#' @return A boolean expression or a logical value.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
#' to_cnf(Or(And(a, b), c))
#' to_dnf(And(Or(a, b), c))
#' @export
eliminate_implications <- function(e, form = NULL) .lg_to_nnf(as_sym(e), simplify = FALSE, form = form)

#' @rdname eliminate_implications
#' @param form For `eliminate_implications()` and `to_nnf()`: `"dnf"` writes
#'   a `Xor` as a disjunction, otherwise as a conjunction.
#' @export
to_nnf <- function(e, simplify = TRUE, form = NULL) .lg_to_nnf(as_sym(e), simplify, form)

# sympy.logic.boolalg.to_nnf: unchanged when already in negation normal form
.lg_to_nnf <- function(e, simplify = TRUE, form = NULL) {
  if (.is_nnf(e, simplify)) return(e)
  .lg_to_nnf_m(e, simplify, form)
}

# the to_nnf method of each boolean class
.lg_to_nnf_m <- function(e, simplify, form) {
  a <- e$args
  switch(e$type,
         And = , Or = .lattice_to_nnf(e$type, a, simplify, form),
         Not = .not_to_nnf(e, simplify, form),
         Implies = .lattice_to_nnf("Or", list(Not(a[[1L]]), a[[2L]]), simplify, form),
         Equivalent = {
           cl <- lapply(seq_along(a)[-length(a)], function(i) Or(Not(a[[i]]), a[[i + 1L]]))
           cl[[length(cl) + 1L]] <- Or(Not(a[[length(a)]]), a[[1L]])
           .lattice_to_nnf("And", cl, simplify, form)
         },
         Xor = if (identical(form, "dnf"))
           .lattice_to_nnf("Or", lapply(.odd_parity_terms(length(a)), function(m)
             do.call(And, lapply(seq_along(a), function(i) if (m[i] == 1L) a[[i]] else Not(a[[i]])))), simplify, form)
         else .lattice_to_nnf("And", .xor_clauses(a, seq(0L, length(a), by = 2L)), simplify, form),
         ITE = .lattice_to_nnf("And", list(Or(Not(a[[1L]]), a[[2L]]), Or(a[[1L]], a[[3L]])), simplify, form),
         e)
}

# Or(*clause) for each way of negating i of the arguments, i in counts
.xor_clauses <- function(a, counts) {
  out <- list()
  for (i in counts) {
    neg <- if (i == 0L) list(integer()) else if (i == length(a)) list(seq_along(a)) else
      lapply(seq_len(choose(length(a), i)), function(j) utils::combn(length(a), i)[, j])
    for (ng in neg) out[[length(out) + 1L]] <- do.call(Or, lapply(seq_along(a), function(k) if (k %in% ng) Not(a[[k]]) else a[[k]]))
  }
  out
}

# _get_odd_parity_terms(n): bit masks (most significant first) of odd weight
.odd_parity_terms <- function(n) {
  out <- list()
  for (v in seq_len(2^n) - 1L) {
    m <- as.integer(rev(as.integer(intToBits(v))[seq_len(n)]))
    if (sum(m) %% 2L == 1L) out[[length(out) + 1L]] <- m
  }
  out
}

# LatticeOp._to_nnf: with simplify, complementary arguments give cls.zero
.lattice_to_nnf <- function(cls, args, simplify, form) {
  argset <- list()
  keys <- character()
  for (arg in args) {
    if (!.is_literal(arg)) arg <- .lg_to_nnf_m(arg, simplify, form)
    parts <- if (simplify && arg$type == cls) arg$args else list(arg)
    for (a in parts) {
      if (simplify && Not(a)$key %in% keys) return(if (cls == "And") S_false else S_true)
      if (!(a$key %in% keys)) { argset[[length(argset) + 1L]] <- a; keys <- c(keys, a$key) }
    }
  }
  do.call(if (cls == "And") And else Or, argset)
}

# Not.to_nnf
.not_to_nnf <- function(e, simplify, form) {
  if (.is_literal(e)) return(e)
  x <- e$args[[1L]]
  a <- x$args
  switch(x$type,
         And = .lattice_to_nnf("Or", lapply(a, Not), simplify, form),
         Or = .lattice_to_nnf("And", lapply(a, Not), simplify, form),
         Implies = .lattice_to_nnf("And", list(a[[1L]], Not(a[[2L]])), simplify, form),
         Equivalent = .lattice_to_nnf("And", list(do.call(Or, a), do.call(Or, lapply(a, Not))), simplify, form),
         Xor = .lattice_to_nnf("And", .xor_clauses(a, seq(1L, length(a), by = 2L)), simplify, form),
         ITE = .lattice_to_nnf("And", list(Or(a[[1L]], Not(a[[3L]])), Or(Not(a[[1L]]), Not(a[[2L]]))), simplify, form),
         stop("Illegal operator ", x$type, " in expression", call. = FALSE))
}

# _distribute((expr, outer, inner))
.distribute <- function(e, outer, inner) {
  fo <- if (outer == "And") And else Or
  fi <- if (inner == "And") And else Or
  if (e$type == inner) {
    k <- which(vapply(e$args, function(a) a$type == outer, NA))
    if (!length(k)) return(e)
    conj <- e$args[[k[1L]]]
    rest <- do.call(fi, e$args[-k[1L]])
    return(do.call(fo, lapply(conj$args, function(c) .distribute(fi(c, rest), outer, inner))))
  }
  if (e$type == outer) return(do.call(fo, lapply(e$args, .distribute, outer = outer, inner = inner)))
  e
}

.BOOL_FUNCS <- c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE")

#' @rdname eliminate_implications
#' @export
to_cnf <- function(e, simplify = FALSE) {
  e <- as_sym(e)
  if (!(e$type %in% .BOOL_FUNCS)) return(e)
  if (simplify) return(simplify_logic(e, "cnf"))
  if (is_cnf(e)) return(e)
  .distribute(eliminate_implications(e, form = "cnf"), "And", "Or")
}

#' @rdname eliminate_implications
#' @export
to_dnf <- function(e, simplify = FALSE) {
  e <- as_sym(e)
  if (!(e$type %in% .BOOL_FUNCS)) return(e)
  if (simplify) return(simplify_logic(e, "dnf"))
  if (is_dnf(e)) return(e)
  .distribute(eliminate_implications(e, form = "dnf"), "Or", "And")
}

# sympy.logic.boolalg.is_literal
.is_literal <- function(e) {
  if (e$type == "Not") return(.is_literal(e$args[[1L]]))
  if (e$type == "BooleanAtom" || !length(e$args)) return(TRUE)
  !(e$type %in% .BOOL_FUNCS) && all(vapply(e$args, function(a) !length(a$args), NA))
}

# is_nnf(expr, simplified)
.is_nnf <- function(e, simplified = TRUE) {
  if (.is_literal(e)) return(TRUE)
  stack <- list(e)
  while (length(stack)) {
    x <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (x$type %in% c("And", "Or")) {
      if (simplified) {
        keys <- vapply(x$args, .k, "")
        for (a in x$args) if (Not(a)$key %in% keys) return(FALSE)
      }
      stack <- c(stack, x$args)
    } else if (!.is_literal(x)) return(FALSE)
  }
  TRUE
}

#' @rdname eliminate_implications
#' @export
is_nnf <- function(e) .is_nnf(as_sym(e))

#' @rdname eliminate_implications
#' @export
is_cnf <- function(e) {
  e <- as_sym(e)
  clause_ok <- function(c) .is_literal(c) || (c$type == "Or" && all(vapply(c$args, .is_literal, TRUE)))
  if (e$type == "And") return(all(vapply(e$args, clause_ok, TRUE)))
  clause_ok(e)
}

#' @rdname eliminate_implications
#' @export
is_dnf <- function(e) {
  e <- as_sym(e)
  term_ok <- function(c) .is_literal(c) || (c$type == "And" && all(vapply(c$args, .is_literal, TRUE)))
  if (e$type == "Or") return(all(vapply(e$args, term_ok, TRUE)))
  term_ok(e)
}

# Evaluate a boolean expression under an assignment (named logical vector
# keyed by atom keys).
.beval <- function(e, env) {
  switch(e$type,
         BooleanAtom = e$name == "true",
         And = all(vapply(e$args, .beval, TRUE, env = env)),
         Or = any(vapply(e$args, .beval, TRUE, env = env)),
         Not = !.beval(e$args[[1L]], env),
         Xor = sum(vapply(e$args, .beval, TRUE, env = env)) %% 2L == 1L,
         Implies = !.beval(e$args[[1L]], env) || .beval(e$args[[2L]], env),
         Equivalent = length(unique(vapply(e$args, .beval, TRUE, env = env))) == 1L,
         ITE = if (.beval(e$args[[1L]], env)) .beval(e$args[[2L]], env) else .beval(e$args[[3L]], env),
         {
           v <- env[[e$key]]
           if (is.null(v)) stop("no value for ", .str(e))
           v
         })
}

#' Truth table of a boolean expression
#'
#' @param e A boolean expression.
#' @param vars Variables, by default all atoms in `e`.
#' @return A data frame with one logical column per variable and a `value`
#'   column.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b")
#' truth_table(Implies(a, b))
#' @export
truth_table <- function(e, vars = NULL) {
  e <- as_sym(e)
  vars <- if (is.null(vars)) .bool_atoms(e) else lapply(vars, as_sym)
  n <- length(vars)
  rows <- if (n) expand.grid(rep(list(c(FALSE, TRUE)), n))[, rev(seq_len(n)), drop = FALSE] else
    data.frame(matrix(nrow = 1L, ncol = 0L))
  vals <- logical(nrow(rows))
  keys <- vapply(vars, .k, "")
  for (r in seq_len(nrow(rows))) {
    env <- if (n) as.list(stats::setNames(unlist(rows[r, ]), keys)) else list()
    vals[r] <- .beval(e, env)
  }
  df <- as.data.frame(rows)
  names(df) <- vapply(vars, .str, "")
  df$value <- vals
  df
}

# Quine-McCluskey after sympy/logic/boolalg.py: _input_to_binlist,
# _simplified_pairs, _rem_redundancy and the SOP and POS conversions.
.ibin <- function(v, bits) as.integer(rev(as.integer(intToBits(v))[seq_len(bits)]))

.input_to_binlist <- function(inputlist, variables) {
  bits <- length(variables)
  names_ <- vapply(variables, function(v) .str(v), "")
  out <- list()
  for (val in as.list(inputlist)) {
    if (is.list(val) && !is.null(names(val))) {
      nonspec <- setdiff(names_, names(val))
      grid <- if (length(nonspec)) expand.grid(rep(list(0:1), length(nonspec)))[, rev(seq_along(nonspec)), drop = FALSE] else data.frame(row = 1)
      for (r in seq_len(nrow(grid))) {
        d <- stats::setNames(as.list(unlist(grid[r, seq_along(nonspec)])), nonspec)
        d[names(val)] <- lapply(val, as.integer)
        out[[length(out) + 1L]] <- vapply(names_, function(nm) as.integer(d[[nm]]), 0L, USE.NAMES = FALSE)
      }
    } else if ((is.numeric(val) || is.logical(val)) && length(val) == 1L) {
      out[[length(out) + 1L]] <- .ibin(as.integer(val), bits)
    } else if (is.numeric(val) || is.logical(val)) {
      if (length(val) != bits)
        stop(.sympy_condition("ValueError", sprintf("Each term must contain %d bits as there are\n%d variables (or be an integer).", bits, bits)))
      out[[length(out) + 1L]] <- as.integer(val)
    } else stop("A term list can only contain lists, ints or dicts.", call. = FALSE)
  }
  out
}

.check_pair <- function(a, b) {
  d <- which(a != b)
  if (length(d) == 1L) d else -1L
}

.simplified_pairs <- function(terms) {
  if (!length(terms)) return(list())
  simplified <- list()
  skeys <- character()
  todo <- rep(TRUE, length(terms))
  ones <- vapply(terms, function(t) sum(t == 1L), 0L)
  for (k in seq_len(length(terms[[1L]])) - 1L) for (i in which(ones == k)) for (j in which(ones == k + 1L)) {
    idx <- .check_pair(terms[[i]], terms[[j]])
    if (idx != -1L) {
      todo[i] <- FALSE; todo[j] <- FALSE
      nt <- terms[[i]]
      nt[idx] <- 3L
      key <- paste(nt, collapse = ",")
      if (!(key %in% skeys)) { simplified[[length(simplified) + 1L]] <- nt; skeys <- c(skeys, key) }
    }
  }
  if (length(simplified)) simplified <- .simplified_pairs(simplified)
  c(simplified, terms[todo])
}

.rem_redundancy <- function(l1, terms) {
  if (!length(terms)) return(list())
  nt <- length(terms); nl <- length(l1)
  dom <- matrix(0L, nt, nl)
  colcount <- integer(nl); rowcount <- integer(nt)
  for (pi in seq_len(nl)) for (ti in seq_len(nt)) if (all(l1[[pi]] == 3L | l1[[pi]] == terms[[ti]])) {
    dom[ti, pi] <- 1L
    colcount[pi] <- colcount[pi] + 1L
    rowcount[ti] <- rowcount[ti] + 1L
  }
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (ri in seq_len(nt)) if (rowcount[ri]) {
      row <- dom[ri, ]
      for (r2 in seq_len(nt)) if (ri != r2 && rowcount[ri] && rowcount[ri] <= rowcount[r2]) {
        row2 <- dom[r2, ]
        if (all(row2 >= row)) {
          rowcount[r2] <- 0L
          changed <- TRUE
          for (pi in which(row2 != 0L)) { dom[r2, pi] <- 0L; colcount[pi] <- colcount[pi] - 1L }
        }
      }
    }
    colcache <- vector("list", nl)
    getcol <- function(ci) { if (is.null(colcache[[ci]])) colcache[[ci]] <<- dom[, ci]; colcache[[ci]] }
    for (ci in seq_len(nl)) if (colcount[ci]) {
      col <- getcol(ci)
      for (c2 in seq_len(nl)) if (ci != c2 && colcount[c2] && colcount[ci] >= colcount[c2]) {
        col2 <- getcol(c2)
        if (all(col >= col2)) {
          colcount[c2] <- 0L
          changed <- TRUE
          for (ti in seq_len(nt)) if (col2[ti] && dom[ti, c2]) { dom[ti, c2] <- 0L; rowcount[ti] <- rowcount[ti] - 1L }
        }
      }
    }
    if (!changed) {
      best <- -1L; mx <- 0L
      for (ci in seq_len(nl)) if (colcount[ci] > mx) { best <- ci; mx <- colcount[ci] }
      if (best != -1L && mx > 1L) {
        bc <- colcache[[best]]
        for (pi in seq_len(nl)) if (pi != best) for (ti in seq_len(nt)) if (bc[ti] && dom[ti, pi]) {
          dom[ti, pi] <- 0L
          changed <- TRUE
          rowcount[ti] <- rowcount[ti] - 1L
          colcount[pi] <- colcount[pi] - 1L
        }
      }
    }
  }
  l1[colcount != 0L]
}

.check_dontcares <- function(minterms, dontcares) {
  mk <- vapply(minterms, paste, "", collapse = ",")
  for (d in dontcares) if (paste(d, collapse = ",") %in% mk)
    stop(.sympy_condition("ValueError", paste0("[", paste(d, collapse = ", "), "] in minterms is also in dontcares")))
}

#' Sum of products and product of sums forms
#'
#' `SOPform(vars, minterms)` and `POSform(vars, minterms)` build minimal
#' disjunctive and conjunctive forms from the rows where the function is
#' true, with SymPy's Quine-McCluskey reduction. A row is an integer (most
#' significant variable first), a vector of bits, or a named list giving
#' the values of some variables (the others take every value).
#'
#' @param vars A list of symbols.
#' @param minterms Rows where the function is true: an integer vector, or a
#'   list of integers, bit vectors and named lists.
#' @param dontcares Rows whose value does not matter, in the same forms.
#' @return A boolean expression.
#' @examples
#' w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' SOPform(list(w, x, y, z), c(1, 3, 7, 11, 15), dontcares = c(0, 2, 5))
#' @export
SOPform <- function(vars, minterms, dontcares = integer()) {
  if (!length(minterms)) return(S_false)
  vars <- lapply(vars, as_sym)
  minterms <- .input_to_binlist(minterms, vars)
  dontcares <- .input_to_binlist(dontcares, vars)
  .check_dontcares(minterms, dontcares)
  .sop_form(vars, minterms, dontcares)
}

.sop_form <- function(vars, minterms, dontcares) {
  essential <- .rem_redundancy(.simplified_pairs(c(minterms, dontcares)), minterms)
  do.call(Or, lapply(essential, function(p)
    do.call(And, lapply(which(p != 3L), function(i) if (p[i] == 1L) vars[[i]] else Not(vars[[i]])))))
}

#' @rdname SOPform
#' @export
POSform <- function(vars, minterms, dontcares = integer()) {
  if (!length(minterms)) return(S_false)
  vars <- lapply(vars, as_sym)
  minterms <- .input_to_binlist(minterms, vars)
  dontcares <- .input_to_binlist(dontcares, vars)
  .check_dontcares(minterms, dontcares)
  used <- vapply(c(minterms, dontcares), paste, "", collapse = ",")
  maxterms <- Filter(function(t) !(paste(t, collapse = ",") %in% used),
                     lapply(seq_len(2^length(vars)) - 1L, .ibin, bits = length(vars)))
  essential <- .rem_redundancy(.simplified_pairs(c(maxterms, dontcares)), maxterms)
  do.call(And, lapply(essential, function(p)
    do.call(Or, lapply(which(p != 3L), function(i) if (p[i] == 0L) vars[[i]] else Not(vars[[i]])))))
}

#' Simplify boolean expressions
#'
#' Minimizes a boolean expression into the shortest equivalent `"cnf"` or
#' `"dnf"` form with the Quine-McCluskey algorithm (up to `max_vars`
#' variables), like SymPy's `simplify_logic()`.
#'
#' @param e A boolean expression.
#' @param form `"cnf"`, `"dnf"`, or `NULL` to pick the shorter.
#' @param max_vars Largest number of variables to minimize.
#' @return A boolean expression.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b")
#' simplify_logic(Or(And(a, b), And(a, Not(b))))
#' @export
simplify_logic <- function(e, form = NULL, max_vars = 10) {
  e <- as_sym(e)
  vars <- .bool_atoms(e)
  if (length(vars) > max_vars) return(e)
  tt <- truth_table(e, vars)
  minterms <- which(tt$value) - 1L
  dnf <- SOPform(vars, minterms)
  cnf <- POSform(vars, minterms)
  if (identical(form, "dnf")) return(dnf)
  if (identical(form, "cnf")) return(cnf)
  if (.cops(cnf) < .cops(dnf)) cnf else dnf
}

# CNF as list of clauses, each a named integer vector (+1 / -1 by atom key).
.cnf_clauses <- function(e) {
  c <- to_cnf(e)
  if (.is_true(c)) return(list())
  if (.is_false(c)) return(list(integer()))
  clauses <- if (c$type == "And") c$args else list(c)
  lapply(clauses, function(cl) {
    lits <- if (cl$type == "Or") cl$args else list(cl)
    v <- vapply(lits, function(l) if (l$type == "Not") -1L else 1L, 0L)
    names(v) <- vapply(lits, function(l) if (l$type == "Not") l$args[[1L]]$key else l$key, "")
    v
  })
}

#' Satisfiability
#'
#' `satisfiable(e)` returns a model (a named list of logical values) making
#' `e` true, or `FALSE`, as SymPy's `satisfiable()`. `algorithm` picks the
#' method: `"dpll2"` (the default; `"pycosat"`, `"minisat22"` and `"z3"`
#' fall back to it as in SymPy when those libraries are absent, except that
#' an explicit `"pycosat"` is an `ImportError`) or `"dpll"` (see
#' [dpll_satisfiable()]). With `all_models = TRUE` the result is the list
#' of all models, or `list(FALSE)`. A tautology without symbols has the
#' model `list(True = TRUE)`. `valid()` tests whether `e` is a tautology.
#'
#' @param e A boolean expression.
#' @param algorithm `NULL` or the name of an algorithm.
#' @param all_models Return all models.
#' @param minimal Ask for minimal models (used only by `"minisat22"`).
#' @param use_lra_theory Use the theory of linear real arithmetic, which is
#'   not available here.
#' @return A named list, a list of them, or `FALSE`.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b")
#' satisfiable(And(a, Not(b)))
#' satisfiable(And(a, Not(a)))
#' satisfiable(Equivalent(a, b), all_models = TRUE)
#' @export
satisfiable <- function(e, algorithm = NULL, all_models = FALSE, minimal = FALSE, use_lra_theory = FALSE) {
  if (use_lra_theory) {
    if (!is.null(algorithm) && algorithm != "dpll2")
      stop(.sympy_condition("ValueError", sprintf("Currently only dpll2 can handle using lra theory. %s is not handled.", algorithm)))
    stop(.sympy_condition("NotImplementedError", "the theory of linear real arithmetic is not available"))
  }
  if (is.null(algorithm)) algorithm <- "dpll2"
  else if (algorithm == "pycosat") stop(.sympy_condition("ImportError", "pycosat module is not present"))
  if (algorithm %in% c("minisat22", "z3")) algorithm <- "dpll2"
  if (algorithm == "dpll") return(dpll_satisfiable(e))
  if (algorithm != "dpll2") stop(.sympy_condition("NotImplementedError", ""))
  e <- .as_bool(e)
  if (.is_false(e)) return(if (all_models) list(FALSE) else FALSE)
  if (.is_true(e)) return(if (all_models) list(list(True = TRUE)) else list(True = TRUE))
  atoms_ <- .bool_atoms(e)
  keymap <- stats::setNames(atoms_, vapply(atoms_, .k, ""))
  pretty <- function(m) {
    out <- lapply(names(keymap), function(k) if (is.null(m[[k]])) TRUE else m[[k]])
    names(out) <- vapply(keymap, .str, "")
    out
  }
  clauses <- .cnf_clauses(e)
  if (!length(clauses)) return(if (all_models) list(pretty(list())) else pretty(list()))
  if (any(vapply(clauses, length, 0L) == 0L)) return(if (all_models) list(FALSE) else FALSE)
  keys <- names(keymap)[names(keymap) %in% unlist(lapply(clauses, names))]
  ints <- lapply(clauses, function(cl) as.integer(cl * match(names(cl), keys)))
  models <- .dpll2_models(ints, length(keys), all_models)
  # an unsatisfiable formula has the single model False, as SymPy's generator
  if (!length(models)) return(if (all_models) list(FALSE) else FALSE)
  models <- lapply(models, function(m) {
    m <- unique(m)
    out <- lapply(seq_along(keys), function(v) v %in% m)
    names(out) <- vapply(keymap[keys], .str, "")
    out
  })
  if (all_models) models else models[[1L]]
}

#' @rdname satisfiable
#' @export
valid <- function(e) identical(satisfiable(Not(as_sym(e))), FALSE)

# StrPrinter._print_Xor: p ^ q at the BitwiseXor precedence
register_str_printer("Xor", function(e)
  paste(vapply(.bool_args(e), function(a) .paren(a, PREC$BitwiseXor), ""), collapse = " ^ "))
register_str_printer("Implies", function(e) paste0("Implies(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))
register_str_printer("Equivalent", function(e) paste0("Equivalent(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
register_str_printer("ITE", function(e) paste0("ITE(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
