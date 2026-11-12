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
#' @param simplify Minimize the result.
#' @return A boolean expression or a logical value.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
#' to_cnf(Or(And(a, b), c))
#' to_dnf(And(Or(a, b), c))
#' @export
eliminate_implications <- function(e) {
  e <- as_sym(e)
  .map_bottom_up(e, function(z) {
    a <- z$args
    switch(z$type,
           Implies = Or(Not(a[[1L]]), a[[2L]]),
           Equivalent = Or(do.call(And, a), do.call(And, lapply(a, Not))),
           Xor = Reduce(function(p, q) Or(And(p, Not(q)), And(Not(p), q)), a),
           ITE = Or(And(a[[1L]], a[[2L]]), And(Not(a[[1L]]), a[[3L]])),
           z)
  })
}

#' @rdname eliminate_implications
#' @export
to_nnf <- function(e) {
  e <- eliminate_implications(as_sym(e))
  .nnf(e)
}

.nnf <- function(e) {
  if (e$type == "Not") {
    a <- e$args[[1L]]
    if (a$type == "And") return(do.call(Or, lapply(a$args, function(z) .nnf(Not(z)))))
    if (a$type == "Or") return(do.call(And, lapply(a$args, function(z) .nnf(Not(z)))))
    if (a$type == "Not") return(.nnf(a$args[[1L]]))
    return(e)
  }
  if (e$type %in% c("And", "Or")) return(sym_rebuild(e, lapply(e$args, .nnf)))
  e
}

.distribute <- function(e, outer, inner) {
  if (e$type == outer) return(sym_rebuild(e, lapply(e$args, .distribute, outer = outer, inner = inner)))
  if (e$type == inner) {
    args <- lapply(e$args, .distribute, outer = outer, inner = inner)
    w <- which(vapply(args, function(z) z$type == outer, TRUE))
    if (length(w)) {
      k <- w[1L]
      rest <- args[-k]
      return(.distribute(do.call(get(outer), lapply(args[[k]]$args, function(z) do.call(get(inner), c(list(z), rest)))),
                         outer, inner))
    }
    return(do.call(get(inner), args))
  }
  e
}

#' @rdname eliminate_implications
#' @export
to_cnf <- function(e, simplify = FALSE) {
  e <- as_sym(e)
  if (simplify) return(simplify_logic(e, "cnf"))
  .distribute(to_nnf(e), "And", "Or")
}

#' @rdname eliminate_implications
#' @export
to_dnf <- function(e, simplify = FALSE) {
  e <- as_sym(e)
  if (simplify) return(simplify_logic(e, "dnf"))
  .distribute(to_nnf(e), "Or", "And")
}

.is_literal <- function(e) e$type != "And" && e$type != "Or" &&
  (e$type != "Not" || !(e$args[[1L]]$type %in% .BOOL_TYPES)) &&
  !(e$type %in% c("Implies", "Equivalent", "Xor", "ITE"))

#' @rdname eliminate_implications
#' @export
is_nnf <- function(e) {
  e <- as_sym(e)
  if (.is_literal(e)) return(TRUE)
  if (e$type %in% c("And", "Or")) return(all(vapply(e$args, is_nnf, TRUE)))
  FALSE
}

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

# Quine-McCluskey prime implicants and a greedy cover.
.qm <- function(minterms, n, dontcares = integer()) {
  if (!length(minterms)) return(list())
  terms <- lapply(unique(c(minterms, dontcares)), function(m) as.integer(rev(as.integer(intToBits(m))[seq_len(n)])))
  primes <- list()
  cur <- terms
  repeat {
    used <- rep(FALSE, length(cur))
    nxt <- list()
    keys <- character()
    if (length(cur) > 1L) for (i in seq_len(length(cur) - 1L)) for (j in (i + 1L):length(cur)) {
      d <- which(cur[[i]] != cur[[j]])
      if (length(d) == 1L && cur[[i]][d] != 2L && cur[[j]][d] != 2L) {
        t <- cur[[i]]
        t[d] <- 2L
        k <- paste(t, collapse = "")
        if (!(k %in% keys)) {
          nxt[[length(nxt) + 1L]] <- t
          keys <- c(keys, k)
        }
        used[i] <- TRUE
        used[j] <- TRUE
      }
    }
    primes <- c(primes, cur[!used])
    if (!length(nxt)) break
    cur <- nxt
  }
  primes <- primes[!duplicated(vapply(primes, paste, "", collapse = ""))]
  covers <- function(p, m) {
    bits <- as.integer(rev(as.integer(intToBits(m))[seq_len(n)]))
    all(p == 2L | p == bits)
  }
  left <- minterms
  chosen <- list()
  # essential primes first, then greedy
  repeat {
    if (!length(left)) break
    ess <- NULL
    for (m in left) {
      cov <- which(vapply(primes, covers, TRUE, m = m))
      if (length(cov) == 1L) {
        ess <- cov
        break
      }
    }
    if (is.null(ess)) {
      cnt <- vapply(primes, function(p) sum(vapply(left, function(m) covers(p, m), TRUE)), 0L)
      lits <- vapply(primes, function(p) sum(p != 2L), 0L)
      ess <- order(-cnt, lits)[1L]
    }
    p <- primes[[ess]]
    chosen[[length(chosen) + 1L]] <- p
    left <- left[!vapply(left, function(m) covers(p, m), TRUE)]
  }
  chosen
}

#' Sum of products and product of sums forms
#'
#' `SOPform(vars, minterms)` and `POSform(vars, minterms)` build minimal
#' disjunctive and conjunctive forms from the rows (as integers, most
#' significant variable first) where the function is true.
#'
#' @param vars A list of symbols.
#' @param minterms Integer row indices where the function is true.
#' @param dontcares Rows whose value does not matter.
#' @return A boolean expression.
#' @examples
#' w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' SOPform(list(w, x, y, z), c(1, 3, 7, 11, 15), dontcares = c(0, 2, 5))
#' @export
SOPform <- function(vars, minterms, dontcares = integer()) {
  vars <- lapply(vars, as_sym)
  n <- length(vars)
  imps <- .qm(as.integer(minterms), n, as.integer(dontcares))
  if (!length(imps)) return(S_false)
  do.call(Or, lapply(imps, function(p) {
    lits <- lapply(which(p != 2L), function(i) if (p[i] == 1L) vars[[i]] else Not(vars[[i]]))
    if (!length(lits)) S_true else do.call(And, lits)
  }))
}

#' @rdname SOPform
#' @export
POSform <- function(vars, minterms, dontcares = integer()) {
  vars <- lapply(vars, as_sym)
  n <- length(vars)
  maxterms <- setdiff(0:(2^n - 1), c(minterms, dontcares))
  imps <- .qm(as.integer(maxterms), n, as.integer(dontcares))
  if (!length(imps)) return(S_true)
  do.call(And, lapply(imps, function(p) {
    lits <- lapply(which(p != 2L), function(i) if (p[i] == 1L) Not(vars[[i]]) else vars[[i]])
    if (!length(lits)) S_false else do.call(Or, lits)
  }))
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

.dpll <- function(clauses, model) {
  repeat {
    clauses <- Filter(function(cl) {
      for (k in names(cl)) {
        mv <- model[[k]]
        if (!is.null(mv) && ((cl[[k]] > 0) == mv)) return(FALSE)
      }
      TRUE
    }, clauses)
    clauses <- lapply(clauses, function(cl) cl[vapply(names(cl), function(k) is.null(model[[k]]), TRUE)])
    if (!length(clauses)) return(model)
    if (any(vapply(clauses, length, 0L) == 0L)) return(NULL)
    units <- Filter(function(cl) length(cl) == 1L, clauses)
    if (!length(units)) break
    for (u in units) model[[names(u)]] <- u[[1L]] > 0
  }
  # pure literals
  allv <- unlist(clauses)
  for (k in unique(names(allv))) {
    s <- unique(allv[names(allv) == k])
    if (length(s) == 1L) model[[k]] <- s > 0
  }
  k <- names(clauses[[1L]])[1L]
  if (!is.null(model[[k]])) return(.dpll(clauses, model))
  for (val in c(TRUE, FALSE)) {
    m2 <- model
    m2[[k]] <- val
    r <- .dpll(clauses, m2)
    if (!is.null(r)) return(r)
  }
  NULL
}

#' Satisfiability
#'
#' `satisfiable(e)` returns a model (a named logical list) making `e` true,
#' or `FALSE`; with `all_models = TRUE` it returns every model. It uses the
#' DPLL algorithm on the CNF form, as SymPy does. `valid()` tests whether
#' `e` is a tautology.
#'
#' @param e A boolean expression.
#' @param all_models Return all models.
#' @return A named list, a list of them, or `FALSE`.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b")
#' satisfiable(And(a, Not(b)))
#' satisfiable(And(a, Not(a)))
#' @export
satisfiable <- function(e, all_models = FALSE) {
  e <- as_sym(e)
  atoms_ <- .bool_atoms(e)
  keymap <- stats::setNames(atoms_, vapply(atoms_, .k, ""))
  pretty <- function(m) {
    out <- lapply(names(keymap), function(k) if (is.null(m[[k]])) TRUE else m[[k]])
    names(out) <- vapply(keymap, .str, "")
    out
  }
  if (all_models) {
    tt <- truth_table(e, atoms_)
    rows <- which(tt$value)
    if (!length(rows)) return(FALSE)
    return(lapply(rows, function(r) as.list(tt[r, seq_along(atoms_), drop = FALSE])))
  }
  m <- .dpll(.cnf_clauses(e), list())
  if (is.null(m)) return(FALSE)
  pretty(m)
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
