# Conjunctive normal forms for the assumptions system, after
# sympy/assumptions/cnf.py, and a small SAT solver for the encoded clauses,
# which answers what sympy/logic/algorithms/dpll2.py answers for ask(),
# satask() and the ReasoningEngine of sympy/assumptions/reasoning_engine.py.
#
# A literal is list(lit, neg) with class "ask_literal"; its key is the key
# of lit, prefixed by "~" when negated. A clause is a named list of literals
# (names are the literal keys), and a CNF is list(clauses) whose clauses are
# named by their sorted literal keys, so the clause set has no duplicates.

.lit_key <- function(l) if (l$neg) paste0("~", l$lit$key) else l$lit$key

# .Literal(lit, is_Not)
.Literal <- function(lit, is_Not = FALSE) {
  if (inherits(lit, "sym_predicate_fn")) lit <- attr(lit, "predicate")
  if (!inherits(lit, "ask_nnf")) lit <- .ask_sympify(lit)
  if (inherits(lit, "ask_nnf")) return(if (is_Not) .nnf_not(lit) else lit)
  if (lit$type == "Not") {
    lit <- lit$args[[1L]]
    is_Not <- TRUE
  }
  structure(list(lit = lit, neg = isTRUE(is_Not)), class = c("ask_literal", "ask_nnf"))
}

#' @export
format.ask_literal <- function(x, ...) paste0("Literal(", .str(x$lit), ", ", if (x$neg) "True" else "False", ")")

#' @export
print.ask_literal <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

.nnf_and <- function(args) structure(list(kind = "and", args = args), class = "ask_nnf")
.nnf_or <- function(args) structure(list(kind = "or", args = args), class = "ask_nnf")

.nnf_not <- function(x) {
  if (inherits(x, "ask_literal")) return(structure(list(lit = x$lit, neg = !x$neg), class = c("ask_literal", "ask_nnf")))
  if (x$kind == "and") .nnf_or(lapply(x$args, .nnf_not)) else .nnf_and(lapply(x$args, .nnf_not))
}

# to_NNF
.to_nnf <- function(expr) {
  if (inherits(expr, "ask_nnf")) return(expr)
  if (inherits(expr, "sym_predicate_fn")) expr <- attr(expr, "predicate")
  expr <- .ask_sympify(expr)
  if (.is_true(expr)) return(.nnf_and(list()))
  if (.is_false(expr)) return(.nnf_or(list()))
  if (is_relational(expr)) expr <- .rel_to_binrel(expr)
  a <- expr$args
  switch(expr$type,
    Not = .nnf_not(.to_nnf(a[[1L]])),
    Or = .nnf_or(lapply(a, .to_nnf)),
    And = .nnf_and(lapply(a, .to_nnf)),
    Nand = .nnf_not(.nnf_and(lapply(a, .to_nnf))),
    Nor = .nnf_not(.nnf_or(lapply(a, .to_nnf))),
    Xor = , Xnor = {
      n <- length(a)
      cnfs <- list()
      for (i in seq(0L, n, by = 2L)) {
        combs <- if (i == 0L) list(integer()) else utils::combn(n, i, simplify = FALSE)
        for (neg in combs) {
          clause <- lapply(seq_len(n), function(j) if (j %in% neg) .nnf_not(.to_nnf(a[[j]])) else .to_nnf(a[[j]]))
          cnfs[[length(cnfs) + 1L]] <- .nnf_or(clause)
        }
      }
      if (expr$type == "Xor") .nnf_and(cnfs) else .nnf_not(.nnf_and(cnfs))
    },
    Implies = .nnf_or(list(.nnf_not(.to_nnf(a[[1L]])), .to_nnf(a[[2L]]))),
    Equivalent = {
      n <- length(a)
      .nnf_and(lapply(seq_len(n), function(i) {
        b <- if (i < n) a[[i + 1L]] else a[[1L]]
        .nnf_or(list(.nnf_not(.to_nnf(a[[i]])), .to_nnf(b)))
      }))
    },
    ITE = {
      L <- .to_nnf(a[[1L]]); M <- .to_nnf(a[[2L]]); R <- .to_nnf(a[[3L]])
      .nnf_and(list(.nnf_or(list(.nnf_not(L), M)), .nnf_or(list(L, R))))
    },
    .Literal(expr))
}

.cnf <- function(clauses = list()) structure(list(clauses = clauses), class = "ask_cnf")
.clause_id <- function(cl) paste(sort(as.character(names(cl)), method = "radix"), collapse = "\r")
.clause <- function(lits) {
  names(lits) <- vapply(lits, .lit_key, "")
  lits[!duplicated(names(lits))]
}
.named_clauses <- function(cls) {
  names(cls) <- vapply(cls, .clause_id, "")
  cls[!duplicated(names(cls))]
}

# distribute_AND_over_OR
.cnf_distribute <- function(expr) {
  if (inherits(expr, "ask_literal")) return(.cnf(.named_clauses(list(.clause(list(expr))))))
  if (expr$kind == "or") return(.cnf_all_or(lapply(expr$args, .cnf_distribute)))
  .cnf_all_and(lapply(expr$args, .cnf_distribute))
}

.cnf_or <- function(a, b) {
  out <- list()
  for (x in .subset2(a, "clauses")) for (y in .subset2(b, "clauses")) {
    z <- c(x, y)
    out[[length(out) + 1L]] <- z[!duplicated(names(z))]
  }
  .cnf(.named_clauses(out))
}

.cnf_and <- function(a, b) {
  cl <- c(.subset2(a, "clauses"), .subset2(b, "clauses"))
  .cnf(cl[!duplicated(names(cl))])
}

.cnf_all_or <- function(cnfs) {
  if (!length(cnfs)) {
    cl <- list(structure(list(), names = character()))
    names(cl) <- ""
    return(.cnf(cl))
  }
  b <- cnfs[[1L]]
  for (r in cnfs[-1L]) b <- .cnf_or(b, r)
  b
}

.cnf_all_and <- function(cnfs) {
  cl <- unlist(lapply(cnfs, `[[`, "clauses"), recursive = FALSE)
  if (is.null(cl)) cl <- list()
  .cnf(cl[!duplicated(names(cl))])
}

# .CNF(), CNF.to_CNF, CNF.from_prop, CNF.add_clauses and CNF.all_predicates
.CNF <- function(clauses = NULL) if (is.null(clauses)) .cnf() else .cnf(clauses)
.CNF_to_CNF <- function(expr) .cnf_distribute(.to_nnf(expr))
.CNF_from_prop <- function(prop) .CNF_to_CNF(prop)
.cnf_add_clauses <- function(cnf, clauses) {
  if (inherits(clauses, "ask_cnf")) clauses <- .subset2(clauses, "clauses")
  else clauses <- .named_clauses(lapply(clauses, function(cl) .clause(unclass(cl))))
  .cnf_and(cnf, .cnf(clauses))
}
.cnf_all_predicates <- function(cnf) {
  lits <- unlist(lapply(.subset2(cnf, "clauses"), function(cl) lapply(cl, `[[`, "lit")), recursive = FALSE)
  if (is.null(lits)) return(list())
  lits[!duplicated(vapply(lits, .k, ""))]
}

#' @export
`$.ask_cnf` <- function(x, name) {
  v <- .subset2(x, name)
  if (identical(name, "clauses") && !is.null(v)) {
    v <- lapply(unname(v), function(cl) structure(unname(cl), pyset = TRUE))
    attr(v, "pyset") <- TRUE
  }
  v
}

#' @export
format.ask_cnf <- function(x, ...) {
  paste(vapply(.subset2(x, "clauses"), function(cl)
    paste0("(", paste(vapply(cl, function(l) paste0(if (l$neg) "~", .str(l$lit)), ""), collapse = " | "), ")"), ""),
    collapse = " & ")
}

# EncodedCNF: symbols are literal syms, encoding maps keys to integers, and
# data is a list of integer vectors (0 stands for the literal False).
.enc_cnf <- function() {
  e <- new.env(parent = emptyenv())
  e$symbols <- list()
  e$encoding <- new.env(parent = emptyenv())
  e$data <- list()
  e
}

.enc_arg <- function(enc, l) {
  k <- l$lit$key
  v <- enc$encoding[[k]]
  if (is.null(v)) {
    n <- length(enc$symbols)
    enc$symbols[[n + 1L]] <- l$lit
    v <- n + 1L
    assign(k, v, envir = enc$encoding)
  }
  if (l$neg) -v else v
}

.enc_clause <- function(enc, cl) {
  unique(vapply(unname(cl), function(l) if (identical(l$lit$key, S_false$key)) 0L else .enc_arg(enc, l), 0L))
}

.enc_from_cnf <- function(enc, cnf) {
  for (l in .cnf_all_predicates(cnf)) .enc_arg(enc, list(lit = l, neg = FALSE))
  enc$data <- lapply(unname(.subset2(cnf, "clauses")), function(cl) .enc_clause(enc, cl))
  invisible(enc)
}

.enc_add_from_cnf <- function(enc, cnf) {
  enc$data <- c(enc$data, lapply(unname(.subset2(cnf, "clauses")), function(cl) .enc_clause(enc, cl)))
  invisible(enc)
}

# ---------------------------------------------------------------------------
# SAT: unit propagation and a DPLL search on integer clauses.

# Unit propagation from the assignment val (1 true, -1 false, 0 free).
# Returns the new assignment, or NULL on a conflict.
.sat_propagate <- function(lits, cid, ncl, val) {
  repeat {
    lv <- sign(lits) * val[abs(lits)]
    sat <- tabulate(cid[lv == 1L], nbins = ncl) > 0L
    free <- lv == 0L
    nfree <- tabulate(cid[free], nbins = ncl)
    open <- !sat
    if (any(open & nfree == 0L)) return(NULL)
    unit <- open & nfree == 1L
    if (!any(unit)) return(val)
    u <- unique(lits[free & unit[cid]])
    if (any(-u %in% u)) return(NULL)
    val[abs(u)] <- sign(u)
  }
}

.sat_flatten <- function(clauses) {
  lens <- lengths(clauses)
  list(lits = as.integer(unlist(clauses)), cid = rep.int(seq_along(clauses), lens), ncl = length(clauses))
}

.sat_search <- function(f, val) {
  val <- .sat_propagate(f$lits, f$cid, f$ncl, val)
  if (is.null(val)) return(NULL)
  lv <- sign(f$lits) * val[abs(f$lits)]
  sat <- tabulate(f$cid[lv == 1L], nbins = f$ncl) > 0L
  if (all(sat)) return(val)
  # branch on the free literal that occurs most often in open clauses
  cand <- f$lits[lv == 0L & !sat[f$cid]]
  tab <- table(cand)
  l <- as.integer(names(tab)[which.max(tab)])
  v1 <- val
  v1[abs(l)] <- sign(l)
  r <- .sat_search(f, v1)
  if (!is.null(r)) return(r)
  v1[abs(l)] <- -sign(l)
  .sat_search(f, v1)
}

# Satisfiability of integer clauses over nvars variables; a literal 0 is
# False. Returns a model (the assignment vector) or NULL.
.sat_solve <- function(clauses, nvars, assume = integer()) {
  clauses <- lapply(clauses, function(cl) cl[cl != 0L])
  if (any(lengths(clauses) == 0L)) return(NULL)
  val <- integer(nvars)
  if (length(assume)) {
    if (any(-assume %in% assume)) return(NULL)
    val[abs(assume)] <- sign(assume)
  }
  if (!length(clauses)) return(val)
  .sat_search(.sat_flatten(clauses), val)
}

# The assignment that unit propagation alone fixes, or NULL on a conflict.
.sat_fixed <- function(clauses, nvars) {
  clauses <- lapply(clauses, function(cl) cl[cl != 0L])
  if (any(lengths(clauses) == 0L)) return(NULL)
  val <- integer(nvars)
  if (!length(clauses)) return(val)
  f <- .sat_flatten(clauses)
  .sat_propagate(f$lits, f$cid, f$ncl, val)
}

.enc_satisfiable <- function(enc) !is.null(.sat_solve(enc$data, length(enc$symbols)))
