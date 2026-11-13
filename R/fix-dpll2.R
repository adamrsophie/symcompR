# The DPLL search of sympy/logic/algorithms/dpll2.py (SATSolver with the
# VSIDS heuristic and no clause learning), used by satisfiable(). SymPy
# numbers the variables in set iteration order, which depends on the
# Python hash seed; here they follow the sorted atom order.

# Iteration order of a small Python set of int literals: hash(-1) is -2 and
# every other int hashes to itself; the slot is the hash modulo the table
# size. ponytail: collisions keep insertion order instead of CPython probing.
.py_int_set <- function(lits) {
  lits <- unique(lits)
  n <- length(lits)
  size <- if (n <= 4L) 8L else if (n <= 18L) 32L else 128L
  h <- ifelse(lits == -1L, -2L, lits)
  lits[order(h %% size, seq_len(n))]
}

# All models (a list of integer vectors of set literals), in SymPy order.
.dpll2_models <- function(clauses, nvars, all_models = TRUE) {
  st <- new.env(parent = emptyenv())
  st$clauses <- lapply(clauses, .py_int_set)
  st$sentinels <- list()
  st$occ <- list()
  key <- function(l) as.character(l)
  get0l <- function(lst, l) {
    v <- lst[[key(l)]]
    if (is.null(v)) integer() else v
  }
  st$queue <- integer()
  for (i in seq_along(st$clauses)) {
    cl <- st$clauses[[i]]
    if (length(cl) == 1L) {
      st$queue <- c(st$queue, cl)
      next
    }
    for (s in unique(c(cl[1L], cl[length(cl)]))) {
      st$sentinels[[key(s)]] <- union(get0l(st$sentinels, s), i)
    }
    for (l in cl) st$occ[[key(l)]] <- length(get0l(st$occ, l)) + 1L
  }
  occ <- function(l) {
    v <- st$occ[[key(l)]]
    if (is.null(v)) 0 else v
  }
  st$scores <- list()
  st$heap <- list(score = numeric(), lit = integer())
  push <- function(sc, l) {
    st$heap$score <- c(st$heap$score, sc)
    st$heap$lit <- c(st$heap$lit, l)
  }
  for (v in seq_len(nvars)) {
    st$scores[[key(v)]] <- -occ(v)
    st$scores[[key(-v)]] <- -occ(-v)
    push(-occ(v), v)
    push(-occ(-v), -v)
  }
  st$varset <- rep(FALSE, nvars)
  st$settings <- integer()
  st$levels <- list(list(decision = 0L, settings = integer(), flipped = FALSE))
  st$unsat <- FALSE
  cur <- function() st$levels[[length(st$levels)]]
  set_cur <- function(lv) st$levels[[length(st$levels)]] <- lv
  heap_pop <- function() {
    o <- order(st$heap$score, st$heap$lit)[1L]
    l <- st$heap$lit[o]
    st$heap$score <- st$heap$score[-o]
    st$heap$lit <- st$heap$lit[-o]
    l
  }
  calculate <- function() {
    repeat {
      if (!length(st$heap$lit)) return(0L)
      o <- order(st$heap$score, st$heap$lit)[1L]
      if (!st$varset[abs(st$heap$lit[o])]) return(heap_pop())
      heap_pop()
    }
  }
  clause_sat <- function(i) any(st$clauses[[i]] %in% st$settings)
  assign_lit <- function(l) {
    st$settings <- c(st$settings, l)
    lv <- cur()
    lv$settings <- c(lv$settings, l)
    set_cur(lv)
    st$varset[abs(l)] <- TRUE
    for (cls in get0l(st$sentinels, -l)) {
      if (clause_sat(cls)) next
      other <- NULL
      for (nl in st$clauses[[cls]]) {
        if (nl == -l) next
        if (cls %in% get0l(st$sentinels, nl)) {
          other <- nl
        } else if (!st$varset[abs(nl)]) {
          st$sentinels[[key(-l)]] <- setdiff(get0l(st$sentinels, -l), cls)
          st$sentinels[[key(nl)]] <- union(get0l(st$sentinels, nl), cls)
          other <- NULL
          break
        }
      }
      if (!is.null(other)) st$queue <- c(st$queue, other)
    }
  }
  undo <- function() {
    for (l in cur()$settings) {
      st$settings <- st$settings[st$settings != l]
      v <- abs(l)
      push(st$scores[[key(v)]], v)
      push(st$scores[[key(-v)]], -v)
      st$varset[v] <- FALSE
    }
    st$levels[[length(st$levels)]] <- NULL
  }
  create_level <- function(l, flipped = FALSE) {
    st$levels[[length(st$levels) + 1L]] <- list(decision = l, settings = integer(), flipped = flipped)
  }
  simplify <- function() {
    while (length(st$queue)) {
      nl <- st$queue[length(st$queue)]
      st$queue <- st$queue[-length(st$queue)]
      if ((-nl) %in% st$settings) {
        st$unsat <- TRUE
      } else {
        assign_lit(nl)
      }
      if (st$unsat) {
        st$queue <- integer()
        return(invisible())
      }
    }
  }
  models <- list()
  simplify()
  if (st$unsat) return(models)
  flip <- FALSE
  repeat {
    if (flip) {
      flip <- FALSE
      l <- cur()$decision
    } else {
      l <- calculate()
      if (l == 0L) {
        models[[length(models) + 1L]] <- st$settings
        if (!all_models) return(models)
        while (cur()$flipped) undo()
        if (length(st$levels) == 1L) return(models)
        fl <- -cur()$decision
        undo()
        create_level(fl, TRUE)
        flip <- TRUE
        next
      }
      create_level(l)
    }
    assign_lit(l)
    simplify()
    if (st$unsat) {
      st$unsat <- FALSE
      while (cur()$flipped) {
        undo()
        if (length(st$levels) == 1L) return(models)
      }
      if (length(st$levels) == 1L) return(models)
      fl <- -cur()$decision
      undo()
      create_level(fl, TRUE)
      flip <- TRUE
    }
  }
}
