# Circuit utilities (sympy/physics/quantum/circuitutils.py) and the
# search for gate identities (sympy/physics/quantum/identitysearch.py).
# Python tuples are R lists, Python sets are lists of distinct elements with
# a "qset" attribute and dictionaries are lists of key-value pairs with a
# "qdict" attribute.

.qcirc_seq <- function(x) {
  if (is.character(x) && length(x) == 1L) return(as.list(strsplit(x, "")[[1L]]))
  if (is_sym(x) && x$type == "Mul") return(.qargs(x))
  if (is_sym(x) && x$type == "Tuple") return(x$args)
  if (is_sym(x)) return(list(x))
  as.list(x)
}

.qtkey <- function(x) {
  if (is_sym(x)) return(.k(x))
  if (is.character(x) && length(x) == 1L) return(paste0("str:", x))
  if (is.numeric(x)) return(.k(as_sym(x)))
  paste0("(", paste(vapply(x, .qtkey, ""), collapse = ","), ")")
}
.qteq <- function(a, b) identical(.qtkey(a), .qtkey(b))

.qset_of <- function(xs) {
  keys <- vapply(xs, .qtkey, "")
  structure(xs[!duplicated(keys)], qset = TRUE)
}

#' Circuit utilities
#'
#' Tools of `sympy.physics.quantum.circuitutils` for circuits given as
#' products of gates or lists of gates: `kmp_table()` is the
#' Knuth-Morris-Pratt table of a list, `find_subcircuit()` finds the
#' 0-based position of a subcircuit (or -1), `replace_subcircuit()`
#' replaces it, `convert_to_symbolic_indices()` and
#' `convert_to_real_indices()` switch between numeric and symbolic qubit
#' labels, and `random_reduce()` and `random_insert()` remove or insert
#' gate identities. A `seed` given as a list of integers is used as the
#' sequence of random choices.
#'
#' @param word A list.
#' @param circuit,subcircuit,replace Circuits.
#' @param start,end,pos Positions (0-based).
#' @param seq A circuit.
#' @param gen A list of symbols to draw new indices from.
#' @param qubit_map A dictionary of symbols and qubit indices.
#' @param gate_ids A list of [GateIdentity()] objects.
#' @param choices A list of circuits.
#' @param seed `NULL`, an integer or a list of integers.
#' @return A list, a position or a circuit.
#' @examples
#' x <- XGate(0); y <- YGate(0)
#' find_subcircuit(x * y * x, y * x)
#' replace_subcircuit(x * y * x, y)
#' @name circuitutils
NULL

#' @rdname circuitutils
#' @export
kmp_table <- function(word) {
  word <- .qcirc_seq(word)
  pos <- 2L; cnd <- 0L
  table <- c(-1L, 0L)
  while (pos < length(word)) {
    if (.qteq(word[[pos]], word[[cnd + 1L]])) {
      cnd <- cnd + 1L
      table <- c(table, cnd)
      pos <- pos + 1L
    } else if (cnd > 0L) {
      cnd <- table[cnd + 1L]
    } else {
      table <- c(table, 0L)
      pos <- pos + 1L
    }
  }
  as.list(table)
}

#' @rdname circuitutils
#' @export
find_subcircuit <- function(circuit, subcircuit, start = 0, end = 0) {
  circuit <- .qcirc_seq(circuit); subcircuit <- .qcirc_seq(subcircuit)
  start <- .qint(start); end <- .qint(end)
  if (!length(subcircuit) || length(subcircuit) > length(circuit)) return(-1L)
  if (end < 1L) end <- length(circuit)
  pos <- start; index <- 0L
  table <- unlist(kmp_table(subcircuit))
  while (pos + index < end) {
    if (.qteq(subcircuit[[index + 1L]], circuit[[pos + index + 1L]])) {
      index <- index + 1L
    } else {
      pos <- pos + index - table[index + 1L]
      index <- if (table[index + 1L] > -1L) table[index + 1L] else 0L
    }
    if (index == length(subcircuit)) return(pos)
  }
  -1L
}

#' @rdname circuitutils
#' @export
replace_subcircuit <- function(circuit, subcircuit, replace = NULL, pos = 0) {
  pos <- max(.qint(pos), 0L)
  circuit <- .qcirc_seq(circuit); subcircuit <- .qcirc_seq(subcircuit)
  replace <- if (is.null(replace)) list() else .qcirc_seq(replace)
  loc <- find_subcircuit(circuit, subcircuit, start = pos)
  if (loc > -1L) {
    left <- circuit[seq_len(loc)]
    right <- circuit[seq.int(loc + length(subcircuit) + 1L, length.out = length(circuit) - loc - length(subcircuit))]
    circuit <- c(left, replace, right)
  }
  circuit
}

# dictionaries
.qd_new <- function() structure(list(), qdict = TRUE)
.qd_find <- function(d, k) {
  kk <- .qtkey(k)
  for (i in seq_along(d)) if (identical(.qtkey(d[[i]][[1L]]), kk)) return(i)
  0L
}
.qd_set <- function(d, k, v) {
  i <- .qd_find(d, k)
  a <- attributes(d)
  if (i) d[[i]] <- list(k, v) else d[[length(d) + 1L]] <- list(k, v)
  attributes(d) <- a
  d
}
.qd_get <- function(d, k) {
  i <- .qd_find(d, k)
  if (!i) .qstop("KeyError", .qtkey(k))
  d[[i]][[2L]]
}
.qis_dict <- function(x) isTRUE(attr(x, "qdict"))

# a generator of numbered symbols: prefix0, prefix1, ...
.qsymbol_gen <- function(prefix = "x", start = 0L) {
  g <- new.env()
  g$i <- start
  g$next_value <- function() { s <- Symbol(paste0(prefix, g$i)); g$i <- g$i + 1L; s }
  class(g) <- "qsymbol_gen"
  g
}
.qgen_wrap <- function(x) {
  if (inherits(x, "qsymbol_gen")) return(x)
  g <- new.env(); g$items <- x; g$pos <- 0L
  g$next_value <- function() {
    g$pos <- g$pos + 1L
    if (g$pos > length(g$items)) .qstop("StopIteration", "")
    g$items[[g$pos]]
  }
  class(g) <- "qsymbol_gen"
  g
}

#' @rdname circuitutils
#' @export
convert_to_symbolic_indices <- function(seq, start = NULL, gen = NULL, qubit_map = NULL) {
  seq <- .qcirc_seq(seq)
  index_gen <- .qsymbol_gen("i", -1L)
  cur_ndx <- index_gen$next_value()
  ndx_map <- .qd_new()
  inverse <- function(m) {
    out <- .qd_new()
    for (p in m) out <- .qd_set(out, p[[2L]], p[[1L]])
    out
  }
  if (!is.null(start)) {
    if (!(is_sym(start) && start$type %in% c("Symbol", "Dummy")))
      .qstop("TypeError", "Expected Symbol for starting index, got ", .qrepr(start), ".")
    cur_ndx <- start
  }
  if (!is.null(gen)) {
    if (!(inherits(gen, "qsymbol_gen") || (is.list(gen) && !is_sym(gen))))
      .qstop("TypeError", "Expected a generator, got ", .qrepr(gen), ".")
    gen <- .qgen_wrap(gen)
    index_gen <- gen
  }
  if (!is.null(qubit_map)) {
    if (!.qis_dict(qubit_map)) .qstop("TypeError", "Expected dict for existing map, got ", .qrepr(qubit_map), ".")
    ndx_map <- qubit_map
  }
  ndx_map <- structure(lapply(ndx_map, function(p) list(p[[1L]], as_sym(p[[2L]]))), qdict = TRUE)
  inv_map <- inverse(ndx_map)
  sym_seq <- list()
  for (item in seq) {
    if (is_sym(item) && .qisa(item, "Gate")) {
      r <- convert_to_symbolic_indices(item$args, qubit_map = ndx_map, start = cur_ndx, gen = index_gen)
    } else if ((is_sym(item) && item$type == "Tuple") || (is.list(item) && !is_sym(item))) {
      r <- convert_to_symbolic_indices(item, qubit_map = ndx_map, start = cur_ndx, gen = index_gen)
    } else {
      r <- NULL
    }
    if (!is.null(r)) {
      sym_item <- r[[1L]]
      for (p in r[[2L]]) ndx_map <- .qd_set(ndx_map, p[[1L]], p[[2L]])
      cur_ndx <- r[[3L]]; index_gen <- r[[4L]]
      inv_map <- inverse(ndx_map)
    } else if (.qd_find(inv_map, as_sym(item))) {
      sym_item <- .qd_get(inv_map, as_sym(item))
    } else {
      if (is.null(gen)) .qstop("TypeError", "'NoneType' object is not an iterator")
      cur_ndx <- gen$next_value()
      ndx_map <- .qd_set(ndx_map, cur_ndx, as_sym(item))
      inv_map <- .qd_set(inv_map, as_sym(item), cur_ndx)
      sym_item <- cur_ndx
    }
    if (is_sym(item) && .qisa(item, "Gate")) sym_item <- .qfunc(item$type, sym_item)
    sym_seq[[length(sym_seq) + 1L]] <- sym_item
  }
  list(sym_seq, ndx_map, cur_ndx, index_gen)
}

#' @rdname circuitutils
#' @export
convert_to_real_indices <- function(seq, qubit_map) {
  seq <- .qcirc_seq(seq)
  if (!.qis_dict(qubit_map)) .qstop("TypeError", "Expected dict for qubit_map, got ", .qrepr(qubit_map), ".")
  qubit_map <- structure(lapply(qubit_map, function(p) list(p[[1L]], as_sym(p[[2L]]))), qdict = TRUE)
  real_seq <- list()
  for (item in seq) {
    if (is_sym(item) && .qisa(item, "Gate")) {
      real_item <- .qfunc(item$type, convert_to_real_indices(item$args, qubit_map))
    } else if ((is_sym(item) && item$type == "Tuple") || (is.list(item) && !is_sym(item))) {
      real_item <- convert_to_real_indices(item, qubit_map)
    } else {
      real_item <- .qd_get(qubit_map, item)
    }
    real_seq[[length(real_seq) + 1L]] <- real_item
  }
  real_seq
}

# sympy.core.random._randrange
.qrandrange <- function(seed = NULL) {
  if (is.null(seed)) return(function(a, b = NULL) { if (is.null(b)) { b <- a; a <- 0L }; a + sample.int(b - a, 1L) - 1L })
  # an integer reseeds CPython's generator (R/utilities-pyrandom.R)
  if (is.numeric(seed) && length(seed) == 1L || is_sym(seed)) return(.sympy_randrange(.qint(seed)))
  st <- new.env()
  st$seq <- vapply(seed, .qint, 0)
  give <- function(a, b = NULL) {
    if (is.null(b)) { b <- a; a <- 0L }
    if (b - a < 1) .qstop("ValueError", "_randrange got empty range")
    repeat {
      if (!length(st$seq)) .qstop("ValueError", "_randrange sequence was too short")
      x <- st$seq[1L]; st$seq <- st$seq[-1L]
      if (a <= x && x < b) return(x)
    }
  }
  give
}

.qtuple_sort_key <- function(t) list(list(10, 0, "tuple"), list(length(t), lapply(t, .sort_key)), .ONE_SORT_KEY, 1)

flatten_ids <- function(ids) {
  acc <- list()
  for (an_id in ids) {
    eq <- .qa(an_id, "equivalent_ids")
    acc <- c(acc, .key_sort(lapply(eq, .qtuple_sort_key), unclass(eq)))
  }
  .key_sort(lapply(acc, .qtuple_sort_key), acc)
}

#' @rdname circuitutils
#' @export
random_reduce <- function(circuit, gate_ids, seed = NULL) {
  if (!length(gate_ids)) return(circuit)
  circuit <- .qcirc_seq(circuit)
  ids <- flatten_ids(gate_ids)
  randrange <- .qrandrange(seed)
  found <- FALSE
  while (length(ids)) {
    i <- randrange(length(ids))
    id <- ids[[i + 1L]]; ids <- ids[-(i + 1L)]
    if (find_subcircuit(circuit, id) != -1L) { found <- TRUE; break }
  }
  if (!found) return(circuit)
  replace_subcircuit(circuit, id)
}

#' @rdname circuitutils
#' @export
random_insert <- function(circuit, choices, seed = NULL) {
  if (!length(choices)) return(circuit)
  circuit <- .qcirc_seq(circuit)
  randrange <- .qrandrange(seed)
  loc <- randrange(length(circuit) + 1L)
  choice <- .qcirc_seq(choices[[randrange(length(choices)) + 1L]])
  c(circuit[seq_len(loc)], choice, circuit[seq.int(loc + 1L, length.out = length(circuit) - loc)])
}

# ---------------------------------------------------------------------------
# identitysearch

#' Gate identities
#'
#' Tools of `sympy.physics.quantum.identitysearch`.
#' `is_scalar_nonsparse_matrix()` (also `is_scalar_matrix()`) tests whether
#' a circuit is a scalar multiple of the identity, `generate_gate_rules()`
#' and `generate_equivalent_ids()` give the equivalent forms of a gate
#' identity, `GateIdentity(...)` holds an identity and
#' `bfs_identity_search()` and `random_identity_search()` look for them.
#' The scipy sparse version `is_scalar_sparse_matrix()` uses the same
#' exact computation here.
#'
#' @param circuit,gate_seq A circuit (a product or list of gates).
#' @param nqubits The number of qubits.
#' @param identity_only Accept only the identity, not its multiples.
#' @param eps Unused tolerance of the sparse version.
#' @param return_as_muls Return products instead of lists.
#' @param ... Gates.
#' @param identity_set A list of `GateIdentity` objects.
#' @param gate_identity A list of gates.
#' @param begin,end The range of gates (0-based, end excluded).
#' @param gate_list A list of gates.
#' @param max_depth The maximum circuit length.
#' @param numgates The number of gates.
#' @return A logical, a list (a set) or a `GateIdentity`.
#' @examples
#' x <- XGate(0)
#' is_scalar_nonsparse_matrix(list(x, x), 1, TRUE)
#' generate_equivalent_ids(list(x, YGate(0), ZGate(0)))
#' @name identitysearch
NULL

#' @rdname identitysearch
#' @export
is_scalar_nonsparse_matrix <- function(circuit, nqubits, identity_only, eps = NULL) {
  matrix <- represent(.mul(.qcirc_seq(circuit)), nqubits = nqubits)
  if (is_sym(matrix) && is_number_type(matrix))
    return(if (isTRUE(identity_only)) .qeq(matrix, S_One) else TRUE)
  tr <- sym_trace(matrix)
  first <- matrix$data[[1L]]
  adjusted <- if (!isTRUE(identity_only)) .mul(list(tr, .pow(first, S_NegativeOne))) else tr
  is_identity <- if (isTRUE(identity_only)) .qeq(first, S_One) else TRUE
  has_correct_trace <- .qeq(adjusted, .pow(S_Two, as_sym(nqubits)))
  diag <- TRUE
  for (i in seq_len(matrix$nrow)) for (j in seq_len(matrix$ncol))
    if (i != j && .qtruthy(.mget(matrix, i, j))) diag <- FALSE
  diag && has_correct_trace && is_identity
}

#' @rdname identitysearch
#' @export
is_scalar_sparse_matrix <- function(circuit, nqubits, identity_only, eps = 1e-11)
  is_scalar_nonsparse_matrix(circuit, nqubits, identity_only, eps)

#' @rdname identitysearch
#' @export
is_scalar_matrix <- is_scalar_nonsparse_matrix

.qget_min_qubits <- function(g) if (g$type == "Pow") .qa(g$args[[1L]], "min_qubits") else .qa(g, "min_qubits")

.qunitary_end <- function(g) is_scalar_matrix(list(Dagger(g), g), .qget_min_qubits(g), TRUE)

ll_op <- function(left, right) {
  if (length(left) && .qunitary_end(left[[1L]])) return(list(left[-1L], c(list(Dagger(left[[1L]])), right)))
  NULL
}
lr_op <- function(left, right) {
  n <- length(left)
  if (n && .qunitary_end(left[[n]])) return(list(left[-n], c(right, list(Dagger(left[[n]])))))
  NULL
}
rl_op <- function(left, right) {
  if (length(right) && .qunitary_end(right[[1L]])) return(list(c(list(Dagger(right[[1L]])), left), right[-1L]))
  NULL
}
rr_op <- function(left, right) {
  n <- length(right)
  if (n && .qunitary_end(right[[n]])) return(list(c(left, list(Dagger(right[[n]]))), right[-n]))
  NULL
}

#' @rdname identitysearch
#' @export
generate_gate_rules <- function(gate_seq, return_as_muls = FALSE) {
  if (is_sym(gate_seq) && is_number_type(gate_seq))
    return(if (return_as_muls) .qset_of(list(list(S_One, S_One))) else .qset_of(list(list(list(), list()))))
  gate_seq <- .qcirc_seq(gate_seq)
  queue <- list(list(gate_seq, list(), 0L))
  rules <- list(list(gate_seq, list()))
  rule_keys <- .qtkey(rules[[1L]])
  max_ops <- length(gate_seq)
  process <- function(new_rule, ops) {
    if (is.null(new_rule)) return()
    k <- .qtkey(new_rule); kr <- .qtkey(list(new_rule[[2L]], new_rule[[1L]]))
    if (!(k %in% rule_keys) && !(kr %in% rule_keys)) {
      rules[[length(rules) + 1L]] <<- new_rule
      rule_keys <<- c(rule_keys, k)
    }
    if (ops + 1L < max_ops) queue[[length(queue) + 1L]] <<- c(new_rule, list(ops + 1L))
  }
  while (length(queue)) {
    q <- queue[[1L]]; queue <- queue[-1L]
    left <- q[[1L]]; right <- q[[2L]]; ops <- q[[3L]]
    process(ll_op(left, right), ops)
    process(lr_op(left, right), ops)
    process(rl_op(left, right), ops)
    process(rr_op(left, right), ops)
  }
  if (return_as_muls) rules <- lapply(rules, function(r) list(.mul(r[[1L]]), .mul(r[[2L]])))
  .qset_of(rules)
}

#' @rdname identitysearch
#' @export
generate_equivalent_ids <- function(gate_seq, return_as_muls = FALSE) {
  if (is_sym(gate_seq) && is_number_type(gate_seq)) return(.qset_of(list(S_One)))
  gate_seq <- .qcirc_seq(gate_seq)
  eq_ids <- list()
  for (rule in generate_gate_rules(gate_seq)) {
    if (!length(rule[[1L]])) eq_ids[[length(eq_ids) + 1L]] <- rule[[2L]]
    else if (!length(rule[[2L]])) eq_ids[[length(eq_ids) + 1L]] <- rule[[1L]]
  }
  if (return_as_muls) eq_ids <- lapply(eq_ids, .mul)
  .qset_of(eq_ids)
}

.qclass("GateIdentity", "Basic",
  `__new__` = .cm(function(cls, args, ...) {
    args <- lapply(args, as_sym)
    new_sym("GateIdentity", args, extra = list(circuit = .mul(args), rules = generate_gate_rules(args),
                                                 eq_ids = generate_equivalent_ids(args)))
  }),
  circuit = .prop(function(self) self$extra$circuit),
  gate_rules = .prop(function(self) self$extra$rules),
  equivalent_ids = .prop(function(self) self$extra$eq_ids),
  sequence = .prop(function(self) self$args),
  `_sympystr` = function(self, p) .str(self$extra$circuit),
  `_sympyrepr` = function(self, p) paste0("GateIdentity(", paste(vapply(self$args, srepr, ""), collapse = ", "), ")"))

#' @rdname identitysearch
#' @export
GateIdentity <- function(...) .qfunc("GateIdentity", list(...))

#' @rdname identitysearch
#' @export
is_degenerate <- function(identity_set, gate_identity) {
  k <- .qtkey(.qcirc_seq(gate_identity))
  for (an_id in identity_set)
    if (k %in% vapply(.qa(an_id, "equivalent_ids"), .qtkey, "")) return(TRUE)
  FALSE
}

#' @rdname identitysearch
#' @export
is_reducible <- function(circuit, nqubits, begin, end) {
  circuit <- .qcirc_seq(circuit)
  current <- list()
  for (ndx in rev(seq.int(.qint(begin), length.out = .qint(end) - .qint(begin)))) {
    current <- c(list(circuit[[ndx + 1L]]), current)
    if (is_scalar_matrix(current, nqubits, FALSE)) return(TRUE)
  }
  FALSE
}

#' @rdname identitysearch
#' @export
bfs_identity_search <- function(gate_list, nqubits, max_depth = NULL, identity_only = FALSE) {
  gate_list <- .qcirc_seq(gate_list)
  if (is.null(max_depth) || max_depth <= 0) max_depth <- length(gate_list)
  queue <- list(list())
  ids <- list()
  while (length(queue)) {
    current <- queue[[1L]]; queue <- queue[-1L]
    for (next_gate in gate_list) {
      new_circuit <- c(current, list(next_gate))
      reducible <- is_reducible(new_circuit, nqubits, 1L, length(new_circuit))
      if (is_scalar_matrix(new_circuit, nqubits, identity_only) && !is_degenerate(ids, new_circuit) && !reducible) {
        gi <- do.call(GateIdentity, new_circuit)
        if (!any(vapply(ids, function(o) identical(.k(o), .k(gi)), NA))) ids[[length(ids) + 1L]] <- gi
      } else if (length(new_circuit) < max_depth && !reducible) {
        queue[[length(queue) + 1L]] <- new_circuit
      }
    }
  }
  structure(ids, qset = TRUE)
}

#' @rdname identitysearch
#' @export
random_identity_search <- function(gate_list, numgates, nqubits) {
  gate_list <- .qcirc_seq(gate_list)
  circuit <- lapply(seq_len(numgates), function(i) gate_list[[sample.int(length(gate_list), 1L)]])
  if (is_scalar_matrix(circuit, nqubits, FALSE)) circuit else NULL
}
