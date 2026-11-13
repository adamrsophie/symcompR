# Prufer sequences of labeled trees, after sympy/combinatorics/prufer.py.
# prufer_from_tree() and tree_from_prufer() in combinatorics.R are
# Prufer.to_prufer() and Prufer.to_tree().

.pr_missing_msg <- function(missing) {
  if (length(missing) == 1L) sprintf("Node %s is missing.", missing)
  else sprintf("Nodes [%s] are missing.", paste(sort(missing), collapse = ", "))
}

#' Prufer sequences
#'
#' `Prufer(x, n)` is SymPy's `Prufer`: a labeled tree on the nodes
#' `0, ..., n - 1`, given either as a list of edges (pairs) or as its Prufer
#' sequence. The fields `$prufer_repr`, `$tree_repr`, `$nodes`, `$rank` and
#' `$size` mirror SymPy's. `prufer_next()` and `prufer_prev()` move through the
#' trees in rank order, `prufer_unrank(rank, n)` builds the tree of a rank,
#' `prufer_edges(...)` turns runs of connected nodes into a list of edges and
#' a node count (relabeled from 0), and `prufer_to_prufer()` and
#' `prufer_to_tree()` convert between the two representations.
#'
#' @param x A list of edges or a Prufer sequence.
#' @param n Number of nodes.
#' @param p A `Prufer`.
#' @param delta Steps to move.
#' @param rank A rank.
#' @param ... Runs of nodes (vectors).
#' @param tree A list of edges.
#' @param prufer A Prufer sequence.
#' @return A `Prufer`, a sequence, or a list of edges.
#' @examples
#' a <- Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)))
#' a$prufer_repr
#' a$rank
#' prufer_next(a)$tree_repr
#' Prufer(c(1, 1))$tree_repr
#' prufer_edges(c(0, 1, 2, 3), c(1, 4, 5), c(1, 4, 6))
#' @export
Prufer <- function(x, n = NULL) {
  tree <- is.list(x) && length(x) > 0L
  if (tree) {
    edges <- lapply(x, as.integer)
    if (!length(edges[[1L]]))
      stop(.sympy_condition("ValueError", "Prufer expects at least one edge in the tree."))
    if (is.null(n)) {
      nodes <- unique(unlist(edges))
      nnodes <- max(nodes) + 1L
      if (nnodes != length(nodes))
        stop(.sympy_condition("ValueError", .pr_missing_msg(setdiff(0:(nnodes - 1L), nodes))))
    } else nnodes <- as.integer(n)
    return(structure(list(tree = edges, prufer = NULL, nodes = nnodes, n_arg = n), class = "Prufer"))
  }
  seq <- as.integer(unlist(x))
  structure(list(tree = NULL, prufer = seq, nodes = length(seq) + 2L, n_arg = n), class = "Prufer")
}

#' @export
`$.Prufer` <- function(x, name) {
  switch(name,
         prufer_repr = {
           pr <- .subset2(x, "prufer")
           if (is.null(pr)) prufer_to_prufer(.subset2(x, "tree"), .subset2(x, "nodes")) else pr
         },
         tree_repr = {
           tr <- .subset2(x, "tree")
           if (is.null(tr)) prufer_to_tree(.subset2(x, "prufer")) else tr
         },
         rank = .pr_rank(x),
         size = .gc_small(gmp::as.bigz(prufer_prev(prufer_prev(x, x$rank))$rank) + 1),
         .subset2(x, name))
}

.pr_rank <- function(p) {
  r <- gmp::as.bigz(0)
  m <- gmp::as.bigz(1)
  pr <- p$prufer_repr
  n <- .subset2(p, "nodes")
  if (n - 3L >= 0L) for (i in (n - 3L):0L) {
    r <- r + m * pr[i + 1L]
    m <- m * n
  }
  .gc_small(r)
}

#' @export
format.Prufer <- function(x, ...) {
  tr <- .subset2(x, "tree")
  s <- if (is.null(tr)) paste0("[", paste(.subset2(x, "prufer"), collapse = ", "), "]")
  else paste0("[", paste(vapply(tr, function(e) paste0("[", paste(e, collapse = ", "), "]"), ""), collapse = ", "), "]")
  n <- .subset2(x, "n_arg")
  paste0("Prufer(", s, if (!is.null(n)) paste0(", ", n), ")")
}

#' @export
print.Prufer <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @rdname Prufer
#' @export
prufer_next <- function(p, delta = 1) prufer_unrank(gmp::as.bigz(p$rank) + delta, .subset2(p, "nodes"))

#' @rdname Prufer
#' @export
prufer_prev <- function(p, delta = 1) prufer_unrank(gmp::as.bigz(p$rank) - delta, .subset2(p, "nodes"))

#' @rdname Prufer
#' @export
prufer_unrank <- function(rank, n) {
  n <- as.integer(.as_int(n))
  rank <- .as_int(rank)
  L <- integer(max(n - 2L, 0L))
  if (n - 3L >= 0L) for (i in (n - 3L):0L) {
    L[i + 1L] <- as.integer(rank %% n)
    rank <- (rank - L[i + 1L]) %/% n
  }
  Prufer(L)
}

#' @rdname Prufer
#' @export
prufer_to_prufer <- function(tree, n) prufer_from_tree(tree, n)

#' @rdname Prufer
#' @export
prufer_to_tree <- function(prufer) tree_from_prufer(as.integer(prufer))

#' @rdname Prufer
#' @export
prufer_edges <- function(...) {
  runs <- lapply(list(...), as.integer)
  e <- list()
  for (r in runs) for (i in seq_len(length(r) - 1L)) e[[length(e) + 1L]] <- sort(r[i:(i + 1L)])
  e <- unique(e)
  got <- unique(unlist(e))
  nmin <- min(vapply(e, `[`, 0L, 1L))
  nmax <- max(vapply(e, `[`, 0L, 2L))
  missing <- setdiff(nmin:nmax, got)
  if (length(missing)) stop(.sympy_condition("ValueError", .pr_missing_msg(missing + nmin)))
  if (nmin != 0L) {
    e <- lapply(e, function(ei) ei - nmin)
    nmax <- nmax - nmin
  }
  e <- e[order(vapply(e, `[`, 0L, 1L), vapply(e, `[`, 0L, 2L))]
  list(e, nmax + 1L)
}
