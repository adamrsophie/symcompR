# Combinatorics, after sympy/combinatorics (partitions.py, subsets.py,
# graycode.py, prufer.py) and sympy/utilities/iterables.py. Permutations and
# permutation groups are in combinatorics-permgroups.R.

# ---------------------------------------------------------------------------
# Iterables

#' Combinatorial generators
#'
#' `permutations(x, k)`, `combinations(x, k)`,
#' `combinations_with_replacement(x, k)`, `subsets(x)`, `set_partitions(n)`,
#' `derangements(x)`, `cartes(...)` (Cartesian product) and
#' `gray_code(n)` return lists. `nT(n, k)` counts partitions of `n` into `k`
#' parts and `stirling(n, k, kind)` gives Stirling numbers.
#'
#' @param x A vector or list.
#' @param k Size.
#' @param n An integer.
#' @param kind 1 or 2 for Stirling numbers of the first or second kind.
#' @param signed For first kind, return signed numbers.
#' @param ... Vectors for `cartes()`.
#' @return A list, or a `sym` integer.
#' @examples
#' length(permutations(1:4))
#' stirling(5, 2)
#' @export
permutations <- function(x, k = length(x)) {
  x <- as.list(x)
  rec <- function(avail, k) {
    if (k == 0L) return(list(list()))
    out <- list()
    for (i in seq_along(avail)) for (rest in rec(avail[-i], k - 1L)) out[[length(out) + 1L]] <- c(avail[i], rest)
    out
  }
  lapply(rec(x, k), unlist)
}

#' @rdname permutations
#' @export
combinations <- function(x, k) {
  if (k > length(x)) return(list())
  if (k == 0L) return(list(x[0]))
  cm <- utils::combn(length(x), k)
  lapply(seq_len(ncol(cm)), function(j) x[cm[, j]])
}

#' @rdname permutations
#' @export
combinations_with_replacement <- function(x, k) {
  n <- length(x)
  out <- list()
  rec <- function(start, cur) {
    if (length(cur) == k) {
      out[[length(out) + 1L]] <<- x[cur]
      return(invisible())
    }
    for (i in seq.int(start, length.out = max(0L, n - start + 1L))) rec(i, c(cur, i))
  }
  if (!n && k > 0) return(list())
  rec(1L, integer())
  out
}

#' @rdname permutations
#' @export
subsets <- function(x, k = NULL) {
  if (!is.null(k)) return(combinations(x, k))
  unlist(lapply(0:length(x), function(j) combinations(x, j)), recursive = FALSE)
}


#' @rdname permutations
#' @export
set_partitions <- function(n) {
  out <- list()
  rec <- function(i, blocks) {
    if (i > n) {
      out[[length(out) + 1L]] <<- blocks
      return(invisible())
    }
    for (b in seq_along(blocks)) {
      nb <- blocks
      nb[[b]] <- c(nb[[b]], i)
      rec(i + 1L, nb)
    }
    rec(i + 1L, c(blocks, list(i)))
  }
  rec(1L, list())
  out
}


#' @rdname permutations
#' @export
derangements <- function(x) {
  ps <- permutations(seq_along(x))
  keep <- Filter(function(p) all(p != seq_along(x)), ps)
  lapply(keep, function(p) x[p])
}


#' @rdname permutations
#' @export
cartes <- function(...) {
  g <- expand.grid(rev(list(...)), stringsAsFactors = FALSE)
  g <- g[, rev(seq_len(ncol(g))), drop = FALSE]
  lapply(seq_len(nrow(g)), function(i) unlist(g[i, ], use.names = FALSE))
}

#' @rdname permutations
#' @export
gray_code <- function(n) {
  lapply(0:(2^n - 1), function(i) {
    g <- bitwXor(i, bitwShiftR(i, 1L))
    as.integer(rev(as.integer(intToBits(g))[seq_len(n)]))
  })
}

#' @rdname permutations
#' @export
nT <- function(n, k) {
  n <- as.integer(n); k <- as.integer(k)
  memo <- new.env()
  f <- function(n, k) {
    if (n == 0L && k == 0L) return(gmp::as.bigz(1))
    if (n <= 0L || k <= 0L) return(gmp::as.bigz(0))
    key <- paste(n, k)
    if (!is.null(memo[[key]])) return(memo[[key]])
    r <- f(n - 1L, k - 1L) + f(n - k, k)
    assign(key, r, envir = memo)
    r
  }
  make_rational(gmp::as.bigq(f(n, k)))
}

#' @rdname permutations
#' @export
stirling <- function(n, k, kind = 2, signed = FALSE) {
  n <- as.integer(n); k <- as.integer(k)
  S <- vector("list", n + 1L)
  for (i in 0:n) S[[i + 1L]] <- gmp::as.bigz(rep(0, k + 1L))
  S[[1L]][1L] <- 1
  if (n >= 1L) for (i in 1:n) for (j in seq_len(min(i, k))) {
    if (kind == 2) S[[i + 1L]][j + 1L] <- j * S[[i]][j + 1L] + S[[i]][j]
    else S[[i + 1L]][j + 1L] <- (i - 1L) * S[[i]][j + 1L] + S[[i]][j]
  }
  v <- S[[n + 1L]][k + 1L]
  if (kind == 1 && signed && (n - k) %% 2L == 1L) v <- -v
  make_rational(gmp::as.bigq(v))
}

# ---------------------------------------------------------------------------
# Prufer sequences

#' Prufer sequences
#'
#' `prufer_from_tree(edges, n)` encodes a labeled tree (0 based edges) and
#' `tree_from_prufer(seq)` decodes it.
#'
#' @param edges A list of length 2 integer vectors.
#' @param n Number of nodes.
#' @param seq A Prufer sequence.
#' @return An integer vector or a list of edges.
#' @export
prufer_from_tree <- function(edges, n) {
  adj <- lapply(seq_len(n), function(i) integer())
  for (e in edges) {
    e <- as.integer(e)
    adj[[e[1] + 1L]] <- c(adj[[e[1] + 1L]], e[2])
    adj[[e[2] + 1L]] <- c(adj[[e[2] + 1L]], e[1])
  }
  deg <- vapply(adj, length, 0L)
  out <- integer()
  removed <- logical(n)
  for (step in seq_len(n - 2L)) {
    leaf <- which(deg == 1L & !removed)[1L]
    nb <- Filter(function(v) !removed[v + 1L], adj[[leaf]])[1L]
    out <- c(out, nb)
    removed[leaf] <- TRUE
    deg[leaf] <- 0L
    deg[nb + 1L] <- deg[nb + 1L] - 1L
  }
  out
}

#' @rdname prufer_from_tree
#' @export
tree_from_prufer <- function(seq) {
  seq <- as.integer(seq)
  n <- length(seq) + 2L
  deg <- rep(1L, n)
  for (v in seq) deg[v + 1L] <- deg[v + 1L] + 1L
  edges <- list()
  for (v in seq) {
    leaf <- which(deg == 1L)[1L] - 1L
    edges[[length(edges) + 1L]] <- sort(c(leaf, v))
    deg[leaf + 1L] <- deg[leaf + 1L] - 1L
    deg[v + 1L] <- deg[v + 1L] - 1L
  }
  last <- which(deg == 1L) - 1L
  edges[[length(edges) + 1L]] <- sort(last)
  edges
}
