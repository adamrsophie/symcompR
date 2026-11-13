# Subsets of a superset, after sympy/combinatorics/subsets.py.

# a comparable key for an element (Python ==)
.subs_keys <- function(x) {
  vapply(as.list(x), function(e) {
    if (is_sym(e)) return(e$key)
    if (is.numeric(e)) return(paste0("n:", format(as.numeric(e), digits = 17)))
    paste0(class(e)[1L], ":", format(e))
  }, "")
}

.subs_elem_str <- function(e) if (is_sym(e)) sstr(e) else format(e)

#' Subsets
#'
#' `Subset(subset, superset)` is SymPy's `Subset`: a subset of an ordered
#' superset, with its ranks in binary, lexicographic and Gray code order. The
#' fields `$subset`, `$superset`, `$size`, `$superset_size`, `$cardinality`,
#' `$rank_binary`, `$rank_lexicographic` and `$rank_gray` mirror SymPy's. The
#' `subset_next_*()`, `subset_prev_*()` and `subset_iterate_*()` functions
#' step through the subsets in each order; `subset_unrank_binary()` and
#' `subset_unrank_gray()` build the subset of a rank, `subset_from_bitlist()`
#' and `bitlist_from_subset()` convert to and from bit strings, and
#' `subset_indices()` gives the positions of a subset in a superset.
#' `ksubsets(superset, k)` lists the subsets of size `k`.
#'
#' @param subset,superset Vectors or lists.
#' @param a A `Subset`.
#' @param k Steps to move, or for `ksubsets()` the subset size.
#' @param rank A rank.
#' @param super_set A vector or list.
#' @param bitlist A bit string.
#' @return A `Subset`, a bit string, an integer, or a list.
#' @examples
#' a <- Subset(c("c", "d"), c("a", "b", "c", "d"))
#' subset_next_binary(a)$subset
#' a$rank_gray
#' subset_unrank_gray(10, c("a", "b", "c", "d"))$subset
#' ksubsets(1:3, 2)
#' @export
Subset <- function(subset, superset) {
  if (length(subset) > length(superset))
    stop(.sympy_condition("ValueError", "Invalid arguments have been provided. The superset must be larger than the subset."))
  sk <- .subs_keys(superset)
  for (e in as.list(subset)) if (!(.subs_keys(list(e)) %in% sk))
    stop(.sympy_condition("ValueError", paste0("The superset provided is invalid as it does not contain the element ", .subs_elem_str(e))))
  structure(list(subset = subset, superset = superset), class = "Subset")
}

#' @export
`$.Subset` <- function(x, name) {
  sub <- .subset2(x, "subset"); sup <- .subset2(x, "superset")
  switch(name,
         size = length(sub),
         superset_size = length(sup),
         cardinality = .gc_small(gmp::as.bigz(2)^length(sup)),
         rank_binary = .gc_bin2int(bitlist_from_subset(sub, sup)),
         rank_lexicographic = .subs_ranklex(subset_indices(sub, sup), length(sup)),
         rank_gray = {
           bits <- bitlist_from_subset(sub, sup)
           GrayCode(nchar(bits), start = bits)$rank
         },
         .subset2(x, name))
}

.subs_ranklex <- function(idx, n) {
  r <- gmp::as.bigz(0)
  i <- 0L
  while (length(idx) && i <= n) {
    if (i %in% idx) { idx <- idx[-match(i, idx)]; r <- r + 1 }
    else r <- r + gmp::as.bigz(2)^(n - i - 1L)
    i <- i + 1L
  }
  .gc_small(r)
}

#' @export
format.Subset <- function(x, ...) {
  paste0("Subset(", paste(format(.subset2(x, "subset")), collapse = ", "), "; ",
         paste(format(.subset2(x, "superset")), collapse = ", "), ")")
}

#' @export
print.Subset <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
Ops.Subset <- function(e1, e2) {
  eq <- inherits(e2, "Subset") && identical(.subs_keys(e1$subset), .subs_keys(e2$subset)) &&
    identical(.subs_keys(e1$superset), .subs_keys(e2$superset))
  switch(.Generic, "==" = eq, "!=" = !eq,
         stop(.sympy_condition("TypeError", paste0("unsupported operand for ", .Generic))))
}

#' @rdname Subset
#' @export
subset_iterate_binary <- function(a, k) {
  bits <- bitlist_from_subset(a$subset, a$superset)
  n <- (gmp::as.bigz(.gc_bin2int(bits)) + k) %% (gmp::as.bigz(2)^a$superset_size)
  subset_from_bitlist(a$superset, .gc_int2bin(n, a$superset_size))
}

#' @rdname Subset
#' @export
subset_next_binary <- function(a) subset_iterate_binary(a, 1)

#' @rdname Subset
#' @export
subset_prev_binary <- function(a) subset_iterate_binary(a, -1)

#' @rdname Subset
#' @export
subset_next_lexicographic <- function(a) {
  i <- a$superset_size - 1L
  indices <- subset_indices(a$subset, a$superset)
  if (i %in% indices) {
    if ((i - 1L) %in% indices) indices <- indices[-match(i - 1L, indices)]
    else {
      indices <- indices[-match(i, indices)]
      i <- i - 1L
      while (i >= 0L && !(i %in% indices)) i <- i - 1L
      if (i >= 0L) {
        indices <- indices[-match(i, indices)]
        indices <- c(indices, i + 1L)
      }
    }
  } else {
    while (!(i %in% indices) && i >= 0L) i <- i - 1L
    indices <- c(indices, i + 1L)
  }
  Subset(a$superset[indices + 1L], a$superset)
}

#' @rdname Subset
#' @export
subset_prev_lexicographic <- function(a) {
  i <- a$superset_size - 1L
  indices <- subset_indices(a$subset, a$superset)
  while (i >= 0L && !(i %in% indices)) i <- i - 1L
  if (i == 0L || (i - 1L) %in% indices) indices <- indices[-match(i, indices)]
  else {
    if (i >= 0L) {
      indices <- indices[-match(i, indices)]
      indices <- c(indices, i - 1L)
    }
    indices <- c(indices, a$superset_size - 1L)
  }
  Subset(a$superset[indices + 1L], a$superset)
}

#' @rdname Subset
#' @export
subset_iterate_graycode <- function(a, k) {
  code <- graycode_unrank(a$superset_size, (gmp::as.bigz(a$rank_gray) + k) %% (gmp::as.bigz(2)^a$superset_size))
  subset_from_bitlist(a$superset, code)
}

#' @rdname Subset
#' @export
subset_next_gray <- function(a) subset_iterate_graycode(a, 1)

#' @rdname Subset
#' @export
subset_prev_gray <- function(a) subset_iterate_graycode(a, -1)

#' @rdname Subset
#' @export
subset_from_bitlist <- function(super_set, bitlist) {
  bits <- strsplit(bitlist, "", fixed = TRUE)[[1L]]
  if (length(super_set) != length(bits))
    stop(.sympy_condition("ValueError", "The sizes of the lists are not equal"))
  Subset(super_set[bits == "1"], super_set)
}

#' @rdname Subset
#' @export
bitlist_from_subset <- function(subset, superset) {
  bits <- rep("0", length(superset))
  if (inherits(subset, "Subset")) subset <- subset$subset
  bits[subset_indices(subset, superset) + 1L] <- "1"
  paste(bits, collapse = "")
}

#' @rdname Subset
#' @export
subset_unrank_binary <- function(rank, superset) {
  subset_from_bitlist(superset, .gc_int2bin(rank, length(superset)))
}

#' @rdname Subset
#' @export
subset_unrank_gray <- function(rank, superset) {
  subset_from_bitlist(superset, graycode_unrank(length(superset), rank))
}

#' @rdname Subset
#' @export
subset_indices <- function(subset, superset) {
  a <- .subs_keys(superset)
  b <- .subs_keys(subset)
  sb <- unique(b)
  d <- integer()
  done <- FALSE
  for (i in seq_along(a)) {
    if (a[i] %in% sb) {
      d[a[i]] <- i - 1L
      sb <- setdiff(sb, a[i])
      if (!length(sb)) { done <- TRUE; break }
    }
  }
  if (!done) return(integer())
  unname(d[b])
}

#' @rdname Subset
#' @export
ksubsets <- function(superset, k) combinations(superset, k)
