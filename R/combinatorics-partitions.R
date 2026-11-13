# Set and integer partitions, after sympy/combinatorics/partitions.py.

# Python's sorted() with a less-than that may fail: CPython runs a binary
# insertion sort on short lists after finding the leading run, so the
# comparisons (and the first one to fail) are the ones SymPy makes.
# ponytail: no merge step, so lists longer than 64 compare differently.
.cp_py_sorted <- function(x, lt) {
  n <- length(x)
  if (n < 2L) return(x)
  run <- 2L
  if (lt(x[[2L]], x[[1L]])) {
    while (run < n && lt(x[[run + 1L]], x[[run]])) run <- run + 1L
    x[seq_len(run)] <- rev(x[seq_len(run)])
  } else {
    while (run < n && !lt(x[[run + 1L]], x[[run]])) run <- run + 1L
  }
  if (run < n) for (s in (run + 1L):n) {
    pivot <- x[[s]]
    l <- 1L; r <- s
    while (l < r) {
      p <- l + ((r - l) %/% 2L)
      if (lt(pivot, x[[p]])) r <- p else l <- p + 1L
    }
    if (l < s) x[(l + 1L):s] <- x[l:(s - 1L)]
    x[[l]] <- pivot
  }
  x
}

# bool(a < b) for SymPy objects: a TypeError when the relation stays symbolic
.py_lt_sym <- function(a, b) {
  r <- Lt(a, b)
  if (.is_true(r)) return(TRUE)
  if (.is_false(r)) return(FALSE)
  stop(.sympy_condition("TypeError", paste0("cannot determine truth value of Relational: ", sstr(r))))
}

# Python list comparison a < b over SymPy elements
.py_list_lt <- function(a, b) {
  for (i in seq_len(min(length(a), length(b)))) {
    if (identical(a[[i]]$key, b[[i]]$key)) next
    return(.py_lt_sym(a[[i]], b[[i]]))
  }
  length(a) < length(b)
}

# sympy.utilities.misc.as_int
.as_int <- function(n) {
  if (is_sym(n)) {
    if (n$type == "Integer") return(gmp::as.bigz(n$value))
    stop(.sympy_condition("ValueError", paste0(sstr(n), " is not an integer")))
  }
  if (inherits(n, "bigz")) return(n)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n != round(n))
    stop(.sympy_condition("ValueError", paste0(if (is.numeric(n)) .pyfloat_str(n) else format(n), " is not an integer")))
  gmp::as.bigz(n)
}

.pyfloat_str <- function(x) {
  s <- format(x, digits = 17)
  for (d in 15:17) { s <- formatC(x, digits = d, format = "g"); if (as.numeric(s) == x) break }
  s <- trimws(s)
  if (!grepl("[.e]", s)) s <- paste0(s, ".0")
  s
}

# Set._infimum_key then ordered(): sets by their numeric infimum, every
# other element at infinity, ties by node count and default_sort_key
.infimum_ordered <- function(els) {
  infk <- vapply(els, function(e) {
    if (e$type == "FiniteSet" && all(vapply(e$args, is_number_type, TRUE)))
      return(min(vapply(e$args, num_value_dbl, 0)))
    if (e$type == "Interval" && is_number_type(e$args[[1L]])) return(num_value_dbl(e$args[[1L]]))
    Inf
  }, 0)
  keys <- lapply(seq_along(els), function(i) list(infk[i], .nodes(els[[i]]), .sort_key(els[[i]])))
  .key_sort(keys, els)
}

.part_block <- function(arg) {
  if (is_sym(arg)) return(list(set = arg, dups = FALSE))
  els <- lapply(as.list(arg), as_sym)
  k <- vapply(els, .k, "")
  if (anyDuplicated(k)) return(list(set = NULL, dups = TRUE))
  list(set = if (length(els)) new_sym("FiniteSet", .infimum_ordered(els)) else S_EmptySet, dups = FALSE)
}

#' Set partitions
#'
#' `Partition(...)` builds the partition of a set into the given disjoint
#' blocks, as SymPy's `Partition`. Each block is a vector, a list or a
#' `FiniteSet`. The fields `$members`, `$size`, `$partition` (the blocks
#' as sorted lists), `$RGS` (the restricted growth string) and `$rank` mirror
#' SymPy's attributes. `p + k` and `p - k` move `k` steps in rank order,
#' `<` and `<=` compare by [partition_sort_key()], `>` and `>=` test for
#' (proper) supersets of blocks as SymPy's sets do, and
#' `partition_from_rgs(rgs, elements)` rebuilds a partition from its
#' restricted growth string.
#'
#' @param ... Blocks.
#' @param p A `Partition`.
#' @param order Ordering for the members (`NULL` keeps their order).
#' @param rgs A restricted growth string (0 based block indices).
#' @param elements The elements, as a vector or list.
#' @return A `Partition`; `partition_sort_key()` returns a sort key.
#' @examples
#' p <- Partition(c(1, 2, 3), 4)
#' p
#' p$RGS
#' p$rank
#' (p + 1)$partition
#' partition_from_rgs(c(0, 1, 2, 0, 1), c(1, 3, 2, 4, 5))
#' @export
Partition <- function(...) {
  blocks <- list()
  dups <- FALSE
  for (arg in list(...)) {
    b <- .part_block(arg)
    if (b$dups) { dups <- TRUE; break }
    blocks[[length(blocks) + 1L]] <- b$set
  }
  if (!all(vapply(blocks, function(b) b$type == "FiniteSet", TRUE)))
    stop(.sympy_condition("ValueError", "Each argument to Partition should be a list, set, or a FiniteSet"))
  allk <- unlist(lapply(blocks, function(b) vapply(b$args, .k, "")))
  if (dups || anyDuplicated(allk))
    stop(.sympy_condition("ValueError", "Partition contained duplicate elements."))
  bk <- vapply(blocks, .k, "")
  blocks <- .infimum_ordered(blocks[!duplicated(bk)])
  members <- .infimum_ordered(unlist(lapply(blocks, `[[`, "args"), recursive = FALSE))
  structure(list(args = blocks, members = members, size = length(members)), class = "Partition")
}

.part_partition <- function(p) {
  ps <- lapply(.subset2(p, "args"), function(b) .key_sort(lapply(b$args, .sort_key), b$args))
  .cp_py_sorted(ps, .py_list_lt)
}

.part_rgs <- function(p) {
  partition <- .part_partition(p)
  keys <- character(); vals <- integer(); els <- list()
  for (i in seq_along(partition)) for (j in partition[[i]]) {
    keys <- c(keys, j$key); vals <- c(vals, i - 1L); els[[length(els) + 1L]] <- j
  }
  o <- .key_sort(lapply(els, .sort_key), seq_along(els))
  vals[unlist(o)]
}

#' @export
`$.Partition` <- function(x, name) {
  switch(name,
         partition = .part_partition(x),
         RGS = .part_rgs(x),
         rank = RGS_rank(.part_rgs(x)),
         .subset2(x, name))
}

.part_str <- function(p) {
  b <- .subset2(p, "args")
  items <- .key_sort(lapply(b, .sort_key), b)
  paste0("Partition(", paste(vapply(items, sstr, ""), collapse = ", "), ")")
}

#' @export
format.Partition <- function(x, ...) .part_str(x)

#' @export
print.Partition <- function(x, ...) {
  cat(.part_str(x), "\n", sep = "")
  invisible(x)
}

# default_sort_key of a Python tuple of SymPy objects
.tuple_sort_key <- function(items) list(list(10, 0, "tuple"), list(length(items), lapply(items, .sort_key)), .ONE_SORT_KEY, 1)

#' @rdname Partition
#' @export
partition_sort_key <- function(p, order = NULL) {
  members <- p$members
  if (!is.null(order)) members <- .key_sort(lapply(members, .sort_key), members)
  list(.sort_key(.int(p$size)), .tuple_sort_key(members), .sort_key(p$rank))
}

#' @rdname Partition
#' @export
partition_from_rgs <- function(rgs, elements) {
  rgs <- as.integer(rgs)
  elements <- as.list(elements)
  if (length(rgs) != length(elements))
    stop(.sympy_condition("ValueError", "mismatch in rgs and element lengths"))
  partition <- rep(list(list()), max(rgs) + 1L)
  for (j in seq_along(rgs)) partition[[rgs[j] + 1L]][[length(partition[[rgs[j] + 1L]]) + 1L]] <- elements[[j]]
  if (!all(lengths(partition) > 0L))
    stop(.sympy_condition("ValueError", "some blocks of the partition were empty."))
  do.call(Partition, partition)
}

.set_superset <- function(a, b) all(vapply(.subset2(b, "args"), .k, "") %in% vapply(.subset2(a, "args"), .k, ""))

#' @export
Ops.Partition <- function(e1, e2) {
  switch(.Generic,
         "+" = .part_add(e1, .as_int(e2)),
         "-" = .part_add(e1, -.as_int(e2)),
         "==" = inherits(e2, "Partition") && setequal(vapply(.subset2(e1, "args"), .k, ""), vapply(.subset2(e2, "args"), .k, "")),
         "!=" = !(inherits(e2, "Partition") && setequal(vapply(.subset2(e1, "args"), .k, ""), vapply(.subset2(e2, "args"), .k, ""))),
         "<" = .key_cmp(partition_sort_key(e1), partition_sort_key(e2)) < 0,
         "<=" = .key_cmp(partition_sort_key(e1), partition_sort_key(e2)) <= 0,
         ">=" = .set_superset(e1, e2),
         ">" = .set_superset(e1, e2) && !setequal(vapply(.subset2(e1, "args"), .k, ""), vapply(.subset2(e2, "args"), .k, "")),
         stop(.sympy_condition("TypeError", paste0("unsupported operand for ", .Generic))))
}

.part_add <- function(p, other) {
  offset <- gmp::as.bigz(.sym_bigz(p$rank)) + other
  size <- p$size
  result <- RGS_unrank(offset %% .sym_bigz(RGS_enum(size)), size)
  partition_from_rgs(result, p$members)
}

.sym_bigz <- function(e) gmp::as.bigz(as.character(e$value))

# ---------------------------------------------------------------------------
# Integer partitions

#' Integer partitions
#'
#' `IntegerPartition(partition, integer)` stores a partition of an integer,
#' as SymPy's `IntegerPartition`. The partition is a vector of positive
#' integers or a named vector of multiplicities (names are the parts); when
#' `integer` is given first, as in SymPy, the parts must add up to it. The
#' fields `$partition`, `$integer` and `$conjugate` mirror SymPy's.
#' `next_lex()` and `prev_lex()` step through the partitions of the integer
#' in lexical order, `as_dict()` gives the multiplicities and `as_ferrers()`
#' draws the Ferrers diagram. Comparisons follow SymPy's (`<` and `<=` compare
#' the parts from smallest to largest).
#'
#' `random_integer_partition(n, seed)` draws a random partition; `seed`
#' follows SymPy's `_randint` (an integer seed or a vector of values to use).
#'
#' @param partition Parts, or the integer when `integer` gives the parts.
#' @param integer The parts, or `NULL`.
#' @param p An `IntegerPartition`.
#' @param char The character for the Ferrers diagram.
#' @param n A positive integer.
#' @param seed A seed for SymPy's `_randint`.
#' @return An `IntegerPartition`; `as_dict()` a named integer vector,
#'   `as_ferrers()` a string and `random_integer_partition()` an integer
#'   vector.
#' @examples
#' a <- IntegerPartition(c(5, 4, 3, 1, 1))
#' a
#' a$conjugate
#' next_lex(IntegerPartition(c(3, 1)))
#' as_dict(a)
#' cat(as_ferrers(IntegerPartition(c(1, 2, 3))))
#' random_integer_partition(10, seed = c(1, 3, 2, 1, 5, 1))
#' @export
IntegerPartition <- function(partition, integer = NULL) {
  if (!is.null(integer)) { tmp <- integer; integer <- partition; partition <- tmp }
  if (inherits(partition, "pydict") || (!is.null(names(partition)) && all(nzchar(names(partition))))) {
    if (inherits(partition, "pydict")) { ks <- partition$keys; vs <- partition$values }
    else { ks <- as.list(as.numeric(names(partition))); vs <- as.list(unname(partition)) }
    kk <- vapply(ks, function(k) as.integer(.as_int(k)), 0L)
    vv <- vapply(vs, function(v) as.integer(.as_int(v)), 0L)
    o <- order(kk, vv, decreasing = TRUE)
    parts <- integer()
    for (i in o) if (vv[i] != 0L) parts <- c(parts, rep(kk[i], vv[i]))
  } else {
    parts <- vapply(as.list(partition), function(k) as.integer(.as_int(k)), 0L)
    parts <- sort(parts, decreasing = TRUE)
  }
  if (is.null(integer)) integer <- sum(parts)
  else {
    integer <- as.integer(.as_int(integer))
    if (sum(parts) != integer) stop(.sympy_condition("ValueError", sprintf("Partition did not add to %d", integer)))
  }
  if (any(parts < 1L)) stop(.sympy_condition("ValueError", "All integer summands must be greater than one"))
  structure(list(partition = parts, integer = integer), class = "IntegerPartition")
}

#' @export
`$.IntegerPartition` <- function(x, name) {
  if (name == "conjugate") return(.ip_conjugate(.subset2(x, "partition")))
  .subset2(x, name)
}

.ip_conjugate <- function(partition) {
  j <- 1L
  temp <- c(partition, 0L)
  k <- temp[1L]
  b <- integer(k)
  while (k > 0L) {
    while (k > temp[j + 1L]) { b[k] <- j; k <- k - 1L }
    j <- j + 1L
  }
  b
}

#' @export
format.IntegerPartition <- function(x, repr = FALSE, ...) {
  p <- .subset2(x, "partition")
  if (repr) return(paste0("IntegerPartition(", .subset2(x, "integer"), ", (", paste(p, collapse = ", "),
                          if (length(p) == 1L) ",", "))"))
  paste0("[", paste(p, collapse = ", "), "]")
}

#' @export
print.IntegerPartition <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

.int_list_lt <- function(a, b) {
  for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(a[i] < b[i])
  length(a) < length(b)
}

#' @export
Ops.IntegerPartition <- function(e1, e2) {
  same <- function() inherits(e2, "IntegerPartition") && .subset2(e1, "integer") == .subset2(e2, "integer") &&
    identical(.subset2(e1, "partition"), .subset2(e2, "partition"))
  ra <- rev(.subset2(e1, "partition")); rb <- rev(.subset2(e2, "partition"))
  switch(.Generic,
         "==" = same(), "!=" = !same(),
         "<" = .int_list_lt(ra, rb), "<=" = !.int_list_lt(rb, ra),
         ">" = .int_list_lt(rb, ra), ">=" = !.int_list_lt(ra, rb),
         stop(.sympy_condition("TypeError", paste0("unsupported operand for ", .Generic))))
}

#' @rdname IntegerPartition
#' @export
as_dict <- function(p) {
  parts <- .subset2(p, "partition")
  keys <- unique(parts)
  stats::setNames(vapply(keys, function(k) sum(parts == k), 0L), keys)
}

#' @rdname IntegerPartition
#' @export
as_ferrers <- function(p, char = "#") paste(strrep(char, .subset2(p, "partition")), collapse = "\n")

# multiplicities d[k] of the parts k = 1..n
.ip_mult <- function(p) {
  d <- integer(.subset2(p, "integer"))
  for (k in .subset2(p, "partition")) d[k] <- d[k] + 1L
  d
}

.ip_from_mult <- function(n, d) {
  ks <- which(d > 0L)
  IntegerPartition(n, stats::setNames(d[ks], ks))
}

#' @rdname IntegerPartition
#' @export
prev_lex <- function(p) {
  n <- .subset2(p, "integer")
  d <- .ip_mult(p)
  keys <- unique(.subset2(p, "partition"))
  m <- length(keys)
  grow <- function(k) if (k > length(d)) d <<- c(d, integer(k - length(d)))
  if (identical(keys, 1L)) return(IntegerPartition(stats::setNames(1L, n)))
  if (keys[m] != 1L) {
    d[keys[m]] <- d[keys[m]] - 1L
    if (keys[m] == 2L) d[1L] <- 2L
    else { d[keys[m] - 1L] <- 1L; d[1L] <- 1L }
  } else {
    d[keys[m - 1L]] <- d[keys[m - 1L]] - 1L
    left <- d[1L] + keys[m - 1L]
    new <- keys[m - 1L]
    d[1L] <- 0L
    while (left) {
      new <- new - 1L
      if (left - new >= 0L) {
        d[new] <- d[new] + left %/% new
        left <- left - d[new] * new
      }
    }
  }
  .ip_from_mult(n, d)
}

#' @rdname IntegerPartition
#' @export
next_lex <- function(p) {
  n <- .subset2(p, "integer")
  d <- c(.ip_mult(p), 0L)
  key <- unique(.subset2(p, "partition"))
  m <- length(key)
  a <- key[m]
  if (a == n) {
    d[] <- 0L
    d[1L] <- n
  } else if (a == 1L) {
    if (d[a] > 1L) {
      d[a + 1L] <- d[a + 1L] + 1L
      d[a] <- d[a] - 2L
    } else {
      b <- key[m - 1L]
      d[b + 1L] <- d[b + 1L] + 1L
      d[1L] <- (d[b] - 1L) * b
      d[b] <- 0L
    }
  } else {
    if (d[a] > 1L) {
      if (m == 1L) {
        d[] <- 0L
        d[a + 1L] <- 1L
        d[1L] <- n - a - 1L
      } else {
        a1 <- a + 1L
        d[a1] <- d[a1] + 1L
        d[1L] <- d[a] * a - a1
        d[a] <- 0L
      }
    } else {
      b <- key[m - 1L]
      b1 <- b + 1L
      d[b1] <- d[b1] + 1L
      need <- d[b] * b + d[a] * a - b1
      d[a] <- 0L; d[b] <- 0L
      d[1L] <- need
    }
  }
  .ip_from_mult(n, d)
}

#' @rdname IntegerPartition
#' @export
random_integer_partition <- function(n, seed = NULL) {
  n <- as.integer(.as_int(n))
  if (n < 1L) stop(.sympy_condition("ValueError", "n must be a positive integer"))
  randint <- .sympy_randint(seed)
  ks <- integer(); ms <- integer()
  while (n > 0L) {
    k <- as.integer(randint(1, n))
    mult <- as.integer(randint(1, n %/% k))
    ks <- c(ks, k); ms <- c(ms, mult)
    n <- n - k * mult
  }
  o <- order(ks, ms, decreasing = TRUE)
  rep(ks[o], ms[o])
}

# ---------------------------------------------------------------------------
# Restricted growth strings

.rgs_generalized <- function(m) {
  d <- matrix(list(gmp::as.bigz(0)), m + 1L, m + 1L)
  for (i in 0:m) d[[1L, i + 1L]] <- gmp::as.bigz(1)
  if (m >= 1L) for (i in 1:m) for (j in 0:(m - 1L)) {
    d[[i + 1L, j + 1L]] <- if (j <= m - i) j * d[[i, j + 1L]] + d[[i, j + 2L]] else gmp::as.bigz(0)
  }
  d
}

#' Restricted growth strings
#'
#' `RGS_generalized(m)` returns the matrix of generalized unrestricted growth
#' string counts, `RGS_enum(m)` the number of restricted growth strings of
#' length `m` (a Bell number), `RGS_unrank(rank, m)` the string of a given
#' rank and `RGS_rank(rgs)` the rank of a string, after
#' `sympy.combinatorics.partitions`.
#'
#' @param m The superset size.
#' @param rank A rank.
#' @param rgs A restricted growth string.
#' @return A `Matrix`, an `Integer` or an integer vector.
#' @examples
#' RGS_generalized(3)
#' RGS_enum(4)
#' RGS_unrank(7, 5)
#' RGS_rank(c(0, 1, 2, 1, 2, 3))
#' @export
RGS_generalized <- function(m) {
  m <- as.integer(m)
  d <- .rgs_generalized(m)
  Matrix(lapply(seq_len(m + 1L), function(i) lapply(seq_len(m + 1L), function(j) .int(d[[i, j]]))))
}

#' @rdname RGS_generalized
#' @export
RGS_enum <- function(m) {
  if (m < 1) return(.int(0))
  if (m == 1) return(.int(1))
  bell(as.integer(m))
}

#' @rdname RGS_generalized
#' @export
RGS_unrank <- function(rank, m) {
  m <- as.integer(m)
  if (m < 1L) stop(.sympy_condition("ValueError", "The superset size must be >= 1"))
  rank <- gmp::as.bigz(if (is_sym(rank)) as.character(rank$value) else rank)
  if (rank < 0 || .sym_bigz(RGS_enum(m)) <= rank) stop(.sympy_condition("ValueError", "Invalid arguments"))
  L <- rep(1L, m + 1L)
  j <- 1L
  D <- .rgs_generalized(m)
  if (m >= 2L) for (i in 2:m) {
    v <- D[[m - i + 1L, j + 1L]]
    cr <- j * v
    if (cr <= rank) {
      L[i + 1L] <- j + 1L
      rank <- rank - cr
      j <- j + 1L
    } else {
      L[i + 1L] <- as.integer(rank %/% v) + 1L
      rank <- rank %% v
    }
  }
  L[-1L] - 1L
}

#' @rdname RGS_generalized
#' @export
RGS_rank <- function(rgs) {
  rgs <- as.integer(rgs)
  n <- length(rgs)
  rank <- gmp::as.bigz(0)
  D <- .rgs_generalized(n)
  if (n >= 2L) for (i in 1:(n - 1L)) {
    k <- n - i - 1L
    m <- max(rgs[1:i])
    rank <- rank + D[[k + 1L, m + 2L]] * rgs[i + 1L]
  }
  .int(rank)
}
