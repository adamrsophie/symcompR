# Multiset partitions by Knuth's Algorithm M (TAOCP 7.1.2.5), after
# sympy/utilities/enumerative.py. The part stack holds components (c, u, v):
# the component number, its unpartitioned amount and its amount in the part.
# Positions are 0 based as in SymPy and shifted by one for R vectors.

.mpt_new <- function(mult) {
  m <- length(mult); n <- sum(mult)
  st <- new.env()
  size <- n * m + 1L
  st$c <- integer(size); st$u <- integer(size); st$v <- integer(size)
  st$f <- integer(n + 1L)
  for (j in seq_len(m)) { st$c[j] <- j - 1L; st$u[j] <- mult[j]; st$v[j] <- mult[j] }
  st$f[1L] <- 0L; st$f[2L] <- m; st$lpart <- 0L
  st$discarded <- 0L; st$k1 <- 0L; st$k2 <- 0L; st$p1 <- 0L
  st
}

# the top part as positions (0 based) of the stack
.mpt_top <- function(st) {
  a <- st$f[st$lpart + 1L]; b <- st$f[st$lpart + 2L]
  if (b > a) a:(b - 1L) else integer()
}

.mpt_decrement <- function(st, part) {
  plen <- length(part)
  for (j in rev(seq_len(plen)) - 1L) {
    pj <- part[j + 1L] + 1L
    if ((j == 0L && st$v[pj] > 1L) || (j > 0L && st$v[pj] > 0L)) {
      st$v[pj] <- st$v[pj] - 1L
      if (j + 1L < plen) for (k in (j + 1L):(plen - 1L)) st$v[part[k + 1L] + 1L] <- st$u[part[k + 1L] + 1L]
      return(TRUE)
    }
  }
  FALSE
}

.mpt_decrement_small <- function(st, part, ub) {
  if (st$lpart >= ub - 1L) { st$p1 <- st$p1 + 1L; return(FALSE) }
  plen <- length(part)
  p0 <- part[1L] + 1L
  for (j in rev(seq_len(plen)) - 1L) {
    pj <- part[j + 1L] + 1L
    if (j == 0L && (st$v[p0] - 1L) * (ub - st$lpart) < st$u[p0]) { st$k1 <- st$k1 + 1L; return(FALSE) }
    if ((j == 0L && st$v[pj] > 1L) || (j > 0L && st$v[pj] > 0L)) {
      st$v[pj] <- st$v[pj] - 1L
      if (j + 1L < plen) for (k in (j + 1L):(plen - 1L)) st$v[part[k + 1L] + 1L] <- st$u[part[k + 1L] + 1L]
      if (plen > 1L && st$v[part[2L] + 1L] == 0L && (st$u[p0] - st$v[p0]) == ((ub - st$lpart - 1L) * st$v[p0])) {
        st$k2 <- st$k2 + 1L
        return(FALSE)
      }
      return(TRUE)
    }
  }
  FALSE
}

.mpt_decrement_large <- function(st, part, amt, lb) {
  if (amt == 1L && !.mpt_decrement(st, part)) return(FALSE)
  min_unalloc <- lb - st$lpart
  if (min_unalloc <= 0L) return(TRUE)
  idx <- part + 1L
  total_mult <- sum(st$u[idx]); total_alloc <- sum(st$v[idx])
  if (total_mult <= min_unalloc) return(FALSE)
  deficit <- min_unalloc - (total_mult - total_alloc)
  if (deficit <= 0L) return(TRUE)
  for (i in rev(seq_along(part)) - 1L) {
    pi <- part[i + 1L] + 1L
    if (i == 0L) {
      if (st$v[pi] > deficit) { st$v[pi] <- st$v[pi] - deficit; return(TRUE) }
      return(FALSE)
    }
    if (st$v[pi] >= deficit) { st$v[pi] <- st$v[pi] - deficit; return(TRUE) }
    deficit <- deficit - st$v[pi]; st$v[pi] <- 0L
  }
  FALSE
}

.mpt_decrement_range <- function(st, part, lb, ub)
  .mpt_decrement_small(st, part, ub) && .mpt_decrement_large(st, .mpt_top(st), 0L, lb)

.mpt_spread <- function(st) {
  j0 <- st$f[st$lpart + 1L]; k <- st$f[st$lpart + 2L]
  base <- k; changed <- FALSE
  # the loop bound is fixed before k moves, as Python's range
  jend <- st$f[st$lpart + 2L]
  if (jend > j0) for (j in j0:(jend - 1L)) {
    st$u[k + 1L] <- st$u[j + 1L] - st$v[j + 1L]
    if (st$u[k + 1L] == 0L) {
      changed <- TRUE
    } else {
      st$c[k + 1L] <- st$c[j + 1L]
      if (changed) st$v[k + 1L] <- st$u[k + 1L]
      else if (st$u[k + 1L] < st$v[j + 1L]) { st$v[k + 1L] <- st$u[k + 1L]; changed <- TRUE }
      else st$v[k + 1L] <- st$v[j + 1L]
      k <- k + 1L
    }
  }
  if (k > base) {
    st$lpart <- st$lpart + 1L
    st$f[st$lpart + 2L] <- k
    return(TRUE)
  }
  FALSE
}

# list_visitor: the partition as lists of components
.mpt_visit <- function(st, components) {
  lapply(0:st$lpart, function(i) {
    a <- st$f[i + 1L]; b <- st$f[i + 2L]
    part <- components[0L]
    if (b > a) for (p in (a:(b - 1L)) + 1L) if (st$v[p] > 0L)
      part <- c(part, rep(components[st$c[p] + 1L], st$v[p]))
    part
  })
}

# factoring_visitor
.mpt_factor_visit <- function(st, primes) {
  vapply(0:st$lpart, function(i) {
    a <- st$f[i + 1L]; b <- st$f[i + 2L]
    fac <- 1
    if (b > a) for (p in (a:(b - 1L)) + 1L) if (st$v[p] > 0L) fac <- fac * primes[st$c[p] + 1L]^st$v[p]
    fac
  }, 0)
}

#' Multiset partition enumeration
#'
#' Knuth's Algorithm M as SymPy's `MultisetPartitionTraverser`:
#' `enum_all(mult)` lists every partition of the multiset with component
#' multiplicities `mult`, `enum_small(mult, ub)` those with fewer than `ub`
#' parts, `enum_large(mult, lb)` those with more than `lb` parts and
#' `enum_range(mult, lb, ub)` those in between; `count_partitions(mult)`
#' counts them. Each partition is a list of parts, each part a vector of
#' component indices (from 0), or of `components` when given.
#' `multiset_partitions_taocp(mult)` is `enum_all()`.
#'
#' @param mult Multiplicities of the components.
#' @param lb,ub Bounds on the number of parts.
#' @param components Labels of the components.
#' @return A list of partitions, or a count.
#' @examples
#' enum_all(c(2, 2))
#' count_partitions(c(2, 2, 1))
#' @name multiset_enumeration
NULL

.mpt_run <- function(mult, mode, lb = 0L, ub = 0L, components = NULL) {
  mult <- as.integer(mult)
  comp <- if (is.null(components)) seq_along(mult) - 1L else components
  out <- list()
  if (mode == "small" && ub <= 0L) return(out)
  if (mode == "large" && lb >= sum(mult)) return(out)
  if (mode == "range" && (ub <= 0L || lb >= sum(mult))) return(out)
  st <- .mpt_new(mult)
  if (mode %in% c("large", "range")) .mpt_decrement_large(st, .mpt_top(st), 0L, lb)
  repeat {
    good <- TRUE
    repeat {
      if (!.mpt_spread(st)) break
      if (mode %in% c("large", "range") && !.mpt_decrement_large(st, .mpt_top(st), 0L, lb)) {
        st$discarded <- st$discarded + 1L; good <- FALSE; break
      }
      if (mode %in% c("small", "range") && st$lpart >= ub) {
        st$discarded <- st$discarded + 1L; good <- FALSE; st$lpart <- ub - 2L; break
      }
    }
    if (good) out[[length(out) + 1L]] <- .mpt_visit(st, comp)
    repeat {
      top <- .mpt_top(st)
      ok <- switch(mode,
        all = .mpt_decrement(st, top),
        small = .mpt_decrement_small(st, top, ub),
        large = .mpt_decrement_large(st, top, 1L, lb),
        range = .mpt_decrement_range(st, top, lb, ub))
      if (ok) break
      if (st$lpart == 0L) return(out)
      st$lpart <- st$lpart - 1L
    }
  }
}

#' @rdname multiset_enumeration
#' @export
enum_all <- function(mult, components = NULL) .mpt_run(mult, "all", components = components)
#' @rdname multiset_enumeration
#' @export
enum_small <- function(mult, ub, components = NULL) .mpt_run(mult, "small", ub = as.integer(ub), components = components)
#' @rdname multiset_enumeration
#' @export
enum_large <- function(mult, lb, components = NULL) .mpt_run(mult, "large", lb = as.integer(lb), components = components)
#' @rdname multiset_enumeration
#' @export
enum_range <- function(mult, lb, ub, components = NULL)
  .mpt_run(mult, "range", lb = as.integer(lb), ub = as.integer(ub), components = components)
#' @rdname multiset_enumeration
#' @export
multiset_partitions_taocp <- function(mult, components = NULL) enum_all(mult, components)
#' @rdname multiset_enumeration
#' @export
count_partitions <- function(mult) length(enum_all(mult))
