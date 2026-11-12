# Permutations and permutation groups, after sympy/combinatorics/permutations.py,
# perm_groups.py, util.py, named_groups.py, group_constructs.py, generators.py
# (rubik) and the permutation group part of homomorphisms.py. Randomized
# algorithms draw from the Mersenne Twister of R/utilities-pyrandom.R with
# SymPy's call sequence, and code that iterates a Python set follows CPython's
# set order (the emulator below), so results match SymPy's exactly.

# CPython's set table (Objects/setobject.c, 3.12), so that code iterating a
# Python set of ints or of int tuples visits the elements in CPython's order.
# A set is a list: slot keys `k` (NA empty, .PYSET_DUMMY deleted), the 64 bit
# hashes split in 32 bit halves `hi` and `lo`, the elements `obj`, and the
# counters `mask`, `fill`, `used` and `finger`. Elements are scalars, or
# integer vectors standing for tuples when `tuples` is TRUE.

.PYSET_DUMMY <- "\001dummy"
.PYSET_U32 <- 4294967296

.pyset_new <- function(tuples = FALSE, size = 8L) {
  list(k = rep(NA_character_, size), hi = numeric(size), lo = numeric(size),
       obj = vector("list", size), mask = size - 1L, fill = 0L, used = 0L,
       finger = 0L, tuples = tuples)
}

.pyset_key <- function(s, x) if (s$tuples) paste0("t", paste(x, collapse = ",")) else as.character(x)

# hash(int) for a nonnegative int below 2^52, and CPython's xxHash based
# tuplehash for a tuple of such ints, as c(hi, lo) 32 bit halves.
.pyset_hash <- function(s, x) {
  if (!s$tuples) return(c(x %/% .PYSET_U32, x %% .PYSET_U32))
  M <- gmp::as.bigz(2)^64
  acc <- gmp::as.bigz("2870177450012600261")
  for (v in x) {
    acc <- (acc + gmp::as.bigz(v) * gmp::as.bigz("14029467366897019727")) %% M
    acc <- (acc * gmp::as.bigz(2)^31) %% M + acc %/% gmp::as.bigz(2)^33
    acc <- (acc * gmp::as.bigz("11400714785074694791")) %% M
  }
  acc <- (acc + .pyset_xor_len(length(x))) %% M
  if (acc == M - 1) acc <- gmp::as.bigz(1546275796)
  c(as.numeric(acc %/% .PYSET_U32), as.numeric(acc %% .PYSET_U32))
}

# len ^ (XXPRIME_5 ^ 3527539) for a small len
.pyset_xor_len <- function(n) {
  C <- gmp::as.bigz("2870177450013471926")
  lo <- as.numeric(C %% 65536)
  C - lo + bitwXor(as.integer(lo), as.integer(n))
}

# perturb >>= 5 on the 64 bit value c(hi, lo)
.pyset_shift <- function(h) c(h[1L] %/% 32, h[2L] %/% 32 + (h[1L] %% 32) * 134217728)

.pyset_insert_clean <- function(s, x, key, h) {
  mask <- s$mask
  i <- h[2L] %% (mask + 1)
  perturb <- h
  repeat {
    if (is.na(s$k[i + 1])) break
    found <- FALSE
    if (i + 9 <= mask) {
      for (j in 1:9) if (is.na(s$k[i + j + 1])) { i <- i + j; found <- TRUE; break }
    }
    if (found) break
    perturb <- .pyset_shift(perturb)
    i <- (i * 5 + 1 + perturb[2L]) %% (mask + 1)
  }
  s$k[i + 1] <- key; s$hi[i + 1] <- h[1L]; s$lo[i + 1] <- h[2L]; s$obj[i + 1] <- list(x)
  s
}

.pyset_resize <- function(s, minused) {
  newsize <- 8
  while (newsize <= minused) newsize <- newsize * 2
  old <- s
  s$k <- rep(NA_character_, newsize); s$hi <- numeric(newsize); s$lo <- numeric(newsize)
  s$obj <- vector("list", newsize); s$mask <- newsize - 1
  s$fill <- s$used
  for (i in seq_along(old$k)) {
    kk <- old$k[i]
    if (!is.na(kk) && kk != .PYSET_DUMMY) s <- .pyset_insert_clean(s, old$obj[[i]], kk, c(old$hi[i], old$lo[i]))
  }
  s
}

# The slot holding x, or the slot where set_add_entry would put it (as a
# negative number minus one).
.pyset_add <- function(s, x) {
  key <- .pyset_key(s, x)
  h <- .pyset_hash(s, x)
  .pyset_add_entry(s, x, key, h)
}

.pyset_add_entry <- function(s, x, key, h) {
  mask <- s$mask
  i <- h[2L] %% (mask + 1)
  perturb <- h
  freeslot <- NA
  repeat {
    probes <- if (i + 9 <= mask) 9L else 0L
    j <- i
    repeat {
      kk <- s$k[j + 1]
      if (is.na(kk)) {
        if (is.na(freeslot)) {
          s$fill <- s$fill + 1L; s$used <- s$used + 1L
          s$k[j + 1] <- key; s$hi[j + 1] <- h[1L]; s$lo[j + 1] <- h[2L]; s$obj[j + 1] <- list(x)
          if (s$fill * 5 < mask * 3) return(s)
          return(.pyset_resize(s, if (s$used > 50000) s$used * 2 else s$used * 4))
        }
        s$used <- s$used + 1L
        s$k[freeslot + 1] <- key; s$hi[freeslot + 1] <- h[1L]; s$lo[freeslot + 1] <- h[2L]
        s$obj[freeslot + 1] <- list(x)
        return(s)
      }
      if (kk == key) return(s)
      if (kk == .PYSET_DUMMY) freeslot <- j  # the last dummy seen, as CPython 3.12
      j <- j + 1
      if (probes == 0L) break
      probes <- probes - 1L
    }
    perturb <- .pyset_shift(perturb)
    i <- (i * 5 + 1 + perturb[2L]) %% (mask + 1)
  }
}

.pyset_find <- function(s, x) {
  key <- .pyset_key(s, x)
  h <- .pyset_hash(s, x)
  mask <- s$mask
  i <- h[2L] %% (mask + 1)
  perturb <- h
  repeat {
    probes <- if (i + 9 <= mask) 9L else 0L
    j <- i
    repeat {
      kk <- s$k[j + 1]
      if (is.na(kk)) return(NA)
      if (kk == key) return(j)
      j <- j + 1
      if (probes == 0L) break
      probes <- probes - 1L
    }
    perturb <- .pyset_shift(perturb)
    i <- (i * 5 + 1 + perturb[2L]) %% (mask + 1)
  }
}

.pyset_contains <- function(s, x) !is.na(.pyset_find(s, x))

.pyset_discard <- function(s, x) {
  j <- .pyset_find(s, x)
  if (is.na(j)) return(s)
  s$k[j + 1] <- .PYSET_DUMMY; s$obj[j + 1] <- list(NULL); s$used <- s$used - 1L
  s
}

# the elements in iteration order
.pyset_items <- function(s) {
  idx <- which(!is.na(s$k) & s$k != .PYSET_DUMMY)
  if (s$tuples) s$obj[idx] else as.numeric(unlist(s$obj[idx]))
}

.pyset_len <- function(s) s$used

# set.pop(): the element at the search finger
.pyset_pop <- function(s) {
  if (s$used == 0L) stop(structure(class = c("KeyError", "error", "condition"), list(message = "'pop from an empty set'", call = NULL)))
  i <- s$finger %% (s$mask + 1)
  while (is.na(s$k[i + 1]) || s$k[i + 1] == .PYSET_DUMMY) {
    i <- i + 1
    if (i > s$mask) i <- 0
  }
  x <- s$obj[[i + 1]]
  s$k[i + 1] <- .PYSET_DUMMY; s$obj[i + 1] <- list(NULL); s$used <- s$used - 1L
  s$finger <- i + 1
  list(set = s, item = x)
}

# set(iterable) for a list or vector of elements
.pyset_from <- function(xs, tuples = FALSE) {
  s <- .pyset_new(tuples)
  for (x in xs) s <- .pyset_add(s, x)
  s
}

# set_merge(so, other): so.update(other) for a set other
.pyset_merge <- function(s, other) {
  if (other$used == 0L) return(s)
  if ((s$fill + other$used) * 5 >= s$mask * 3) s <- .pyset_resize(s, (s$used + other$used) * 2)
  act <- which(!is.na(other$k) & other$k != .PYSET_DUMMY)
  if (s$fill == 0L && s$mask == other$mask && other$fill == other$used) {
    s$k <- other$k; s$hi <- other$hi; s$lo <- other$lo; s$obj <- other$obj
    s$fill <- other$fill; s$used <- other$used
    return(s)
  }
  if (s$fill == 0L) {
    s$fill <- other$used; s$used <- other$used
    for (i in act) s <- .pyset_insert_clean(s, other$obj[[i]], other$k[i], c(other$hi[i], other$lo[i]))
    return(s)
  }
  for (i in act) s <- .pyset_add_entry(s, other$obj[[i]], other$k[i], c(other$hi[i], other$lo[i]))
  s
}

.pyset_copy <- function(s) .pyset_merge(.pyset_new(s$tuples), s)

# s.union(*others)
.pyset_union <- function(s, ...) {
  r <- .pyset_copy(s)
  for (o in list(...)) r <- .pyset_merge(r, o)
  r
}

# s.difference(other): a new set of the elements of s, in the order of s, not in other
.pyset_difference <- function(s, other) {
  keep <- Filter(function(x) !.pyset_contains(other, x), if (s$tuples) .pyset_items(s) else as.list(.pyset_items(s)))
  .pyset_from(keep, s$tuples)
}

# Iteration order of set(xs) for distinct nonnegative ints xs: sorted when
# every element is below the final table size, else by simulation.
.pyset_int_order <- function(xs) {
  n <- length(unique(xs))
  size <- 8; fill <- 0
  # table size after n insertions without deletions
  for (u in seq_len(n)) {
    if (u * 5 >= (size - 1) * 3) { size <- 8; while (size <= u * 4) size <- size * 2 }
  }
  if (!length(xs) || max(xs) < size) return(sort(unique(as.numeric(xs))))
  .pyset_items(.pyset_from(xs))
}

# ---------------------------------------------------------------------------
# Permutations and cycles, after sympy/combinatorics/permutations.py. A
# Permutation holds its 0 based array form; arrays are composed as SymPy's
# _af_rmul: `a[b + 1L]` is [a[i] for i in b].

.perr <- function(cls, msg) stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))

.perm_af <- function(a) structure(list(array = as.integer(a)), class = "Permutation")
.pkey <- function(p) paste(p$array, collapse = ",")
.pkeys <- function(ps) vapply(ps, .pkey, "")
.perm_in <- function(p, ps) .pkey(p) %in% .pkeys(ps)
.af_invert <- function(a) { inv <- integer(length(a)); inv[a + 1L] <- seq_along(a) - 1L; inv }
.af_rmul <- function(a, b) a[b + 1L]
.af_rmuln <- function(...) {
  l <- list(...)
  if (!length(l)) .perr("ValueError", "String must not be empty")
  r <- l[[length(l)]]
  for (k in rev(seq_len(length(l) - 1L))) r <- l[[k]][r + 1L]
  r
}
.af_pow <- function(a, n) {
  n <- as.numeric(n)
  if (n < 0) { a <- .af_invert(a); n <- -n }
  b <- seq_along(a) - 1L
  while (n > 0) {
    if (n %% 2 == 1) b <- b[a + 1L]
    n <- n %/% 2
    if (n > 0) a <- a[a + 1L]
  }
  b
}
.af_parity <- function(pi) {
  n <- length(pi); seen <- logical(n); c <- 0L
  for (j in seq_len(n)) if (!seen[j]) {
    c <- c + 1L; seen[j] <- TRUE; i <- j
    while (pi[i] != j - 1L) { i <- pi[i] + 1L; seen[i] <- TRUE }
  }
  (n - c) %% 2L
}
.af_commutes_with <- function(a, b) {
  n <- length(a)
  if (n < 2L) return(TRUE)
  k <- seq_len(n - 1L)
  !any(a[b[k] + 1L] != b[a[k] + 1L])
}
# Python's Random._randbelow(n) for 0 < n < 2^31, drawing from the generator
# of R/utilities-pyrandom.R without big integers
.randbelow <- function(n) {
  k <- 0; m <- n
  while (m > 0) { k <- k + 1; m <- m %/% 2 }
  repeat {
    r <- .py_genrand() %/% 2^(32 - k)
    if (r < n) return(as.integer(r))
  }
}
.af_is_identity <- function(a) all(a == seq_along(a) - 1L)
.is_identity <- function(p) .af_is_identity(p$array)
.perm_check <- function(p) if (!inherits(p, "Permutation")) .perr("TypeError", "expected a Permutation")

#' Permutations
#'
#' `Permutation()` builds a permutation of `0..n-1`, following SymPy's
#' constructor: an integer vector of length two or more is the array form
#' (`Permutation(c(2, 0, 1))`), a list gives disjoint or composed cycles
#' (`Permutation(list(c(0, 1), c(2, 3)))`), several integers form one cycle
#' (`Permutation(1, 2)`), a single integer `n` is the identity on `0..n`
#' (`Permutation(3)`), and no argument is the empty permutation. `size`
#' pads with fixed points and must be named. `*` composes (`p * q` applies
#' `p` first, then `q`), `^` with a number is a power, `p ^ q` with a
#' permutation is the conjugate `q^-1 p q`, `i ^ p` is the image of the
#' integer `i`, and `p + k` or `p - k` move `k` steps in lexicographic
#' rank.
#'
#' The helpers mirror the methods of SymPy's `Permutation`:
#' `array_form()`, `cyclic_form()`, `full_cyclic_form()`, `perm_size()`,
#' `perm_list()`, `perm_support()`, `perm_inverse()`, `apply_perm()`,
#' `perm_call()`, `perm_rmul()`, `perm_prod()`, `perm_mul_inv()`,
#' `perm_commutes_with()`, `perm_commutator()`, `perm_transpositions()`,
#' `perm_atoms()`, `perm_order()`, `perm_parity()`, `perm_sign()`,
#' `perm_signature()`, `perm_is_even()`, `perm_is_odd()`,
#' `perm_is_identity()`, `perm_is_singleton()`, `perm_is_empty()`,
#' `perm_cycle_structure()`, `perm_cycles()`, `perm_length()`,
#' `perm_max()`, `perm_min()`, `perm_ascents()`, `perm_descents()`,
#' `perm_runs()`, `perm_index()`, `perm_inversions()`,
#' `perm_inversion_vector()`, `perm_cardinality()`, `perm_resize()`.
#'
#' @param ... For `Permutation()`, the array form, a list of cycles,
#'   integers forming one cycle, a `Permutation` or a `Cycle`. For
#'   `perm_call()`, an index, a vector to permute, or integers of a cycle to
#'   multiply by. For `perm_rmul()`, permutations.
#' @param size Total size; pads with fixed points.
#' @param p,q Permutations.
#' @param i An element or vector of elements.
#' @param perms A list of permutations.
#' @param n A size.
#' @return A `Permutation`, an integer, logical, vector or list.
#' @examples
#' p <- Permutation(c(1, 2, 0))
#' q <- Permutation(list(c(0, 1)), size = 3)
#' p * q
#' perm_order(p)
#' cyclic_form(p * q)
#' Permutation(1, 2, size = 4)
#' @export
Permutation <- function(..., size = NULL) {
  args <- list(...)
  if (!is.null(size)) size <- as.integer(size)
  if (!length(args)) return(.perm_af(seq_len(if (is.null(size)) 0L else size) - 1L))
  if (length(args) > 1L) return(.perm_af(cycle_as_list(do.call(Cycle, args), size)))
  a <- args[[1L]]
  if (inherits(a, "Permutation")) {
    if (is.null(size) || size == length(a$array)) return(a)
    return(Permutation(a$array, size = size))
  }
  if (inherits(a, "Cycle")) return(.perm_af(cycle_as_list(a, size)))
  if (!is.list(a) && length(a) == 1L) {
    a <- as.integer(a)
    if (!is.null(size) && a + 1L > size) .perr("ValueError", sprintf("size is too small when max is %d", a))
    return(.perm_af(seq_len(a + 1L) - 1L))
  }
  if (is.list(a) && length(a)) {
    c <- Cycle()
    for (ci in a) c <- .cycle_call(c, as.integer(ci))
    aform <- cycle_as_list(c)
  } else {
    aform <- as.integer(unlist(a))
    if (anyDuplicated(aform)) .perr("ValueError", "there were repeated elements.")
    if (length(aform) && !setequal(aform, seq_along(aform) - 1L))
      .perr("ValueError", sprintf("Integers 0 through %d must be present.", max(aform)))
    if (!is.null(size) && length(aform) && max(aform) + 1L > size)
      .perr("ValueError", sprintf("max element should not exceed %d", size - 1L))
  }
  if (!is.null(size) && size > length(aform)) aform <- c(aform, seq(length(aform), size - 1L))
  .perm_af(aform)
}

# SymPy's str (and, with repr = TRUE, repr) of a permutation: cyclic
# notation with the largest element shown as a singleton when it is fixed.
.perm_str <- function(p, repr = FALSE) {
  n <- length(p$array)
  if (!n) return(if (repr) "Permutation()" else "()")
  cf <- cyclic_form(p)
  sep <- if (repr) ", " else " "
  s <- paste0(vapply(cf, function(c) paste0("(", paste(c, collapse = sep), ")"), ""), collapse = "")
  if (!((n - 1L) %in% unlist(cf))) s <- paste0("(", n - 1L, ")", s)
  if (repr) paste0("Permutation", s) else s
}

#' @export
format.Permutation <- function(x, repr = FALSE, ...) .perm_str(x, repr)

#' @export
print.Permutation <- function(x, ...) {
  cat(.perm_str(x), "\n", sep = "")
  invisible(x)
}

#' @rdname Permutation
#' @export
array_form <- function(p) p$array

#' @rdname Permutation
#' @export
perm_size <- function(p) length(p$array)

#' @rdname Permutation
#' @export
perm_inverse <- function(p) .perm_af(.af_invert(p$array))

#' @rdname Permutation
#' @export
apply_perm <- function(p, i) p$array[as.integer(i) + 1L]

#' @rdname Permutation
#' @export
cyclic_form <- function(p) {
  a <- p$array
  seen <- logical(length(a))
  out <- list()
  for (s in seq_along(a) - 1L) {
    if (seen[s + 1L]) next
    cyc <- s
    seen[s + 1L] <- TRUE
    j <- a[s + 1L]
    while (!seen[j + 1L]) {
      cyc <- c(cyc, j)
      seen[j + 1L] <- TRUE
      j <- a[j + 1L]
    }
    if (length(cyc) > 1L) out[[length(out) + 1L]] <- cyc
  }
  out
}

#' @rdname Permutation
#' @export
full_cyclic_form <- function(p) {
  cf <- cyclic_form(p)
  need <- setdiff(seq_along(p$array) - 1L, unlist(cf))
  rv <- c(cf, lapply(need, function(i) i))
  rv[order(vapply(rv, `[`, 0L, 1L))]
}

#' @rdname Permutation
#' @export
perm_list <- function(p, size = NULL) {
  rv <- p$array
  n <- length(rv)
  if (!is.null(size)) {
    if (size > n) {
      rv <- c(rv, seq(n, size - 1L))
    } else {
      i <- n - 1L
      while (length(rv) && rv[length(rv)] == i) { rv <- rv[-length(rv)]; i <- i - 1L }
    }
  }
  as.integer(rv)
}

#' @rdname Permutation
#' @export
perm_support <- function(p) which(p$array != seq_along(p$array) - 1L) - 1L

#' @rdname Permutation
#' @export
perm_rmul <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) && !inherits(args[[1L]], "Permutation")) args <- args[[1L]]
  for (a in args) .perm_check(a)
  rv <- args[[1L]]
  for (k in seq_along(args)[-1L]) rv <- args[[k]] * rv
  rv
}

#' @rdname Permutation
#' @export
perm_prod <- function(perms) {
  if (!length(perms)) return(Permutation(integer()))
  do.call(perm_rmul, rev(perms))
}

#' @rdname Permutation
#' @export
perm_mul_inv <- function(p, q) .perm_af(.af_rmul(.af_invert(p$array), q$array))

#' @rdname Permutation
#' @export
perm_commutes_with <- function(p, q) .af_commutes_with(p$array, q$array)

.perm_mul <- function(p, q) {
  a <- p$array; b <- q$array
  if (!length(b)) return(.perm_af(a))
  if (length(b) < length(a)) b <- c(b, seq(length(b), length(a) - 1L))
  .perm_af(c(b[a + 1L], b[setdiff(seq_along(b), seq_along(a))]))
}

.perm_conj <- function(p, h) {
  n <- length(p$array)
  if (n != length(h$array)) .perr("ValueError", "The permutations must be of equal size.")
  a <- integer(n)
  hh <- h$array
  a[hh + 1L] <- hh[p$array + 1L]
  .perm_af(a)
}

.perm_image <- function(p, i) {
  if (i < 0 || i >= length(p$array))
    .perr("TypeError", sprintf("%s should be an integer between 0 and %d", format(i), length(p$array) - 1L))
  p$array[i + 1L]
}

# One Ops method for permutations, groups and cycles, so that R dispatches
# mixed products such as p * G to the same function.
.perm_group_ops <- function(e1, e2) {
  g <- .Generic
  if (g == "*") {
    if (inherits(e1, "Permutation") && inherits(e2, "Permutation")) return(.perm_mul(e1, e2))
    if (inherits(e1, "Permutation") && inherits(e2, "PermutationGroup")) return(Coset(e1, e2, dir = "-"))
    if (inherits(e1, "PermutationGroup") && inherits(e2, "Permutation")) return(Coset(e2, e1, dir = "+"))
    if (inherits(e1, "PermutationGroup") && inherits(e2, "PermutationGroup")) return(.group_mul(e1, e2))
    .perr("TypeError", "unsupported operand types for *")
  }
  if (g == "^") {
    if (inherits(e1, "Permutation") && inherits(e2, "Permutation")) return(.perm_conj(e1, e2))
    if (inherits(e1, "Permutation")) return(.perm_af(.af_pow(e1$array, as.numeric(e2))))
    if (inherits(e2, "Permutation")) {
      if (is.numeric(e1) && e1 == round(e1)) return(.perm_image(e2, e1))
      .perr("NotImplementedError", sprintf("i^p = p(i) when i is an integer, not %s.", format(e1)))
    }
  }
  if (g %in% c("+", "-") && inherits(e1, "Permutation") && !missing(e2)) {
    if (!is.numeric(e2) && !inherits(e2, "bigz"))
      .perr("TypeError", sprintf("unsupported operand type(s) for +: 'int' and '%s'", class(e2)[1L]))
    if (g == "-") e2 <- -e2
    card <- gmp::factorialZ(length(e1$array))
    rk <- (gmp::as.bigz(perm_rank(e1)) + gmp::as.bigz(e2)) %% card
    return(unrank_lex(length(e1$array), rk))
  }
  if (g %in% c("==", "!=")) {
    eq <- FALSE
    if (inherits(e1, "Permutation") && inherits(e2, "Permutation")) eq <- identical(e1$array, e2$array)
    if (inherits(e1, "PermutationGroup") && inherits(e2, "PermutationGroup"))
      eq <- identical(.pkeys(e1$gens), .pkeys(e2$gens))
    if (inherits(e1, "Cycle") && inherits(e2, "Cycle"))
      eq <- setequal(paste(e1$keys, e1$vals), paste(e2$keys, e2$vals))
    return(if (g == "==") eq else !eq)
  }
  .perr("TypeError", paste0("unsupported operator for permutations: ", g))
}

#' @export
Ops.Permutation <- .perm_group_ops

#' @rdname Permutation
#' @export
perm_call <- function(p, ...) {
  i <- list(...)
  if (length(i) == 1L) {
    i <- i[[1L]]
    if (!is.list(i) && length(i) == 1L) return(.perm_image(p, as.numeric(i)))
    if (length(i) != length(p$array))
      .perr("TypeError", sprintf("%s should have the length %d.", .py_list_str(i), length(p$array)))
    return(i[p$array + 1L])
  }
  p * Permutation(do.call(Cycle, i), size = length(p$array))
}

.py_list_str <- function(x) paste0("[", paste(unlist(x), collapse = ", "), "]")

#' @rdname Permutation
#' @export
perm_transpositions <- function(p) {
  res <- list()
  for (x in cyclic_form(p)) {
    nx <- length(x)
    if (nx == 2L) res[[length(res) + 1L]] <- x
    else if (nx > 2L) for (y in rev(x[-1L])) res[[length(res) + 1L]] <- c(x[1L], y)
  }
  res
}

#' @rdname Permutation
#' @param x A character string (its characters) or a vector.
#' @param key An optional function mapping each element to its sort key.
#' @export
perm_from_sequence <- function(x, key = NULL) {
  if (is.character(x) && length(x) == 1L) x <- strsplit(x, "")[[1L]]
  k <- if (is.null(key)) x else vapply(x, key, x[1L], USE.NAMES = FALSE)
  o <- order(k, seq_along(k), method = "radix") - 1L
  perm_inverse(Permutation(o))
}

#' @rdname Permutation
#' @export
perm_atoms <- function(p) sort(unique(p$array))

#' @rdname Permutation
#' @export
perm_next_lex <- function(p) {
  perm <- p$array
  n <- length(perm)
  i <- n - 2L
  # Python's perm[-1] is the last element
  at <- function(k) if (k < 0L) perm[n + k + 1L] else perm[k + 1L]
  while (at(i + 1L) < at(i)) i <- i - 1L
  if (i == -1L) return(NULL)
  j <- n - 1L
  while (perm[j + 1L] < perm[i + 1L]) j <- j - 1L
  t <- perm[j + 1L]; perm[j + 1L] <- perm[i + 1L]; perm[i + 1L] <- t
  i <- i + 1L; j <- n - 1L
  while (i < j) {
    t <- perm[j + 1L]; perm[j + 1L] <- perm[i + 1L]; perm[i + 1L] <- t
    i <- i + 1L; j <- j - 1L
  }
  .perm_af(perm)
}

#' @rdname Permutation
#' @param r A rank.
#' @export
perm_unrank_nonlex <- function(n, r) {
  n <- as.integer(n)
  a <- seq_len(n) - 1L
  r <- gmp::as.bigz(r) %% gmp::factorialZ(n)
  while (n > 0L) {
    k <- as.integer(r %% n)
    t <- a[n]; a[n] <- a[k + 1L]; a[k + 1L] <- t
    r <- r %/% n
    n <- n - 1L
  }
  .perm_af(a)
}

#' @rdname Permutation
#' @export
perm_rank_nonlex <- function(p) {
  perm <- p$array
  inv <- .af_invert(perm)
  if (!length(inv)) return(0)
  n <- length(perm)
  ss <- integer()
  while (n > 1L) {
    s <- perm[n]; t <- inv[n]
    perm[n] <- perm[t + 1L]; perm[t + 1L] <- s
    inv[n] <- inv[s + 1L]; inv[s + 1L] <- t
    ss <- c(ss, s)
    n <- n - 1L
  }
  r <- gmp::as.bigz(0)
  m <- length(p$array) - length(ss) + 1L
  for (k in rev(seq_along(ss))) {
    r <- ss[k] + (m + (length(ss) - k)) * r
  }
  .nt_out(r)
}

#' @rdname Permutation
#' @export
perm_next_nonlex <- function(p) {
  r <- gmp::as.bigz(perm_rank_nonlex(p))
  if (r == gmp::factorialZ(length(p$array)) - 1) return(NULL)
  perm_unrank_nonlex(length(p$array), r + 1)
}

#' @rdname Permutation
#' @export
perm_rank <- function(p) {
  rho <- p$array
  size <- length(rho)
  rank <- gmp::as.bigz(0)
  n <- size - 1L
  psize <- gmp::factorialZ(max(n, 0L))
  for (j in seq_len(max(size - 1L, 0L))) {
    rank <- rank + rho[j] * psize
    if (j < size) {
      k <- (j + 1L):size
      rho[k] <- rho[k] - (rho[k] > rho[j])
    }
    psize <- psize %/% n
    n <- n - 1L
  }
  .nt_out(rank)
}

#' @rdname Permutation
#' @export
perm_cardinality <- function(p) .nt_out(gmp::factorialZ(length(p$array)))

#' @rdname Permutation
#' @export
perm_parity <- function(p) .af_parity(p$array)

#' @rdname Permutation
#' @export
perm_is_even <- function(p) perm_parity(p) == 0L

#' @rdname Permutation
#' @export
perm_is_odd <- function(p) perm_parity(p) == 1L

#' @rdname Permutation
#' @export
perm_is_singleton <- function(p) length(p$array) == 1L

#' @rdname Permutation
#' @export
perm_is_empty <- function(p) length(p$array) == 0L

#' @rdname Permutation
#' @export
perm_is_identity <- function(p) .is_identity(p)

#' @rdname Permutation
#' @export
perm_ascents <- function(p) {
  a <- p$array
  if (length(a) < 2L) return(integer())
  which(a[-length(a)] < a[-1L]) - 1L
}

#' @rdname Permutation
#' @export
perm_descents <- function(p) {
  a <- p$array
  if (length(a) < 2L) return(integer())
  which(a[-length(a)] > a[-1L]) - 1L
}

#' @rdname Permutation
#' @export
perm_max <- function(p) {
  a <- p$array
  if (!length(a)) return(0L)
  m <- a[a != seq_along(a) - 1L]
  if (!length(m)) .perr("ValueError", "max() iterable argument is empty")
  max(m)
}

#' @rdname Permutation
#' @export
perm_min <- function(p) {
  a <- p$array
  if (!length(a)) return(0L)
  m <- a[a != seq_along(a) - 1L]
  if (!length(m)) .perr("ValueError", "min() iterable argument is empty")
  min(m)
}

#' @rdname Permutation
#' @export
perm_inversions <- function(p) {
  a <- p$array
  n <- length(a)
  if (n < 2L) return(0L)
  sum(vapply(seq_len(n - 1L), function(i) sum(a[i] > a[(i + 1L):n]), 0L))
}

#' @rdname Permutation
#' @export
perm_commutator <- function(p, q) {
  a <- p$array; b <- q$array
  n <- length(a)
  if (length(b) != n) .perr("ValueError", "The permutations must be of equal size.")
  inva <- .af_invert(a); invb <- .af_invert(b)
  .perm_af(a[b[inva[invb + 1L] + 1L] + 1L])
}

#' @rdname Permutation
#' @export
perm_signature <- function(p) if (perm_is_even(p)) 1L else -1L

#' @rdname Permutation
#' @export
perm_sign <- function(p) perm_signature(p)

#' @rdname Permutation
#' @export
perm_order <- function(p) {
  lens <- vapply(cyclic_form(p), length, 0L)
  if (!length(lens)) return(1)
  r <- gmp::as.bigz(1)
  for (l in lens) r <- gmp::lcm.bigz(r, l)
  .nt_out(r)
}

#' @rdname Permutation
#' @export
perm_length <- function(p) length(perm_support(p))

#' @rdname Permutation
#' @export
perm_cycle_structure <- function(p) {
  lens <- vapply(cyclic_form(p), length, 0L)
  u <- unique(lens)
  rv <- stats::setNames(vapply(u, function(l) sum(lens == l), 0L), u)
  singletons <- length(p$array) - sum(lens)
  if (singletons) rv <- c(rv, `1` = singletons)
  rv
}

#' @rdname Permutation
#' @export
perm_cycles <- function(p) length(full_cyclic_form(p))

#' @rdname Permutation
#' @export
perm_index <- function(p) {
  a <- p$array
  if (length(a) < 2L) return(0L)
  sum(which(a[-length(a)] > a[-1L]) - 1L)
}

#' @rdname Permutation
#' @export
perm_runs <- function(p) runs(p$array)

#' @rdname Permutation
#' @export
perm_inversion_vector <- function(p) {
  a <- p$array
  n <- length(a)
  vapply(seq_len(max(n - 1L, 0L)), function(i) sum(a[(i + 1L):n] < a[i]), 0L)
}

#' @rdname Permutation
#' @export
perm_rank_trotterjohnson <- function(p) {
  perm <- p$array
  if (!length(perm) || .is_identity(p)) return(0)
  if (identical(perm, c(1L, 0L))) return(1)
  n <- length(perm)
  rank <- gmp::as.bigz(0)
  for (j in seq_len(n - 1L)) {
    k <- 1L; i <- 0L
    while (perm[i + 1L] != j) {
      if (perm[i + 1L] < j) k <- k + 1L
      i <- i + 1L
    }
    j1 <- j + 1L
    rank <- if (rank %% 2 == 0) j1 * rank + j1 - k else j1 * rank + k - 1L
  }
  .nt_out(rank)
}

#' @rdname Permutation
#' @export
perm_unrank_trotterjohnson <- function(size, r) {
  size <- as.integer(size)
  perm <- integer(size)
  r2 <- gmp::as.bigz(0)
  n <- gmp::factorialZ(size)
  pj <- gmp::as.bigz(1)
  rank <- gmp::as.bigz(r)
  for (j in seq_len(max(size - 1L, 0L)) + 1L) {
    pj <- pj * j
    r1 <- (rank * pj) %/% n
    k <- as.integer(r1 - j * r2)
    if (r2 %% 2 == 0) {
      if (k > 0L) for (i in (j - 1L):(j - k)) perm[i + 1L] <- perm[i]
      perm[j - k] <- j - 1L
    } else {
      if (j - 1L > k) for (i in (j - 1L):(k + 1L)) perm[i + 1L] <- perm[i]
      perm[k + 1L] <- j - 1L
    }
    r2 <- r1
  }
  .perm_af(perm)
}

#' @rdname Permutation
#' @export
perm_next_trotterjohnson <- function(p) {
  pi <- p$array
  n <- length(pi)
  st <- 0L
  rho <- pi
  done <- FALSE
  m <- n - 1L
  while (m > 0L && !done) {
    d <- match(m, rho) - 1L
    if (d < m) for (i in d:(m - 1L)) rho[i + 1L] <- rho[i + 2L]
    par <- .af_parity(rho[seq_len(m)])
    if (par == 1L) {
      if (d == m) m <- m - 1L
      else {
        t <- pi[st + d + 1L]; pi[st + d + 1L] <- pi[st + d + 2L]; pi[st + d + 2L] <- t
        done <- TRUE
      }
    } else if (d == 0L) {
      m <- m - 1L
      st <- st + 1L
    } else {
      t <- pi[st + d + 1L]; pi[st + d + 1L] <- pi[st + d]; pi[st + d] <- t
      done <- TRUE
    }
  }
  if (m == 0L) return(NULL)
  .perm_af(pi)
}

.perm_bool_mat <- function(m) {
  n <- nrow(m)
  .new_mat(lapply(as.vector(t(m)), function(v) if (v) S_One else S_Zero), n, n)
}

#' @rdname Permutation
#' @export
perm_precedence_matrix <- function(p) {
  n <- length(p$array); perm <- p$array
  m <- matrix(0L, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) if (j > i) m[perm[i] + 1L, perm[j] + 1L] <- 1L
  .perm_bool_mat(m)
}

#' @rdname Permutation
#' @export
perm_precedence_distance <- function(p, q) {
  n <- length(p$array)
  if (n != length(q$array)) .perr("ValueError", "The permutations must be of equal size.")
  a <- matrix(0L, n, n); b <- matrix(0L, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) if (j > i) {
    a[p$array[i] + 1L, p$array[j] + 1L] <- 1L
    b[q$array[i] + 1L, q$array[j] + 1L] <- 1L
  }
  diag(a) <- 0L
  (n * (n - 1L)) %/% 2L - sum(a * b)
}

.perm_adj <- function(p) {
  n <- length(p$array); perm <- p$array
  m <- matrix(0L, n, n)
  for (i in seq_len(max(n - 1L, 0L))) m[perm[i] + 1L, perm[i + 1L] + 1L] <- 1L
  m
}

#' @rdname Permutation
#' @export
perm_adjacency_matrix <- function(p) .perm_bool_mat(.perm_adj(p))

#' @rdname Permutation
#' @export
perm_adjacency_distance <- function(p, q) {
  n <- length(p$array)
  if (n != length(q$array)) .perr("ValueError", "The permutations must be of the same size.")
  a <- .perm_adj(p) * .perm_adj(q)
  diag(a) <- 0L
  n - sum(a) - 1L
}

#' @rdname Permutation
#' @export
perm_positional_distance <- function(p, q) {
  if (length(p$array) != length(q$array)) .perr("ValueError", "The permutations must be of the same size.")
  sum(abs(p$array - q$array))
}

#' @rdname Permutation
#' @param m,s Josephus parameters: every `m`th of `n` is removed until `s`
#'   remain.
#' @export
perm_josephus <- function(m, n, s = 1) {
  m <- m - 1
  Q <- seq_len(n) - 1L
  perm <- integer()
  while (length(Q) > max(s, 1)) {
    for (dp in seq_len(m)) Q <- c(Q[-1L], Q[1L])
    perm <- c(perm, Q[1L])
    Q <- Q[-1L]
  }
  Permutation(c(perm, Q))
}

#' @rdname Permutation
#' @param inversion An inversion vector.
#' @export
perm_from_inversion_vector <- function(inversion) {
  size <- length(inversion)
  N <- seq_len(size + 1L) - 1L
  perm <- integer()
  for (k in seq_len(size)) {
    v <- inversion[k]
    if (v < 0 || v >= length(N)) .perr("ValueError", "The inversion vector is not valid.")
    val <- N[v + 1L]
    perm <- c(perm, val)
    N <- N[N != val]
  }
  .perm_af(c(perm, N))
}

#' @rdname Permutation
#' @export
perm_random <- function(n) {
  x <- seq_len(n) - 1L
  for (i in rev(seq_len(max(n - 1L, 0L)))) {
    j <- .randbelow(i + 1L)
    t <- x[i + 1L]; x[i + 1L] <- x[j + 1L]; x[j + 1L] <- t
  }
  .perm_af(x)
}

#' @rdname Permutation
#' @export
unrank_lex <- function(n, r) {
  size <- as.integer(n)
  perm <- integer(size)
  rank <- gmp::as.bigz(r)
  psize <- gmp::as.bigz(1)
  for (i in seq_len(size) - 1L) {
    new_psize <- psize * (i + 1L)
    d <- as.integer(rank %% new_psize %/% psize)
    rank <- rank - d * psize
    perm[size - i] <- d
    if (i > 0L) {
      k <- (size - i + 1L):size
      perm[k] <- perm[k] + (perm[k] > d - 1L)
    }
    psize <- new_psize
  }
  .perm_af(perm)
}

#' @rdname Permutation
#' @export
perm_resize <- function(p, n) {
  aform <- p$array
  l <- length(aform)
  if (n > l) return(.perm_af(c(aform, seq(l, n - 1L))))
  if (n < l) {
    new <- list()
    for (cyc in full_cyclic_form(p)) {
      if (min(cyc) <= n - 1L) {
        if (max(cyc) > n - 1L)
          .perr("ValueError", sprintf("The permutation cannot be resized to %d because the cycle %s may break.",
                                      as.integer(n), .py_tuple_str(cyc)))
        new[[length(new) + 1L]] <- cyc
      }
    }
    return(Permutation(new))
  }
  p
}

.py_tuple_str <- function(x) if (length(x) == 1L) paste0("(", x, ",)") else paste0("(", paste(x, collapse = ", "), ")")

# ---------------------------------------------------------------------------
# Cycles: SymPy's Cycle is a dict mapping each moved point to its image,
# with missing keys mapping to themselves; keys keep insertion order.

#' Cycles
#'
#' `Cycle(...)` is SymPy's `Cycle`: `Cycle(1, 2)` (or `Cycle(c(1, 2))`) maps
#' 1 to 2 and 2 to 1. `cycle_mul(c, ...)` is SymPy's call `c(...)`: the
#' product of `c` and the cycle given by `...`, applying `c` first, so
#' `cycle_mul(Cycle(1, 2), 2, 3)` is SymPy's `Cycle(1, 2)(2, 3)`.
#' `cycle_as_list()` is SymPy's `Cycle.list()` and `cycle_size()` its
#' `size`. A `Cycle` converts to a permutation with [Permutation()].
#'
#' @param ... Integers of one cycle, a `Cycle` or a `Permutation`.
#' @param c A `Cycle`.
#' @param size Optional size of the list.
#' @return A `Cycle`, an integer vector or an integer.
#' @examples
#' Cycle(1, 2)
#' cycle_mul(Cycle(1, 2), 2, 3)
#' cycle_as_list(Cycle(1, 2), 4)
#' @export
Cycle <- function(...) {
  args <- list(...)
  c <- structure(list(keys = integer(), vals = integer()), class = "Cycle")
  if (!length(args)) return(c)
  if (length(args) == 1L) {
    a <- args[[1L]]
    if (inherits(a, "Permutation")) {
      for (cy in cyclic_form(a)) c <- .cyc_update(c, .cycle_call(c, cy))
      return(c)
    }
    if (inherits(a, "Cycle")) return(a)
    if (length(a) > 1L) args <- as.list(a)
  }
  if (any(lengths(args) != 1L)) .perr("ValueError", "a cycle takes integers")
  v <- vapply(args, function(z) as.integer(z), 0L)
  if (any(v < 0L)) .perr("ValueError", "negative integers are not allowed in a cycle.")
  if (anyDuplicated(v)) .perr("ValueError", "All elements must be unique in a cycle.")
  c$keys <- v
  c$vals <- c(v[-1L], v[1L])
  c
}

.cyc_get <- function(c, ks) {
  m <- match(ks, c$keys)
  out <- as.integer(ks)
  out[!is.na(m)] <- c$vals[m[!is.na(m)]]
  out
}

.cyc_set <- function(c, k, v) {
  i <- match(k, c$keys)
  if (is.na(i)) { c$keys <- c(c$keys, as.integer(k)); c$vals <- c(c$vals, as.integer(v)) } else c$vals[i] <- as.integer(v)
  c
}

.cyc_update <- function(c, other) {
  for (i in seq_along(other$keys)) c <- .cyc_set(c, other$keys[i], other$vals[i])
  c
}

.cycle_call <- function(c, other) {
  rv <- if (inherits(other, c("Cycle", "Permutation"))) Cycle(other) else do.call(Cycle, as.list(other))
  v <- .cyc_get(rv, .cyc_get(c, c$keys))
  for (i in seq_along(c$keys)) rv <- .cyc_set(rv, c$keys[i], v[i])
  rv
}

#' @rdname Cycle
#' @export
cycle_mul <- function(c, ...) {
  a <- list(...)
  if (length(a) == 1L && (inherits(a[[1L]], c("Cycle", "Permutation")) || length(a[[1L]]) > 1L)) a <- a[[1L]]
  .cycle_call(c, if (is.list(a) && !inherits(a, c("Cycle", "Permutation"))) unlist(a) else a)
}

#' @rdname Cycle
#' @export
cycle_size <- function(c) if (!length(c$keys)) 0L else max(c$keys) + 1L

#' @rdname Cycle
#' @export
cycle_as_list <- function(c, size = NULL) {
  if (!length(c$keys) && is.null(size)) .perr("ValueError", "must give size for empty Cycle")
  if (!is.null(size)) {
    big <- max(c(c$keys[c$vals != c$keys], -1L))
    size <- max(size, big + 1L)
  } else size <- cycle_size(c)
  .cyc_get(c, seq_len(size) - 1L)
}

.cycle_str <- function(c, repr = FALSE) {
  if (!length(c$keys)) return(if (repr) "Cycle()" else "()")
  cycles <- cyclic_form(Permutation(c))
  s <- paste0(vapply(cycles, .py_tuple_str, ""), collapse = "")
  big <- cycle_size(c) - 1L
  if (!(big %in% unlist(cycles))) s <- paste0(s, "(", big, ")")
  if (repr) paste0("Cycle", s) else gsub(",", "", s, fixed = TRUE)
}

#' @export
format.Cycle <- function(x, repr = FALSE, ...) .cycle_str(x, repr)

#' @export
print.Cycle <- function(x, ...) {
  cat(.cycle_str(x), "\n", sep = "")
  invisible(x)
}

#' @export
Ops.Cycle <- .perm_group_ops

# ---------------------------------------------------------------------------
# Helpers of sympy/combinatorics/util.py and the module level functions of
# perm_groups.py (_orbit, _orbits, _orbit_transversal, _stabilizer). A
# transversal is a named list of permutations (or of arrays) whose names
# are the orbit points in SymPy's dict order; an orbit returned as a Python
# set is an integer vector in CPython's set iteration order.

.tr_get <- function(tr, b) tr[[as.character(b)]]
.tr_pts <- function(tr) as.integer(names(tr))
.ident <- function(degree) .perm_af(seq_len(degree) - 1L)

.base_ordering <- function(base, degree) {
  base_len <- length(base)
  ordering <- integer(degree)
  for (i in seq_len(base_len)) ordering[base[i] + 1L] <- i - 1L
  current <- base_len
  for (i in seq_len(degree) - 1L) if (!(i %in% base)) {
    ordering[i + 1L] <- current
    current <- current + 1L
  }
  ordering
}

.check_cycles_alt_sym <- function(perm) {
  n <- length(perm$array); af <- perm$array
  total_len <- 0L
  used <- integer()
  for (i in seq_len(n %/% 2L) - 1L) {
    if (!(i %in% used) && i < n %/% 2L - total_len) {
      current_len <- 1L
      used <- c(used, i)
      j <- i
      while (af[j + 1L] != i) {
        current_len <- current_len + 1L
        j <- af[j + 1L]
        used <- c(used, j)
      }
      total_len <- total_len + current_len
      if (current_len > n %/% 2L && current_len < n - 2L && isTRUE(isprime(current_len))) return(TRUE)
    }
  }
  FALSE
}

.distribute_gens_by_base <- function(base, gens) {
  base_len <- length(base)
  degree <- length(gens[[1L]]$array)
  stabs <- rep(list(list()), base_len)
  if (!base_len) return(stabs)
  max_stab_index <- 0L
  for (gen in gens) {
    j <- 0L
    while (j < base_len - 1L && gen$array[base[j + 1L] + 1L] == base[j + 1L]) j <- j + 1L
    if (j > max_stab_index) max_stab_index <- j
    for (k in 0:j) stabs[[k + 1L]][[length(stabs[[k + 1L]]) + 1L]] <- gen
  }
  if (max_stab_index + 1L < base_len)
    for (i in (max_stab_index + 1L):(base_len - 1L)) stabs[[i + 1L]] <- c(stabs[[i + 1L]], list(.ident(degree)))
  stabs
}

# _orbit_transversal(..., pairs=True, af=True, slp=True): the orbit points
# in discovery order, the transversal arrays, and for each point the
# indices (first occurrence) of the generators of its word.
.orbit_tr <- function(degree, generators, alpha) {
  gens <- lapply(generators, `[[`, "array")
  gk <- vapply(gens, paste, "", collapse = ",")
  gidx <- match(gk, gk) - 1L
  pts <- as.integer(alpha)
  trs <- list(seq_len(degree) - 1L)
  slp <- list(integer())
  used <- logical(degree); used[alpha + 1L] <- TRUE
  k <- 1L
  while (k <= length(pts)) {
    x <- pts[k]; px <- trs[[k]]; pslp <- slp[[k]]
    for (gi in seq_along(gens)) {
      temp <- gens[[gi]][x + 1L]
      if (!used[temp + 1L]) {
        pts <- c(pts, temp)
        trs[[length(trs) + 1L]] <- gens[[gi]][px + 1L]
        slp[[length(slp) + 1L]] <- c(gidx[gi], pslp)
        used[temp + 1L] <- TRUE
      }
    }
    k <- k + 1L
  }
  list(pts = pts, tr = trs, slp = slp)
}

.orbits_transversals_from_bsgs <- function(base, strong_gens_distr, transversals_only = FALSE, slp = FALSE) {
  base_len <- length(base)
  degree <- length(strong_gens_distr[[1L]][[1L]]$array)
  transversals <- vector("list", base_len)
  slps <- vector("list", base_len)
  basic_orbits <- vector("list", base_len)
  for (i in seq_len(base_len)) {
    ot <- .orbit_tr(degree, strong_gens_distr[[i]], base[i])
    transversals[[i]] <- stats::setNames(lapply(ot$tr, .perm_af), ot$pts)
    slps[[i]] <- stats::setNames(ot$slp, ot$pts)
    basic_orbits[[i]] <- ot$pts
  }
  if (transversals_only) return(transversals)
  if (!slp) return(list(basic_orbits = basic_orbits, transversals = transversals))
  list(basic_orbits = basic_orbits, transversals = transversals, slps = slps)
}

.handle_precomputed_bsgs <- function(base, strong_gens, transversals = NULL, basic_orbits = NULL,
                                     strong_gens_distr = NULL) {
  if (is.null(strong_gens_distr)) strong_gens_distr <- .distribute_gens_by_base(base, strong_gens)
  if (is.null(transversals)) {
    if (is.null(basic_orbits)) {
      ot <- .orbits_transversals_from_bsgs(base, strong_gens_distr)
      basic_orbits <- ot$basic_orbits; transversals <- ot$transversals
    } else {
      transversals <- .orbits_transversals_from_bsgs(base, strong_gens_distr, transversals_only = TRUE)
    }
  } else if (is.null(basic_orbits)) {
    basic_orbits <- lapply(transversals, .tr_pts)
  }
  list(transversals = transversals, basic_orbits = basic_orbits, strong_gens_distr = strong_gens_distr)
}

.remove_gens <- function(base, strong_gens, basic_orbits = NULL, strong_gens_distr = NULL) {
  base_len <- length(base)
  degree <- length(strong_gens[[1L]]$array)
  if (is.null(strong_gens_distr)) strong_gens_distr <- .distribute_gens_by_base(base, strong_gens)
  if (is.null(basic_orbits))
    basic_orbits <- lapply(seq_len(base_len), function(i) .orbit(degree, strong_gens_distr[[i]], base[i]))
  strong_gens_distr[[base_len + 1L]] <- list()
  res <- strong_gens
  for (i in rev(seq_len(base_len))) {
    gens_copy <- strong_gens_distr[[i]]
    for (gen in strong_gens_distr[[i]]) {
      if (!.perm_in(gen, strong_gens_distr[[i + 1L]])) {
        temp_gens <- gens_copy
        temp_gens[[match(.pkey(gen), .pkeys(temp_gens))]] <- NULL
        if (!length(temp_gens)) next
        temp_orbit <- .orbit(degree, temp_gens, base[i])
        if (setequal(temp_orbit, basic_orbits[[i]])) {
          gens_copy[[match(.pkey(gen), .pkeys(gens_copy))]] <- NULL
          k <- match(.pkey(gen), .pkeys(res))
          if (is.na(k)) .perr("ValueError", "list.remove(x): x not in list")
          res[[k]] <- NULL
        }
      }
    }
  }
  res
}

.strip <- function(g, base, orbits, transversals) {
  h <- g$array
  base_len <- length(base)
  for (i in seq_len(base_len)) {
    beta <- h[base[i] + 1L]
    if (beta == base[i]) next
    if (!(beta %in% orbits[[i]])) return(list(h = .perm_af(h), j = i))
    u <- .tr_get(transversals[[i]], beta)$array
    h <- .af_invert(u)[h + 1L]
  }
  list(h = .perm_af(h), j = base_len + 1L)
}

# _strip_af: h is an array (or NULL for SymPy's False); transversals hold
# arrays; an slp entry is c(level, index, inverse flag)
.strip_af <- function(h, base, orbits, transversals, j, slp = list(), slps = list()) {
  base_len <- length(base)
  if (j + 1L <= base_len - 1L) for (i in (j + 1L):(base_len - 1L)) {
    beta <- h[base[i + 1L] + 1L]
    if (beta == base[i + 1L]) next
    if (!(beta %in% orbits[[i + 1L]])) return(list(h = h, j = i + 1L, slp = slp))
    u <- .tr_get(transversals[[i + 1L]], beta)
    if (identical(h, u)) return(list(h = NULL, j = base_len + 1L, slp = slp))
    h <- .af_invert(u)[h + 1L]
    if (length(slp)) {
      u_slp <- rev(.tr_get(slps[[i + 1L]], beta))
      slp <- c(lapply(u_slp, function(g) c(i, g, 1L)), slp)
    }
  }
  list(h = h, j = base_len + 1L, slp = slp)
}

.strong_gens_from_distr <- function(strong_gens_distr) {
  if (length(strong_gens_distr) == 1L) return(strong_gens_distr[[1L]])
  result <- strong_gens_distr[[1L]]
  for (gen in strong_gens_distr[[2L]]) if (!.perm_in(gen, result)) result[[length(result) + 1L]] <- gen
  result
}

# BFS orbit of the points alpha (a vector), in discovery order
.orbit_bfs <- function(degree, generators, alpha) {
  gens <- lapply(generators, `[[`, "array")
  orb <- as.integer(alpha)
  used <- logical(degree); used[orb + 1L] <- TRUE
  k <- 1L
  while (k <= length(orb)) {
    b <- orb[k]
    for (gen in gens) {
      temp <- gen[b + 1L]
      if (!used[temp + 1L]) { orb <- c(orb, temp); used[temp + 1L] <- TRUE }
    }
    k <- k + 1L
  }
  orb
}

# _orbit: an integer vector in set order for a point or action "union", a
# list of integer vectors (tuples) for "tuples" and "sets"
.orbit <- function(degree, generators, alpha, action = "tuples") {
  if (length(alpha) == 1L || action == "union")
    return(as.integer(.pyset_int_order(.orbit_bfs(degree, generators, unlist(alpha)))))
  gens <- lapply(generators, `[[`, "array")
  if (action == "tuples") {
    orb <- list(as.integer(alpha)); used <- paste(alpha, collapse = ",")
    k <- 1L
    while (k <= length(orb)) {
      b <- orb[[k]]
      for (gen in gens) {
        temp <- gen[b + 1L]; key <- paste(temp, collapse = ",")
        if (!(key %in% used)) { orb[[length(orb) + 1L]] <- temp; used <- c(used, key) }
      }
      k <- k + 1L
    }
    return(lapply(.pyset_items(.pyset_from(orb, tuples = TRUE)), as.integer))
  }
  if (action == "sets") {
    fs <- function(v) as.integer(.pyset_int_order(v))
    orb <- list(fs(unlist(alpha))); used <- paste(sort(unlist(alpha)), collapse = ",")
    k <- 1L
    while (k <= length(orb)) {
      b <- orb[[k]]
      for (gen in gens) {
        temp <- gen[b + 1L]; key <- paste(sort(temp), collapse = ",")
        if (!(key %in% used)) { orb[[length(orb) + 1L]] <- fs(temp); used <- c(used, key) }
      }
      k <- k + 1L
    }
    return(lapply(.pyset_items(.pyset_from(orb, tuples = TRUE)), as.integer))
  }
  NULL
}

.orbits <- function(degree, generators) {
  orbs <- list()
  sorted_I <- seq_len(degree) - 1L
  while (length(sorted_I)) {
    orb <- .orbit(degree, generators, sorted_I[1L])
    orbs[[length(orbs) + 1L]] <- orb
    sorted_I <- sorted_I[!(sorted_I %in% orb)]
  }
  orbs
}

.stabilizer <- function(degree, generators, alpha) {
  orb <- as.integer(alpha)
  table <- list(); table_inv <- list()
  table[[as.character(alpha)]] <- seq_len(degree) - 1L
  table_inv[[as.character(alpha)]] <- seq_len(degree) - 1L
  used <- logical(degree); used[alpha + 1L] <- TRUE
  gens <- lapply(generators, `[[`, "array")
  stab_gens <- list(); stab_keys <- character()
  k <- 1L
  while (k <= length(orb)) {
    b <- orb[k]
    tb <- table[[as.character(b)]]
    for (gen in gens) {
      temp <- gen[b + 1L]
      if (!used[temp + 1L]) {
        gen_temp <- gen[tb + 1L]
        orb <- c(orb, temp)
        table[[as.character(temp)]] <- gen_temp
        table_inv[[as.character(temp)]] <- .af_invert(gen_temp)
        used[temp + 1L] <- TRUE
      } else {
        sg <- table_inv[[as.character(temp)]][gen[tb + 1L] + 1L]
        key <- paste(sg, collapse = ",")
        if (!(key %in% stab_keys)) { stab_gens[[length(stab_gens) + 1L]] <- sg; stab_keys <- c(stab_keys, key) }
      }
    }
    k <- k + 1L
  }
  lapply(stab_gens, .perm_af)
}

# ---------------------------------------------------------------------------
# Permutation groups, after sympy/combinatorics/perm_groups.py. A group is
# an environment, so that the base and strong generating set, the random
# product replacement state and the cached properties persist between calls
# the way SymPy's attributes do.

.new_group <- function(gens) {
  G <- new.env(parent = emptyenv())
  G$gens <- gens
  G$degree <- length(gens[[1L]]$array)
  G$r <- length(gens)
  G$order <- NULL; G$elements <- list(); G$center <- NULL
  for (f in c("is_abelian", "is_transitive", "is_sym", "is_alt", "is_primitive", "is_nilpotent",
              "is_solvable", "is_trivial", "transitivity_degree", "max_div", "is_perfect",
              "is_cyclic", "is_dihedral")) assign(f, NULL, envir = G)
  G$base <- integer(); G$strong_gens <- list(); G$strong_gens_slp <- list()
  G$basic_orbits <- list(); G$transversals <- list(); G$transversal_slp <- list()
  G$random_gens <- list()
  structure(G, class = "PermutationGroup")
}

.group_check <- function(G) if (!inherits(G, "PermutationGroup")) .perr("TypeError", "expected a PermutationGroup")

#' Permutation groups
#'
#' `PermutationGroup()` is the group generated by permutations, given as
#' arguments or as one list, as SymPy's `PermutationGroup`: generators are
#' resized to a common degree, duplicates (with `dups = TRUE`) and
#' identities (when there are other generators) are dropped. The group
#' keeps SymPy's caches, so the base and strong generating set computed by
#' the Schreier-Sims algorithm, the state of the product replacement
#' generator and the known properties persist on the object.
#'
#' The functions below are SymPy's methods. Properties: `group_generators()`,
#' `group_degree()`, `group_identity()`, `group_order()`, `group_elements()`,
#' `is_abelian()`, `is_cyclic()`, `is_dihedral()`, `is_trivial()`,
#' `is_solvable()`, `is_nilpotent()`, `is_perfect_group()`,
#' `is_polycyclic()`, `is_symmetric_group()`, `is_alternating_group()`,
#' `is_alt_sym()`, `is_transitive()`, `is_primitive()`,
#' `transitivity_degree()`, `max_div()`, `is_elementary()`.
#' Base and strong generating sets: `group_base()`, `strong_gens()`,
#' `basic_orbits()`, `basic_transversals()`, `basic_stabilizers()`,
#' `schreier_sims()`, `schreier_sims_incremental()`,
#' `schreier_sims_random()`, `schreier_vector()`, `baseswap()`.
#' Membership and cosets: `group_contains()`, `coset_factor()`,
#' `coset_rank()`, `coset_unrank()`, `coset_transversal()`, `coset_table()`,
#' `generator_product()`, `is_subgroup()`, `group_index()`, `subgroup()`,
#' `group_equals()`. Orbits: `orbit()`, `orbits()`, `orbit_rep()`,
#' `orbit_transversal()`, `minimal_block()`, `minimal_blocks()`.
#' Subgroups: `stabilizer()`, `pointwise_stabilizer()`, `centralizer()`,
#' `group_center()`, `normal_closure()`, `group_commutator()`,
#' `derived_subgroup()`, `derived_series()`, `lower_central_series()`,
#' `composition_series()`, `sylow_subgroup()`, `subgroup_search()`,
#' `is_normal()`, `quotient_group()`, `abelian_invariants()`.
#' Elements: `generate_dimino()`, `generate_schreier_sims()`,
#' `group_generate()`, `conjugacy_class()`, `conjugacy_classes()`,
#' `molien()`. Random elements, drawn from the generator that
#' [sympy_seed()] seeds: `random_pr()`, `random_stab()`, `group_random()`,
#' `make_perm()`.
#'
#' @param ... Permutations or one list of permutations.
#' @param dups Remove duplicate generators.
#' @param degree Degree to resize generators of different sizes to.
#' @param G,H,K,N Permutation groups.
#' @return A `PermutationGroup`, list, vector, integer or logical.
#' @examples
#' G <- DihedralGroup(4)
#' group_order(G)
#' group_base(SymmetricGroup(4))
#' orbits(PermutationGroup(Permutation(0, 1), Permutation(2, 3)))
#' @export
PermutationGroup <- function(..., dups = TRUE, degree = NULL) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) && !inherits(args[[1L]], c("Permutation", "Cycle")))
    args <- args[[1L]]
  if (!length(args)) args <- list(Permutation())
  args <- lapply(args, function(a) if (inherits(a, "Cycle")) Permutation(a) else a)
  sizes <- vapply(args, function(a) length(a$array), 0L)
  if (!is.null(degree) && length(unique(sizes)) <= 1L)
    .perr("TypeError", "Basic.__new__() got an unexpected keyword argument 'degree'")
  if (length(unique(sizes)) > 1L) {
    if (is.null(degree)) degree <- max(sizes)
    args <- lapply(args, function(a) if (length(a$array) != degree) Permutation(a, size = degree) else a)
  }
  if (dups) args <- args[!duplicated(.pkeys(args))]
  if (length(args) > 1L) {
    ni <- Filter(Negate(.is_identity), args)
    args <- if (length(ni)) ni else args[1L]
  }
  .new_group(args)
}

.group_str <- function(G) paste0("PermutationGroup([\n", paste0("    ", vapply(G$gens, .perm_str, ""), collapse = ",\n"), "])")

#' @export
format.PermutationGroup <- function(x, ...) .group_str(x)

#' @export
print.PermutationGroup <- function(x, ...) {
  cat(.group_str(x), "\n", sep = "")
  invisible(x)
}

#' @export
Ops.PermutationGroup <- .perm_group_ops

#' @rdname PermutationGroup
#' @export
group_generators <- function(G) G$gens

#' @rdname PermutationGroup
#' @export
group_degree <- function(G) G$degree

#' @rdname PermutationGroup
#' @export
group_identity <- function(G) .ident(G$degree)

# G * H: the direct product of SymPy's PermutationGroup.__mul__
.group_mul <- function(G, H) {
  n1 <- G$degree; n2 <- H$degree
  g1 <- lapply(G$gens, function(p) .perm_af(c(p$array, n1 + seq_len(n2) - 1L)))
  g2 <- lapply(H$gens, function(p) .perm_af(c(seq_len(n1) - 1L, p$array + n1)))
  PermutationGroup(c(g1, g2))
}

.random_pr_init <- function(G, r, n, random_prec_n = NULL) {
  deg <- G$degree
  random_gens <- lapply(G$gens, `[[`, "array")
  k <- length(random_gens)
  if (k < r) for (i in k:(r - 1L)) random_gens[[i + 1L]] <- random_gens[[i - k + 1L]]
  random_gens[[length(random_gens) + 1L]] <- seq_len(deg) - 1L
  G$random_gens <- random_gens
  if (is.null(random_prec_n)) {
    for (i in seq_len(n)) random_pr(G)
  } else {
    for (i in seq_len(n)) random_pr(G, .random_prec = random_prec_n[[i]])
  }
  invisible(NULL)
}

#' @rdname PermutationGroup
#' @param gen_count,iterations Size and warm up of the product replacement
#'   state.
#' @param .random_prec Precomputed choices, as SymPy's `_random_prec`.
#' @export
random_pr <- function(G, gen_count = 11, iterations = 50, .random_prec = NULL) {
  if (!length(G$random_gens)) .random_pr_init(G, gen_count, iterations)
  rg <- G$random_gens
  r <- length(rg) - 1L
  if (is.null(.random_prec)) {
    s <- .randbelow(r)
    t <- .randbelow(r - 1L)
    if (t == s) t <- r - 1L
    x <- c(1L, 2L)[.randbelow(2L) + 1L]
    e <- c(-1L, 1L)[.randbelow(2L) + 1L]
  } else {
    s <- .random_prec$s; t <- .random_prec$t
    if (t == s) t <- r - 1L
    x <- .random_prec$x; e <- .random_prec$e
  }
  if (x == 1L) {
    rg[[s + 1L]] <- rg[[s + 1L]][.af_pow(rg[[t + 1L]], e) + 1L]
    rg[[r + 1L]] <- rg[[r + 1L]][rg[[s + 1L]] + 1L]
  } else {
    rg[[s + 1L]] <- .af_pow(rg[[t + 1L]], e)[rg[[s + 1L]] + 1L]
    rg[[r + 1L]] <- rg[[s + 1L]][rg[[r + 1L]] + 1L]
  }
  G$random_gens <- rg
  .perm_af(rg[[r + 1L]])
}

.union_find_rep <- function(num, parents) {
  rep <- num; parent <- parents[num + 1L]
  while (parent != rep) { rep <- parent; parent <- parents[rep + 1L] }
  temp <- num; parent <- parents[num + 1L]
  while (parent != rep) { parents[temp + 1L] <- rep; temp <- parent; parent <- parents[temp + 1L] }
  list(rep = rep, parents = parents)
}

#' @rdname PermutationGroup
#' @export
group_base <- function(G) {
  if (!length(G$base)) schreier_sims(G)
  G$base
}

#' @rdname PermutationGroup
#' @param base A base (0 based points).
#' @param strong_gens A strong generating set.
#' @param pos Position of the base points to swap.
#' @param randomized Use the randomized version.
#' @param transversals,basic_orbits,strong_gens_distr Precomputed data.
#' @export
baseswap <- function(G, base, strong_gens, pos, randomized = FALSE, transversals = NULL,
                     basic_orbits = NULL, strong_gens_distr = NULL) {
  hp <- .handle_precomputed_bsgs(base, strong_gens, transversals, basic_orbits, strong_gens_distr)
  transversals <- hp$transversals; basic_orbits <- hp$basic_orbits; strong_gens_distr <- hp$strong_gens_distr
  base_len <- length(base)
  degree <- G$degree
  i0 <- pos + 1L
  size <- (length(basic_orbits[[i0]]) * length(basic_orbits[[i0 + 1L]])) %/%
    length(.orbit_bfs(degree, strong_gens_distr[[i0]], base[i0 + 1L]))
  T <- if (pos + 2L > base_len - 1L) list() else strong_gens_distr[[pos + 3L]]
  if (isTRUE(randomized)) {
    stab_pos <- PermutationGroup(strong_gens_distr[[i0]])
    sv <- schreier_vector(stab_pos, base[i0 + 1L])
    while (length(.orbit_bfs(degree, T, base[i0])) != size) {
      T[[length(T) + 1L]] <- random_stab(stab_pos, base[i0 + 1L], schreier_vector = sv)
    }
  } else {
    Gamma <- .pyset_from(basic_orbits[[i0]])
    Gamma <- .pyset_discard(Gamma, base[i0])
    Gamma <- .pyset_discard(Gamma, base[i0 + 1L])
    while (length(.orbit_bfs(degree, T, base[i0])) != size) {
      gamma <- .pyset_items(Gamma)[1L]
      x <- .tr_get(transversals[[i0]], gamma)
      temp <- match(base[i0 + 1L], x$array) - 1L
      if (!(temp %in% basic_orbits[[i0 + 1L]])) {
        Gamma <- .pyset_difference(Gamma, .pyset_from(.orbit(degree, T, gamma)))
      } else {
        y <- .tr_get(transversals[[i0 + 1L]], temp)
        el <- .perm_af(x$array[y$array + 1L])
        if (!(el$array[base[i0] + 1L] %in% .orbit_bfs(degree, T, base[i0]))) {
          T[[length(T) + 1L]] <- el
          Gamma <- .pyset_difference(Gamma, .pyset_from(.orbit(degree, T, base[i0])))
        }
      }
    }
  }
  sgnd <- strong_gens_distr
  sgnd[[i0 + 1L]] <- T
  base_new <- base
  base_new[c(i0, i0 + 1L)] <- base[c(i0 + 1L, i0)]
  strong_gens_new <- .strong_gens_from_distr(sgnd)
  for (gen in T) if (!.perm_in(gen, strong_gens_new)) strong_gens_new[[length(strong_gens_new) + 1L]] <- gen
  list(base_new, strong_gens_new)
}

#' @rdname PermutationGroup
#' @export
basic_orbits <- function(G) {
  if (!length(G$basic_orbits)) schreier_sims(G)
  G$basic_orbits
}

#' @rdname PermutationGroup
#' @export
basic_stabilizers <- function(G) {
  if (!length(G$transversals)) schreier_sims(G)
  if (!length(G$base)) return(list())
  lapply(.distribute_gens_by_base(G$base, G$strong_gens), PermutationGroup)
}

#' @rdname PermutationGroup
#' @export
basic_transversals <- function(G) {
  if (!length(G$transversals)) schreier_sims(G)
  G$transversals
}

#' @rdname PermutationGroup
#' @export
composition_series <- function(G) {
  der <- derived_series(G)
  if (!all(vapply(der[[length(der)]]$gens, .is_identity, TRUE)))
    .perr("NotImplementedError", "Group should be solvable")
  series <- list()
  for (i in seq_len(length(der) - 1L)) {
    H <- der[[i + 1L]]
    up_seg <- list()
    for (g in der[[i]]$gens) {
      K <- PermutationGroup(c(list(g), H$gens))
      ord <- gmp::as.bigz(group_order(K)) %/% gmp::as.bigz(group_order(H))
      down_seg <- list()
      fi <- factorint(.nt_out(ord))
      for (pk in seq_along(fi)) {
        p <- as.numeric(names(fi)[pk])
        for (rep_ in seq_len(fi[[pk]])) {
          down_seg[[length(down_seg) + 1L]] <- PermutationGroup(c(list(g), H$gens))
          g <- g^p
        }
      }
      up_seg <- c(down_seg, up_seg)
      H <- K
    }
    up_seg[[1L]] <- der[[i]]
    series <- c(series, up_seg)
  }
  c(series, list(der[[length(der)]]))
}

#' @rdname PermutationGroup
#' @export
coset_transversal <- function(G, H) {
  if (!is_subgroup(H, G)) .perr("ValueError", "The argument must be a subgroup")
  if (group_order(H) == 1) return(group_elements(G))
  .schreier_sims_base(G, group_base(H))
  base <- group_base(G)
  base_ordering <- .base_ordering(base, G$degree)
  identity <- .ident(G$degree)
  transversals <- basic_transversals(G)
  for (l in seq_along(transversals)) {
    t <- unname(transversals[[l]])
    key <- vapply(t, function(x) base_ordering[x$array[base[l] + 1L] + 1L], 0L)
    transversals[[l]] <- t[order(key, method = "radix")]
  }
  orbits <- basic_orbits(H)
  h_stabs <- basic_stabilizers(H)
  g_stabs <- basic_stabilizers(G)
  nz <- min(length(g_stabs), length(h_stabs))
  indices <- Map(function(x, y) gmp::as.bigz(group_order(x)) %/% gmp::as.bigz(group_order(y)),
                 g_stabs[seq_len(nz)], h_stabs[seq_len(nz)])
  T <- if (length(g_stabs) > length(h_stabs)) group_elements(g_stabs[[length(h_stabs) + 1L]]) else list(identity)
  l <- length(h_stabs) - 1L
  t_len <- length(T)
  while (l > -1L) {
    T_next <- list()
    for (u in transversals[[l + 1L]]) {
      if (u == identity) next
      b <- base_ordering[u$array[base[l + 1L] + 1L] + 1L]
      for (t in T) {
        p <- t * u
        if (all(vapply(orbits[[l + 1L]], function(h) base_ordering[p$array[h + 1L] + 1L] >= b, TRUE)))
          T_next[[length(T_next) + 1L]] <- p
        if (t_len + length(T_next) == indices[[l + 1L]]) break
      }
      if (t_len + length(T_next) == indices[[l + 1L]]) break
    }
    T <- c(T, T_next)
    t_len <- t_len + length(T_next)
    l <- l - 1L
  }
  T[[match(.pkey(identity), .pkeys(T))]] <- NULL
  c(list(identity), T)
}

.coset_representative <- function(G, g, H) {
  if (group_order(H) == 1) return(g)
  hb <- group_base(H)
  if (!identical(group_base(G)[seq_along(hb)], hb)) .schreier_sims_base(G, hb)
  orbits <- basic_orbits(H)
  h_transversals <- lapply(basic_transversals(H), unname)
  transversals <- lapply(basic_transversals(G), unname)
  base <- group_base(G)
  base_ordering <- .base_ordering(base, G$degree)
  step <- function(l, x) {
    ob <- orbits[[l + 1L]]
    gamma <- ob[which.min(vapply(ob, function(y) base_ordering[x$array[y + 1L] + 1L], 0L))]
    imgs <- vapply(h_transversals[[l + 1L]], function(h) h$array[base[l + 1L] + 1L], 0L)
    i <- match(gamma, imgs)
    x <- h_transversals[[l + 1L]][[i]] * x
    if (l < length(orbits) - 1L) {
      for (u in transversals[[l + 1L]]) if (u$array[base[l + 1L] + 1L] == x$array[base[l + 1L] + 1L]) break
      x <- step(l + 1L, x * u^-1) * u
    }
    x
  }
  step(0L, g)
}

#' @rdname PermutationGroup
#' @export
coset_table <- function(G, H) {
  if (!is_subgroup(H, G)) .perr("ValueError", "The argument must be a subgroup")
  T <- coset_transversal(G, H)
  n <- length(T)
  Tk <- .pkeys(T)
  A <- unlist(lapply(G$gens, function(gen) list(gen, gen^-1)), recursive = FALSE)
  table <- lapply(seq_len(n), function(i)
    vapply(A, function(x) match(.pkey(.coset_representative(G, T[[i]] * x, H)), Tk) - 1L, 0L))
  nA <- length(A)
  gamma <- 1L
  for (alpha in seq_len(n) - 1L) for (a in seq_len(nA) - 1L) {
    beta <- table[[alpha + 1L]][a + 1L]
    if (beta >= gamma) {
      if (beta > gamma) {
        for (x in seq_len(nA)) {
          z <- table[[gamma + 1L]][x]
          table[[gamma + 1L]][x] <- table[[beta + 1L]][x]
          table[[beta + 1L]][x] <- z
          for (i in seq_len(n)) {
            if (table[[i]][x] == beta) table[[i]][x] <- gamma
            else if (table[[i]][x] == gamma) table[[i]][x] <- beta
          }
        }
      }
      gamma <- gamma + 1L
    }
    if (gamma >= n - 1L) return(table)
  }
  NULL
}

#' @rdname PermutationGroup
#' @export
group_center <- function(G) {
  if (is.null(G$center) || !length(G$center$gens)) G$center <- centralizer(G, G)
  G$center
}

#' @rdname PermutationGroup
#' @param other A group, a list of permutations or a permutation.
#' @export
centralizer <- function(G, other) {
  if (inherits(other, "PermutationGroup")) {
    if (is_trivial(other) || is_trivial(G)) return(G)
    degree <- G$degree
    identity <- .ident(degree)
    orbits <- orbits(other)
    num_orbits <- length(orbits)
    orbits <- orbits[order(-lengths(orbits), method = "radix")]
    long_base <- integer()
    orbit_reps <- integer(num_orbits)
    orbit_reps_indices <- integer(num_orbits)
    orbit_descr <- integer(degree)
    for (i in seq_len(num_orbits)) {
      orbit <- orbits[[i]]
      orbit_reps[i] <- orbit[1L]
      orbit_reps_indices[i] <- length(long_base)
      orbit_descr[orbit + 1L] <- i - 1L
      long_base <- c(long_base, orbit)
    }
    ssi <- schreier_sims_incremental(G, base = long_base)
    base <- ssi[[1L]]; strong_gens <- ssi[[2L]]
    sgd <- .distribute_gens_by_base(base, strong_gens)
    i <- 0L
    for (i in seq_along(base) - 1L) {
      s <- sgd[[i + 1L]]
      if (length(s) == 1L && .is_identity(s[[1L]])) break
    }
    base <- base[seq_len(i)]
    base_len <- i
    for (j in seq_len(num_orbits) - 1L) if (base[base_len] %in% orbits[[j + 1L]]) break
    rel_orbits <- orbits[seq_len(j + 1L)]
    transversals <- lapply(seq_along(rel_orbits), function(j) {
      ot <- .orbit_tr(other$degree, other$gens, orbit_reps[j])
      stats::setNames(lapply(ot$tr, .perm_af), ot$pts)
    })
    tests <- vector("list", base_len)
    for (l in seq_len(base_len) - 1L) {
      if (base[l + 1L] %in% orbit_reps) tests[[l + 1L]] <- function(x) TRUE
      else tests[[l + 1L]] <- local({
        ll <- l
        function(computed_words) {
          g <- computed_words[[ll + 1L]]
          rep_orb_index <- orbit_descr[base[ll + 1L] + 1L]
          rep <- orbit_reps[rep_orb_index + 1L]
          im <- g$array[base[ll + 1L] + 1L]
          im_rep <- g$array[rep + 1L]
          tr_el <- .tr_get(transversals[[rep_orb_index + 1L]], base[ll + 1L])
          if (im_rep >= length(tr_el$array)) .perr("IndexError", "list index out of range")
          im == tr_el$array[im_rep + 1L]
        }
      })
    }
    ogens <- lapply(other$gens, `[[`, "array")
    prop <- function(g) {
      a <- g$array
      if (any(vapply(ogens, function(gen) any(a >= length(gen)), TRUE))) .perr("IndexError", "list index out of range")
      identical(lapply(ogens, function(gen) a[gen + 1L]), lapply(ogens, function(gen) gen[a + 1L]))
    }
    return(subgroup_search(G, prop, base = base, strong_gens = strong_gens, tests = tests))
  }
  if (inherits(other, "Permutation")) return(centralizer(G, PermutationGroup(list(other))))
  if (is.list(other)) return(centralizer(G, PermutationGroup(other)))
  NULL
}

#' @rdname PermutationGroup
#' @export
group_commutator <- function(G, H, K) {
  commutators <- list()
  for (ggen in H$gens) for (hgen in K$gens) {
    comm <- .perm_af(.af_rmuln(hgen$array, ggen$array, .af_invert(hgen$array), .af_invert(ggen$array)))
    if (!.perm_in(comm, commutators)) commutators[[length(commutators) + 1L]] <- comm
  }
  normal_closure(G, commutators)
}

#' @rdname PermutationGroup
#' @param g A permutation (or its array form).
#' @param factor_index Return the base images instead of the factors.
#' @export
coset_factor <- function(G, g, factor_index = FALSE) {
  if (inherits(g, "Permutation")) g <- g$array
  if (inherits(g, "Cycle")) g <- cycle_as_list(g)
  g <- as.integer(g)
  if (length(g) != G$degree) .perr("ValueError", "g should be the same size as permutations of G")
  I <- seq_len(G$degree) - 1L
  bo <- basic_orbits(G)
  transversals <- G$transversals
  factors <- integer()
  base <- group_base(G)
  h <- g
  for (i in seq_along(base)) {
    beta <- h[base[i] + 1L]
    if (beta == base[i]) { factors <- c(factors, beta); next }
    if (!(beta %in% bo[[i]])) return(list())
    u <- .tr_get(transversals[[i]], beta)$array
    h <- .af_invert(u)[h + 1L]
    factors <- c(factors, beta)
  }
  if (!identical(h, I)) return(list())
  if (factor_index) return(as.list(factors))
  tr <- basic_transversals(G)
  lapply(seq_along(base), function(i) .tr_get(tr[[i]], factors[i]))
}

#' @rdname PermutationGroup
#' @param original Express in the original generators.
#' @export
generator_product <- function(G, g, original = FALSE) {
  product <- list()
  if (.is_identity(g)) return(list())
  sg <- strong_gens(G)
  if (.perm_in(g, sg)) {
    if (!original || .perm_in(g, G$gens)) return(list(g))
    for (s in G$strong_gens_slp[[.pkey(g)]]) product <- c(product, generator_product(G, s, original = TRUE))
    return(product)
  } else if (.perm_in(g^-1, sg)) {
    g <- g^-1
    if (!original || .perm_in(g, G$gens)) return(list(g^-1))
    for (s in G$strong_gens_slp[[.pkey(g)]]) product <- c(product, generator_product(G, s, original = TRUE))
    return(lapply(rev(product), function(q) q^-1))
  }
  f <- coset_factor(G, g, TRUE)
  for (i in seq_along(f)) {
    slp <- .tr_get(G$transversal_slp[[i]], f[[i]])
    for (s in slp) {
      if (!original) product[[length(product) + 1L]] <- sg[[s + 1L]]
      else product <- c(product, generator_product(G, sg[[s + 1L]], original = TRUE))
    }
  }
  product
}

#' @rdname PermutationGroup
#' @export
coset_rank <- function(G, g) {
  factors <- coset_factor(G, g, TRUE)
  if (!length(factors)) return(NULL)
  rank <- gmp::as.bigz(0); b <- gmp::as.bigz(1)
  for (i in seq_along(G$base)) {
    j <- match(factors[[i]], G$basic_orbits[[i]]) - 1L
    rank <- rank + b * j
    b <- b * length(G$transversals[[i]])
  }
  .nt_out(rank)
}

#' @rdname PermutationGroup
#' @param rank A rank.
#' @param af Return array forms.
#' @export
coset_unrank <- function(G, rank, af = FALSE) {
  rank <- gmp::as.bigz(rank)
  if (rank < 0 || rank >= gmp::as.bigz(group_order(G))) return(NULL)
  base <- group_base(G)
  transversals <- basic_transversals(G)
  bo <- basic_orbits(G)
  m <- length(base)
  if (!m) {
    h <- seq_len(G$degree) - 1L
    return(if (af) h else .perm_af(h))
  }
  v <- integer(m)
  for (i in seq_len(m)) {
    nt <- length(transversals[[i]])
    c <- as.integer(rank %% nt)
    rank <- rank %/% nt
    v[i] <- bo[[i]][c + 1L]
  }
  a <- lapply(seq_len(m), function(i) .tr_get(transversals[[i]], v[i])$array)
  h <- do.call(.af_rmuln, a)
  if (af) h else .perm_af(h)
}

#' @rdname PermutationGroup
#' @export
group_elements <- function(G) {
  if (!length(G$elements)) G$elements <- group_generate(G)
  G$elements
}

#' @rdname PermutationGroup
#' @export
derived_series <- function(G) {
  res <- list(G)
  current <- G
  nxt <- derived_subgroup(G)
  while (!is_subgroup(current, nxt)) {
    res[[length(res) + 1L]] <- nxt
    current <- nxt
    nxt <- derived_subgroup(nxt)
  }
  res
}

#' @rdname PermutationGroup
#' @export
derived_subgroup <- function(G) {
  r <- G$r
  gens <- lapply(G$gens, `[[`, "array")
  degree <- G$degree
  s <- .pyset_new(tuples = TRUE)
  for (i in seq_len(r)) for (j in seq_len(r)) {
    p1 <- gens[[i]]; p2 <- gens[[j]]
    c <- integer(degree)
    c[p2[p1 + 1L] + 1L] <- p1[p2 + 1L]
    s <- .pyset_add(s, c)
  }
  cms <- lapply(.pyset_items(s), .perm_af)
  normal_closure(G, cms)
}

#' @rdname PermutationGroup
#' @param method `"coset"` (Schreier-Sims) or `"dimino"`.
#' @export
group_generate <- function(G, method = "coset", af = FALSE) {
  if (method == "coset") return(generate_schreier_sims(G, af))
  if (method == "dimino") return(generate_dimino(G, af))
  .perr("NotImplementedError", sprintf("No generation defined for %s", method))
}

#' @rdname PermutationGroup
#' @export
generate_dimino <- function(G, af = FALSE) {
  idn <- seq_len(G$degree) - 1L
  element_list <- list(idn)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  assign(paste(idn, collapse = ","), TRUE, envir = seen)
  out <- list(idn)
  gens <- lapply(G$gens, `[[`, "array")
  for (i in seq_along(gens)) {
    D <- element_list
    N <- list(idn)
    while (length(N)) {
      A <- N
      N <- list()
      for (a in A) for (g in gens[seq_len(i)]) {
        ag <- a[g + 1L]
        if (is.null(seen[[paste(ag, collapse = ",")]])) {
          for (d in D) {
            ap <- d[ag + 1L]
            out[[length(out) + 1L]] <- ap
            element_list[[length(element_list) + 1L]] <- ap
            assign(paste(ap, collapse = ","), TRUE, envir = seen)
            N[[length(N) + 1L]] <- ap
          }
        }
      }
    }
  }
  G$order <- length(element_list)
  if (af) out else lapply(out, .perm_af)
}

#' @rdname PermutationGroup
#' @export
generate_schreier_sims <- function(G, af = FALSE) {
  n <- G$degree
  u <- basic_transversals(G)
  bo <- G$basic_orbits
  wrap <- if (af) function(a) a else .perm_af
  if (!length(u)) return(lapply(G$gens, function(x) wrap(x$array)))
  if (length(u) == 1L) return(lapply(bo[[1L]], function(i) wrap(.tr_get(u[[1L]], i)$array)))
  u <- rev(u)
  bo <- rev(bo)
  stg <- list(seq_len(n) - 1L)
  posmax <- lengths(u)
  n1 <- length(posmax) - 1L
  pos <- integer(n1)
  h <- 0L
  out <- list()
  repeat {
    if (pos[h + 1L] >= posmax[h + 1L]) {
      if (h == 0L) return(out)
      pos[h + 1L] <- 0L
      h <- h - 1L
      stg[[length(stg)]] <- NULL
      next
    }
    p <- .tr_get(u[[h + 1L]], bo[[h + 1L]][pos[h + 1L] + 1L])$array[stg[[length(stg)]] + 1L]
    pos[h + 1L] <- pos[h + 1L] + 1L
    stg[[length(stg) + 1L]] <- p
    h <- h + 1L
    if (h == n1) {
      last <- stg[[length(stg)]]
      for (i in bo[[length(bo)]]) out[[length(out) + 1L]] <- wrap(.tr_get(u[[length(u)]], i)$array[last + 1L])
      stg[[length(stg)]] <- NULL
      h <- h - 1L
    }
  }
}

#' @rdname PermutationGroup
#' @param strict Require equal degrees.
#' @export
group_contains <- function(G, g, strict = TRUE) {
  if (!inherits(g, "Permutation")) return(FALSE)
  if (length(g$array) != G$degree) {
    if (strict) return(FALSE)
    g <- Permutation(g, size = G$degree)
  }
  if (.perm_in(g, G$gens)) return(TRUE)
  length(coset_factor(G, g$array, TRUE)) > 0L
}

#' @rdname PermutationGroup
#' @export
is_perfect_group <- function(G) {
  if (is.null(G$is_perfect)) G$is_perfect <- group_equals(G, derived_subgroup(G))
  G$is_perfect
}

#' @rdname PermutationGroup
#' @export
is_abelian <- function(G) {
  if (inherits(G, c("FpGroup", "FreeGroup"))) return(fp_is_abelian(G))
  if (!is.null(G$is_abelian)) return(G$is_abelian)
  G$is_abelian <- TRUE
  gens <- lapply(G$gens, `[[`, "array")
  for (x in gens) for (y in gens) {
    if (.py_list_le(y, x)) next
    if (!.af_commutes_with(x, y)) {
      G$is_abelian <- FALSE
      return(FALSE)
    }
  }
  TRUE
}

# Python's list comparison y <= x
.py_list_le <- function(y, x) {
  n <- min(length(x), length(y))
  d <- which(y[seq_len(n)] != x[seq_len(n)])
  if (length(d)) return(y[d[1L]] < x[d[1L]])
  length(y) <= length(x)
}

#' @rdname PermutationGroup
#' @export
abelian_invariants <- function(G) {
  if (is_trivial(G)) return(integer())
  gns <- G$gens
  inv <- numeric()
  H <- derived_subgroup(G)
  Hgens <- H$gens
  Gc <- G
  for (p in primefactors(group_order(G))) {
    ranks <- numeric()
    repeat {
      pows <- list()
      for (g in gns) {
        elm <- g^p
        if (!group_contains(H, elm)) pows[[length(pows) + 1L]] <- elm
      }
      K <- if (length(pows)) PermutationGroup(c(Hgens, pows)) else H
      r <- gmp::as.bigz(group_order(Gc)) %/% gmp::as.bigz(group_order(K))
      Gc <- K
      gns <- pows
      if (r == 1) break
      ranks <- c(ranks, multiplicity(p, .nt_out(r)))
    }
    if (length(ranks)) {
      pw <- rep(1, ranks[1L])
      for (i in ranks) for (j in seq_len(i)) pw[j] <- pw[j] * p
      inv <- c(inv, pw)
    }
  }
  sort(inv)
}

#' @rdname PermutationGroup
#' @param p A prime.
#' @export
is_elementary <- function(G, p) is_abelian(G) && all(vapply(G$gens, function(g) perm_order(g) == p, TRUE))

.eval_is_alt_sym_naive <- function(G, only_sym = FALSE, only_alt = FALSE) {
  if (only_sym && only_alt) .perr("ValueError", sprintf("Both %s and %s cannot be set to True", "True", "True"))
  n <- G$degree
  sym_order <- gmp::factorialZ(n)
  ord <- gmp::as.bigz(group_order(G))
  if (ord == sym_order) { G$is_sym <- TRUE; G$is_alt <- FALSE; return(!only_alt) }
  if (2 * ord == sym_order) { G$is_sym <- FALSE; G$is_alt <- TRUE; return(!only_sym) }
  FALSE
}

.eval_is_alt_sym_monte_carlo <- function(G, eps = 0.05, perms = NULL) {
  if (is.null(perms)) {
    n <- G$degree
    c_n <- if (n < 17) 0.34 else 0.57
    d_n <- c_n * log(2) / log(n)
    N_eps <- as.integer(-log(eps) / d_n)
    for (i in seq_len(N_eps)) if (.check_cycles_alt_sym(random_pr(G))) return(TRUE)
    return(FALSE)
  }
  for (perm in perms) if (.check_cycles_alt_sym(perm)) return(TRUE)
  FALSE
}

#' @rdname PermutationGroup
#' @param eps Error bound of the Monte Carlo test.
#' @export
is_alt_sym <- function(G, eps = 0.05, .random_prec = NULL) {
  if (!is.null(.random_prec)) {
    N_eps <- .random_prec$N_eps
    return(.eval_is_alt_sym_monte_carlo(G, perms = .random_prec$perms[seq_len(N_eps)]))
  }
  if (isTRUE(G$is_sym) || isTRUE(G$is_alt)) return(TRUE)
  if (identical(G$is_sym, FALSE) && identical(G$is_alt, FALSE)) return(FALSE)
  n <- G$degree
  if (n < 8) return(.eval_is_alt_sym_naive(G))
  if (is_transitive(G)) return(.eval_is_alt_sym_monte_carlo(G, eps = eps))
  G$is_sym <- FALSE; G$is_alt <- FALSE
  FALSE
}

#' @rdname PermutationGroup
#' @export
is_nilpotent <- function(G) {
  if (is.null(G$is_nilpotent)) {
    lcs <- lower_central_series(G)
    if (all(vapply(lcs[[length(lcs)]]$gens, .is_identity, TRUE))) {
      G$is_solvable <- TRUE; G$is_nilpotent <- TRUE
      return(TRUE)
    }
    G$is_nilpotent <- FALSE
    return(FALSE)
  }
  G$is_nilpotent
}

.group_copy <- function(G) PermutationGroup(G$gens)

#' @rdname PermutationGroup
#' @param gr The group to test normality in.
#' @export
is_normal <- function(G, gr, strict = TRUE) {
  if (!is_subgroup(G, gr, strict = strict)) return(FALSE)
  d_self <- G$degree; d_gr <- gr$degree
  if (is_trivial(G) && (d_self == d_gr || !strict)) return(TRUE)
  if (isTRUE(gr$is_abelian)) return(TRUE)
  new_self <- .group_copy(G)
  if (!strict && d_self != d_gr) {
    if (d_self < d_gr) new_self <- PermutationGroup(c(new_self$gens, list(Permutation(d_gr - 1L))))
    else gr <- PermutationGroup(c(gr$gens, list(Permutation(d_self - 1L))))
  }
  gens2 <- lapply(new_self$gens, `[[`, "array")
  gens1 <- lapply(gr$gens, `[[`, "array")
  for (g1 in gens1) for (g2 in gens2) {
    p <- .af_rmuln(g1, g2, .af_invert(g1))
    if (!length(coset_factor(new_self, p, TRUE))) return(FALSE)
  }
  TRUE
}

#' @rdname PermutationGroup
#' @export
is_primitive <- function(G, randomized = TRUE) {
  if (!is.null(G$is_primitive)) return(G$is_primitive)
  if (identical(is_transitive(G), FALSE)) return(FALSE)
  if (randomized) {
    v <- schreier_vector(G, 0L)
    stab <- PermutationGroup(lapply(seq_along(G$gens), function(i) random_stab(G, 0L, v)))
  } else stab <- stabilizer(G, 0L)
  for (orb in orbits(stab)) {
    x <- orb[1L]
    if (x != 0L && any(minimal_block(G, c(0L, x)) != 0L)) {
      G$is_primitive <- FALSE
      return(FALSE)
    }
  }
  G$is_primitive <- TRUE
  TRUE
}

#' @rdname PermutationGroup
#' @export
minimal_blocks <- function(G, randomized = TRUE) {
  number_blocks <- function(blocks) {
    u <- unique(blocks)
    match(blocks, u) - 1L
  }
  if (!is_transitive(G)) return(FALSE)
  blocks <- list(); num_blocks <- list(); rep_blocks <- list()
  if (randomized) {
    v <- schreier_vector(G, 0L)
    stab <- PermutationGroup(lapply(seq_along(G$gens), function(i) random_stab(G, 0L, v)))
  } else stab <- stabilizer(G, 0L)
  for (orb in orbits(stab)) {
    x <- orb[1L]
    if (x != 0L) {
      block <- minimal_block(G, c(0L, x))
      num_block <- number_blocks(block)
      rep <- which(num_block == 0L) - 1L
      minimal <- TRUE
      mask <- logical(length(blocks))
      for (i in seq_along(rep_blocks)) {
        r <- rep_blocks[[i]]
        if (length(r) > length(rep) && all(rep %in% r)) mask[i] <- TRUE
        else if (length(r) < length(rep) && all(r %in% rep)) { minimal <- FALSE; break }
      }
      blocks <- blocks[!mask]; num_blocks <- num_blocks[!mask]; rep_blocks <- rep_blocks[!mask]
      if (minimal && !any(vapply(num_blocks, identical, TRUE, num_block))) {
        blocks[[length(blocks) + 1L]] <- block
        num_blocks[[length(num_blocks) + 1L]] <- num_block
        rep_blocks[[length(rep_blocks) + 1L]] <- rep
      }
    }
  }
  blocks
}

#' @rdname PermutationGroup
#' @export
is_solvable <- function(G) {
  if (is.null(G$is_solvable)) {
    if (gmp::as.bigz(group_order(G)) %% 2 != 0) return(TRUE)
    ds <- derived_series(G)
    if (all(vapply(ds[[length(ds)]]$gens, .is_identity, TRUE))) {
      G$is_solvable <- TRUE
      return(TRUE)
    }
    G$is_solvable <- FALSE
    return(FALSE)
  }
  G$is_solvable
}

#' @rdname PermutationGroup
#' @export
is_subgroup <- function(G, H, strict = TRUE) {
  if (inherits(H, "SymmetricPermutationGroup")) return(G$degree == H$degree)
  if (!inherits(H, "PermutationGroup")) return(FALSE)
  if (G == H || identical(G$gens[[1L]]$array, integer())) return(TRUE)
  if (gmp::as.bigz(group_order(H)) %% gmp::as.bigz(group_order(G)) != 0) return(FALSE)
  if (G$degree == H$degree || (G$degree < H$degree && !strict)) gens <- G$gens else return(FALSE)
  all(vapply(gens, function(g) group_contains(H, g, strict = strict), TRUE))
}

#' @rdname PermutationGroup
#' @export
is_polycyclic <- function(G) is_solvable(G)

#' @rdname PermutationGroup
#' @export
is_transitive <- function(G, strict = TRUE) {
  if (isTRUE(G$is_transitive)) return(TRUE)
  if (strict) {
    if (!is.null(G$is_transitive)) return(G$is_transitive)
    ans <- length(orbit(G, 0L)) == G$degree
    G$is_transitive <- ans
    return(ans)
  }
  got_orb <- FALSE
  for (x in orbits(G)) if (length(x) > 1L) {
    if (got_orb) return(FALSE)
    got_orb <- TRUE
  }
  got_orb
}

#' @rdname PermutationGroup
#' @export
is_trivial <- function(G) {
  if (is.null(G$is_trivial)) G$is_trivial <- length(G$gens) == 1L && .is_identity(G$gens[[1L]])
  G$is_trivial
}

#' @rdname PermutationGroup
#' @export
lower_central_series <- function(G) {
  res <- list(G)
  current <- G
  nxt <- group_commutator(G, G, current)
  while (!is_subgroup(current, nxt)) {
    res[[length(res) + 1L]] <- nxt
    current <- nxt
    nxt <- group_commutator(G, G, current)
  }
  res
}

#' @rdname PermutationGroup
#' @export
max_div <- function(G) {
  if (!is.null(G$max_div)) return(G$max_div)
  n <- G$degree
  if (n == 1L) return(1L)
  for (x in 2:n) if (n %% x == 0L) {
    d <- n %/% x
    G$max_div <- d
    return(d)
  }
}

#' @rdname PermutationGroup
#' @param points Points of the block or of the stabilizer.
#' @export
minimal_block <- function(G, points) {
  if (!is_transitive(G)) return(FALSE)
  n <- G$degree
  gens <- G$gens
  parents <- seq_len(n) - 1L
  ranks <- rep(1L, n)
  not_rep <- integer()
  k <- length(points)
  if (k > max_div(G)) return(integer(n))
  for (i in seq_len(k - 1L)) {
    parents[points[i + 1L] + 1L] <- points[1L]
    not_rep <- c(not_rep, points[i + 1L])
  }
  ranks[points[1L] + 1L] <- k
  i <- 0L
  len_not_rep <- k - 1L
  while (i < len_not_rep) {
    gamma <- not_rep[i + 1L]
    i <- i + 1L
    for (gen in gens) {
      ur <- .union_find_rep(gamma, parents); parents <- ur$parents; delta <- ur$rep
      r1 <- .union_find_rep(gen$array[gamma + 1L], parents); parents <- r1$parents
      r2 <- .union_find_rep(gen$array[delta + 1L], parents); parents <- r2$parents
      temp <- 0L
      if (r1$rep != r2$rep) {
        if (ranks[r1$rep + 1L] >= ranks[r2$rep + 1L]) { new1 <- r1$rep; new2 <- r2$rep } else { new1 <- r2$rep; new2 <- r1$rep }
        total_rank <- ranks[new1 + 1L] + ranks[new2 + 1L]
        if (total_rank > max_div(G)) return(integer(n))
        parents[new2 + 1L] <- new1
        ranks[new1 + 1L] <- total_rank
        not_rep <- c(not_rep, new2)
        temp <- 1L
      }
      len_not_rep <- len_not_rep + temp
    }
  }
  for (i in seq_len(n) - 1L) parents <- .union_find_rep(i, parents)$parents
  u <- unique(parents)
  first <- match(u, parents) - 1L
  first[match(parents, u)]
}

#' @rdname PermutationGroup
#' @param x A permutation.
#' @export
conjugacy_class <- function(G, x) {
  new_class <- list(x); keys <- .pkey(x)
  last_iteration <- list(x)
  while (length(last_iteration)) {
    this_iteration <- list(); tkeys <- character()
    for (y in last_iteration) for (s in G$gens) {
      conjugated <- s * y * perm_inverse(s)
      k <- .pkey(conjugated)
      if (!(k %in% keys) && !(k %in% tkeys)) {
        this_iteration[[length(this_iteration) + 1L]] <- conjugated
        tkeys <- c(tkeys, k)
      }
    }
    for (y in last_iteration) if (!(.pkey(y) %in% keys)) { new_class[[length(new_class) + 1L]] <- y; keys <- c(keys, .pkey(y)) }
    last_iteration <- this_iteration
  }
  new_class
}

#' @rdname PermutationGroup
#' @export
conjugacy_classes <- function(G) {
  identity <- .ident(G$degree)
  known <- .pkey(identity)
  classes <- list(list(identity))
  for (x in group_generate(G)) {
    if (!(.pkey(x) %in% known)) {
      nc <- conjugacy_class(G, x)
      classes[[length(classes) + 1L]] <- nc
      known <- c(known, .pkeys(nc))
    }
  }
  classes
}

#' @rdname PermutationGroup
#' @param t The variable of the Molien series (default the symbol `t`).
#' @export
molien <- function(G, t = NULL) {
  t <- if (is.null(t)) Symbol("t") else as_sym(t)
  if (G$degree == 0L) return(S_One)
  ord <- group_order(G)
  terms <- list()
  for (cc in conjugacy_classes(G)) {
    cs <- perm_cycle_structure(cc[[1L]])
    factors <- list(.int(length(cc)))
    for (k in seq_along(cs))
      factors[[length(factors) + 1L]] <- .pow(.sub(S_One, .pow(t, .int(as.integer(names(cs)[k])))), .int(-cs[[k]]))
    terms[[length(terms) + 1L]] <- .mul(factors)
  }
  .div(.add(terms), as_sym(ord))
}

#' @rdname PermutationGroup
#' @param k Number of random elements per round.
#' @export
normal_closure <- function(G, other, k = 10) {
  if (inherits(other, "PermutationGroup")) {
    degree <- G$degree
    identity <- .ident(degree)
    if (all(vapply(other$gens, .is_identity, TRUE))) return(other)
    Z <- PermutationGroup(other$gens)
    ssi <- schreier_sims_incremental(Z)
    base <- ssi[[1L]]; strong_gens <- ssi[[2L]]
    sgd <- .distribute_gens_by_base(base, strong_gens)
    ot <- .orbits_transversals_from_bsgs(base, sgd)
    basic_orbits <- ot$basic_orbits; basic_transversals <- ot$transversals
    .random_pr_init(G, r = 10, n = 20)
    loop <- TRUE
    while (loop) {
      .random_pr_init(Z, r = 10, n = 10)
      for (i in seq_len(k)) {
        g <- random_pr(G)
        h <- random_pr(Z)
        conj <- h ^ g
        res <- .strip(conj, base, basic_orbits, basic_transversals)
        if (!.is_identity(res$h) || res$j != length(base) + 1L) {
          Z <- PermutationGroup(c(Z$gens, list(conj)))
          strong_gens[[length(strong_gens) + 1L]] <- conj
          ssi <- schreier_sims_incremental(Z, base, strong_gens)
          base <- ssi[[1L]]; strong_gens <- ssi[[2L]]
          sgd <- .distribute_gens_by_base(base, strong_gens)
          ot <- .orbits_transversals_from_bsgs(base, sgd)
          basic_orbits <- ot$basic_orbits; basic_transversals <- ot$transversals
        }
      }
      loop <- FALSE
      for (g in G$gens) {
        for (h in Z$gens) {
          conj <- h ^ g
          res <- .strip(conj, base, basic_orbits, basic_transversals)
          if (!.is_identity(res$h) || res$j != length(base) + 1L) { loop <- TRUE; break }
        }
        if (loop) break
      }
    }
    return(Z)
  }
  if (inherits(other, "Permutation")) return(normal_closure(G, PermutationGroup(list(other)), k))
  if (is.list(other)) return(normal_closure(G, PermutationGroup(other), k))
  NULL
}

#' @rdname PermutationGroup
#' @param alpha A point, or a vector of points for the action on tuples or
#'   sets.
#' @param action `"tuples"`, `"sets"` or `"union"` for a vector `alpha`.
#' @export
orbit <- function(G, alpha, action = "tuples") .orbit(G$degree, G$gens, alpha, action)

#' @rdname PermutationGroup
#' @param beta A point.
#' @param schreier_vector A Schreier vector of `alpha`.
#' @export
orbit_rep <- function(G, alpha, beta, schreier_vector = NULL) {
  if (is.null(schreier_vector)) schreier_vector <- schreier_vector(G, alpha)
  if (is.na(schreier_vector[beta + 1L])) return(FALSE)
  k <- schreier_vector[beta + 1L]
  gens <- lapply(G$gens, `[[`, "array")
  a <- list()
  while (k != -1L) {
    a[[length(a) + 1L]] <- gens[[k + 1L]]
    beta <- match(beta, gens[[k + 1L]]) - 1L
    k <- schreier_vector[beta + 1L]
  }
  if (length(a)) .perm_af(do.call(.af_rmuln, a)) else .ident(G$degree)
}

#' @rdname PermutationGroup
#' @param pairs Return (point, permutation) pairs.
#' @export
orbit_transversal <- function(G, alpha, pairs = FALSE) {
  ot <- .orbit_tr(G$degree, G$gens, alpha)
  tr <- lapply(ot$tr, .perm_af)
  if (!pairs) return(tr)
  Map(function(x, y) list(x, y), ot$pts, tr)
}

#' @rdname PermutationGroup
#' @export
orbits <- function(G) .orbits(G$degree, G$gens)

#' @rdname PermutationGroup
#' @export
group_order <- function(G) {
  if (inherits(G, "WeylGroup")) return(.weyl_group_order(G))
  if (inherits(G, c("FpGroup", "FreeGroup"))) return(fp_order(G))
  if (!is.null(G$order)) return(G$order)
  if (isTRUE(G$is_sym)) {
    G$order <- .nt_out(gmp::factorialZ(G$degree))
    return(G$order)
  }
  if (isTRUE(G$is_alt)) {
    G$order <- .nt_out(gmp::factorialZ(G$degree) %/% 2)
    return(G$order)
  }
  m <- gmp::as.bigz(1)
  for (t in basic_transversals(G)) m <- m * length(t)
  G$order <- .nt_out(m)
  G$order
}

#' @rdname PermutationGroup
#' @export
group_index <- function(G, H) {
  if (is_subgroup(H, G)) return(.nt_out(gmp::as.bigz(group_order(G)) %/% gmp::as.bigz(group_order(H))))
  NULL
}

#' @rdname PermutationGroup
#' @export
is_symmetric_group <- function(G) {
  if (!is.null(G$is_sym)) return(G$is_sym)
  n <- G$degree
  if (n >= 8) {
    if (is_transitive(G)) {
      if (.eval_is_alt_sym_monte_carlo(G)) {
        if (any(vapply(G$gens, perm_is_odd, TRUE))) { G$is_sym <- TRUE; G$is_alt <- FALSE; return(TRUE) }
        G$is_sym <- FALSE; G$is_alt <- TRUE
        return(FALSE)
      }
      return(.eval_is_alt_sym_naive(G, only_sym = TRUE))
    }
    G$is_sym <- FALSE; G$is_alt <- FALSE
    return(FALSE)
  }
  .eval_is_alt_sym_naive(G, only_sym = TRUE)
}

#' @rdname PermutationGroup
#' @export
is_alternating_group <- function(G) {
  if (!is.null(G$is_alt)) return(G$is_alt)
  n <- G$degree
  if (n >= 8) {
    if (is_transitive(G)) {
      if (.eval_is_alt_sym_monte_carlo(G)) {
        if (all(vapply(G$gens, perm_is_even, TRUE))) { G$is_sym <- FALSE; G$is_alt <- TRUE; return(TRUE) }
        G$is_sym <- TRUE; G$is_alt <- FALSE
        return(FALSE)
      }
      return(.eval_is_alt_sym_naive(G, only_alt = TRUE))
    }
    G$is_sym <- FALSE; G$is_alt <- FALSE
    return(FALSE)
  }
  .eval_is_alt_sym_naive(G, only_alt = TRUE)
}

.distinct_primes_lemma <- function(primes) {
  primes <- sort(primes)
  l <- length(primes)
  for (i in seq_len(l)) for (j in seq_len(l)) if (j > i && primes[j] %% primes[i] == 1) return(NULL)
  TRUE
}

#' @rdname PermutationGroup
#' @export
is_cyclic <- function(G) {
  if (!is.null(G$is_cyclic)) return(G$is_cyclic)
  if (length(G$gens) == 1L) { G$is_cyclic <- TRUE; G$is_abelian <- TRUE; return(TRUE) }
  if (identical(G$is_abelian, FALSE)) { G$is_cyclic <- FALSE; return(FALSE) }
  ord <- group_order(G)
  if (ord < 6) {
    G$is_abelian <- TRUE
    if (ord != 4) { G$is_cyclic <- TRUE; return(TRUE) }
  }
  factors <- factorint(ord)
  if (all(factors == 1)) {
    if (isTRUE(G$is_abelian)) { G$is_cyclic <- TRUE; return(TRUE) }
    if (isTRUE(.distinct_primes_lemma(as.numeric(names(factors))))) {
      G$is_cyclic <- TRUE; G$is_abelian <- TRUE
      return(TRUE)
    }
  }
  if (!is_abelian(G)) { G$is_cyclic <- FALSE; return(FALSE) }
  id <- .ident(G$degree)
  ps <- as.numeric(names(factors))[factors > 1]
  G$is_cyclic <- all(vapply(ps, function(p) any(vapply(G$gens, function(g) !(g^(ord %/% p) == id), TRUE)), TRUE))
  G$is_cyclic
}

#' @rdname PermutationGroup
#' @export
is_dihedral <- function(G) {
  if (!is.null(G$is_dihedral)) return(G$is_dihedral)
  ord <- group_order(G)
  if (ord %% 2 == 1) { G$is_dihedral <- FALSE; return(FALSE) }
  if (ord == 2) { G$is_dihedral <- TRUE; return(TRUE) }
  if (ord == 4) { G$is_dihedral <- !is_cyclic(G); return(G$is_dihedral) }
  if (is_abelian(G)) { G$is_dihedral <- FALSE; return(FALSE) }
  n <- ord %/% 2
  gens <- G$gens
  if (length(gens) == 2L) {
    x <- gens[[1L]]; y <- gens[[2L]]
    a <- perm_order(x); b <- perm_order(y)
    if (a < b) { t <- x; x <- y; y <- t; t <- a; a <- b; b <- t }
    if (a == 2 && b == 2) { G$is_dihedral <- TRUE; return(TRUE) }
    if (a == n && b == 2 && y * x * y == perm_inverse(x)) { G$is_dihedral <- TRUE; return(TRUE) }
  }
  order_2 <- list(); order_n <- list()
  for (p in group_elements(G)) {
    k <- perm_order(p)
    if (k == 2) order_2[[length(order_2) + 1L]] <- p
    else if (k == n) order_n[[length(order_n) + 1L]] <- p
  }
  if (length(order_2) != n + 1 - n %% 2) { G$is_dihedral <- FALSE; return(FALSE) }
  if (!length(order_n)) { G$is_dihedral <- FALSE; return(FALSE) }
  x <- order_n[[1L]]
  y <- order_2[[1L]]
  if (n %% 2 == 0 && y == x^(n %/% 2)) y <- order_2[[2L]]
  G$is_dihedral <- y * x * y == perm_inverse(x)
  G$is_dihedral
}

#' @rdname PermutationGroup
#' @param incremental Use the incremental Schreier-Sims algorithm.
#' @export
pointwise_stabilizer <- function(G, points, incremental = TRUE) {
  points <- as.integer(points)
  if (incremental) {
    ssi <- schreier_sims_incremental(G, base = points)
    stab_gens <- Filter(function(gen) identical(gen$array[points + 1L], points), ssi[[2L]])
    if (!length(stab_gens)) stab_gens <- list(.ident(G$degree))
    return(PermutationGroup(stab_gens))
  }
  gens <- G$gens
  for (x in points) gens <- .stabilizer(G$degree, gens, x)
  PermutationGroup(gens)
}

#' @rdname PermutationGroup
#' @param n For `make_perm()`, the number of generators to multiply, or a
#'   vector of generator indices to use in turn.
#' @param seed `NULL` (the global generator), an integer seed, or a vector of
#'   values to use in turn, as SymPy's `_randrange`.
#' @export
make_perm <- function(G, n, seed = NULL) {
  if (length(n) > 1L || is.list(n)) {
    if (!is.null(seed)) .perr("ValueError", "If n is a sequence, seed should be None")
    seed <- structure(as.list(n), sequence = TRUE); n <- length(n)
  }
  randomrange <- .sympy_randrange(seed)
  result <- .ident(G$degree)
  m <- length(G$gens)
  for (i in seq_len(n)) {
    p <- G$gens[[as.integer(randomrange(m)) + 1L]]
    result <- .perm_af(result$array[p$array + 1L])
  }
  result
}

#' @rdname PermutationGroup
#' @export
group_random <- function(G, af = FALSE) {
  rank <- .py_randrange(group_order(G))
  coset_unrank(G, rank, af)
}

#' @rdname PermutationGroup
#' @export
random_stab <- function(G, alpha, schreier_vector = NULL, .random_prec = NULL) {
  if (is.null(schreier_vector)) schreier_vector <- schreier_vector(G, alpha)
  rand <- if (is.null(.random_prec)) random_pr(G) else .random_prec$rand
  beta <- rand$array[alpha + 1L]
  h <- orbit_rep(G, alpha, beta, schreier_vector)
  .perm_af(.af_invert(h$array)[rand$array + 1L])
}

#' @rdname PermutationGroup
#' @export
schreier_sims <- function(G) {
  if (length(G$transversals)) return(invisible(NULL))
  .schreier_sims_base(G)
  invisible(NULL)
}

.schreier_sims_base <- function(G, base = NULL) {
  ssi <- schreier_sims_incremental(G, base = base, slp_dict = TRUE)
  base <- ssi[[1L]]; strong_gens <- ssi[[2L]]
  G$base <- base
  G$strong_gens <- strong_gens
  G$strong_gens_slp <- ssi[[3L]]
  if (!length(base)) {
    G$transversals <- list(); G$basic_orbits <- list()
    return(invisible(NULL))
  }
  sgd <- .distribute_gens_by_base(base, strong_gens)
  ot <- .orbits_transversals_from_bsgs(base, sgd, slp = TRUE)
  sk <- .pkeys(strong_gens)
  slps <- lapply(seq_along(ot$slps), function(i) {
    gk <- .pkeys(sgd[[i]])
    lapply(ot$slps[[i]], function(s) match(gk[s + 1L], sk) - 1L)
  })
  G$transversals <- ot$transversals
  G$basic_orbits <- lapply(ot$basic_orbits, sort)
  G$transversal_slp <- slps
  invisible(NULL)
}

#' @rdname PermutationGroup
#' @param gens Generators to use instead of the group's.
#' @param slp_dict Also return the straight line programs of the strong
#'   generators.
#' @export
schreier_sims_incremental <- function(G, base = NULL, gens = NULL, slp_dict = FALSE) {
  if (is.null(base)) base <- integer()
  base <- as.integer(base)
  if (is.null(gens)) gens <- G$gens
  degree <- G$degree
  id_af <- seq_len(degree) - 1L
  if (length(gens) == 1L && .is_identity(gens[[1L]])) {
    if (slp_dict) return(list(base, gens, stats::setNames(list(list(gens[[1L]])), .pkey(gens[[1L]]))))
    return(list(base, gens))
  }
  .base <- base
  .gens <- Filter(Negate(.is_identity), gens)
  for (gen in .gens) {
    if (all(gen$array[.base + 1L] == .base)) {
      new <- which(gen$array != id_af)[1L] - 1L
      .base <- c(.base, new)
    }
  }
  sgd <- .distribute_gens_by_base(.base, .gens)
  strong_gens_slp <- list()
  orbs <- list(); transversals <- list(); slps <- list()
  base_len <- length(.base)
  set_level <- function(i) {
    ot <- .orbit_tr(degree, sgd[[i + 1L]], .base[i + 1L])
    transversals[[i + 1L]] <<- stats::setNames(ot$tr, ot$pts)
    slps[[i + 1L]] <<- stats::setNames(ot$slp, ot$pts)
    orbs[[i + 1L]] <<- ot$pts
  }
  for (i in seq_len(base_len) - 1L) set_level(i)
  i <- base_len - 1L
  while (i >= 0L) {
    continue_i <- FALSE
    tri <- transversals[[i + 1L]]
    for (bidx in seq_along(tri)) {
      beta <- as.integer(names(tri)[bidx]); u_beta <- tri[[bidx]]
      for (jj in seq_along(sgd[[i + 1L]])) {
        gen <- sgd[[i + 1L]][[jj]]$array
        gb <- gen[beta + 1L]
        u1 <- .tr_get(transversals[[i + 1L]], gb)
        g1 <- gen[u_beta + 1L]
        slp <- c(list(c(i, jj - 1L, 0L)), lapply(.tr_get(slps[[i + 1L]], beta), function(g) c(i, g, 0L)))
        if (!identical(g1, u1)) {
          y <- TRUE
          u1_inv <- .af_invert(u1)
          schreier_gen <- u1_inv[g1 + 1L]
          u1_inv_slp <- rev(.tr_get(slps[[i + 1L]], gb))
          slp <- c(lapply(u1_inv_slp, function(g) c(i, g, 1L)), slp)
          st <- .strip_af(schreier_gen, .base, orbs, transversals, i, slp = slp, slps = slps)
          h <- st$h; j <- st$j; slp <- st$slp
          if (j <= base_len) {
            y <- FALSE
          } else if (!is.null(h)) {
            y <- FALSE
            moved <- which(h != id_af)[1L] - 1L
            .base <- c(.base, moved)
            base_len <- base_len + 1L
            sgd[[length(sgd) + 1L]] <- list()
          }
          if (!y) {
            hp <- .perm_af(h)
            strong_gens_slp[[length(strong_gens_slp) + 1L]] <- list(hp, slp)
            if (i + 1L <= j - 1L) for (l in (i + 1L):(j - 1L)) {
              sgd[[l + 1L]][[length(sgd[[l + 1L]]) + 1L]] <- hp
              set_level(l)
            }
            i <- j - 1L
            continue_i <- TRUE
          }
        }
        if (continue_i) break
      }
      if (continue_i) break
    }
    if (continue_i) next
    i <- i - 1L
  }
  strong_gens <- .gens
  if (slp_dict) {
    sgs <- list()
    for (ks in strong_gens_slp) {
      k <- ks[[1L]]; slp <- ks[[2L]]
      strong_gens[[length(strong_gens) + 1L]] <- k
      words <- lapply(slp, function(s) {
        q <- sgd[[s[1L] + 1L]][[s[2L] + 1L]]
        if (s[3L] == 1L) q^-1 else q
      })
      sgs[[.pkey(k)]] <- words
    }
    for (g in .gens) sgs[[.pkey(g)]] <- list(g)
    return(list(.base, strong_gens, sgs))
  }
  c(list(.base), list(c(strong_gens, lapply(strong_gens_slp, `[[`, 1L))))
}

#' @rdname PermutationGroup
#' @param consec_succ Number of consecutive successful sifts to stop at.
#' @export
schreier_sims_random <- function(G, base = NULL, gens = NULL, consec_succ = 10, .random_prec = NULL) {
  if (is.null(base)) base <- integer()
  base <- as.integer(base)
  if (is.null(gens)) gens <- G$gens
  base_len <- length(base)
  n <- G$degree
  for (gen in gens) {
    if (all(gen$array[base + 1L] == base)) {
      new <- 0L
      while (gen$array[new + 1L] == new) new <- new + 1L
      base <- c(base, new)
      base_len <- base_len + 1L
    }
  }
  sgd <- .distribute_gens_by_base(base, gens)
  transversals <- list(); orbs <- list()
  for (i in seq_len(base_len)) {
    ot <- .orbit_tr(n, sgd[[i]], base[i])
    transversals[[i]] <- stats::setNames(lapply(ot$tr, .perm_af), ot$pts)
    orbs[[i]] <- ot$pts
  }
  rp <- .random_prec$g
  c <- 0L
  while (c < consec_succ) {
    if (is.null(.random_prec)) g <- random_pr(G)
    else { g <- rp[[length(rp)]]; rp[[length(rp)]] <- NULL }
    st <- .strip(g, base, orbs, transversals)
    h <- st$h; j <- st$j
    y <- TRUE
    if (j <= base_len) {
      y <- FALSE
    } else if (!.is_identity(h)) {
      y <- FALSE
      moved <- 0L
      while (h$array[moved + 1L] == moved) moved <- moved + 1L
      base <- c(base, moved)
      base_len <- base_len + 1L
      sgd[[length(sgd) + 1L]] <- list()
    }
    if (!y) {
      if (j >= 2L) for (l0 in 1:(j - 1L)) {
        sgd[[l0 + 1L]][[length(sgd[[l0 + 1L]]) + 1L]] <- h
        ot <- .orbit_tr(n, sgd[[l0 + 1L]], base[l0 + 1L])
        transversals[[l0 + 1L]] <- stats::setNames(lapply(ot$tr, .perm_af), ot$pts)
        orbs[[l0 + 1L]] <- ot$pts
      }
      c <- 0L
    } else c <- c + 1L
  }
  strong_gens <- sgd[[1L]]
  for (gen in sgd[[2L]]) if (!.perm_in(gen, strong_gens)) strong_gens[[length(strong_gens) + 1L]] <- gen
  list(base, strong_gens)
}

#' @rdname PermutationGroup
#' @export
schreier_vector <- function(G, alpha) {
  n <- G$degree
  v <- rep(NA_integer_, n)
  v[alpha + 1L] <- -1L
  orb <- as.integer(alpha)
  used <- logical(n); used[alpha + 1L] <- TRUE
  gens <- lapply(G$gens, `[[`, "array")
  k <- 1L
  while (k <= length(orb)) {
    b <- orb[k]
    for (i in seq_along(gens)) {
      temp <- gens[[i]][b + 1L]
      if (!used[temp + 1L]) { orb <- c(orb, temp); used[temp + 1L] <- TRUE; v[temp + 1L] <- i - 1L }
    }
    k <- k + 1L
  }
  v
}

#' @rdname PermutationGroup
#' @export
stabilizer <- function(G, alpha) PermutationGroup(.stabilizer(G$degree, G$gens, alpha))

#' @rdname PermutationGroup
#' @export
strong_gens <- function(G) {
  if (!length(G$strong_gens)) schreier_sims(G)
  G$strong_gens
}

#' @rdname PermutationGroup
#' @export
subgroup <- function(G, gens) {
  if (!all(vapply(gens, function(g) group_contains(G, g), TRUE)))
    .perr("ValueError", "The group does not contain the supplied generators")
  PermutationGroup(gens)
}

#' @rdname PermutationGroup
#' @param prop A predicate on permutations defining the subgroup.
#' @param tests Optional list of level tests.
#' @param init_subgroup A known subgroup of the result.
#' @export
subgroup_search <- function(G, prop, base = NULL, strong_gens = NULL, tests = NULL, init_subgroup = NULL) {
  if (is.null(base)) {
    ssi <- schreier_sims_incremental(G)
    base <- ssi[[1L]]; strong_gens <- ssi[[2L]]
  }
  base <- as.integer(base)
  base_len <- length(base)
  degree <- G$degree
  identity <- .ident(degree)
  base_ordering <- c(.base_ordering(base, degree), degree, -1L)
  bo <- function(pt) base_ordering[pt + 1L]
  get_reps <- function(orbits) vapply(orbits, function(orbit) orbit[which.min(bo(orbit))], 0L)
  sgd <- .distribute_gens_by_base(base, strong_gens)
  ot <- .orbits_transversals_from_bsgs(base, sgd)
  basic_orbits <- ot$basic_orbits; transversals <- ot$transversals
  if (is.null(init_subgroup)) init_subgroup <- PermutationGroup(list(identity))
  if (is.null(tests)) tests <- rep(list(function(x) TRUE), base_len)
  res <- init_subgroup
  f <- base_len - 1L
  l <- base_len - 1L
  res_base <- base
  ssi <- schreier_sims_incremental(res, base = res_base)
  res_base <- ssi[[1L]]; res_strong_gens <- ssi[[2L]]
  res_sgd <- .distribute_gens_by_base(res_base, res_strong_gens)
  res_generators <- res$gens
  res_init <- lapply(seq_len(base_len), function(i) .orbit_bfs(degree, res_sgd[[i]], res_base[i]))
  orbit_reps <- vector("list", base_len)
  orbits <- .orbits(degree, res_sgd[[f + 1L]])
  orbit_reps[[f + 1L]] <- get_reps(orbits)
  k <- match(base[f + 1L], orbit_reps[[f + 1L]])
  if (is.na(k)) .perr("ValueError", "list.remove(x): x not in list")
  orbit_reps[[f + 1L]] <- orbit_reps[[f + 1L]][-k]
  c <- integer(base_len)
  u <- rep(list(identity), base_len)
  sorted_orbits <- lapply(seq_len(base_len), function(i) {
    o <- basic_orbits[[i]]
    o[order(bo(o), method = "radix")]
  })
  mu <- rep(NA_integer_, base_len)
  nu <- rep(NA_integer_, base_len)
  update_nu <- function(l) {
    temp_index <- length(basic_orbits[[l + 1L]]) + 1L - length(res_init[[l + 1L]])
    if (temp_index >= length(sorted_orbits[[l + 1L]])) nu[l + 1L] <<- base_ordering[degree + 1L]
    else nu[l + 1L] <<- sorted_orbits[[l + 1L]][temp_index + 1L]
  }
  mu[l + 1L] <- degree + 1L
  update_nu(l)
  computed_words <- rep(list(identity), base_len)
  img <- function(w, pt) w$array[pt + 1L]
  repeat {
    while (l < base_len - 1L) {
      cw <- img(computed_words[[l + 1L]], base[l + 1L])
      if (!(cw %in% orbit_reps[[l + 1L]] && bo(mu[l + 1L]) < bo(cw) && bo(cw) < bo(nu[l + 1L]) &&
            tests[[l + 1L]](computed_words))) break
      new_point <- cw
      res_base[l + 1L] <- new_point
      new_stab_gens <- .stabilizer(degree, res_sgd[[l + 1L]], new_point)
      res_sgd[[l + 2L]] <- new_stab_gens
      orbits <- .orbits(degree, new_stab_gens)
      orbit_reps[[l + 2L]] <- get_reps(orbits)
      l <- l + 1L
      temp_orbit <- vapply(basic_orbits[[l + 1L]], function(point) img(computed_words[[l]], point), 0L)
      sorted_orbits[[l + 1L]] <- temp_orbit[order(bo(temp_orbit), method = "radix")]
      new_mu <- degree + 1L
      for (i in seq_len(l) - 1L) {
        if (base[l + 1L] %in% res_init[[i + 1L]]) {
          candidate <- img(computed_words[[i + 1L]], base[i + 1L])
          if (bo(candidate) > bo(new_mu)) new_mu <- candidate
        }
      }
      mu[l + 1L] <- new_mu
      update_nu(l)
      c[l + 1L] <- 0L
      temp_point <- sorted_orbits[[l + 1L]][c[l + 1L] + 1L]
      gamma <- match(temp_point, computed_words[[l]]$array) - 1L
      u[[l + 1L]] <- .tr_get(transversals[[l + 1L]], gamma)
      computed_words[[l + 1L]] <- .perm_af(computed_words[[l]]$array[u[[l + 1L]]$array + 1L])
    }
    g <- computed_words[[l + 1L]]
    temp_point <- img(g, base[l + 1L])
    if (l == base_len - 1L && bo(mu[l + 1L]) < bo(temp_point) && bo(temp_point) < bo(nu[l + 1L]) &&
        temp_point %in% orbit_reps[[l + 1L]] && tests[[l + 1L]](computed_words) && isTRUE(prop(g))) {
      res_generators[[length(res_generators) + 1L]] <- g
      res_base <- base
      res_strong_gens[[length(res_strong_gens) + 1L]] <- g
      res_sgd <- .distribute_gens_by_base(res_base, res_strong_gens)
      res_init <- lapply(seq_len(base_len), function(i) .orbit_bfs(degree, res_sgd[[i]], res_base[i]))
      orbit_reps[[f + 1L]] <- get_reps(orbits)
      l <- f
    }
    while (l >= 0L && c[l + 1L] == length(basic_orbits[[l + 1L]]) - 1L) l <- l - 1L
    if (l == -1L) return(PermutationGroup(res_generators))
    if (l < f) {
      f <- l
      c[l + 1L] <- 0L
      temp_orbits <- .orbits(degree, res_sgd[[f + 1L]])
      orbit_reps[[f + 1L]] <- get_reps(temp_orbits)
      mu[l + 1L] <- degree + 1L
      temp_index <- length(basic_orbits[[l + 1L]]) + 1L - length(res_init[[l + 1L]])
      if (temp_index >= length(sorted_orbits[[l + 1L]])) nu[l + 1L] <- base_ordering[degree + 1L]
      else nu[l + 1L] <- sorted_orbits[[l + 1L]][temp_index + 1L]
    }
    c[l + 1L] <- c[l + 1L] + 1L
    if (l == 0L) gamma <- sorted_orbits[[l + 1L]][c[l + 1L] + 1L]
    else gamma <- match(sorted_orbits[[l + 1L]][c[l + 1L] + 1L], computed_words[[l]]$array) - 1L
    u[[l + 1L]] <- .tr_get(transversals[[l + 1L]], gamma)
    if (l == 0L) computed_words[[l + 1L]] <- u[[l + 1L]]
    else computed_words[[l + 1L]] <- .perm_af(computed_words[[l]]$array[u[[l + 1L]]$array + 1L])
  }
}

#' @rdname PermutationGroup
#' @export
transitivity_degree <- function(G) {
  if (is.null(G$transitivity_degree)) {
    n <- G$degree
    Gc <- G
    for (i in seq_len(n) - 1L) {
      orb <- orbit(Gc, i)
      if (length(orb) != n - i) { G$transitivity_degree <- i; return(i) }
      Gc <- stabilizer(Gc, i)
    }
    G$transitivity_degree <- n
    return(n)
  }
  G$transitivity_degree
}

.p_elements_group <- function(G, p) {
  gens <- G$gens
  gens <- gens[order(-vapply(gens, function(x) as.numeric(perm_order(x)), 0), method = "radix")]
  gens_p <- lapply(gens, function(g) g^(as.numeric(perm_order(g)) %/% p))
  gens_r <- list()
  for (i in seq_along(gens)) {
    x <- gens[[i]]
    x_order <- as.numeric(perm_order(x))
    x_p <- x^(x_order %/% p)
    P <- if (i > 1L) PermutationGroup(gens_p[seq_len(i - 1L)]) else PermutationGroup(group_identity(G))
    if (!group_contains(P, x^(x_order %/% p))) {
      gens_r[[length(gens_r) + 1L]] <- x^(x_order %/% p)
    } else {
      for (s in generator_product(P, x_p, original = TRUE)) x <- x * s^-1
      x_order <- x_order / p
      gens <- gens[-i]; gens_p <- gens_p[-i]
      j0 <- i - 2L
      while (j0 < length(gens) && as.numeric(perm_order(gens[[j0 + 1L]])) >= x_order) j0 <- j0 + 1L
      gens <- append(gens, list(x), after = j0)
      gens_p <- append(gens_p, list(x), after = j0)
    }
  }
  PermutationGroup(gens_r)
}

.sylow_alt_sym <- function(G, p) {
  n <- G$degree
  gens <- list()
  identity <- Permutation(n - 1L)
  alt <- p == 2 && all(vapply(G$gens, perm_is_even, TRUE))
  coeffs <- numeric()
  m <- n
  while (m > 0) { coeffs <- c(coeffs, m %% p); m <- m %/% p }
  power <- length(coeffs) - 1L
  for (i in seq_len(power)) {
    if (i == 1L && alt) next
    pi1 <- p^(i - 1L); pii <- p^i
    gen <- Permutation(as.integer((seq_len(pii) - 1L + pi1) %% pii))
    gens[[length(gens) + 1L]] <- identity * gen
    if (alt) {
      gen <- Permutation(0, 1) * gen * Permutation(0, 1) * gen
      gens[[length(gens) + 1L]] <- gen
    }
  }
  start <- 0
  while (power > 0L) {
    a <- coeffs[power + 1L]
    for (rep_ in seq_len(a)) {
      shift <- Permutation()
      if (start > 0) {
        for (i in seq_len(p^power) - 1L) shift <- perm_call(shift, i, start + i)
        if (alt) {
          gen <- Permutation(0, 1) * shift * Permutation(0, 1) * shift
          gens[[length(gens) + 1L]] <- gen
          j <- 2L * (power - 1L)
        } else j <- power
        for (i in seq_len(min(j, length(gens)))) {
          if (alt && (i - 1L) %% 2L == 1L) next
          gens[[length(gens) + 1L]] <- shift * gens[[i]] * shift
        }
      }
      start <- start + p^power
    }
    power <- power - 1L
  }
  gens
}

#' @rdname PermutationGroup
#' @export
sylow_subgroup <- function(G, p) {
  if (!isTRUE(isprime(p))) .perr("ValueError", "p must be a prime")
  is_p_group <- function(G) {
    m <- gmp::as.bigz(group_order(G)); n <- 0L
    while (m %% p == 0) {
      m <- m %/% p; n <- n + 1L
      if (m == 1) return(list(TRUE, n))
    }
    list(FALSE, n)
  }
  sylow_reduce <- function(mu, nu) {
    Q <- sylow_subgroup(.hom_image(mu), p)
    Q <- .hom_invert_subgroup(mu, Q)
    nu <- .hom_restrict_to(nu, Q)
    R <- sylow_subgroup(.hom_image(nu), p)
    .hom_invert_subgroup(nu, R)
  }
  ord <- gmp::as.bigz(group_order(G))
  if (ord %% p != 0) return(PermutationGroup(list(group_identity(G))))
  pg <- is_p_group(G)
  n <- pg[[2L]]
  if (pg[[1L]]) return(G)
  if (is_alt_sym(G)) return(PermutationGroup(.sylow_alt_sym(G, p)))
  orbs <- orbits(G)
  non_p_orbits <- Filter(function(o) length(o) %% p != 0 && length(o) != 1L, orbs)
  if (length(non_p_orbits)) {
    o <- non_p_orbits[[1L]]
    H <- stabilizer(G, o[length(o)])
    return(sylow_subgroup(H, p))
  }
  if (!is_transitive(G)) {
    orbs <- orbs[order(lengths(orbs), method = "radix")]
    omega1 <- orbs[[length(orbs)]]
    orbs <- orbs[-length(orbs)]
    sets <- lapply(orbs, function(o) .pyset_from(.orbit_bfs(G$degree, G$gens, min(o))))
    omega2 <- .pyset_items(do.call(.pyset_union, sets))
    mu <- orbit_homomorphism(G, omega1)
    nu <- orbit_homomorphism(G, omega2)
    return(sylow_reduce(mu, nu))
  }
  blocks <- minimal_blocks(G)
  if (length(blocks) > 1L) {
    mu <- block_homomorphism(G, blocks[[1L]])
    nu <- block_homomorphism(G, blocks[[2L]])
    return(sylow_reduce(mu, nu))
  } else if (length(blocks) == 1L) {
    block <- blocks[[1L]]
    if (any(block != 0L)) {
      mu <- block_homomorphism(G, block)
      if (!is_p_group(.hom_image(mu))[[1L]]) {
        S <- sylow_subgroup(.hom_image(mu), p)
        return(sylow_subgroup(.hom_invert_subgroup(mu, S), p))
      }
    }
  }
  g <- group_random(G)
  g_order <- as.numeric(perm_order(g))
  while (g_order %% p != 0 || g_order == 0) {
    g <- group_random(G)
    g_order <- as.numeric(perm_order(g))
  }
  g <- g^(g_order %/% p)
  if (ord %% p^2 != 0) return(PermutationGroup(g))
  C <- centralizer(G, g)
  while (gmp::as.bigz(group_order(C)) %% gmp::as.bigz(p)^n != 0) {
    S <- sylow_subgroup(C, p)
    s_order <- gmp::as.bigz(group_order(S))
    Z <- group_center(S)
    P <- .p_elements_group(Z, p)
    h <- group_random(P)
    C_h <- centralizer(G, h)
    while (gmp::as.bigz(group_order(C_h)) %% p * s_order != 0) {
      h <- group_random(P)
      C_h <- centralizer(G, h)
    }
    C <- C_h
  }
  sylow_subgroup(C, p)
}

#' @rdname PermutationGroup
#' @export
quotient_group <- function(G, N) {
  if (!is_normal(N, G)) .perr("ValueError", "N must be a normal subgroup of G")
  t <- coset_table(G, N)
  gen <- lapply(seq_along(t[[1L]]), function(j) Permutation(vapply(t, function(row) row[j], 0L)))
  PermutationGroup(gen)
}

#' @rdname PermutationGroup
#' @export
group_equals <- function(G, H) {
  if (!inherits(H, "PermutationGroup")) return(FALSE)
  if (setequal(.pkeys(G$gens), .pkeys(H$gens))) return(TRUE)
  for (g in G$gens) if (!group_contains(H, g)) return(FALSE)
  for (g in H$gens) if (!group_contains(G, g)) return(FALSE)
  TRUE
}

# ---------------------------------------------------------------------------
# Named groups (named_groups.py), direct products (group_constructs.py),
# Rubik generators (generators.py), cosets, and the permutation group part
# of homomorphisms.py that sylow_subgroup() needs.

#' Named permutation groups
#'
#' SymPy's named groups, with the same generators and the properties SymPy
#' records on them: `SymmetricGroup(n)`, `AlternatingGroup(n)`,
#' `CyclicGroup(n)`, `DihedralGroup(n)`, `AbelianGroup(...)` (the direct
#' product of cyclic groups of the given orders) and `RubikGroup(n)`.
#' `DirectProduct(...)` is the direct product of permutation groups;
#' `rubik(n)` and `rubik_cube_generators()` give the generators of the
#' Rubik's cube groups.
#'
#' @param n Degree, or the cube dimension.
#' @param ... Cyclic orders for `AbelianGroup()`, groups for
#'   `DirectProduct()`.
#' @return A `PermutationGroup`, or a list of permutations.
#' @examples
#' group_order(AlternatingGroup(5))
#' group_generators(SymmetricGroup(4))
#' group_order(AbelianGroup(2, 3, 4))
#' @export
SymmetricGroup <- function(n) {
  n <- as.integer(n)
  if (n == 1L) G <- PermutationGroup(list(Permutation(0L)))
  else if (n == 2L) G <- PermutationGroup(list(Permutation(c(1L, 0L))))
  else {
    a <- seq_len(n) - 1L
    a[1:2] <- c(1L, 0L)
    G <- PermutationGroup(list(.perm_af(c(seq_len(n - 1L), 0L)), .perm_af(a)))
  }
  .set_symmetric_group_properties(G, n, n)
  G$is_sym <- TRUE
  G
}

.set_symmetric_group_properties <- function(G, n, degree) {
  G$is_abelian <- n < 3
  G$is_nilpotent <- n < 3
  G$is_solvable <- n < 5
  G$degree <- degree
  G$is_transitive <- TRUE
  G$is_dihedral <- n %in% c(2, 3)
}

#' @rdname SymmetricGroup
#' @export
AlternatingGroup <- function(n) {
  n <- as.integer(n)
  if (n == 1L) G <- PermutationGroup(list(Permutation(0L)))
  else if (n == 2L) G <- PermutationGroup(list(Permutation(c(0L, 1L))))
  else {
    gen1 <- seq_len(n) - 1L
    gen1[1:3] <- c(1L, 2L, 0L)
    gen2 <- if (n %% 2L) c(seq_len(n - 1L), 0L) else c(0L, seq_len(n - 2L) + 1L, 1L)
    gens <- if (identical(gen1, gen2)) list(gen1) else list(gen1, gen2)
    G <- PermutationGroup(lapply(gens, .perm_af), dups = FALSE)
  }
  .set_alternating_group_properties(G, n, n)
  G$is_alt <- TRUE
  G
}

.set_alternating_group_properties <- function(G, n, degree) {
  if (n < 3) G$order <- 1
  G$is_abelian <- n < 4
  G$is_nilpotent <- n < 4
  G$is_solvable <- n < 5
  G$degree <- degree
  G$is_transitive <- n != 2
  G$is_dihedral <- FALSE
}

#' @rdname SymmetricGroup
#' @export
CyclicGroup <- function(n) {
  n <- as.integer(n)
  G <- PermutationGroup(list(.perm_af(c(seq_len(n - 1L), 0L))))
  G$is_abelian <- TRUE; G$is_nilpotent <- TRUE; G$is_solvable <- TRUE
  G$degree <- n; G$is_transitive <- TRUE; G$order <- n
  G$is_dihedral <- n == 2L
  G
}

#' @rdname SymmetricGroup
#' @export
DihedralGroup <- function(n) {
  n <- as.integer(n)
  if (n == 1L) return(PermutationGroup(list(Permutation(c(1L, 0L)))))
  if (n == 2L) return(PermutationGroup(list(Permutation(c(1L, 0L, 3L, 2L)), Permutation(c(2L, 3L, 0L, 1L)),
                                            Permutation(c(3L, 2L, 1L, 0L)))))
  G <- PermutationGroup(list(.perm_af(c(seq_len(n - 1L), 0L)), .perm_af(rev(seq_len(n) - 1L))))
  G$is_nilpotent <- bitwAnd(n, n - 1L) == 0L
  G$is_dihedral <- TRUE; G$is_abelian <- FALSE; G$is_solvable <- TRUE
  G$degree <- n; G$is_transitive <- TRUE; G$order <- 2 * n
  G
}

#' @rdname SymmetricGroup
#' @export
AbelianGroup <- function(...) {
  orders <- as.integer(unlist(list(...)))
  groups <- lapply(orders, CyclicGroup)
  G <- do.call(DirectProduct, groups)
  G$is_abelian <- TRUE
  G$degree <- sum(orders)
  G$order <- prod(orders)
  G
}

#' @rdname SymmetricGroup
#' @export
DirectProduct <- function(...) {
  groups <- list(...)
  if (length(groups) == 1L && is.list(groups[[1L]]) && !inherits(groups[[1L]], "PermutationGroup"))
    groups <- groups[[1L]]
  degrees <- vapply(groups, function(g) g$degree, 0L)
  gens_count <- vapply(groups, function(g) length(g$gens), 0L)
  total_degree <- sum(degrees)
  array_gens <- rep(list(seq_len(total_degree) - 1L), sum(gens_count))
  current_gen <- 0L; current_deg <- 0L
  for (i in seq_along(groups)) {
    for (j in seq_len(gens_count[i])) {
      gen <- groups[[i]]$gens[[j]]$array
      array_gens[[current_gen + j]][current_deg + seq_len(degrees[i])] <- gen + current_deg
    }
    current_gen <- current_gen + gens_count[i]
    current_deg <- current_deg + degrees[i]
  }
  keys <- vapply(array_gens, paste, "", collapse = ",")
  PermutationGroup(lapply(array_gens[!duplicated(keys)], .perm_af), dups = FALSE)
}

#' @rdname SymmetricGroup
#' @export
RubikGroup <- function(n) {
  if (n <= 1) .perr("ValueError", "Invalid cube. n has to be greater than 1")
  PermutationGroup(rubik(n))
}

#' @rdname SymmetricGroup
#' @export
rubik_cube_generators <- function() {
  a <- list(
    list(c(1, 3, 8, 6), c(2, 5, 7, 4), c(9, 33, 25, 17), c(10, 34, 26, 18), c(11, 35, 27, 19)),
    list(c(9, 11, 16, 14), c(10, 13, 15, 12), c(1, 17, 41, 40), c(4, 20, 44, 37), c(6, 22, 46, 35)),
    list(c(17, 19, 24, 22), c(18, 21, 23, 20), c(6, 25, 43, 16), c(7, 28, 42, 13), c(8, 30, 41, 11)),
    list(c(25, 27, 32, 30), c(26, 29, 31, 28), c(3, 38, 43, 19), c(5, 36, 45, 21), c(8, 33, 48, 24)),
    list(c(33, 35, 40, 38), c(34, 37, 39, 36), c(3, 9, 46, 32), c(2, 12, 47, 29), c(1, 14, 48, 27)),
    list(c(41, 43, 48, 46), c(42, 45, 47, 44), c(14, 22, 30, 38), c(15, 23, 31, 39), c(16, 24, 32, 40)))
  lapply(a, function(x) Permutation(lapply(x, function(xi) xi - 1L), size = 48))
}

#' @rdname SymmetricGroup
#' @export
rubik <- function(n) {
  n <- as.integer(n)
  if (n < 2L) .perr("ValueError", "dimension of cube must be > 1")
  nm <- c("U", "F", "R", "B", "L", "D")
  faces <- list()
  count <- 0L
  for (f in nm) {
    faces[[f]] <- matrix(count + seq_len(n * n) - 1L, n, n, byrow = TRUE)
    count <- count + n * n
  }
  getr <- function(f, i) faces[[f]][, n - i + 1L]
  getl <- function(f, i) faces[[f]][, i]
  getu <- function(f, i) faces[[f]][i, ]
  getd <- function(f, i) faces[[f]][n - i + 1L, ]
  setr <- function(f, i, s) faces[[f]][, n - i + 1L] <<- s
  setl <- function(f, i, s) faces[[f]][, i] <<- s
  setu <- function(f, i, s) faces[[f]][i, ] <<- s
  setd <- function(f, i, s) faces[[f]][n - i + 1L, ] <<- s
  cw <- function(Fc, r = 1L) for (k in seq_len(r)) faces[[Fc]] <<- t(faces[[Fc]][n:1, , drop = FALSE])
  ccw <- function(Fc) cw(Fc, 3L)
  fcw <- function(i, r = 1L) {
    for (k in seq_len(r)) {
      if (i == 0L) cw("F")
      i <- i + 1L
      temp <- getr("L", i)
      setr("L", i, getu("D", i))
      setu("D", i, rev(getl("R", i)))
      setl("R", i, getd("U", i))
      setd("U", i, rev(temp))
      i <- i - 1L
    }
  }
  fccw <- function(i) fcw(i, 3L)
  FCW <- function(r = 1L) {
    for (k in seq_len(r)) {
      cw("F"); ccw("B"); cw("U")
      t <- faces[["U"]]
      cw("L"); faces[["U"]] <<- faces[["L"]]
      cw("D"); faces[["L"]] <<- faces[["D"]]
      cw("R"); faces[["D"]] <<- faces[["R"]]
      faces[["R"]] <<- t
    }
  }
  FCCW <- function() FCW(3L)
  UCW <- function(r = 1L) {
    for (k in seq_len(r)) {
      cw("U"); ccw("D")
      t <- faces[["F"]]
      faces[["F"]] <<- faces[["R"]]
      faces[["R"]] <<- faces[["B"]]
      faces[["B"]] <<- faces[["L"]]
      faces[["L"]] <<- t
    }
  }
  UCCW <- function() UCW(3L)
  g <- list()
  perm <- function() g[[length(g) + 1L]] <<- Permutation(unlist(lapply(nm, function(f) as.vector(t(faces[[f]])))))
  for (i in seq_len(n - 1L) - 1L) { fcw(i); perm(); fccw(i) }
  UCW()
  for (i in seq_len(n - 1L) - 1L) { fcw(i); UCCW(); perm(); UCW(); fccw(i) }
  UCCW()
  FCW(); UCCW(); FCCW()
  for (i in seq_len(n - 1L) - 1L) {
    fcw(i); FCW(); UCW(); FCCW(); perm(); FCW(); UCCW(); FCCW(); fccw(i)
  }
  FCW(); UCW(); FCCW()
  g
}

#' Cosets and the symmetric group of a degree
#'
#' `Coset(g, H, G, dir)` is SymPy's `Coset`: the right coset `H*g` for
#' `dir = "+"` and the left coset `g*H` for `dir = "-"`; `G` defaults to the
#' full symmetric group, `SymmetricPermutationGroup(n)`. `G * g` and `g * G`
#' build cosets too. `is_left_coset()`, `is_right_coset()` and
#' `coset_as_list()` are its methods.
#'
#' @param g A permutation.
#' @param H,G Permutation groups.
#' @param dir `"+"` or `"-"`.
#' @param n The degree.
#' @param x A `Coset`.
#' @return A `Coset`, a `SymmetricPermutationGroup`, a logical or a list.
#' @examples
#' G <- PermutationGroup(Permutation(1, 2), Permutation(0, 1))
#' cst <- G * Permutation(1, 2)
#' is_right_coset(cst)
#' length(coset_as_list(cst))
#' @export
Coset <- function(g, H, G = NULL, dir = "+") {
  if (!inherits(g, "Permutation")) .perr("NotImplementedError", "")
  if (!inherits(H, "PermutationGroup")) .perr("NotImplementedError", "")
  if (!is.null(G)) {
    if (!inherits(G, c("PermutationGroup", "SymmetricPermutationGroup"))) .perr("NotImplementedError", "")
    if (!is_subgroup(H, G)) .perr("ValueError", sprintf("%s must be a subgroup of %s.", .group_str(H), format(G)))
    inG <- if (inherits(G, "SymmetricPermutationGroup")) length(g$array) == G$degree else group_contains(G, g)
    if (!inG) .perr("ValueError", sprintf("%s must be an element of %s.", .perm_str(g), format(G)))
  } else {
    if (length(g$array) != H$degree)
      .perr("ValueError", sprintf("The size of the permutation %s and the degree of the permutation group %s should be matching ",
                                  .perm_str(g), .group_str(H)))
    G <- SymmetricPermutationGroup(length(g$array))
  }
  if (!(dir %in% c("+", "-"))) .perr("ValueError", sprintf("dir must be one of '+' or '-' not %s", dir))
  structure(list(g = g, H = H, G = G, dir = dir), class = "Coset")
}

#' @rdname Coset
#' @export
is_left_coset <- function(x) x$dir == "-"

#' @rdname Coset
#' @export
is_right_coset <- function(x) x$dir == "+"

#' @rdname Coset
#' @export
coset_as_list <- function(x) {
  els <- group_elements(x$H)
  if (x$dir == "+") lapply(els, function(h) h * x$g) else lapply(els, function(h) x$g * h)
}

#' @export
format.Coset <- function(x, ...) paste0("Coset(", .perm_str(x$g), ", ", .group_str(x$H), ", ", format(x$G), ", ", x$dir, ")")

#' @export
print.Coset <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @rdname Coset
#' @export
SymmetricPermutationGroup <- function(n) structure(list(degree = as.integer(n)), class = "SymmetricPermutationGroup")

#' @export
format.SymmetricPermutationGroup <- function(x, ...) paste0("SymmetricPermutationGroup(", x$degree, ")")

#' @export
print.SymmetricPermutationGroup <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Homomorphisms between permutation groups (homomorphisms.py)

.hom_new <- function(domain, codomain, from, to) {
  H <- new.env(parent = emptyenv())
  H$domain <- domain; H$codomain <- codomain
  H$from <- from; H$to <- to; H$keys <- .pkeys(from)
  H$inverses <- NULL; H$kernel <- NULL; H$image <- NULL
  structure(H, class = "GroupHomomorphism")
}

#' @export
format.GroupHomomorphism <- function(x, ...) "GroupHomomorphism"

.hom_image <- function(H) {
  if (is.null(H$image)) {
    vals <- H$to[!duplicated(.pkeys(H$to))]
    H$image <- subgroup(H$codomain, vals)
  }
  H$image
}

.hom_invs <- function(H) {
  image <- .hom_image(H)
  inv_from <- list(); inv_to <- list()
  for (k in seq_along(H$from)) {
    v <- H$to[[k]]
    if (!(.pkey(v) %in% .pkeys(inv_from) || .is_identity(v))) {
      inv_from[[length(inv_from) + 1L]] <- v
      inv_to[[length(inv_to) + 1L]] <- H$from[[k]]
    }
  }
  lookup <- function(s) {
    i <- match(.pkey(s), .pkeys(inv_from))
    if (!is.na(i)) return(inv_to[[i]])
    inv_to[[match(.pkey(s^-1), .pkeys(inv_from))]]^-1
  }
  for (g in strong_gens(image)) {
    if (.pkey(g) %in% .pkeys(inv_from) || .is_identity(g)) next
    parts <- rev(image$strong_gens_slp[[.pkey(g)]])
    w <- perm_prod(lapply(parts, lookup))
    inv_from[[length(inv_from) + 1L]] <- g
    inv_to[[length(inv_to) + 1L]] <- w
  }
  list(from = inv_from, to = inv_to)
}

.hom_invert <- function(H, g) {
  if (inherits(g, "Permutation")) {
    if (is.null(H$inverses)) H$inverses <- .hom_invs(H)
    image <- .hom_image(H)
    w <- group_identity(H$domain)
    if (!group_contains(image, g)) .perr("ValueError", "Given element is not in the image of the homomorphism.")
    gens <- rev(generator_product(image, g))
    ik <- .pkeys(H$inverses$from)
    lookup <- function(s) {
      i <- match(.pkey(s), ik)
      if (!is.na(i)) return(H$inverses$to[[i]])
      H$inverses$to[[match(.pkey(s^-1), ik)]]^-1
    }
    return(perm_prod(c(list(w), lapply(Filter(Negate(.is_identity), gens), lookup))))
  }
  if (is.list(g)) return(lapply(g, function(e) .hom_invert(H, e)))
  NULL
}

.hom_apply <- function(H, elem) {
  if (!group_contains(H$domain, elem)) {
    if (is.list(elem) && !inherits(elem, "Permutation")) return(lapply(elem, function(e) .hom_apply(H, e)))
    .perr("ValueError", "The supplied element does not belong to the domain")
  }
  if (.is_identity(elem)) return(group_identity(H$codomain))
  value <- group_identity(H$codomain)
  for (g in generator_product(H$domain, elem, original = TRUE)) {
    i <- match(.pkey(g), H$keys)
    if (!is.na(i)) value <- H$to[[i]] * value
    else value <- H$to[[match(.pkey(g^-1), H$keys)]]^-1 * value
  }
  value
}

.hom_kernel <- function(H) {
  if (is.null(H$kernel)) H$kernel <- subgroup_search(H$domain, function(g) .is_identity(.hom_apply(H, g)))
  H$kernel
}

.hom_restrict_to <- function(H, K) {
  if (!inherits(K, "PermutationGroup") || !is_subgroup(K, H$domain))
    .perr("ValueError", "Given H is not a subgroup of the domain")
  .hom_new(K, H$codomain, K$gens, lapply(K$gens, function(g) .hom_apply(H, g)))
}

.hom_invert_subgroup <- function(H, K) {
  if (!is_subgroup(K, .hom_image(H))) .perr("ValueError", "Given H is not a subgroup of the image")
  gens <- list()
  P <- PermutationGroup(group_identity(.hom_image(H)))
  for (h in K$gens) {
    h_i <- .hom_invert(H, h)
    if (!group_contains(P, h_i)) {
      gens[[length(gens) + 1L]] <- h_i
      P <- PermutationGroup(gens)
    }
    for (k in .hom_kernel(H)$gens) {
      if (!group_contains(P, k * h_i)) {
        gens[[length(gens) + 1L]] <- k * h_i
        P <- PermutationGroup(gens)
      }
    }
  }
  P
}

#' Homomorphisms of permutation groups
#'
#' `orbit_homomorphism(G, omega)` is the action of `G` on its invariant set
#' `omega` and `block_homomorphism(G, blocks)` the action on a block system
#' (as returned by [minimal_block()]), both as in SymPy. `hom_image()`,
#' `hom_kernel()`, `hom_apply()` and `hom_invert()` are the homomorphism's
#' `image()`, `kernel()`, call and `invert()`.
#'
#' @param G A permutation group.
#' @param omega A vector of points.
#' @param blocks A block system.
#' @param hom A homomorphism.
#' @param g A permutation.
#' @return A `GroupHomomorphism`, a group or a permutation.
#' @examples
#' D <- DihedralGroup(6)
#' h <- block_homomorphism(D, minimal_block(D, c(0, 3)))
#' group_order(hom_image(h))
#' @export
orbit_homomorphism <- function(G, omega) {
  omega <- as.integer(omega)
  codomain <- SymmetricGroup(length(omega))
  identity <- group_identity(codomain)
  images <- lapply(G$gens, function(g) identity * Permutation(match(g$array[omega + 1L], omega) - 1L))
  .schreier_sims_base(G, omega)
  H <- .hom_new(G, codomain, G$gens, images)
  bs <- basic_stabilizers(G)
  H$kernel <- if (length(bs) > length(omega)) bs[[length(omega) + 1L]] else PermutationGroup(list(group_identity(G)))
  H
}

#' @rdname orbit_homomorphism
#' @export
block_homomorphism <- function(G, blocks) {
  n <- length(blocks)
  m <- 0L
  p <- integer()
  b <- integer(n)
  for (i in seq_len(n) - 1L) if (blocks[i + 1L] == i) { p <- c(p, i); b[i + 1L] <- m; m <- m + 1L }
  for (i in seq_len(n)) b[i] <- b[blocks[i] + 1L]
  codomain <- SymmetricGroup(m)
  images <- lapply(G$gens, function(g) Permutation(b[g$array[p + 1L] + 1L]))
  .hom_new(G, codomain, G$gens, images)
}

#' @rdname orbit_homomorphism
#' @export
hom_image <- function(hom) if (inherits(hom, "FpGroupHomomorphism")) .fhom_image(hom) else .hom_image(hom)

#' @rdname orbit_homomorphism
#' @export
hom_kernel <- function(hom) if (inherits(hom, "FpGroupHomomorphism")) .fhom_kernel(hom) else .hom_kernel(hom)

#' @rdname orbit_homomorphism
#' @export
hom_apply <- function(hom, g) if (inherits(hom, "FpGroupHomomorphism")) .fhom_apply(hom, g) else .hom_apply(hom, g)

#' @rdname orbit_homomorphism
#' @export
hom_invert <- function(hom, g) if (inherits(hom, "FpGroupHomomorphism")) .fhom_invert(hom, g) else .hom_invert(hom, g)

#' Seed SymPy's random generator
#'
#' `sympy_seed(a)` is `sympy.core.random.seed(a)`: it seeds the Mersenne
#' Twister generator that the randomized group algorithms ([random_pr()],
#' [group_random()], [normal_closure()], [sylow_subgroup()]) and
#' [perm_random()] draw from, so that they reproduce SymPy's values.
#'
#' @param a An integer seed.
#' @return `NULL`, invisibly.
#' @examples
#' sympy_seed(3)
#' perm_random(5)
#' @export
sympy_seed <- function(a) .py_seed(a)
