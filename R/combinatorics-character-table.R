# Character tables of finite permutation groups, after
# sympy/combinatorics/character_table.py: Dixon's algorithm over a prime
# field for any group and the Murnaghan-Nakayama rule for S_n. Elements of
# the prime field F_p are doubles in 0..p - 1; p stays small enough that
# every product and every matrix product is exact.

.ct_gcd <- function(a, b) {
  a <- abs(a); b <- abs(b)
  while (b != 0) { t <- a %% b; a <- b; b <- t }
  a
}
.ct_lcm <- function(a, b) a / .ct_gcd(a, b) * b

# the inverse of a modulo p, by the extended Euclidean algorithm
.ct_inv <- function(a, p) {
  r0 <- p; r1 <- a %% p; s0 <- 0; s1 <- 1
  if (r1 == 0) stop(.sympy_condition("ZeroDivisionError", "inverse of zero"))
  while (r1 != 0) {
    q <- r0 %/% r1
    tmp <- r0 - q * r1; r0 <- r1; r1 <- tmp
    tmp <- s0 - q * s1; s0 <- s1; s1 <- tmp
  }
  s0 %% p
}

.ct_pow <- function(b, e, p) {
  r <- 1; b <- b %% p
  while (e > 0) {
    if (e %% 2 == 1) r <- (r * b) %% p
    b <- (b * b) %% p
    e <- e %/% 2
  }
  r
}

# reduced row echelon form modulo p, with the pivot columns
.ct_rref <- function(M, p) {
  M <- M %% p
  nr <- nrow(M); nc <- ncol(M)
  piv <- integer(); r <- 1L
  for (j in seq_len(nc)) {
    if (r > nr) break
    k <- which(M[r:nr, j] != 0)
    if (!length(k)) next
    k <- k[1L] + r - 1L
    if (k != r) { tmp <- M[r, ]; M[r, ] <- M[k, ]; M[k, ] <- tmp }
    M[r, ] <- (M[r, ] * .ct_inv(M[r, j], p)) %% p
    for (i in seq_len(nr)) if (i != r && M[i, j] != 0) M[i, ] <- (M[i, ] - M[i, j] * M[r, ]) %% p
    piv <- c(piv, j)
    r <- r + 1L
  }
  list(M = M, piv = piv)
}

# a basis (as rows) of the vectors v with M v = 0 modulo p
.ct_nullspace <- function(M, p) {
  rr <- .ct_rref(M, p)
  nc <- ncol(M)
  free <- setdiff(seq_len(nc), rr$piv)
  out <- matrix(0, length(free), nc)
  for (k in seq_along(free)) {
    v <- numeric(nc); v[free[k]] <- 1
    for (i in seq_along(rr$piv)) v[rr$piv[i]] <- (-rr$M[i, free[k]]) %% p
    out[k, ] <- v
  }
  out
}

.ct_matinv <- function(V, p) {
  n <- nrow(V)
  .ct_rref(cbind(V, diag(n)), p)$M[, n + seq_len(n), drop = FALSE]
}

# permutations as zero based array forms
.ct_af <- function(g) as.integer(g$array)
.ct_key <- function(af) paste(af, collapse = ",")
# Permutation.rmul_with_af(a, b): [a[i] for i in b]
.ct_rmul <- function(a, b) a[b + 1L]
.ct_order <- function(af) {
  n <- length(af); seen <- logical(n); o <- 1
  for (i in seq_len(n)) if (!seen[i]) {
    l <- 0; j <- i
    while (!seen[j]) { seen[j] <- TRUE; j <- af[j] + 1L; l <- l + 1 }
    o <- .ct_lcm(o, l)
  }
  o
}

.ct_lookup <- function(cc) {
  env <- new.env(hash = TRUE, parent = emptyenv())
  for (i in rev(seq_along(cc))) for (g in cc[[i]]) assign(.ct_key(.ct_af(g)), i, envir = env)
  env
}
.ct_class_of <- function(env, af) {
  k <- .ct_key(af)
  if (exists(k, envir = env, inherits = FALSE)) get(k, envir = env) else NA_integer_
}

# the class multiplication matrix of class ind (_compute_cmmatrices)
.ct_cmmatrix <- function(cc, ind, reps, lookup, p) {
  n <- length(cc)
  m <- matrix(0, n, n)
  for (g in cc[[ind]]) {
    gaf <- .ct_af(g)
    for (t in seq_len(n)) {
      i <- .ct_class_of(lookup, .ct_rmul(reps[[t]], gaf))
      if (!is.na(i)) m[i, t] <- m[i, t] + 1
    }
  }
  m %% p
}

# _dixon_prime: a prime p > 2*sqrt(order) with p = 1 (mod exponent)
.ct_dixon_prime <- function(order, exponent) {
  if (exponent == 1) return(3)
  p <- 2 * floor(sqrt(order))
  repeat {
    p <- as.numeric(nextprime(p + 1))
    if (p %% exponent == 1) break
  }
  p
}

# _simultaneous_diagonalize: the common left eigenvectors of the class
# matrices, one row for each character
.ct_simdiag <- function(n, nextmat, p) {
  spaces <- list(diag(n))
  ind <- 0L
  while (length(spaces) < n && ind < n) {
    ind <- ind + 1L
    N <- nextmat(ind)
    new <- list()
    for (S in spaces) {
      k <- nrow(S)
      if (k <= 1L) { new[[length(new) + 1L]] <- S; next }
      rr <- .ct_rref(S, p)
      T <- ((rr$M %*% N) %% p)[, rr$piv, drop = FALSE]
      for (lam in 0:(p - 1)) {
        Y <- .ct_nullspace((t(T) - lam * diag(k)) %% p, p)
        if (nrow(Y)) new[[length(new) + 1L]] <- (Y %*% rr$M) %% p
      }
    }
    spaces <- new
  }
  if (length(spaces) != n)
    stop(.sympy_condition("ValueError", "Failed to compute the common eigenspace decomposition."))
  do.call(rbind, spaces)
}

.ct_invmap <- function(cc, reps, lookup) {
  n <- length(cc); inv <- rep(-1L, n)
  for (i in seq_len(n)) {
    if (inv[i] != -1L) next
    g <- reps[[i]]
    ig <- integer(length(g)); ig[g + 1L] <- seq_along(g) - 1L
    j <- .ct_class_of(lookup, ig)
    if (!is.na(j)) { inv[i] <- j; inv[j] <- i }
  }
  inv
}

.ct_powermap <- function(cc, reps, lookup, e) {
  n <- length(cc)
  pm <- matrix(-1L, n, e)
  id <- seq_along(reps[[1L]]) - 1L
  for (i in seq_len(n)) {
    gk <- id
    for (k in seq_len(e)) {
      t <- .ct_class_of(lookup, gk)
      if (!is.na(t)) pm[i, k] <- t
      gk <- .ct_rmul(gk, reps[[i]])
    }
  }
  pm
}

# _normalize_fp: scale each eigenvector to the values of its character
.ct_normalize <- function(cc, X, p, inv_map) {
  n <- length(cc)
  sizes <- lengths(cc) %% p
  order <- sum(lengths(cc)) %% p
  lapply(seq_len(nrow(X)), function(r) {
    row <- (X[r, ] * .ct_inv(X[r, 1L], p)) %% p
    dot <- 0
    for (k in seq_len(n)) dot <- (dot + ((sizes[k] * row[k]) %% p) * row[inv_map[k]]) %% p
    chi1_sq <- (order * .ct_inv(dot, p)) %% p
    root <- sqrt_mod(chi1_sq, p)
    if (is.null(root))
      stop(.sympy_condition("ValueError", sprintf("Failed to compute sqrt(%s) on GF(%s)", chi1_sq, p)))
    (row * as.numeric(root)) %% p
  })
}

# _get_global_conductor: the least k such that all values lie in the k-th
# cyclotomic field
.ct_conductor <- function(rows, pm, m) {
  n <- length(rows)
  gal_m <- Filter(function(a) .ct_gcd(a, m) == 1, 0:(m - 1))
  for (q in rev(sort(as.numeric(primefactors(m))))) {
    while (m %% q == 0) {
      m <- m %/% q
      r <- 1 %% m
      fixed <- TRUE
      for (a in gal_m) {
        if (a %% m != r) next
        for (i in seq_len(n)) {
          row <- rows[[i]]
          for (j in seq_len(n)) if (row[pm[j, a + 1L]] != row[j]) { fixed <- FALSE; break }
          if (!fixed) break
        }
        if (!fixed) break
      }
      if (!fixed) { m <- m * q; break }
    }
  }
  m
}

# Python list comparison of two numeric vectors
.ct_cmp_vec <- function(a, b) {
  for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(if (a[i] < b[i]) -1L else 1L)
  as.integer(sign(length(a) - length(b)))
}
.ct_cmp_row <- function(a, b) {
  for (i in seq_along(a)) {
    c <- .ct_cmp_vec(a[[i]], b[[i]])
    if (c != 0L) return(c)
  }
  0L
}

# _sort_characters: sort the rows, then swap the trivial character into
# the first row
.ct_sort_characters <- function(rows) {
  o <- integer()
  for (i in seq_along(rows)) {
    j <- length(o)
    while (j > 0L && .ct_cmp_row(rows[[o[j]]], rows[[i]]) > 0L) j <- j - 1L
    o <- append(o, i, after = j)
  }
  rows <- rows[o]
  for (i in seq_along(rows)) if (all(vapply(rows[[i]], function(v) identical(as.numeric(v), 1), TRUE))) {
    tmp <- rows[[1L]]; rows[[1L]] <- rows[[i]]; rows[[i]] <- tmp
    break
  }
  rows
}

# _lift_to_minimal_field: the integer table, or the values in the k-th
# cyclotomic field as their dense coefficient vectors (highest power first,
# without leading zeros)
.ct_lift <- function(rows, pm, k, e, p) {
  n <- length(rows)
  half <- p %/% 2
  sgn <- function(v) ifelse(v <= half, v, v - p)
  if (k == 1) return(.ct_sort_characters(lapply(rows, function(r) as.list(sgn(r)))))
  x <- .ct_pow(as.numeric(primitive_root(p)), (p - 1) / k, p)
  gal <- Filter(function(a) .ct_gcd(a, k) == 1, 0:(k - 1))
  phi <- length(gal)
  V <- matrix(0, phi, phi)
  for (ai in seq_len(phi)) {
    xa <- .ct_pow(x, gal[ai], p)
    for (i in seq_len(phi)) V[ai, i] <- .ct_pow(xa, i - 1, p)
  }
  Vinv <- .ct_matinv(V, p)
  lifted <- vapply(gal, function(a) { A <- a; while (.ct_gcd(A, e) != 1) A <- A + k; A }, 0)
  dM <- lapply(seq_len(n), function(i) {
    row <- rows[[i]]
    lapply(seq_len(n), function(j) {
      b <- vapply(lifted, function(A) row[pm[j, (A %% e) + 1L]], 0)
      cf <- rev(sgn(as.vector((Vinv %*% b) %% p)))
      nz <- which(cf != 0)
      if (length(nz)) cf[nz[1L]:length(cf)] else numeric()
    })
  })
  .ct_sort_characters(dM)
}

.ct_new <- function(rows, k, reps) structure(list(rows = rows, zeta_order = as.integer(k), reps = reps),
                                             class = "CharacterTable")

#' Character tables of finite groups
#'
#' `character_table(G)` computes the character table of a permutation
#' group, as SymPy's `PermutationGroup.character_table()` and
#' `CharacterTable.from_perm_group()`: rows are the irreducible
#' characters, with the trivial character first, and columns the
#' conjugacy classes in the order of [conjugacy_classes()]. Symmetric groups
#' use `character_table_sn()`, the Murnaghan-Nakayama rule, whose rows are
#' indexed by the partitions of `n` in reverse lexicographic order and whose
#' columns by the partitions in lexicographic order. Other groups use
#' `dixon_character_table()`, Dixon's algorithm on a list of conjugacy
#' classes over a prime field, with the values lifted to the smallest
#' cyclotomic field `QQ(zeta_k)` that holds them.
#'
#' The table is a `CharacterTable` object: `ct_tolist()` gives its rows as
#' lists of expressions in the root of unity `zeta_k` (an
#' [AlgebraicNumber()]), `ct_as_matrix()` the matrix, `ct_shape()` its
#' dimensions, `ct_zeta_order()` the order `k` of the root of unity (1 for
#' an integer table) and `ct_conjugacy_class_reps()` the first element of
#' each conjugacy class.
#'
#' @param G A permutation group.
#' @param conjugacy_classes A list of conjugacy classes, each a list of
#'   permutations; one of them must hold the identity.
#' @param n The degree of the symmetric group, a positive integer.
#' @param tbl A `CharacterTable`.
#' @return A `CharacterTable`; the accessors return a list of lists of
#'   expressions, a matrix, an integer vector, an integer or a list of
#'   permutations.
#' @examples
#' ct_tolist(character_table_sn(3))
#' ct_as_matrix(character_table(AlternatingGroup(4)))
#' @export
character_table <- function(G) {
  if (isTRUE(is_symmetric_group(G))) return(character_table_sn(G$degree))
  dixon_character_table(conjugacy_classes(G))
}

#' @rdname character_table
#' @export
dixon_character_table <- function(conjugacy_classes) {
  cc <- lapply(conjugacy_classes, as.list)
  identity <- NA_integer_
  for (i in seq_along(cc)) {
    af <- .ct_af(cc[[i]][[1L]])
    if (identical(af, seq_along(af) - 1L)) { identity <- i; break }
  }
  if (is.na(identity)) stop(.sympy_condition("ValueError", "No identity class found"))
  if (identity != 1L) { tmp <- cc[[identity]]; cc[[identity]] <- cc[[1L]]; cc[[1L]] <- tmp }
  reps <- lapply(cc, function(c) .ct_af(c[[1L]]))
  order <- sum(lengths(cc))
  exponent <- Reduce(.ct_lcm, lapply(reps, .ct_order), 1)
  p <- .ct_dixon_prime(order, exponent)
  lookup <- .ct_lookup(cc)
  n <- length(cc)
  X <- .ct_simdiag(n, function(ind) .ct_cmmatrix(cc, ind, reps, lookup, p), p)
  rows <- .ct_normalize(cc, X, p, .ct_invmap(cc, reps, lookup))
  pm <- .ct_powermap(cc, reps, lookup, exponent)
  k <- .ct_conductor(rows, pm, exponent)
  out <- .ct_lift(rows, pm, k, exponent, p)
  # restore the place of the identity class
  if (identity != 1L) out <- lapply(out, function(r) { r[c(1L, identity)] <- r[c(identity, 1L)]; r })
  .ct_new(out, k, lapply(conjugacy_classes, function(c) as.list(c)[[1L]]))
}

# partitions of n with decreasing parts, in lexicographic order
.ct_partitions <- function(n, mx = n) {
  if (n == 0L) return(list(integer()))
  out <- list()
  for (f in seq_len(min(n, mx))) for (r in .ct_partitions(n - f, f)) out[[length(out) + 1L]] <- c(f, r)
  out
}

#' @rdname character_table
#' @export
character_table_sn <- function(n) {
  if (n <= 0) stop(.sympy_condition("ValueError", "n must be a positive integer"))
  n <- as.integer(n)
  partitions <- .ct_partitions(n)
  cache <- new.env(hash = TRUE, parent = emptyenv())
  beta_of <- function(lam) vapply(0:(n - 1L), function(i) (if (i < length(lam)) lam[i + 1L] else 0L) + n - 1L - i, 0L)
  beta_to_partition <- function(beta) {
    idx <- 0:(n - 1L)
    keep <- beta > n - 1L - idx
    as.integer(beta[keep] - (n - 1L - idx[keep]))
  }
  character <- function(lam, mu) {
    key <- paste(paste(lam, collapse = ","), paste(mu, collapse = ","), sep = "|")
    v <- cache[[key]]
    if (!is.null(v)) return(v)
    if (!length(mu)) {
      v <- as.integer(!length(lam))
    } else {
      beta <- beta_of(lam)
      cl <- mu[1L]
      v <- 0L
      for (i in seq_len(n)) {
        nb <- beta[i] - cl
        if (nb < 0L || nb %in% beta) next
        height <- if (i < n) sum(beta[(i + 1L):n] > nb) else 0L
        newb <- beta[-i]
        at <- 0L
        while (at < length(newb) && newb[at + 1L] > nb) at <- at + 1L
        newb <- append(newb, nb, after = at)
        v <- v + (if (height %% 2L == 1L) -1L else 1L) * character(beta_to_partition(newb), mu[-1L])
      }
    }
    assign(key, v, envir = cache)
    v
  }
  representative <- function(partition) {
    af <- seq_len(n) - 1L
    start <- 0L
    for (cl in partition) {
      for (index in seq_len(cl) - 1L) af[start + index + 1L] <- start + (index + 1L) %% cl
      start <- start + cl
    }
    .perm_af(af)
  }
  rows <- lapply(rev(partitions), function(r) lapply(partitions, function(cc) character(r, cc)))
  .ct_new(rows, 1L, lapply(partitions, representative))
}

.ct_entry <- function(v, z) {
  if (is.null(z)) return(.int(v))
  if (!length(v)) return(S_Zero)
  d <- length(v)
  terms <- list()
  for (i in seq_len(d)) if (v[i] != 0) {
    e <- d - i
    terms[[length(terms) + 1L]] <- if (e == 0) .int(v[i]) else .mul(list(.int(v[i]), if (e == 1) z else .pow(z, .int(e))))
  }
  .add(terms)
}

#' @rdname character_table
#' @export
ct_tolist <- function(tbl) {
  k <- tbl$zeta_order
  z <- if (k == 1L) NULL else
    AlgebraicNumber(list(cyclotomic_poly(k, Dummy("x")), exp(.mul(list(S_Two, S_Pi, sym_I, .rat(1, k))))),
                    alias = paste0("zeta", k))
  lapply(tbl$rows, function(r) lapply(r, .ct_entry, z = z))
}

#' @rdname character_table
#' @export
ct_as_matrix <- function(tbl) Matrix(ct_tolist(tbl))

#' @rdname character_table
#' @export
ct_shape <- function(tbl) c(length(tbl$rows), length(tbl$rows))

#' @rdname character_table
#' @export
ct_zeta_order <- function(tbl) tbl$zeta_order

#' @rdname character_table
#' @export
ct_conjugacy_class_reps <- function(tbl) tbl$reps

#' @export
format.CharacterTable <- function(x, ...) sstr(ct_as_matrix(x))

#' @export
print.CharacterTable <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
