# Tensor canonicalization (Butler-Portugal algorithm), after
# sympy/combinatorics/tensor_can.py, with the BSGS helpers it needs from
# sympy/combinatorics/util.py, perm_groups.py and testutil.py ported here
# under a `.tc_` prefix.
#
# Permutations are 0 based array forms. Internally they are integer vectors;
# exported functions take and return `Permutation` objects (R/combinatorics.R)
# where SymPy uses `Permutation`, and plain integer vectors where SymPy uses
# lists. A Python tuple is returned as an unnamed list.

# ---------------------------------------------------------------------------
# array form helpers (sympy/combinatorics/permutations.py)

.tc_af_new <- function(a) structure(list(array = as.integer(a)), class = "Permutation")
.tc_af <- function(p) if (inherits(p, "Permutation")) p$array else as.integer(p)
.tc_invert <- function(a) {
  inv <- integer(length(a))
  inv[a + 1L] <- seq_along(a) - 1L
  inv
}
.tc_is_id <- function(a) all(a == seq_along(a) - 1L)
.tc_key <- function(a) paste(a, collapse = ",")
.tc_value_error <- function(msg) stop(.sympy_condition("ValueError", msg))

# Python list comparison a < b for integer vectors
.tc_list_lt <- function(a, b) {
  n <- min(length(a), length(b))
  d <- which(a[seq_len(n)] != b[seq_len(n)])
  if (length(d)) a[d[1L]] < b[d[1L]] else length(a) < length(b)
}

# order of a list of equal length integer vectors, lexicographic and stable
# like Python's sort
.tc_lex_order <- function(l) {
  if (length(l) < 2L) return(seq_along(l))
  m <- do.call(rbind, l)
  do.call(order, lapply(seq_len(ncol(m)), function(k) m[, k]))
}

# ---------------------------------------------------------------------------
# orbits and transversals (perm_groups.py: _orbit, _orbit_transversal)

# The orbit as a sorted vector. SymPy returns a set; its callers here only
# take minima or build lists whose order does not change any result.
.tc_orbit <- function(degree, gens, alpha) {
  orb <- alpha
  used <- logical(degree)
  used[alpha + 1L] <- TRUE
  k <- 1L
  while (k <= length(orb)) {
    b <- orb[k]
    for (gen in gens) {
      t <- gen[b + 1L]
      if (!used[t + 1L]) {
        orb <- c(orb, t)
        used[t + 1L] <- TRUE
      }
    }
    k <- k + 1L
  }
  sort(orb)
}

# list(x = orbit points in order of discovery, p = transversal elements)
.tc_orbit_transversal <- function(degree, gens, alpha) {
  xs <- alpha
  ps <- list(seq_len(degree) - 1L)
  used <- logical(degree)
  used[alpha + 1L] <- TRUE
  k <- 1L
  while (k <= length(xs)) {
    x <- xs[k]
    px <- ps[[k]]
    for (gen in gens) {
      t <- gen[x + 1L]
      if (!used[t + 1L]) {
        xs <- c(xs, t)
        ps[[length(ps) + 1L]] <- gen[px + 1L]
        used[t + 1L] <- TRUE
      }
    }
    k <- k + 1L
  }
  list(x = xs, p = ps)
}

# ---------------------------------------------------------------------------
# BSGS helpers (util.py)

.tc_distribute_gens_by_base <- function(base, gens) {
  base_len <- length(base)
  degree <- length(gens[[1L]])
  stabs <- rep(list(list()), base_len)
  max_stab_index <- 0L
  for (gen in gens) {
    j <- 0L
    while (j < base_len - 1L && gen[base[j + 1L] + 1L] == base[j + 1L]) j <- j + 1L
    if (j > max_stab_index) max_stab_index <- j
    for (k in 0:j) stabs[[k + 1L]] <- c(stabs[[k + 1L]], list(gen))
  }
  if (max_stab_index + 1L < base_len) {
    for (i in (max_stab_index + 1L):(base_len - 1L)) {
      stabs[[i + 1L]] <- c(stabs[[i + 1L]], list(seq_len(degree) - 1L))
    }
  }
  stabs
}

.tc_orbits_transversals_from_bsgs <- function(base, strong_gens_distr) {
  degree <- length(strong_gens_distr[[1L]][[1L]])
  lapply(seq_along(base), function(i) .tc_orbit_transversal(degree, strong_gens_distr[[i]], base[i]))
}

# _strip_af: returns list(h, j) with h NULL where SymPy returns False
.tc_strip_af <- function(h, base, transversals, j) {
  base_len <- length(base)
  i <- j + 1L
  while (i < base_len) {
    beta <- h[base[i + 1L] + 1L]
    if (beta != base[i + 1L]) {
      tr <- transversals[[i + 1L]]
      pos <- match(beta, tr$x)
      if (is.na(pos)) return(list(h = h, j = i + 1L))
      u <- tr$p[[pos]]
      if (identical(h, u)) return(list(h = NULL, j = base_len + 1L))
      h <- .tc_invert(u)[h + 1L]
    }
    i <- i + 1L
  }
  list(h = h, j = base_len + 1L)
}

# the generators PermutationGroup(gens) keeps: equal sizes, no duplicates,
# no identities unless the group is trivial
.tc_group_gens <- function(gens) {
  if (!length(gens)) gens <- list(integer())
  sizes <- lengths(gens)
  if (length(unique(sizes)) > 1L) {
    degree <- max(sizes)
    gens <- lapply(gens, function(a) if (length(a) < degree) c(a, length(a):(degree - 1L)) else a)
  }
  gens <- gens[!duplicated(vapply(gens, .tc_key, ""))]
  if (length(gens) > 1L) {
    non_identity <- Filter(function(g) !.tc_is_id(g), gens)
    gens <- if (length(non_identity)) non_identity else gens[1L]
  }
  gens
}

# PermutationGroup(gens).schreier_sims_incremental(base); gens are the
# group's generators as returned by .tc_group_gens
.tc_schreier_sims_incremental <- function(gens, base = integer()) {
  base <- as.integer(base)
  degree <- length(gens[[1L]])
  id_af <- seq_len(degree) - 1L
  if (length(gens) == 1L && .tc_is_id(gens[[1L]])) return(list(base = base, gens = gens))
  .base <- base
  .gens <- Filter(function(x) !.tc_is_id(x), gens)
  for (gen in .gens) {
    if (all(gen[.base + 1L] == .base)) .base <- c(.base, which(gen != id_af)[1L] - 1L)
  }
  distr <- .tc_distribute_gens_by_base(.base, .gens)
  base_len <- length(.base)
  trs <- lapply(seq_len(base_len), function(i) .tc_orbit_transversal(degree, distr[[i]], .base[i]))
  new_gens <- list()
  i <- base_len - 1L
  while (i >= 0L) {
    continue_i <- FALSE
    tr <- trs[[i + 1L]]
    for (bi in seq_along(tr$x)) {
      beta <- tr$x[bi]
      u_beta <- tr$p[[bi]]
      for (gen in distr[[i + 1L]]) {
        gb <- gen[beta + 1L]
        u1 <- tr$p[[match(gb, tr$x)]]
        g1 <- gen[u_beta + 1L]
        if (!identical(g1, u1)) {
          y <- TRUE
          schreier_gen <- .tc_invert(u1)[g1 + 1L]
          st <- .tc_strip_af(schreier_gen, .base, trs, i)
          h <- st$h
          j <- st$j
          if (j <= base_len) {
            y <- FALSE
          } else if (!is.null(h)) {
            y <- FALSE
            .base <- c(.base, which(h != id_af)[1L] - 1L)
            base_len <- base_len + 1L
            distr[[base_len]] <- list()
          }
          if (!y) {
            new_gens <- c(new_gens, list(h))
            for (l in (i + 1L):(j - 1L)) {
              distr[[l + 1L]] <- c(distr[[l + 1L]], list(h))
              trs[[l + 1L]] <- .tc_orbit_transversal(degree, distr[[l + 1L]], .base[l + 1L])
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
  list(base = .base, gens = c(.gens, new_gens))
}

# PermutationGroup(gens).pointwise_stabilizer(points).generators
.tc_pointwise_stabilizer <- function(gens, points) {
  points <- as.integer(points)
  r <- .tc_schreier_sims_incremental(gens, points)
  stab <- Filter(function(gen) identical(gen[points + 1L], points), r$gens)
  if (!length(stab)) stab <- list(seq_len(length(gens[[1L]])) - 1L)
  .tc_group_gens(stab)
}

# all elements of the group generated by gens (PermutationGroup.generate)
.tc_group_elements <- function(gens) {
  els <- list(seq_len(length(gens[[1L]])) - 1L)
  seen <- new.env(hash = TRUE)
  assign(.tc_key(els[[1L]]), TRUE, envir = seen)
  k <- 1L
  while (k <= length(els)) {
    for (gen in gens) {
      e <- gen[els[[k]] + 1L]
      ke <- .tc_key(e)
      if (!exists(ke, envir = seen, inherits = FALSE)) {
        assign(ke, TRUE, envir = seen)
        els[[length(els) + 1L]] <- e
      }
    }
    k <- k + 1L
  }
  els
}

# ---------------------------------------------------------------------------
# tensor_can.py

.tc_dummy_sgs <- function(dummies, sym, n) {
  dummies <- as.integer(dummies)
  n <- as.integer(n)
  if (length(dummies) > n) .tc_value_error("List too large")
  L <- length(dummies)
  res <- list()
  # exchange of contravariant and covariant indices
  if (!is.null(sym)) {
    for (j in dummies[seq_len(L) %% 2L == 1L]) {
      a <- 0:(n + 1L)
      if (isTRUE(sym == 1)) {
        a[n + 1L] <- n + 1L
        a[n + 2L] <- n
      }
      a[c(j + 1L, j + 2L)] <- a[c(j + 2L, j + 1L)]
      res[[length(res) + 1L]] <- a
    }
  }
  # rename dummy indices
  for (j in dummies[seq_len(max(L - 3L, 0L))][seq_len(max(L - 3L, 0L)) %% 2L == 1L]) {
    a <- 0:(n + 1L)
    a[(j + 1L):(j + 4L)] <- a[c(j + 3L, j + 4L, j + 1L, j + 2L)]
    res[[length(res) + 1L]] <- a
  }
  res
}

#' Strong generators of the dummy index symmetry
#'
#' Port of SymPy's `dummy_sgs()`. `dummies[2k]` and `dummies[2k+1]` are
#' paired (contravariant, covariant) dummy indices in consecutive positions.
#' The generators exchange the two indices of a pair (with a sign change for
#' an antisymmetric metric) and rename consecutive pairs.
#'
#' @param dummies Integer vector of dummy indices (0 based).
#' @param sym Metric symmetry: `NULL` (none), `0` (symmetric) or `1`
#'   (antisymmetric).
#' @param n Number of indices.
#' @return A list of array forms (integer vectors of length `n + 2`).
#' @examples
#' dummy_sgs(2:7, 0, 8)
#' @export
dummy_sgs <- function(dummies, sym, n) .tc_dummy_sgs(dummies, sym, n)

# list of minima of the orbits of indices in the group of dummies
.tc_min_dummies <- function(dummies, sym, indices) {
  num_types <- length(sym)
  m <- lapply(dummies, function(dx) if (length(dx)) min(dx) else NULL)
  res <- indices
  # SymPy repeats this loop num_types times with the same result
  if (num_types) {
    for (c in seq_along(indices)) {
      i <- indices[c]
      for (j in seq_len(num_types)) {
        if (i %in% dummies[[j]]) {
          res[c] <- m[[j]]
          break
        }
      }
    }
  }
  res
}

# the representative h satisfying s[h[b]] == j, or NULL
.tc_trace_S <- function(s, j, b, S_cosets) {
  for (h in S_cosets[[b + 1L]]) if (s[h[b + 1L] + 1L] == j) return(h)
  NULL
}

# the representative h satisfying h[gj] == p_i, or NULL
.tc_trace_D <- function(gj, p_i, Dxtrav) {
  for (h in Dxtrav) if (h[gj + 1L] == p_i) return(h)
  NULL
}

# remove p0 and its paired index from dumx
.tc_dumx_remove <- function(dumx, p0) {
  for (k in seq_along(dumx)) {
    dx <- dumx[[k]]
    pos <- match(p0, dx)
    if (is.na(pos)) next
    p0_paired <- if ((pos - 1L) %% 2L == 0L) dx[pos + 1L] else dx[pos - 1L]
    dx <- dx[-match(p0, dx)]
    dx <- dx[-match(p0_paired, dx)]
    dumx[[k]] <- dx
  }
  dumx
}

#' Coset representatives from transversals
#'
#' Port of SymPy's `transversal2coset()`: for each point `i < size`, the
#' sorted transversal elements when `i` is a base point, else the identity;
#' trailing identity entries are dropped.
#'
#' @param size Permutation size.
#' @param base Base points (0 based).
#' @param transversal Transversals as returned by [get_transversals()].
#' @return A list of lists of array forms.
#' @examples
#' b <- riemann_bsgs[[1L]]
#' transversal2coset(6, b, get_transversals(b, riemann_bsgs[[2L]]))
#' @export
transversal2coset <- function(size, base, transversal) {
  size <- as.integer(size)
  idl <- list(seq_len(size) - 1L)
  a <- list()
  j <- 1L
  for (i in seq_len(size) - 1L) {
    if (i %in% base) {
      vals <- unname(lapply(transversal[[j]], as.integer))
      a[[i + 1L]] <- vals[.tc_lex_order(vals)]
      j <- j + 1L
    } else {
      a[[i + 1L]] <- idl
    }
  }
  j <- length(a)
  while (identical(a[[j]], idl)) j <- j - 1L
  a[seq_len(j)]
}

#' Butler-Portugal canonical representative of a double coset
#'
#' Port of SymPy's `double_coset_can_rep()`, the Butler-Portugal algorithm
#' for tensor canonicalization with dummy indices. The tensor is the
#' permutation `g` (size `n + 2`, the last two points carry the sign); the
#' slot symmetry group `S` has minimal BSGS `b_S`, `sgens`; the dummy
#' symmetry group `D` is built from `dummies` and `sym`. The result is the
#' minimum of `D*g*S` under lexicographic order, or `0` when the tensor
#' vanishes.
#'
#' @param dummies List of integer vectors of dummy indices, one per index
#'   type, ordered `d0, -d0, d1, -d1, ...`.
#' @param sym List of metric symmetries, one per index type: `NULL`, `0` or
#'   `1`.
#' @param b_S Base of a minimal slot symmetry BSGS.
#' @param sgens Strong generators of the slot symmetry (permutations).
#' @param S_transversals Transversals of the slot BSGS, from
#'   [get_transversals()].
#' @param g The permutation representing the tensor.
#' @return An integer vector (array form), or `0`.
#' @examples
#' gens <- list(Permutation(c(2, 1, 0, 3, 4, 5, 7, 6)), Permutation(c(4, 1, 2, 3, 0, 5, 7, 6)))
#' base <- c(0, 2)
#' tr <- get_transversals(base, gens)
#' double_coset_can_rep(list(0:5), list(0), base, gens, tr, Permutation(c(4, 2, 0, 1, 3, 5, 6, 7)))
#' double_coset_can_rep(list(0:5), list(0), base, gens, tr, Permutation(c(4, 1, 3, 0, 5, 2, 6, 7)))
#' @export
double_coset_can_rep <- function(dummies, sym, b_S, sgens, S_transversals, g) {
  g <- .tc_af(g)
  size <- length(g)
  num_dummies <- size - 2L
  indices <- seq_len(num_dummies) - 1L
  if (!is.list(sym)) sym <- as.list(sym)
  all_metrics_with_sym <- !any(vapply(sym, is.null, NA))
  num_types <- length(sym)
  dumx <- lapply(dummies, as.integer)
  b_S <- as.integer(b_S)
  sgensx <- lapply(sgens, .tc_af)
  if (length(b_S)) S_transversals <- transversal2coset(size, b_S, S_transversals)
  # strong generating set for D
  dsgs <- function() {
    res <- list()
    for (k in seq_len(num_types)) res <- c(res, .tc_dummy_sgs(dumx[[k]], sym[[k]], num_dummies))
    res
  }
  dsgsx <- dsgs()
  idn <- seq_len(size) - 1L
  # TAB = list of entries (s, d, h) where h = d*g*s
  TAB <- list(list(s = idn, d = idn, h = g))
  for (i in seq_len(size - 2L) - 1L) {
    b <- i
    testb <- (b %in% b_S) && length(sgensx) > 0L
    deltab <- if (testb) .tc_orbit(size, sgensx, b) else b
    # p1 = min(IMAGES) = min(Union D_p*h*deltab for h in TAB)
    md <- if (all_metrics_with_sym) {
      .tc_min_dummies(dumx, sym, indices)
    } else {
      vapply(seq_len(size - 2L) - 1L, function(ii) min(.tc_orbit(size, dsgsx, ii)), 0L)
    }
    p_i <- min(vapply(TAB, function(t) min(md[t$h[deltab + 1L] + 1L]), 0L))
    Dxtrav <- if (length(dsgsx)) lapply(.tc_orbit_transversal(size, dsgsx, p_i)$p, .tc_invert) else NULL
    # compute the orbit of p_i
    deltap <- p_i
    for (ii in seq_len(num_types)) {
      if (p_i %in% dumx[[ii]]) {
        if (!is.null(sym[[ii]])) {
          deltap <- dumx[[ii]]
        } else {
          p_i_index <- (match(p_i, dumx[[ii]]) - 1L) %% 2L
          dxi <- dumx[[ii]]
          deltap <- dxi[seq_along(dxi) > p_i_index & (seq_along(dxi) - 1L - p_i_index) %% 2L == 0L]
        }
        break
      }
    }
    TAB1 <- list()
    while (length(TAB)) {
      ent <- TAB[[length(TAB)]]
      TAB[[length(TAB)]] <- NULL
      s <- ent$s
      d <- ent$d
      h <- ent$h
      mdh <- md[h[deltab + 1L] + 1L]
      if (min(mdh) != p_i) next
      deltab1 <- deltab[mdh == p_i]
      # NEXT = s*deltab1 intersection (d*g)**-1*deltap
      dg <- d[g + 1L]
      dginv <- .tc_invert(dg)
      sdeltab <- s[deltab1 + 1L]
      gdeltap <- dginv[deltap + 1L]
      NEXT <- sdeltab[sdeltab %in% gdeltap]
      for (j in NEXT) {
        if (testb) {
          s1 <- .tc_trace_S(s, j, b, S_transversals)
          if (is.null(s1)) next
          s1 <- s[s1 + 1L]
        } else {
          s1 <- s
        }
        if (!is.null(Dxtrav)) {
          d1 <- .tc_trace_D(dg[j + 1L], p_i, Dxtrav)
          if (is.null(d1)) next
        } else {
          if (p_i != dg[j + 1L]) next
          d1 <- idn
        }
        stopifnot(d1[dg[j + 1L] + 1L] == p_i)
        d1 <- d1[d + 1L]
        h1 <- d1[g[s1 + 1L] + 1L]
        TAB1[[length(TAB1) + 1L]] <- list(s = s1, d = d1, h = h1)
      }
    }
    # if TAB contains equal permutations, keep only one of them;
    # if TAB contains equal permutations up to the sign, return 0
    TAB1 <- TAB1[.tc_lex_order(lapply(TAB1, `[[`, "h"))]
    prev <- integer(size)
    body <- seq_len(size - 2L)
    while (length(TAB1)) {
      ent <- TAB1[[length(TAB1)]]
      TAB1[[length(TAB1)]] <- NULL
      h <- ent$h
      if (identical(h[body], prev[body])) {
        if (h[size] != prev[size]) return(0)
      } else {
        TAB[[length(TAB) + 1L]] <- ent
      }
      prev <- h
    }
    # stabilize the SGS
    sgensx <- Filter(function(h) h[b + 1L] == b, sgensx)
    if (b %in% b_S) b_S <- b_S[-match(b, b_S)]
    dumx <- .tc_dumx_remove(dumx, p_i)
    dsgsx <- dsgs()
  }
  TAB[[1L]]$h
}

#' Canonicalize a tensor with respect to its free indices
#'
#' Port of SymPy's `canonical_free()`: uses the slot symmetry with BSGS
#' `base`, `gens` to bring the free indices (the first `num_free` indices)
#' to the lowest slots in lexicographic order.
#'
#' @param base,gens BSGS of the slot permutation group.
#' @param g The permutation representing the tensor.
#' @param num_free Number of free indices.
#' @return An integer vector (array form).
#' @examples
#' gens <- list(Permutation(c(1, 0, 2, 3, 5, 4)), Permutation(c(2, 3, 0, 1, 4, 5)),
#'              Permutation(c(0, 1, 3, 2, 5, 4)))
#' canonical_free(c(0, 2), gens, Permutation(c(2, 1, 0, 3, 4, 5)), 4)
#' @export
canonical_free <- function(base, gens, g, num_free) {
  g <- .tc_af(g)
  size <- length(g)
  if (!length(base)) return(g)
  transversals <- .tc_get_transversals(as.integer(base), lapply(gens, .tc_af))
  # SymPy then appends the other points to `base`; that does not change h
  h <- g
  free_pts <- seq_len(num_free) - 1L
  for (transv in transversals) {
    h_i <- rep(size, num_free)
    s <- NULL
    for (sk in transv$p) {
      h1 <- h[sk + 1L]
      hi <- match(free_pts, h1) - 1L
      if (.tc_list_lt(hi, h_i)) {
        h_i <- hi
        s <- sk
      }
    }
    if (!is.null(s)) h <- h[s + 1L]
  }
  h
}

.tc_get_map_slots <- function(size, fixed_slots) {
  res <- seq_len(size) - 1L
  pos <- 0L
  for (i in seq_len(size) - 1L) {
    if (i %in% fixed_slots) next
    res[i + 1L] <- pos
    pos <- pos + 1L
  }
  res
}

.tc_lift_sgens <- function(size, fixed_slots, free, s) {
  a <- integer()
  j <- k <- 1L
  fd <- free[order(fixed_slots, free)]
  num_free <- length(free)
  for (i in seq_len(size) - 1L) {
    if (i %in% fixed_slots) {
      a <- c(a, fd[k])
      k <- k + 1L
    } else {
      a <- c(a, s[j] + num_free)
      j <- j + 1L
    }
  }
  a
}

.tc_msym_ok <- function(m) is.null(m) || (length(m) == 1L && !is.na(m) && m %in% c(0, 1))

#' Canonicalize a tensor formed by component tensors
#'
#' Port of SymPy's `canonicalize()`. The tensor is the permutation `g`; the
#' free indices come first, then the dummy indices in the order
#' `d0, -d0, d1, -d1, ...`. Each argument in `...` describes a type of
#' component tensor as `list(base_i, gens_i, n_i, sym_i)`: a BSGS of its
#' slot symmetry, the number of such tensors, and their commutation
#' symmetry (`NULL` none, `0` commuting, `1` anticommuting). The algorithm
#' first minimizes the free indices with [canonical_free()], then the dummy
#' indices with [double_coset_can_rep()]. `riemann_bsgs` is the BSGS of the
#' Riemann tensor symmetry.
#'
#' @param g The permutation representing the tensor.
#' @param dummies Integer vector of dummy indices, or a list of them, one per
#'   index type.
#' @param msym Metric symmetry (`NULL`, `0` or `1`), or a list of them, one
#'   per index type.
#' @param ... Tensor types, each `list(base_i, gens_i, n_i, sym_i)`.
#' @return The canonical array form (integer vector), or `0` when the tensor
#'   vanishes.
#' @examples
#' b2a <- get_symmetric_group_sgs(2, 1)
#' g <- Permutation(c(1, 3, 0, 5, 4, 2, 6, 7))
#' canonicalize(g, 0:5, 0, list(b2a[[1]], b2a[[2]], 1, 0), list(b2a[[1]], b2a[[2]], 2, 0))
#' canonicalize(g, 0:5, 0, list(b2a[[1]], b2a[[2]], 1, 0), list(b2a[[1]], b2a[[2]], 2, 1))
#' canonicalize(Permutation(c(0, 2, 3, 1, 4, 5)), 2:3, 0,
#'              list(riemann_bsgs[[1]], riemann_bsgs[[2]], 1, 0))
#' @export
canonicalize <- function(g, dummies, msym, ...) {
  v <- list(...)
  if (!is.list(msym) && length(msym) > 1L) msym <- as.list(msym)
  if (!is.list(msym)) {
    if (!.tc_msym_ok(msym)) .tc_value_error("msym must be 0, 1 or None")
    num_types <- 1L
  } else {
    num_types <- length(msym)
    if (!all(vapply(msym, .tc_msym_ok, NA))) .tc_value_error("msym entries must be 0, 1 or None")
    if (length(dummies) != num_types) .tc_value_error("dummies and msym must have the same number of elements")
  }
  g_af <- .tc_af(g)
  size <- length(g_af)
  v1 <- list()
  for (t in v) {
    base_i <- as.integer(t[[1L]])
    gens_i <- lapply(t[[2L]], .tc_af)
    # check that the BSGS is minimal; this property is used in
    # double_coset_can_rep; if it is not minimal use canonicalize_naive
    if (!.tc_is_minimal_bsgs(base_i, gens_i)) {
      mbsgs <- .tc_get_minimal_bsgs(gens_i)
      if (is.null(mbsgs)) return(.tc_canonicalize_naive(g, dummies, msym, ...))
      base_i <- mbsgs$base
      gens_i <- mbsgs$gens
    }
    v1[[length(v1) + 1L]] <- list(base_i, gens_i, rep(list(integer()), t[[3L]]), t[[4L]])
  }
  if (num_types == 1L && !is.list(msym)) {
    dummies <- list(dummies)
    msym <- list(msym)
  }
  dummies <- lapply(dummies, as.integer)
  flat_dummies <- unlist(dummies)
  if (is.null(flat_dummies)) flat_dummies <- integer()
  nfd <- length(flat_dummies)
  if (nfd) {
    rng <- if (flat_dummies[nfd] >= flat_dummies[1L]) flat_dummies[1L]:flat_dummies[nfd] else integer()
    if (!identical(flat_dummies, rng)) .tc_value_error("dummies is not valid")
  }
  # slot symmetry of the tensor
  gp <- .tc_gens_products(v1)
  if (size != gp$size) .tc_value_error(sprintf("g has size %d, generators have size %d", size, gp$size))
  free <- (seq_len(size - 2L) - 1L)[!(seq_len(size - 2L) - 1L) %in% flat_dummies]
  num_free <- length(free)
  # g1 minimal tensor under slot symmetry
  g1 <- canonical_free(gp$base, gp$gens, g_af, num_free)
  if (!nfd) return(g1)
  # save the sign of g1
  sign <- if (g1[size] == size - 1L) 0L else 1L
  # the free indices are kept fixed; free_i lists the slots of each
  # component tensor occupied by free indices
  start <- 0L
  for (i in seq_along(v)) {
    t <- v[[i]]
    free_i <- list()
    len_tens <- length(.tc_af(t[[2L]][[1L]])) - 2L
    for (j in seq_len(t[[3L]])) {
      h <- g1[start + seq_len(len_tens)]
      fr <- integer()
      for (k in free) if (k %in% h) fr <- c(fr, match(k, h) - 1L)
      free_i[[length(free_i) + 1L]] <- fr
      start <- start + len_tens
    }
    v1[[i]] <- list(as.integer(t[[1L]]), lapply(t[[2L]], .tc_af), free_i, t[[4L]])
  }
  # BSGS of the tensor with fixed free indices
  gp <- .tc_gens_products(v1)
  size <- gp$size
  # reduce the permutations getting rid of the free indices
  pos_free <- match(seq_len(num_free) - 1L, g1) - 1L
  size_red <- size - num_free
  g1_red <- g1[g1 %in% flat_dummies] - num_free
  g1_red <- c(g1_red, if (sign) c(size_red - 1L, size_red - 2L) else c(size_red - 2L, size_red - 1L))
  map_slots <- .tc_get_map_slots(size, pos_free)
  sbase_red <- map_slots[gp$base[!gp$base %in% pos_free] + 1L]
  sgens_red <- lapply(gp$gens, function(y) map_slots[y[!y %in% pos_free] + 1L])
  dummies_red <- lapply(dummies, function(y) y - num_free)
  transv_red <- get_transversals(sbase_red, sgens_red)
  g2 <- double_coset_can_rep(dummies_red, msym, sbase_red, sgens_red, transv_red, g1_red)
  if (identical(g2, 0)) return(0)
  # lift to the case with the free indices
  .tc_lift_sgens(size, pos_free, free, g2)
}

# canonicalize_naive (sympy/combinatorics/testutil.py): the minimum of the
# whole double coset D*g*S, used when no minimal BSGS is found
.tc_canonicalize_naive <- function(g, dummies, sym, ...) {
  v1 <- lapply(list(...), function(t) list(as.integer(t[[1L]]), lapply(t[[2L]], .tc_af),
                                           rep(list(integer()), t[[3L]]), t[[4L]]))
  gp <- .tc_gens_products(v1)
  size <- gp$size
  if (!is.list(sym) && length(sym) > 1L) sym <- as.list(sym)
  if (!is.list(sym)) {
    .tc_dummy_sgs(dummies, sym, size - 2L)
    if (is.null(sym)) stop(.sympy_condition("TypeError", "object of type 'NoneType' has no len()"))
    dummies <- list(dummies)
    sym <- list(sym)
  }
  dgens <- list()
  for (i in seq_along(sym)) dgens <- c(dgens, .tc_dummy_sgs(dummies[[i]], sym[[i]], size - 2L))
  # SymPy's D is then the group of the empty permutation and indexing fails
  if (!length(dgens)) stop(.sympy_condition("IndexError", "list index out of range"))
  dlist <- .tc_group_elements(.tc_group_gens(dgens))
  g <- .tc_af(g)
  st <- list()
  for (s in .tc_group_elements(.tc_group_gens(gp$gens))) {
    h <- g[s + 1L]
    for (d in dlist) st[[length(st) + 1L]] <- d[h + 1L]
  }
  st <- st[!duplicated(vapply(st, .tc_key, ""))]
  a <- st[.tc_lex_order(st)]
  prev <- integer(size)
  body <- seq_len(size - 2L)
  for (h in a) {
    if (identical(h[body], prev[body]) && h[size] != prev[size]) return(0)
    prev <- h
  }
  a[[1L]]
}

#' Direct products of permutation generators
#'
#' Port of SymPy's `perm_af_direct_product()`. With `signed = TRUE` the last
#' two points of each array form carry the sign.
#'
#' @param gens1,gens2 Lists of array forms.
#' @param signed Whether the permutations are signed.
#' @return A list of array forms.
#' @examples
#' perm_af_direct_product(list(c(1, 0, 2, 3), c(0, 1, 3, 2)), list(c(1, 0)), FALSE)
#' @export
perm_af_direct_product <- function(gens1, gens2, signed = TRUE) {
  gens1 <- lapply(gens1, .tc_af)
  gens2 <- lapply(gens2, .tc_af)
  signed <- isTRUE(as.logical(signed))
  s <- if (signed) 2L else 0L
  n1 <- length(gens1[[1L]]) - s
  n2 <- length(gens2[[1L]]) - s
  start <- seq_len(n1) - 1L
  end <- n1 + seq_len(n2) - 1L
  if (signed) {
    gens1 <- lapply(gens1, function(gen) {
      m <- length(gen)
      c(gen[seq_len(m - 2L)], end, gen[m - 1L] + n2, gen[m] + n2)
    })
  } else {
    gens1 <- lapply(gens1, function(gen) c(gen, end))
  }
  gens2 <- lapply(gens2, function(gen) c(start, gen + n1))
  c(gens1, gens2)
}

.tc_bsgs_direct_product <- function(base1, gens1, base2, gens2, signed = TRUE) {
  s <- if (isTRUE(as.logical(signed))) 2L else 0L
  n1 <- length(gens1[[1L]]) - s
  base <- c(as.integer(base1), as.integer(base2) + n1)
  gens <- perm_af_direct_product(gens1, gens2, signed)
  id_af <- seq_len(length(gens[[1L]])) - 1L
  gens <- Filter(function(h) !identical(h, id_af), gens)
  if (!length(gens)) gens <- list(id_af)
  list(base = base, gens = gens)
}

#' BSGS of tensor symmetries
#'
#' Ports of SymPy's BSGS constructors for tensor slot symmetries.
#' `get_symmetric_group_sgs(n, antisym)` is the minimal BSGS of a
#' (anti)symmetric tensor of rank `n`; `bsgs_direct_product()` is the direct
#' product of two BSGS; `get_transversals()` gives the basic transversals of
#' a BSGS; `get_minimal_bsgs()` recomputes a BSGS with minimal base (or
#' `NULL`); `tensor_gens()` is the BSGS of `n` tensors of the same type, with
#' fixed slots `list_free_indices` and commutation symmetry `sym`;
#' `gens_products()` combines tensors of different types, each given as
#' `list(base_i, gens_i, free_i, sym_i)`. `riemann_bsgs` is the BSGS of the
#' Riemann tensor.
#'
#' @param n Rank of the tensor.
#' @param antisym `FALSE` for a symmetric, `TRUE` for an antisymmetric tensor.
#' @param base,base1,base2 Bases (0 based integer vectors).
#' @param gens,gens1,gens2 Strong generators (lists of permutations).
#' @param signed Whether the permutations are signed.
#' @param list_free_indices List of the slots occupied by fixed indices, one
#'   integer vector per tensor.
#' @param sym Commutation symmetry: `NULL`, `0` or `1`.
#' @param ... Tensor types, each `list(base_i, gens_i, free_i, sym_i)`.
#' @return `get_symmetric_group_sgs()`, `bsgs_direct_product()` and
#'   `get_minimal_bsgs()` return `list(base, gens)`; `tensor_gens()` and
#'   `gens_products()` return `list(size, base, gens)`; `get_transversals()`
#'   returns a list of named lists (point to array form).
#' @examples
#' get_symmetric_group_sgs(3)
#' b1 <- get_symmetric_group_sgs(1); b2 <- get_symmetric_group_sgs(2)
#' bsgs_direct_product(b1[[1]], b1[[2]], b2[[1]], b2[[2]])
#' get_transversals(riemann_bsgs[[1]], riemann_bsgs[[2]])
#' b3 <- get_symmetric_group_sgs(3)
#' tensor_gens(b3[[1]], b3[[2]], list(1, 0))
#' gens_products(list(b2[[1]], b2[[2]], list(integer(), integer()), 0))
#' get_minimal_bsgs(c(2, 0), list(Permutation(c(1, 0, 2, 3, 5, 4)), Permutation(c(2, 3, 0, 1, 4, 5))))
#' @name tensor-bsgs
NULL

#' @rdname tensor-bsgs
#' @export
bsgs_direct_product <- function(base1, gens1, base2, gens2, signed = TRUE) {
  r <- .tc_bsgs_direct_product(base1, lapply(gens1, .tc_af), base2, lapply(gens2, .tc_af), signed)
  list(r$base, lapply(r$gens, .tc_af_new))
}

#' @rdname tensor-bsgs
#' @export
get_symmetric_group_sgs <- function(n, antisym = FALSE) {
  n <- as.integer(n)
  if (n %in% c(0L, 1L)) return(list(integer(), list(.tc_af_new(seq_len(n + 2L) - 1L))))
  gens <- lapply(seq_len(n - 1L) - 1L, function(i) {
    a <- seq_len(n) - 1L
    a[c(i + 1L, i + 2L)] <- c(i + 1L, i)
    c(a, if (isTRUE(antisym == 0)) c(n, n + 1L) else c(n + 1L, n))
  })
  list(seq_len(n - 1L) - 1L, lapply(gens, .tc_af_new))
}

#' @rdname tensor-bsgs
#' @format `riemann_bsgs`: `list(base, gens)` for the Riemann tensor.
#' @export
riemann_bsgs <- list(c(0L, 2L), list(.tc_af_new(c(1L, 0L, 2L, 3L, 5L, 4L)),
                                     .tc_af_new(c(0L, 1L, 3L, 2L, 5L, 4L)),
                                     .tc_af_new(c(2L, 3L, 0L, 1L, 4L, 5L))))

.tc_get_transversals <- function(base, gens) {
  if (!length(base)) return(list())
  stabs <- .tc_distribute_gens_by_base(base, gens)
  .tc_orbits_transversals_from_bsgs(base, stabs)
}

#' @rdname tensor-bsgs
#' @export
get_transversals <- function(base, gens) {
  trs <- .tc_get_transversals(as.integer(base), lapply(gens, .tc_af))
  lapply(trs, function(t) stats::setNames(t$p, t$x))
}

.tc_is_minimal_bsgs <- function(base, gens) {
  base1 <- integer()
  sgs1 <- gens
  size <- length(gens[[1L]])
  for (i in seq_len(size) - 1L) {
    if (!all(vapply(sgs1, function(h) h[i + 1L] == i, NA))) {
      base1 <- c(base1, i)
      sgs1 <- Filter(function(h) h[i + 1L] == i, sgs1)
    }
  }
  identical(base1, as.integer(base))
}

.tc_get_minimal_bsgs <- function(gens) {
  r <- .tc_schreier_sims_incremental(.tc_group_gens(gens))
  if (!.tc_is_minimal_bsgs(r$base, r$gens)) return(NULL)
  r
}

#' @rdname tensor-bsgs
#' @export
get_minimal_bsgs <- function(base, gens) {
  r <- .tc_get_minimal_bsgs(lapply(gens, .tc_af))
  if (is.null(r)) return(NULL)
  list(r$base, lapply(r$gens, .tc_af_new))
}

.tc_tensor_gens <- function(base, gens, list_free_indices, sym = 0) {
  lfi <- lapply(list_free_indices, as.integer)
  base <- as.integer(base)
  # if not base there is no slot symmetry for the component tensors;
  # with fewer than two tensors without free indices there is no
  # commutation symmetry, so there is no resulting slot symmetry
  if (!length(base) && sum(lengths(lfi) == 0L) < 2L) {
    size <- length(lfi) * (length(gens[[1L]]) - 2L) + 2L
    return(list(size = size, base = integer(), gens = list(seq_len(size) - 1L)))
  }
  G <- if (any(lengths(lfi) > 0L)) .tc_group_gens(gens) else NULL
  get_bsgs <- function(free_indices) {
    if (!length(free_indices)) return(list(base = base, gens = gens))
    H <- .tc_pointwise_stabilizer(G, free_indices)
    .tc_schreier_sims_incremental(H)
  }
  no_free <- list()
  size <- length(gens[[1L]])
  id_af <- seq_len(size) - 1L
  num_indices <- size - 2L
  if (!length(lfi[[1L]])) no_free[[1L]] <- seq_len(num_indices) - 1L
  r <- get_bsgs(lfi[[1L]])
  res_base <- r$base
  res_gens <- r$gens
  for (i in seq_along(lfi)[-1L]) {
    r1 <- get_bsgs(lfi[[i]])
    dp <- .tc_bsgs_direct_product(res_base, res_gens, r1$base, r1$gens, TRUE)
    res_base <- dp$base
    res_gens <- dp$gens
    if (!length(lfi[[i]])) no_free[[length(no_free) + 1L]] <- size - 2L + seq_len(num_indices) - 1L
    size <- size + num_indices
  }
  nr <- size - 2L
  # SymPy compares with the identity of a single tensor here
  res_gens <- Filter(function(h) !identical(h, id_af), res_gens)
  # if sym there are no commuting tensors stop here
  if (is.null(sym) || !length(no_free)) {
    if (!length(res_gens)) res_gens <- list(id_af)
    return(list(size = size, base = res_base, gens = res_gens))
  }
  # the slot symmetry group is S = P*C, with P the direct product of the
  # component symmetries and C the group (anti)commuting the component
  # tensors with no free indices; its base is the ordered union of the bases
  base_comm <- integer()
  for (i in seq_len(length(no_free) - 1L)) {
    ind1 <- no_free[[i]]
    ind2 <- no_free[[i + 1L]]
    a <- c(seq_len(ind1[1L]) - 1L, ind2, ind1)
    base_comm <- c(base_comm, ind1[1L])
    last <- ind2[length(ind2)]
    if (last + 1L < nr) a <- c(a, (last + 1L):(nr - 1L))
    a <- c(a, if (isTRUE(sym == 0)) c(nr, nr + 1L) else c(nr + 1L, nr))
    res_gens[[length(res_gens) + 1L]] <- a
  }
  for (i in base_comm) if (!(i %in% res_base)) res_base <- c(res_base, i)
  res_base <- sort(res_base)
  if (!length(res_gens)) res_gens <- list(id_af)
  list(size = size, base = res_base, gens = res_gens)
}

#' @rdname tensor-bsgs
#' @export
tensor_gens <- function(base, gens, list_free_indices, sym = 0) {
  r <- .tc_tensor_gens(base, lapply(gens, .tc_af), list_free_indices, sym)
  list(as.numeric(r$size), r$base, lapply(r$gens, .tc_af_new))
}

.tc_gens_products <- function(v) {
  tg <- function(t) .tc_tensor_gens(t[[1L]], lapply(t[[2L]], .tc_af), t[[3L]], t[[4L]])
  r <- tg(v[[1L]])
  res_base <- r$base
  res_gens <- r$gens
  for (i in seq_along(v)[-1L]) {
    r <- tg(v[[i]])
    dp <- .tc_bsgs_direct_product(res_base, res_gens, r$base, r$gens, TRUE)
    res_base <- dp$base
    res_gens <- dp$gens
  }
  # SymPy then filters `h != id_af` comparing a Permutation with a list,
  # which removes nothing
  list(size = length(res_gens[[1L]]), base = res_base, gens = res_gens)
}

#' @rdname tensor-bsgs
#' @export
gens_products <- function(...) {
  r <- .tc_gens_products(list(...))
  list(as.numeric(r$size), r$base, lapply(r$gens, .tc_af_new))
}
