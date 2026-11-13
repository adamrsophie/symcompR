# Sparse domain matrices, after sympy/polys/matrices/sdm.py. An SDM is
# list(d = rows, shape, domain) where the rows are a named list keyed by the
# 0-based row index (as a string) in insertion order, and each row a named
# list of its nonzero elements keyed by the 0-based column index, also in
# insertion order, as Python's dict of dicts.
# ponytail: Python iterates small integer sets in ascending order; sets here
# are iterated sorted, which differs from CPython only for sparse rows with
# column indices above the hash table size.

.sdm_new <- function(d, shape, K) {
  if (is.null(d)) d <- list()
  if (length(d) && is.null(names(d))) names(d) <- character(length(d))
  structure(list(d = d, shape = as.integer(shape), domain = K), class = "SDM")
}
.sk <- function(i) as.character(i)
.ski <- function(x) as.integer(names(x))
.sset <- function(x) as.character(sort(as.integer(x)))

.sdm_check_bounds <- function(d, m, n) {
  rk <- as.integer(names(d))
  if (any(rk < 0L | rk >= m)) stop(.dm_err("DMBadInputError", "Row out of range"))
  ck <- unlist(lapply(d, function(r) as.integer(names(r))))
  if (any(ck < 0L | ck >= n)) stop(.dm_err("DMBadInputError", "Column out of range"))
}

#' @rdname DDM
#' @export
SDM <- function(elemsdict, shape, domain) {
  shape <- .dm_shape(shape)
  d <- .dm_as_dod(elemsdict)
  .sdm_check_bounds(d, shape[1L], shape[2L])
  .sdm_new(d, shape, domain)
}
# a dict of dicts from a nested named list or a list of (key, value) pairs
.dm_as_dod <- function(x) {
  pairs <- function(z) {
    if (!length(z)) return(list())
    if (!is.null(names(z)) && all(nzchar(names(z)))) return(z)
    out <- list()
    for (p in z) out[[.sk(as.integer(as.character(p[[1L]])))]] <- p[[2L]]
    out
  }
  d <- pairs(x)
  lapply(d, pairs)
}

# binop_dict and unop_dict
.sdm_binop <- function(A, B, fab, fa, fb) {
  Anz <- names(A); Bnz <- names(B)
  C <- list()
  for (i in .sset(intersect(Anz, Bnz))) {
    Ai <- A[[i]]; Bi <- B[[i]]
    Ci <- list()
    Anzi <- names(Ai); Bnzi <- names(Bi)
    for (j in .sset(intersect(Anzi, Bnzi))) { Cij <- fab(Ai[[j]], Bi[[j]]); if (.de_nz(Cij)) Ci[[j]] <- Cij }
    for (j in .sset(setdiff(Anzi, Bnzi))) { Cij <- fa(Ai[[j]]); if (.de_nz(Cij)) Ci[[j]] <- Cij }
    for (j in .sset(setdiff(Bnzi, Anzi))) { Cij <- fb(Bi[[j]]); if (.de_nz(Cij)) Ci[[j]] <- Cij }
    if (length(Ci)) C[[i]] <- Ci
  }
  for (i in .sset(setdiff(Anz, Bnz))) {
    Ci <- list()
    for (j in names(A[[i]])) { Cij <- fa(A[[i]][[j]]); if (.de_nz(Cij)) Ci[[j]] <- Cij }
    if (length(Ci)) C[[i]] <- Ci
  }
  for (i in .sset(setdiff(Bnz, Anz))) {
    Ci <- list()
    for (j in names(B[[i]])) { Cij <- fb(B[[i]][[j]]); if (.de_nz(Cij)) Ci[[j]] <- Cij }
    if (length(Ci)) C[[i]] <- Ci
  }
  C
}
.sdm_unop <- function(A, f) {
  B <- list()
  for (i in names(A)) {
    Bi <- list()
    for (j in names(A[[i]])) { Bij <- f(A[[i]][[j]]); if (.de_nz(Bij)) Bi[[j]] <- Bij }
    if (length(Bi)) B[[i]] <- Bi
  }
  B
}
.sdm_transpose <- function(M) {
  MT <- list()
  for (i in names(M)) for (j in names(M[[i]])) {
    if (is.null(MT[[j]])) MT[[j]] <- list()
    MT[[j]][[i]] <- M[[i]][[j]]
  }
  MT
}
.sdm_dotvec <- function(A, B, K) {
  s <- .de_zero(K)
  for (j in .sset(intersect(names(A), names(B)))) s <- s + A[[j]] * B[[j]]
  s
}
.sdm_matvecmul <- function(A, B, K) {
  C <- list()
  for (i in names(A)) { Ci <- .sdm_dotvec(A[[i]], B, K); if (.de_nz(Ci)) C[[i]] <- Ci }
  C
}
.sdm_matmul <- function(A, B, K) {
  C <- list()
  for (i in names(A)) {
    Ai <- A[[i]]; Ci <- list()
    for (k in .sset(intersect(names(Ai), names(B)))) {
      Aik <- Ai[[k]]
      for (j in names(B[[k]])) {
        Cij <- Ci[[j]]
        if (!is.null(Cij)) {
          Cij <- Cij + Aik * B[[k]][[j]]
          if (.de_nz(Cij)) Ci[[j]] <- Cij else Ci[[j]] <- NULL
        } else {
          Cij <- Aik * B[[k]][[j]]
          if (.de_nz(Cij)) Ci[[j]] <- Cij
        }
      }
    }
    if (length(Ci)) C[[i]] <- Ci
  }
  C
}
.de_recip <- function(a) if (inherits(a, "bigz") || inherits(a, "bigq")) 1 / gmp::as.bigq(a) else a^-1L

# sdm_irref: list(rref, pivots, nonzero_columns)
.sdm_irref <- function(A) {
  rows <- unname(A)
  mins <- vapply(rows, function(r) min(.ski(r)), 0L)
  Arows <- rows[order(mins, method = "radix")]
  pmap <- list(); reduced <- integer(); nonreduced <- integer(); nzc <- list()
  while (length(Arows)) {
    Ai <- Arows[[length(Arows)]]; Arows[[length(Arows)]] <- NULL
    Ai <- Ai[!(.ski(Ai) %in% reduced)]
    for (j in sort(intersect(nonreduced, .ski(Ai)))) {
      jk <- .sk(j); Aj <- pmap[[jk]]; Aij <- Ai[[jk]]
      Ainz <- names(Ai); Ajnz <- names(Aj)
      for (k in .sset(setdiff(Ajnz, Ainz))) Ai[[k]] <- -Aij * Aj[[k]]
      Ai[[jk]] <- NULL
      Ainz <- setdiff(Ainz, jk)
      for (k in .sset(intersect(Ajnz, Ainz))) {
        Aik <- Ai[[k]] - Aij * Aj[[k]]
        if (.de_nz(Aik)) Ai[[k]] <- Aik else Ai[[k]] <- NULL
      }
    }
    if (!length(Ai)) next
    j <- min(.ski(Ai)); jk <- .sk(j)
    Aijinv <- .de_recip(Ai[[jk]])
    for (l in names(Ai)) Ai[[l]] <- Ai[[l]] * Aijinv
    pmap[[jk]] <- Ai
    Ainz <- names(Ai)
    ks <- nzc[[jk]]; nzc[[jk]] <- NULL
    for (k in sort(ks)) {
      kk <- .sk(k); Ak <- pmap[[kk]]; Akj <- Ak[[jk]]; Aknz <- names(Ak)
      for (l in .sset(setdiff(Ainz, Aknz))) { Ak[[l]] <- -Akj * Ai[[l]]; nzc[[l]] <- union(nzc[[l]], k) }
      Ak[[jk]] <- NULL
      Aknz <- setdiff(Aknz, jk)
      for (l in .sset(intersect(Ainz, Aknz))) {
        Akl <- Ak[[l]] - Akj * Ai[[l]]
        if (.de_nz(Akl)) Ak[[l]] <- Akl
        else { Ak[[l]] <- NULL; if (l != jk) nzc[[l]] <- setdiff(nzc[[l]], k) }
      }
      pmap[[kk]] <- Ak
      if (length(Ak) == 1L) { reduced <- c(reduced, k); nonreduced <- setdiff(nonreduced, k) }
    }
    if (length(Ai) == 1L) reduced <- c(reduced, j)
    else {
      nonreduced <- c(nonreduced, j)
      for (l in names(Ai)) if (l != jk) nzc[[l]] <- union(nzc[[l]], j)
    }
  }
  pivots <- sort(c(reduced, nonreduced))
  nzc2 <- lapply(nzc, function(s) match(s, pivots) - 1L)
  rref <- stats::setNames(lapply(pivots, function(p) pmap[[.sk(p)]]), .sk(seq_along(pivots) - 1L))
  list(rref = rref, pivots = pivots, nzc = nzc2)
}

# sdm_rref_den: list(rref, den, pivots)
.sdm_rref_den <- function(A, K) {
  if (!length(A)) return(list(rref = list(), den = .de_one(K), pivots = integer()))
  if (length(A) == 1L) {
    Ai <- A[[1L]]; j <- min(.ski(Ai))
    return(list(rref = list(`0` = Ai), den = Ai[[.sk(j)]], pivots = j))
  }
  exquo <- function(a, b) if (.de_kind(K) %in% c("RR", "CC")) .de_quo(K, a, b) else .de_exquo(K, a, b)
  rows_in_order <- unname(A[.sset(names(A))])
  ctr_red <- list(); ctr_unred <- list()
  rr <- list(); denom <- NULL; divisor <- NULL
  mins <- vapply(rows_in_order, function(r) min(.ski(r)), 0L)
  for (Ai in rows_in_order[order(mins, method = "radix")]) {
    Ai <- Ai[!(names(Ai) %in% names(ctr_red))]
    Ai_cancel <- list()
    for (j in .sset(intersect(names(ctr_unred), names(Ai)))) {
      Aij <- Ai[[j]]; Ai[[j]] <- NULL
      Aj <- rr[[ctr_unred[[j]] + 1L]]
      for (k in names(Aj)) {
        Aik_c <- Ai_cancel[[k]]
        if (is.null(Aik_c)) Ai_cancel[[k]] <- Aij * Aj[[k]]
        else {
          Aik_c <- Aik_c + Aij * Aj[[k]]
          if (.de_nz(Aik_c)) Ai_cancel[[k]] <- Aik_c else Ai_cancel[[k]] <- NULL
        }
      }
    }
    Ai_nz <- names(Ai); Ac_nz <- names(Ai_cancel)
    d <- if (is.null(denom)) .de_one(K) else denom
    for (k in .sset(setdiff(Ac_nz, Ai_nz))) Ai[[k]] <- -Ai_cancel[[k]]
    for (k in .sset(setdiff(Ai_nz, Ac_nz))) Ai[[k]] <- Ai[[k]] * d
    for (k in .sset(intersect(Ac_nz, Ai_nz))) {
      Aik <- Ai[[k]] * d - Ai_cancel[[k]]
      if (.de_nz(Aik)) Ai[[k]] <- Aik else Ai[[k]] <- NULL
    }
    if (!length(Ai)) next
    j <- .sk(min(.ski(Ai)))
    Aij <- Ai[[j]]; Ai[[j]] <- NULL
    for (pk in names(ctr_unred)) {
      k <- ctr_unred[[pk]] + 1L
      Ak <- rr[[k]]
      if (is.null(Ak[[j]])) {
        for (l in names(Ak)) {
          Akl <- Ak[[l]] * Aij
          if (!is.null(divisor)) Akl <- exquo(Akl, divisor)
          Ak[[l]] <- Akl
        }
        rr[[k]] <- Ak
        next
      }
      Akj <- Ak[[j]]; Ak[[j]] <- NULL
      Ai_nz <- names(Ai); Ak_nz <- names(Ak)
      for (l in .sset(setdiff(Ai_nz, Ak_nz))) {
        Ak[[l]] <- -Akj * Ai[[l]]
        if (!is.null(divisor)) Ak[[l]] <- exquo(Ak[[l]], divisor)
      }
      for (l in .sset(setdiff(Ak_nz, Ai_nz))) {
        Ak[[l]] <- Aij * Ak[[l]]
        if (!is.null(divisor)) Ak[[l]] <- exquo(Ak[[l]], divisor)
      }
      for (l in .sset(intersect(Ai_nz, Ak_nz))) {
        Akl <- Aij * Ak[[l]] - Akj * Ai[[l]]
        if (.de_nz(Akl)) {
          if (!is.null(divisor)) Akl <- exquo(Akl, divisor)
          Ak[[l]] <- Akl
        } else Ak[[l]] <- NULL
      }
      rr[[k]] <- if (length(Ak)) Ak else list()
      if (!length(Ak)) { ctr_unred[[pk]] <- NULL; ctr_red[[pk]] <- k - 1L }
    }
    i <- length(rr)
    rr[[i + 1L]] <- if (length(Ai)) Ai else list()
    if (length(Ai)) ctr_unred[[j]] <- i else ctr_red[[j]] <- i
    if (!.de_is_one(K, Aij)) denom <- if (is.null(denom)) Aij else denom * Aij
    if (!is.null(divisor)) denom <- exquo(denom, divisor)
    divisor <- denom
  }
  if (is.null(denom)) denom <- .de_one(K)
  ctr <- c(ctr_red, ctr_unred)
  row_to_col <- integer(length(rr))
  for (cc in names(ctr)) row_to_col[ctr[[cc]] + 1L] <- as.integer(cc)
  o <- order(row_to_col)
  pivots <- row_to_col[o]
  out <- list()
  for (i in seq_along(o)) {
    Ai <- rr[[o[i]]]
    Ai[[.sk(pivots[i])]] <- denom
    out[[.sk(i - 1L)]] <- Ai
  }
  list(rref = out, den = denom, pivots = pivots)
}

.sdm_nullspace_from_rref_core <- function(A, one, ncols, pivots, nzc) {
  nonpivots <- setdiff(seq_len(ncols) - 1L, pivots)
  Kb <- list()
  for (j in nonpivots) {
    Kj <- list(); Kj[[.sk(j)]] <- one
    for (i in sort(nzc[[.sk(j)]])) Kj[[.sk(pivots[i + 1L])]] <- -A[[.sk(i)]][[.sk(j)]]
    Kb[[length(Kb) + 1L]] <- Kj
  }
  list(K = Kb, nonpivots = nonpivots)
}
.sdm_particular_from_rref <- function(A, ncols, pivots) {
  P <- list()
  for (i in seq_along(pivots)) {
    j <- pivots[i]; Ai <- A[[.sk(i - 1L)]]
    Ain <- Ai[[.sk(ncols - 1L)]]
    if (!is.null(Ain)) P[[.sk(j)]] <- Ain / Ai[[.sk(j)]]
  }
  P
}

# sdm_berk: the characteristic polynomial as a dict of coefficients
.sdm_berk <- function(M, n, K) {
  zero <- .de_zero(K); one <- .de_one(K)
  if (n == 0L) return(list(`0` = one))
  a <- zero; R <- list(); C <- list(); A <- list()
  for (i in names(M)) for (j in names(M[[i]])) {
    ii <- as.integer(i); jj <- as.integer(j); Mij <- M[[i]][[j]]
    if (ii && jj) { k <- .sk(ii - 1L); if (is.null(A[[k]])) A[[k]] <- list(); A[[k]][[.sk(jj - 1L)]] <- Mij }
    else if (ii) C[[.sk(ii - 1L)]] <- Mij
    else if (jj) R[[.sk(jj - 1L)]] <- Mij
    else a <- Mij
  }
  AnC <- C
  RC <- .sdm_dotvec(R, C, K)
  Tvals <- list(one, -a, -RC)
  if (n >= 3L) for (i in 3:n) {
    AnC <- .sdm_matvecmul(A, AnC, K)
    if (!length(AnC)) break
    Tvals[[length(Tvals) + 1L]] <- -.sdm_dotvec(R, AnC, K)
  }
  while (length(Tvals) && .de_is_zero(Tvals[[length(Tvals)]])) Tvals[[length(Tvals)]] <- NULL
  q <- .sdm_berk(A, n - 1L, K)
  Tvals <- rev(Tvals)
  Tq <- list()
  qk <- .ski(q)
  lo <- min(qk); hi <- min(max(qk) + length(Tvals), n + 1L)
  if (lo < hi) for (i in lo:(hi - 1L)) {
    Ti <- stats::setNames(Tvals, .sk(seq_along(Tvals) - 1L + i - length(Tvals) + 1L))
    Tqi <- .sdm_dotvec(Ti, q, K)
    if (.de_nz(Tqi)) Tq[[.sk(i)]] <- Tqi
  }
  Tq
}

# _strongly_connected_components (non recursive Tarjan)
.py_scc <- function(V, Gmap) {
  lowlink <- list(); indices <- list(); stack <- character(); callstack <- list(); components <- list()
  start <- function(v) {
    idx <- length(stack)
    indices[[v]] <<- idx; lowlink[[v]] <<- idx
    stack <<- c(stack, v)
    callstack[[length(callstack) + 1L]] <<- list(v = v, nb = Gmap[[v]], pos = 0L)
  }
  finish <- function(v1) {
    if (lowlink[[v1]] == indices[[v1]]) {
      comp <- character()
      repeat { top <- stack[length(stack)]; stack <<- stack[-length(stack)]; comp <- c(comp, top); if (top == v1) break }
      components[[length(components) + 1L]] <<- rev(comp)
    }
    v2 <- callstack[[length(callstack)]]$v
    callstack[[length(callstack)]] <<- NULL
    if (length(callstack)) {
      w <- callstack[[length(callstack)]]$v
      lowlink[[w]] <<- min(lowlink[[w]], lowlink[[v2]])
    }
  }
  for (v in V) {
    if (!is.null(indices[[v]])) next
    start(v)
    while (length(callstack)) {
      top <- callstack[[length(callstack)]]
      if (top$pos >= length(top$nb)) { finish(top$v); next }
      callstack[[length(callstack)]]$pos <- top$pos + 1L
      v2 <- top$nb[[top$pos + 1L]]
      if (is.null(indices[[v2]])) start(v2)
      else if (v2 %in% stack) lowlink[[top$v]] <- min(lowlink[[top$v]], indices[[v2]])
    }
  }
  lapply(components, as.integer)
}

.sdm_to_list <- function(a) {
  s <- .subset2(a, "shape"); K <- .subset2(a, "domain"); z <- .de_zero(K)
  rows <- rep(list(rep(list(z), s[2L])), s[1L])
  d <- .subset2(a, "d")
  for (i in names(d)) for (j in names(d[[i]])) rows[[as.integer(i) + 1L]][[as.integer(j) + 1L]] <- d[[i]][[j]]
  rows
}
.sdm_to_ddm <- function(a) .ddm_new(.sdm_to_list(a), .subset2(a, "shape"), .subset2(a, "domain"))
.sdm_from_list <- function(rows, shape, K) {
  d <- list()
  for (i in seq_len(shape[1L])) {
    nz <- which(vapply(rows[[i]], .de_nz, NA))
    if (length(nz)) d[[.sk(i - 1L)]] <- stats::setNames(rows[[i]][nz], .sk(nz - 1L))
  }
  .sdm_new(d, shape, K)
}
.sdm_to_dod_pd <- function(a) {
  d <- .subset2(a, "d")
  .pydict(lapply(names(d), as.integer), lapply(d, function(r) .pydict(lapply(names(r), as.integer), unname(r))))
}
.sdm_dok_pd <- function(a) {
  d <- .subset2(a, "d"); keys <- list(); vals <- list()
  for (i in names(d)) for (j in names(d[[i]])) {
    keys[[length(keys) + 1L]] <- .pytuple(list(as.integer(i), as.integer(j)))
    vals[[length(vals) + 1L]] <- d[[i]][[j]]
  }
  .pydict(keys, vals)
}
# a Python range()[slice] as 0-based indices; s is list(start, stop, step)
.py_slice_idx <- function(sl, n) {
  if (is.numeric(sl) && !is.list(sl)) return(as.integer(sl))
  st <- sl[[3L]]; if (is.null(st)) st <- 1L
  norm <- function(v, dflt) { if (is.null(v)) return(dflt); v <- as.integer(v); if (v < 0L) v <- max(v + n, if (st > 0) 0L else -1L); if (v > n) v <- if (st > 0) n else n - 1L; v }
  a <- norm(sl[[1L]], if (st > 0) 0L else n - 1L)
  b <- norm(sl[[2L]], if (st > 0) n else -1L)
  if (st > 0) { if (a >= b) integer() else seq.int(a, b - 1L, by = st) }
  else { if (a <= b) integer() else seq.int(a, b + 1L, by = st) }
}

.sdm_extract <- function(a, rws, cls) {
  s <- .subset2(a, "shape"); K <- .subset2(a, "domain"); d <- .subset2(a, "d")
  m <- s[1L]; n <- s[2L]
  if (!(length(d) && length(rws) && length(cls))) return(.sdm_new(list(), c(length(rws), length(cls)), K))
  if (!(-m <= min(rws) && max(rws) < m)) stop(.py_err("IndexError", "Row index out of range"))
  if (!(-n <= min(cls) && max(cls) < n)) stop(.py_err("IndexError", "Column index out of range"))
  rowmap <- list(); colmap <- list()
  for (i2 in seq_along(rws)) { k <- .sk(rws[i2] %% m); rowmap[[k]] <- c(rowmap[[k]], i2 - 1L) }
  for (j2 in seq_along(cls)) { k <- .sk(cls[j2] %% n); colmap[[k]] <- c(colmap[[k]], j2 - 1L) }
  out <- list()
  for (i1 in .sset(intersect(names(rowmap), names(d)))) {
    row1 <- d[[i1]]; row2 <- list()
    for (j1 in .sset(intersect(names(colmap), names(row1)))) for (j2 in colmap[[j1]]) row2[[.sk(j2)]] <- row1[[j1]]
    if (length(row2)) for (i2 in rowmap[[i1]]) out[[.sk(i2)]] <- row2
  }
  .sdm_new(out, c(length(rws), length(cls)), K)
}
.sdm_extract_slice <- function(a, ri, ci) {
  d <- .subset2(a, "d"); out <- list()
  for (i in names(d)) {
    ii <- as.integer(i)
    if (ii %in% ri) {
      row <- list()
      for (j in names(d[[i]])) if (as.integer(j) %in% ci) row[[.sk(match(as.integer(j), ci) - 1L)]] <- d[[i]][[j]]
      if (length(row)) out[[.sk(match(ii, ri) - 1L)]] <- row
    }
  }
  .sdm_new(out, c(length(ri), length(ci)), .subset2(a, "domain"))
}

.sdm_method <- function(a, name) {
  d <- .subset2(a, "d"); s <- .subset2(a, "shape"); K <- .subset2(a, "domain")
  m <- s[1L]; n <- s[2L]
  new <- function(dd, shape = s, KK = K) .sdm_new(dd, shape, KK)
  switch(name,
    shape = .pytuple(as.list(s)), rows = m, cols = n, domain = K, fmt = "sparse", is_DDM = FALSE, is_DFM = FALSE,
    getitem = function(i, j) {
      i <- as.integer(as.character(i)); j <- as.integer(as.character(j))
      if (!(-m <= i && i < m && -n <= j && j < n)) stop(.py_err("IndexError", "index out of range"))
      v <- d[[.sk(i %% m)]][[.sk(j %% n)]]
      if (is.null(v)) .de_zero(K) else v
    },
    setitem = function(i, j, value) {
      i <- as.integer(as.character(i)); j <- as.integer(as.character(j))
      if (!(-m <= i && i < m && -n <= j && j < n)) stop(.py_err("IndexError", "index out of range"))
      ik <- .sk(i %% m); jk <- .sk(j %% n)
      if (.de_nz(value)) { if (is.null(d[[ik]])) d[[ik]] <- list(); d[[ik]][[jk]] <- value }
      else if (!is.null(d[[ik]])) { d[[ik]][[jk]] <- NULL; if (!length(d[[ik]])) d[[ik]] <- NULL }
      new(d)
    },
    extract = function(rws, cls) .sdm_extract(a, as.integer(unlist(rws)), as.integer(unlist(cls))),
    extract_slice = function(sl1, sl2) .sdm_extract_slice(a, .py_slice_idx(sl1, m), .py_slice_idx(sl2, n)),
    to_list = function() .sdm_to_list(a),
    to_list_flat = function() unlist(.sdm_to_list(a), recursive = FALSE),
    to_flat_nz = function() {
      dk <- .sdm_dok_pd(a)
      .pytuple(list(dk$values, .pytuple(list(.pytuple(dk$keys), .pytuple(as.list(s))))))
    },
    to_dod = function() .sdm_to_dod_pd(a),
    to_dok = function() .sdm_dok_pd(a),
    iter_values = function() unname(unlist(lapply(unname(d), unname), recursive = FALSE)),
    iter_items = function() { dk <- .sdm_dok_pd(a); Map(function(k, v) .pytuple(list(k, v)), dk$keys, dk$values) },
    to_ddm = function() .sdm_to_ddm(a), to_sdm = function() a,
    to_dfm = function() stop(.py_err("NotImplementedError", "DFM requires python-flint")),
    to_dfm_or_ddm = function() .sdm_to_ddm(a),
    copy = function() a,
    convert_to = function(K2) {
      if (.de_same(K2, K)) return(a)
      new(.sdm_unop(d, function(e) .de_convert(K2, e, K)), s, K2)
    },
    transpose = function() new(.sdm_transpose(d), rev(s)),
    add = function(b) new(.sdm_binop(d, .subset2(b, "d"), `+`, identity, identity)),
    sub = function(b) new(.sdm_binop(d, .subset2(b, "d"), `-`, identity, function(e) -e)),
    neg = function() new(.sdm_unop(d, function(e) -e)),
    mul = function(b) new(.sdm_unop(d, function(e) e * b)),
    rmul = function(b) new(.sdm_unop(d, function(e) b * e)),
    matmul = function(b) {
      if (!.de_same(K, .subset2(b, "domain"))) stop(.dm_err("DMDomainError"))
      bs <- .subset2(b, "shape")
      if (n != bs[1L]) stop(.dm_err("DMShapeError"))
      new(.sdm_matmul(d, .subset2(b, "d"), K), c(m, bs[2L]))
    },
    mul_elementwise = function(b) {
      if (!.de_same(K, .subset2(b, "domain"))) stop(.dm_err("DMDomainError"))
      if (!identical(s, .subset2(b, "shape"))) stop(.dm_err("DMShapeError"))
      z <- function(e) .de_zero(K)
      new(.sdm_binop(d, .subset2(b, "d"), `*`, z, z))
    },
    nnz = function() sum(vapply(d, length, 0L)),
    scc = function() {
      V <- .sk(seq_len(m) - 1L)
      Emap <- stats::setNames(lapply(V, function(v) if (is.null(d[[v]])) character() else names(d[[v]])), V)
      .py_scc(V, Emap)
    },
    rref = function() {
      r <- .sdm_irref(d)
      .pytuple(list(new(r$rref), as.list(r$pivots)))
    },
    rref_den = function() {
      r <- .sdm_rref_den(d, K)
      .pytuple(list(new(r$rref), r$den, as.list(r$pivots)))
    },
    inv = function() .ddm_to_sdm(.ddm_method(.sdm_to_ddm(a), "inv")()),
    det = function() .ddm_method(.sdm_to_ddm(a), "det")(),
    lu = function() { r <- .ddm_method(.sdm_to_ddm(a), "lu")(); .pytuple(list(.ddm_to_sdm(r[[1L]]), .ddm_to_sdm(r[[2L]]), r[[3L]])) },
    qr = function() { r <- .ddm_qr(.sdm_to_ddm(a)); .pytuple(list(.ddm_to_sdm(r[[1L]]), .ddm_to_sdm(r[[2L]]))) },
    lu_solve = function(b) .ddm_to_sdm(.ddm_method(.sdm_to_ddm(a), "lu_solve")(.sdm_to_ddm(b))),
    fflu = function() .pytuple(lapply(.ddm_fflu(.sdm_to_ddm(a)), .ddm_to_sdm)),
    nullspace = function() {
      r <- .sdm_irref(d)
      ns <- .sdm_nullspace_from_rref_core(r$rref, .de_one(K), n, r$pivots, r$nzc)
      Kd <- stats::setNames(ns$K, .sk(seq_along(ns$K) - 1L))
      .pytuple(list(new(Kd, c(length(ns$K), n)), as.list(ns$nonpivots)))
    },
    nullspace_from_rref = function(pivots = NULL) {
      pivots <- if (is.null(pivots)) sort(vapply(d, function(r) min(.ski(r)), 0L)) else as.integer(unlist(pivots))
      if (!length(pivots)) return(.pytuple(list(.sdm_eye(c(n, n), K), as.list(seq_len(n) - 1L))))
      if (length(pivots) == n) return(.pytuple(list(new(list(), c(0L, n)), list())))
      pivot_val <- d[["0"]][[.sk(pivots[1L])]]
      nzcols <- list()
      for (i in names(d)) for (j in names(d[[i]])) nzcols[[j]] <- c(nzcols[[j]], list(list(as.integer(i), d[[i]][[j]])))
      basis <- list(); nonpivots <- list()
      for (j in seq_len(n) - 1L) {
        if (j %in% pivots) next
        nonpivots[[length(nonpivots) + 1L]] <- j
        vec <- list(); vec[[.sk(j)]] <- pivot_val
        for (pr in nzcols[[.sk(j)]]) vec[[.sk(pivots[pr[[1L]] + 1L])]] <- -pr[[2L]]
        basis[[length(basis) + 1L]] <- vec
      }
      .pytuple(list(new(stats::setNames(basis, .sk(seq_along(basis) - 1L)), c(length(basis), n)), nonpivots))
    },
    particular = function() {
      r <- .sdm_irref(d)
      P <- .sdm_particular_from_rref(r$rref, n, r$pivots)
      new(if (length(P)) list(`0` = P) else list(), c(1L, n - 1L))
    },
    hstack = function(...) {
      Anew <- d; cols <- n
      for (Bk in list(...)) {
        for (i in names(.subset2(Bk, "d"))) {
          if (is.null(Anew[[i]])) Anew[[i]] <- list()
          Bki <- .subset2(Bk, "d")[[i]]
          for (j in names(Bki)) Anew[[i]][[.sk(as.integer(j) + cols)]] <- Bki[[j]]
        }
        cols <- cols + .subset2(Bk, "shape")[2L]
      }
      new(Anew, c(m, cols))
    },
    vstack = function(...) {
      Anew <- d; rws <- m
      for (Bk in list(...)) {
        for (i in names(.subset2(Bk, "d"))) Anew[[.sk(as.integer(i) + rws)]] <- .subset2(Bk, "d")[[i]]
        rws <- rws + .subset2(Bk, "shape")[1L]
      }
      new(Anew, c(rws, n))
    },
    applyfunc = function(func, domain) new(lapply(d, function(r) lapply(r, func)), s, domain),
    charpoly = function() {
      pd <- .sdm_berk(d, m, K)
      pl <- rep(list(.de_zero(K)), m + 1L)
      for (i in names(pd)) pl[[as.integer(i) + 1L]] <- pd[[i]]
      pl
    },
    is_zero_matrix = function() !length(d),
    is_upper = function() all(unlist(lapply(names(d), function(i) as.integer(i) <= .ski(d[[i]])))),
    is_lower = function() all(unlist(lapply(names(d), function(i) as.integer(i) >= .ski(d[[i]])))),
    is_diagonal = function() all(unlist(lapply(names(d), function(i) as.integer(i) == .ski(d[[i]])))),
    diagonal = function() {
      z <- .de_zero(K); out <- list()
      for (i in names(d)) if (as.integer(i) < n) out[[length(out) + 1L]] <- if (is.null(d[[i]][[i]])) z else d[[i]][[i]]
      out
    },
    lll = function(delta = gmp::as.bigq(3, 4)) .ddm_to_sdm(.ddm_lll(.sdm_to_ddm(a), delta)),
    lll_transform = function(delta = gmp::as.bigq(3, 4)) {
      r <- .ddm_lll_transform(.sdm_to_ddm(a), delta)
      .pytuple(list(.ddm_to_sdm(r[[1L]]), .ddm_to_sdm(r[[2L]])))
    },
    zeros = function(shape, domain) .sdm_new(list(), .dm_shape(shape), domain),
    eye = function(size, domain) .sdm_eye(if (is.list(size)) .dm_shape(size) else as.integer(size), domain),
    stop(.py_err("AttributeError", paste0("'SDM' object has no attribute '", name, "'"))))
}
.sdm_eye <- function(shape, K) {
  if (length(shape) == 1L) shape <- c(shape, shape)
  one <- .de_one(K); d <- list()
  for (i in seq_len(min(shape)) - 1L) { r <- list(); r[[.sk(i)]] <- one; d[[.sk(i)]] <- r }
  .sdm_new(d, shape, K)
}
.sdm_diag <- function(diagonal, K, shape) {
  d <- list()
  for (i in seq_along(diagonal)) if (.de_nz(diagonal[[i]])) { r <- list(); r[[.sk(i - 1L)]] <- diagonal[[i]]; d[[.sk(i - 1L)]] <- r }
  .sdm_new(d, shape, K)
}

#' @rdname DDM
#' @export
`$.SDM` <- function(x, name) .sdm_method(x, name)
#' @rdname DDM
#' @export
format.SDM <- function(x, ...) .dm_repr(x)
#' @rdname DDM
#' @export
print.SDM <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname DDM
#' @export
Ops.SDM <- function(e1, e2) .dm_ops(.Generic, e1, e2, "SDM")

# ---------------------------------------------------------------------------
# lll.py

.ddm_lll_core <- function(x, delta, return_transform) {
  if (gmp::as.bigq(1, 4) >= delta || delta >= 1) stop(.dm_err("DMValueError", "delta must lie in range (0.25, 1)"))
  s <- .subset2(x, "shape")
  if (s[1L] > s[2L]) stop(.dm_err("DMShapeError", "input matrix must have shape (m, n) with m <= n"))
  if (.de_kind(.subset2(x, "domain")) != "ZZ") stop(.dm_err("DMDomainError", "input matrix domain must be ZZ"))
  m <- s[1L]; n <- s[2L]
  y <- .subset2(x, "rows")
  q0 <- gmp::as.bigq(0); half <- gmp::as.bigq(1, 2)
  ystar <- rep(list(rep(list(q0), n)), m)
  mu <- rep(list(rep(list(q0), m)), m)
  gstar <- rep(list(q0), m)
  Tm <- if (return_transform) .ddm_eye(c(m, m), ZZ)[["rows"]] else NULL
  lin_err <- "input matrix contains linearly dependent rows"
  closest <- function(v) gmp::as.bigz(floor(v + half))
  lovasz <- function(k) gstar[[k]] >= (delta - mu[[k]][[k - 1L]]^2) * gstar[[k - 1L]]
  mu_small <- function(k, j) abs(mu[[k]][[j]]) <= half
  dot_rows <- function(X, Y, r0, r1) { t <- 0; for (z in seq_len(n)) t <- t + X[[r0]][[z]] * Y[[r1]][[z]]; t }
  reduce_row <- function(r0, r1) {
    r <- closest(mu[[r0]][[r1]])
    y[[r0]] <<- lapply(seq_len(n), function(z) y[[r0]][[z]] - r * y[[r1]][[z]])
    for (z in seq_len(r1 - 1L)) mu[[r0]][[z]] <<- mu[[r0]][[z]] - r * mu[[r1]][[z]]
    mu[[r0]][[r1]] <<- mu[[r0]][[r1]] - r
    if (return_transform) Tm[[r0]] <<- lapply(seq_len(m), function(z) Tm[[r0]][[z]] - r * Tm[[r1]][[z]])
  }
  for (i in seq_len(m)) {
    ystar[[i]] <- lapply(y[[i]], gmp::as.bigq)
    for (j in seq_len(i - 1L)) {
      rd <- dot_rows(y, ystar, i, j)
      if (gstar[[j]] == 0) stop(.dm_err("DMRankError", lin_err))
      mu[[i]][[j]] <- rd / gstar[[j]]
      ystar[[i]] <- lapply(seq_len(n), function(z) ystar[[i]][[z]] - mu[[i]][[j]] * ystar[[j]][[z]])
    }
    gstar[[i]] <- gmp::as.bigq(dot_rows(ystar, ystar, i, i))
  }
  k <- 2L
  while (k <= m) {
    if (!mu_small(k, k - 1L)) reduce_row(k, k - 1L)
    if (lovasz(k)) {
      for (l in rev(seq_len(k - 2L))) if (!mu_small(k, l)) reduce_row(k, l)
      k <- k + 1L
    } else {
      nu <- mu[[k]][[k - 1L]]
      alpha <- gstar[[k]] + nu^2 * gstar[[k - 1L]]
      if (alpha == 0) stop(.dm_err("DMRankError", lin_err))
      beta <- gstar[[k - 1L]] / alpha
      mu[[k]][[k - 1L]] <- nu * beta
      gstar[[k]] <- gstar[[k]] * beta
      gstar[[k - 1L]] <- alpha
      t <- y[[k]]; y[[k]] <- y[[k - 1L]]; y[[k - 1L]] <- t
      if (k > 2L) for (z in seq_len(k - 2L)) { t <- mu[[k]][[z]]; mu[[k]][[z]] <- mu[[k - 1L]][[z]]; mu[[k - 1L]][[z]] <- t }
      if (k < m) for (i in (k + 1L):m) {
        xi <- mu[[i]][[k]]
        mu[[i]][[k]] <- mu[[i]][[k - 1L]] - nu * xi
        mu[[i]][[k - 1L]] <- mu[[k]][[k - 1L]] * mu[[i]][[k]] + xi
      }
      if (return_transform) { t <- Tm[[k]]; Tm[[k]] <- Tm[[k - 1L]]; Tm[[k - 1L]] <- t }
      k <- max(k - 1L, 2L)
    }
  }
  list(y = .ddm_new(y, s, ZZ), T = if (return_transform) .ddm_new(Tm, c(m, m), ZZ))
}
.ddm_lll <- function(x, delta) .ddm_lll_core(x, gmp::as.bigq(delta), FALSE)$y
.ddm_lll_transform <- function(x, delta) { r <- .ddm_lll_core(x, gmp::as.bigq(delta), TRUE); .pytuple(list(r$y, r$T)) }
