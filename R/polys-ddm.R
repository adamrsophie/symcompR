# Dense domain matrices, after sympy/polys/matrices/ddm.py (DDM) and
# sympy/polys/matrices/dense.py (the in place algorithms on lists of rows).
# A DDM is list(rows = list of rows, each a list of domain elements,
# shape = c(m, n), domain = K). The algorithms below take and return the
# list of rows; R's copy semantics replace Python's in place updates.

.dm_err <- function(name, message = "") {
  parents <- switch(name,
    DMError = character(),
    DMBadInputError = c("DMError", "ValueError"),
    DMDomainError = "DMError", DMNotAField = c("DMDomainError", "DMError"),
    DMFormatError = "DMError", DMNonInvertibleMatrixError = "DMError",
    DMRankError = "DMError", DMShapeError = "DMError",
    DMNonSquareMatrixError = c("DMShapeError", "DMError"), DMValueError = "DMError",
    character())
  structure(list(message = paste(as.character(message), collapse = ""), call = NULL),
            class = c(name, parents, "sympy_error", "error", "condition"))
}
.py_err <- function(name, message = "") .sympy_condition(name, message)

# ---------------------------------------------------------------------------
# dense.py

.ddm_transpose <- function(a, n) {
  m <- length(a)
  if (!m) return(rep(list(list()), n))
  lapply(seq_len(n), function(j) lapply(a, `[[`, j))
}
.ddm_iadd <- function(a, b) Map(function(ai, bi) Map(`+`, ai, bi), a, b)
.ddm_isub <- function(a, b) Map(function(ai, bi) Map(`-`, ai, bi), a, b)
.ddm_ineg <- function(a) lapply(a, function(ai) lapply(ai, function(z) -z))
.ddm_imul <- function(a, b) lapply(a, function(ai) lapply(ai, function(z) z * b))
.ddm_irmul <- function(a, b) lapply(a, function(ai) lapply(ai, function(z) b * z))
# a + b %*% c with a given (the start values of the sums)
.ddm_imatmul <- function(a, b, c) {
  if (!length(c)) return(a)
  nc <- length(c[[1L]])
  for (i in seq_along(b)) {
    bi <- b[[i]]
    for (j in seq_len(nc)) {
      s <- a[[i]][[j]]
      for (k in seq_along(bi)) s <- s + bi[[k]] * c[[k]][[j]]
      a[[i]][[j]] <- s
    }
  }
  a
}

# ddm_irref: list(rows, pivots) with 0-based pivots
.ddm_irref <- function(a, partial_pivot = FALSE) {
  m <- length(a)
  if (!m) return(list(rows = a, pivots = integer()))
  n <- length(a[[1L]])
  i <- 1L
  pivots <- integer()
  for (j in seq_len(n)) {
    if (partial_pivot) {
      cand <- i:m
      ip <- cand[which.max(vapply(cand, function(r) base::Mod(as.complex(a[[r]][[j]])), 0))]
      t <- a[[i]]; a[[i]] <- a[[ip]]; a[[ip]] <- t
    }
    aij <- a[[i]][[j]]
    if (.de_is_zero(aij)) {
      found <- FALSE
      if (i < m) for (ip in (i + 1L):m) {
        aij <- a[[ip]][[j]]
        if (.de_nz(aij)) { t <- a[[i]]; a[[i]] <- a[[ip]]; a[[ip]] <- t; found <- TRUE; break }
      }
      if (!found) next
    }
    ai <- a[[i]]
    aijinv <- aij^-1L
    for (l in j:n) ai[[l]] <- ai[[l]] * aijinv
    a[[i]] <- ai
    for (k in seq_len(m)) {
      ak <- a[[k]]
      if (k == i || .de_is_zero(ak[[j]])) next
      akj <- ak[[j]]
      ak[[j]] <- ak[[j]] - akj
      if (j < n) for (l in (j + 1L):n) ak[[l]] <- ak[[l]] - akj * ai[[l]]
      a[[k]] <- ak
    }
    pivots <- c(pivots, j - 1L)
    i <- i + 1L
    if (i > m) break
  }
  list(rows = a, pivots = pivots)
}

# ddm_irref_den: list(rows, den, pivots)
.ddm_irref_den <- function(a, K) {
  m <- length(a)
  if (!m) return(list(rows = a, den = .de_one(K), pivots = integer()))
  n <- length(a[[1L]])
  d <- NULL
  pivots <- integer(); no_pivots <- integer()
  i <- 1L
  for (j in seq_len(n)) {
    aij <- a[[i]][[j]]
    if (.de_is_zero(aij)) {
      found <- FALSE
      if (i < m) for (ip in (i + 1L):m) {
        aij <- a[[ip]][[j]]
        if (.de_nz(aij)) { t <- a[[i]]; a[[i]] <- a[[ip]]; a[[ip]] <- t; found <- TRUE; break }
      }
      if (!found) { no_pivots <- c(no_pivots, j); next }
    }
    if (length(pivots)) {
      pivot_val <- aij * a[[1L]][[pivots[1L]]]
      if (!is.null(d)) pivot_val <- .de_exquo(K, pivot_val, d)
      for (ip in seq_along(pivots)) a[[ip]][[pivots[ip]]] <- pivot_val
    }
    for (jnp in no_pivots) {
      for (ip in seq_len(i - 1L)) {
        aijp <- a[[ip]][[jnp]]
        if (.de_nz(aijp)) {
          aijp <- aijp * aij
          if (!is.null(d)) aijp <- .de_exquo(K, aijp, d)
          a[[ip]][[jnp]] <- aijp
        }
      }
    }
    ai <- a[[i]]
    for (jp in seq_len(m)) {
      if (jp == i) next
      aj <- a[[jp]]
      if (j < n) for (kp in (j + 1L):n) {
        ajk <- aij * aj[[kp]] - aj[[j]] * ai[[kp]]
        if (!is.null(d)) ajk <- .de_exquo(K, ajk, d)
        aj[[kp]] <- ajk
      }
      aj[[j]] <- .de_zero(K)
      a[[jp]] <- aj
    }
    pivots <- c(pivots, j)
    i <- i + 1L
    if (i > m) break
    d <- if (!.de_is_one(K, aij)) aij else NULL
  }
  den <- if (!length(pivots)) .de_one(K) else a[[1L]][[pivots[1L]]]
  list(rows = a, den = den, pivots = pivots - 1L)
}

.ddm_idet <- function(a, K) {
  m <- length(a)
  if (!m) return(.de_one(K))
  n <- length(a[[1L]])
  uf <- .de_one(K)
  for (k in seq_len(n - 1L)) {
    if (.de_is_zero(a[[k]][[k]])) {
      found <- FALSE
      if (k < n) for (i in (k + 1L):n) {
        if (.de_nz(a[[i]][[k]])) { t <- a[[k]]; a[[k]] <- a[[i]]; a[[i]] <- t; uf <- -uf; found <- TRUE; break }
      }
      if (!found) return(.de_zero(K))
    }
    akkm1 <- if (k > 1L) a[[k - 1L]][[k - 1L]] else .de_one(K)
    for (i in (k + 1L):n) for (j in (k + 1L):n)
      a[[i]][[j]] <- .de_exquo(K, a[[i]][[j]] * a[[k]][[k]] - a[[i]][[k]] * a[[k]][[j]], akkm1)
  }
  uf * a[[m]][[length(a[[m]])]]
}

.ddm_iinv <- function(a, K) {
  if (!.de_is_field(K)) stop(.dm_err("DMDomainError", "Not a field"))
  m <- length(a)
  if (!m) return(a)
  n <- length(a[[1L]])
  if (m != n) stop(.dm_err("DMNonSquareMatrixError"))
  aug <- lapply(seq_len(n), function(i) c(a[[i]], lapply(seq_len(n), function(j) if (i == j) .de_one(K) else .de_zero(K))))
  r <- .ddm_irref(aug)
  if (!identical(as.integer(r$pivots), seq_len(n) - 1L)) stop(.dm_err("DMNonInvertibleMatrixError", "Matrix det == 0; not invertible."))
  lapply(r$rows, function(row) row[(n + 1L):(2L * n)])
}

# ddm_ilu: list(rows, swaps) with 0-based swap pairs
.ddm_ilu <- function(a) {
  m <- length(a)
  if (!m) return(list(rows = a, swaps = list()))
  n <- length(a[[1L]])
  swaps <- list()
  for (i in seq_len(min(m, n))) {
    if (.de_is_zero(a[[i]][[i]])) {
      found <- FALSE
      if (i < m) for (ip in (i + 1L):m) {
        if (.de_nz(a[[ip]][[i]])) {
          swaps[[length(swaps) + 1L]] <- c(i - 1L, ip - 1L)
          t <- a[[i]]; a[[i]] <- a[[ip]]; a[[ip]] <- t; found <- TRUE; break
        }
      }
      if (!found) next
    }
    if (i < m) for (j in (i + 1L):m) {
      l_ji <- a[[j]][[i]] / a[[i]][[i]]
      a[[j]][[i]] <- l_ji
      if (i < n) for (k in (i + 1L):n) a[[j]][[k]] <- a[[j]][[k]] - l_ji * a[[i]][[k]]
    }
  }
  list(rows = a, swaps = swaps)
}
.ddm_ilu_split <- function(L, U, K) {
  m <- length(U)
  if (!m) return(list(L = L, U = U, swaps = list()))
  n <- length(U[[1L]])
  r <- .ddm_ilu(U); U <- r$rows
  if (m > 1L) for (i in 2:m) {
    j <- min(i - 1L, n)
    for (k in seq_len(j)) { L[[i]][[k]] <- U[[i]][[k]]; U[[i]][[k]] <- .de_zero(K) }
  }
  list(L = L, U = U, swaps = r$swaps)
}
.ddm_ilu_solve <- function(x, L, U, swaps, b) {
  m <- length(U)
  if (!m) return(x)
  n <- length(U[[1L]])
  m2 <- length(b)
  if (!m2) stop(.dm_err("DMShapeError", "Shape mismtch"))
  o <- length(b[[1L]])
  if (m != m2) stop(.dm_err("DMShapeError", "Shape mismtch"))
  if (m < n) stop(.py_err("NotImplementedError", "Underdetermined"))
  for (sw in swaps) { t <- b[[sw[1L] + 1L]]; b[[sw[1L] + 1L]] <- b[[sw[2L] + 1L]]; b[[sw[2L] + 1L]] <- t }
  y <- rep(list(vector("list", o)), m)
  for (k in seq_len(o)) for (i in seq_len(m)) {
    rhs <- b[[i]][[k]]
    for (j in seq_len(i - 1L)) rhs <- rhs - L[[i]][[j]] * y[[j]][[k]]
    y[[i]][[k]] <- rhs
  }
  if (m > n) for (i in (n + 1L):m) for (j in seq_len(o)) if (.de_nz(y[[i]][[j]])) stop(.dm_err("DMNonInvertibleMatrixError"))
  for (k in seq_len(o)) for (i in rev(seq_len(n))) {
    if (.de_is_zero(U[[i]][[i]])) stop(.dm_err("DMNonInvertibleMatrixError"))
    rhs <- y[[i]][[k]]
    if (i < n) for (j in (i + 1L):n) rhs <- rhs - U[[i]][[j]] * x[[j]][[k]]
    x[[i]][[k]] <- rhs / U[[i]][[i]]
  }
  x
}

# ddm_berk: the characteristic polynomial as a column of n + 1 elements
.ddm_berk <- function(M, K) {
  m <- length(M)
  if (!m) return(list(list(.de_one(K))))
  n <- length(M[[1L]])
  if (m != n) stop(.dm_err("DMShapeError", "Not square"))
  if (n == 1L) return(list(list(.de_one(K)), list(-M[[1L]][[1L]])))
  a <- M[[1L]][[1L]]
  R <- list(M[[1L]][-1L])
  C <- lapply(M[-1L], function(row) list(row[[1L]]))
  A <- lapply(M[-1L], function(row) row[-1L])
  q <- .ddm_berk(A, K)
  T <- rep(list(rep(list(.de_zero(K)), n)), n + 1L)
  for (i in seq_len(n)) { T[[i]][[i]] <- .de_one(K); T[[i + 1L]][[i]] <- -a }
  AnC <- NULL
  for (i in 2:n) {
    if (i == 2L) AnC <- C
    else {
      C <- AnC
      AnC <- lapply(C, function(row) list(.de_zero(K)))
      AnC <- .ddm_imatmul(AnC, A, C)
    }
    RAnC <- .ddm_imatmul(list(list(.de_zero(K))), R, AnC)
    for (j in 0:(n - i)) T[[i + j + 1L]][[j + 1L]] <- -RAnC[[1L]][[1L]]
  }
  .ddm_imatmul(lapply(seq_len(n + 1L), function(i) list(.de_zero(K))), T, q)
}

# ---------------------------------------------------------------------------
# DDM class

.ddm_new <- function(rows, shape, K) structure(list(rows = rows, shape = as.integer(shape), domain = K), class = "DDM")

#' Dense and sparse domain matrices
#'
#' `DDM(rows, shape, domain)` is a dense matrix over a polynomial domain
#' (SymPy's `sympy.polys.matrices.ddm.DDM`), `SDM(elemsdict, shape,
#' domain)` a sparse one (`sympy.polys.matrices.sdm.SDM`) given as a nested
#' list `list(i = list(j = element))` with 0-based names, or as a list of
#' `list(i, list(list(j, element), ...))` pairs. They are the internal
#' representations of [DomainMatrix()] (`A$rep`), with the same methods as
#' in SymPy, called with `$`: `M$rref()`, `M$det()`, `M$matmul(N)`, ...
#' Indices passed to methods are 0-based as in SymPy.
#'
#' @param rows A list of rows, each a list of domain elements.
#' @param elemsdict The nonzero elements by row and column.
#' @param shape The number of rows and columns.
#' @param domain The domain of the elements.
#' @param x A `DDM` or `SDM`.
#' @param name A field or method name.
#' @param e1,e2 Matrices or elements.
#' @param ... Unused.
#' @return A `DDM` or `SDM` object.
#' @examples
#' A <- DDM(list(list(ZZ(1), ZZ(2)), list(ZZ(3), ZZ(4))), c(2, 2), ZZ)
#' A$det()
#' A$charpoly()
#' B <- SDM(list(`0` = list(`1` = ZZ(2))), c(2, 2), ZZ)
#' B$to_ddm()
#' @name DDM
#' @export
DDM <- function(rows, shape, domain) {
  shape <- .dm_shape(shape)
  if (!is.list(rows) || !all(vapply(rows, function(r) is.list(r) && !is.object(r), NA)))
    stop(.dm_err("DMBadInputError", "rowslist must be a list of lists"))
  if (length(rows) != shape[1L] || any(vapply(rows, length, 0L) != shape[2L]))
    stop(.dm_err("DMBadInputError", "Inconsistent row-list/shape"))
  .ddm_new(lapply(rows, function(r) unname(r)), shape, domain)
}
.dm_shape <- function(shape) {
  s <- as.integer(unlist(lapply(shape, function(z) if (inherits(z, "bigz")) as.integer(as.character(z)) else z)))
  if (length(s) != 2L) stop(.dm_err("DMBadInputError", "shape must be a pair"))
  s
}
.ddm_zeros <- function(shape, K) { z <- .de_zero(K); .ddm_new(rep(list(rep(list(z), shape[2L])), shape[1L]), shape, K) }
.ddm_eye <- function(shape, K) {
  if (length(shape) == 1L) shape <- c(shape, shape)
  M <- .ddm_zeros(shape, K)
  for (i in seq_len(min(shape))) M[["rows"]][[i]][[i]] <- .de_one(K)
  M
}
.ddm_check <- function(a, op, b, ashape, bshape) {
  if (!.de_same(a$domain, b$domain))
    stop(.dm_err("DMDomainError", paste0("Domain mismatch: ", format(a$domain), " ", op, " ", format(b$domain))))
  if (!identical(as.integer(ashape), as.integer(bshape)))
    stop(.dm_err("DMShapeError", paste0("Shape mismatch: ", .shape_str(.subset2(a, "shape")), " ", op, " ", .shape_str(.subset2(b, "shape")))))
}
.shape_str <- function(s) paste0("(", s[1L], ", ", s[2L], ")")

.ddm_to_sdm <- function(a) {
  rows <- .subset2(a, "rows")
  d <- list()
  for (i in seq_along(rows)) {
    nz <- which(vapply(rows[[i]], .de_nz, NA))
    if (length(nz)) d[[as.character(i - 1L)]] <- stats::setNames(rows[[i]][nz], as.character(nz - 1L))
  }
  .sdm_new(d, .subset2(a, "shape"), .subset2(a, "domain"))
}
.ddm_convert_to <- function(a, K) {
  K0 <- .subset2(a, "domain")
  if (.de_same(K, K0)) return(a)
  .ddm_new(lapply(.subset2(a, "rows"), function(r) lapply(r, function(z) .de_convert(K, z, K0))), .subset2(a, "shape"), K)
}
.ddm_nullspace_from_rref <- function(a, pivots = NULL) {
  rows <- .subset2(a, "rows"); s <- .subset2(a, "shape"); m <- s[1L]; n <- s[2L]; K <- .subset2(a, "domain")
  if (is.null(pivots)) {
    pivots <- integer(); last <- -1L
    for (i in seq_len(m)) {
      ai <- rows[[i]]
      if (last + 1L < n) for (j in (last + 1L):(n - 1L)) if (.de_nz(ai[[j + 1L]])) { last <- j; pivots <- c(pivots, j); break }
    }
  }
  if (!length(pivots)) return(.pytuple(list(.ddm_eye(c(n, n), K), as.list(seq_len(n) - 1L))))
  pivot_val <- rows[[1L]][[pivots[1L] + 1L]]
  basis <- list(); nonpivots <- list()
  for (i in seq_len(n) - 1L) {
    if (i %in% pivots) next
    nonpivots[[length(nonpivots) + 1L]] <- i
    vec <- lapply(seq_len(n) - 1L, function(j) if (i == j) pivot_val else .de_zero(K))
    for (ii in seq_along(pivots)) { jj <- pivots[ii] + 1L; vec[[jj]] <- vec[[jj]] - rows[[ii]][[i + 1L]] }
    basis[[length(basis) + 1L]] <- vec
  }
  .pytuple(list(.ddm_new(basis, c(length(basis), n), K), nonpivots))
}
.ddm_fflu_core <- function(rows, K, nr, nc) {
  LU <- rows; perm <- seq_len(nr) - 1L; rank <- 0L
  for (j in seq_len(min(nr, nc))) {
    if (all(vapply(seq_len(nr), function(i) .de_is_zero(LU[[i]][[j]]), NA))) next
    pivot_row <- -1L
    if (rank < nr) for (i in (rank + 1L):nr) if (.de_nz(LU[[i]][[j]])) { pivot_row <- i; break }
    if (pivot_row == -1L) next
    if (pivot_row != rank + 1L) {
      t <- LU[[rank + 1L]]; LU[[rank + 1L]] <- LU[[pivot_row]]; LU[[pivot_row]] <- t
      t <- perm[rank + 1L]; perm[rank + 1L] <- perm[pivot_row]; perm[pivot_row] <- t
    }
    pivot <- LU[[rank + 1L]][[j]]
    if (rank + 1L < nr) for (i in (rank + 2L):nr) {
      multiplier <- LU[[i]][[j]]
      denominator <- if (rank > 0L) LU[[rank]][[rank]] else .de_one(K)
      if (j < nc) for (k in (j + 1L):nc) LU[[i]][[k]] <- .de_exquo(K, pivot * LU[[i]][[k]] - LU[[rank + 1L]][[k]] * multiplier, denominator)
      LU[[i]][[j]] <- multiplier
    }
    rank <- rank + 1L
  }
  list(LU = LU, perm = perm)
}
.ddm_fflu <- function(a) {
  s <- .subset2(a, "shape"); rows <- s[1L]; cols <- s[2L]; K <- .subset2(a, "domain")
  r <- .ddm_fflu_core(.subset2(a, "rows"), K, rows, cols)
  U <- r$LU
  P <- .ddm_zeros(c(rows, rows), K)[["rows"]]
  for (i in seq_along(r$perm)) P[[i]][[r$perm[i] + 1L]] <- .de_one(K)
  L <- .ddm_zeros(c(rows, rows), K)[["rows"]]
  i <- 1L; j <- 1L
  while (i <= rows && j <= cols) {
    if (.de_nz(U[[i]][[j]])) {
      L[[i]][[i]] <- U[[i]][[j]]
      if (i < rows) for (l in (i + 1L):rows) { L[[l]][[i]] <- U[[l]][[j]]; U[[l]][[j]] <- .de_zero(K) }
      i <- i + 1L
    }
    j <- j + 1L
  }
  if (i <= rows) for (ii in i:rows) L[[ii]][[ii]] <- .de_one(K)
  D <- .ddm_zeros(c(rows, rows), K)[["rows"]]
  if (rows >= 1L) D[[1L]][[1L]] <- L[[1L]][[1L]]
  if (rows >= 2L) for (ii in 2:rows) D[[ii]][[ii]] <- L[[ii - 1L]][[ii - 1L]] * L[[ii]][[ii]]
  .pytuple(list(.ddm_new(P, c(rows, rows), K), .ddm_new(L, c(rows, rows), K), .ddm_new(D, c(rows, rows), K),
                .ddm_new(U, s, K)))
}
.ddm_qr <- function(a) {
  s <- .subset2(a, "shape"); rows <- s[1L]; cols <- s[2L]; K <- .subset2(a, "domain")
  Q <- .subset2(a, "rows")
  R <- .ddm_zeros(c(min(rows, cols), cols), K)[["rows"]]
  if (!.de_is_field(K)) stop(.dm_err("DMDomainError", "QR decomposition requires a field (e.g. QQ)."))
  dot_cols <- function(i, j) { t <- .de_zero(K); for (k in seq_len(rows)) t <- t + Q[[k]][[i]] * Q[[k]][[j]]; t }
  for (j in seq_len(cols)) {
    for (i in seq_len(min(j - 1L, rows))) {
      dot_ii <- dot_cols(i, i)
      if (.de_nz(dot_ii)) {
        R[[i]][[j]] <- dot_cols(i, j) / dot_ii
        for (k in seq_len(rows)) Q[[k]][[j]] <- Q[[k]][[j]] - R[[i]][[j]] * Q[[k]][[i]]
      }
    }
    if (j <= rows) {
      dot_jj <- dot_cols(j, j)
      if (.de_nz(dot_jj)) R[[j]][[j]] <- .de_one(K)
    }
  }
  k <- min(rows, cols)
  .pytuple(list(.ddm_new(lapply(Q, function(r) r[seq_len(k)]), c(rows, k), K), .ddm_new(R, c(k, cols), K)))
}

.ddm_method <- function(a, name) {
  rows <- .subset2(a, "rows"); s <- .subset2(a, "shape"); K <- .subset2(a, "domain")
  m <- s[1L]; n <- s[2L]
  switch(name,
    shape = .pytuple(as.list(s)), rows = m, cols = n, domain = K, fmt = "dense", is_DDM = TRUE, is_DFM = FALSE,
    getitem = function(i, j) rows[[.dm_idx(i, m)]][[.dm_idx(j, n)]],
    setitem = function(i, j, value) { a[["rows"]][[.dm_idx(i, m)]][[.dm_idx(j, n)]] <- value; a },
    extract = function(rws, cls) .ddm_new(lapply(rws, function(i) lapply(cls, function(j) rows[[.dm_idx(i, m)]][[.dm_idx(j, n)]])), c(length(rws), length(cls)), K),
    to_list = function() rows,
    to_list_flat = function() unlist(rows, recursive = FALSE),
    flat = function() unlist(rows, recursive = FALSE),
    to_flat_nz = function() .ddm_to_sdm(a)$to_flat_nz(),
    to_dod = function() .sdm_to_dod_pd(.ddm_to_sdm(a)),
    to_dok = function() .sdm_method(.ddm_to_sdm(a), "to_dok")(),
    iter_values = function() Filter(.de_nz, unlist(rows, recursive = FALSE)),
    to_ddm = function() a, to_sdm = function() .ddm_to_sdm(a),
    to_dfm = function() stop(.py_err("NotImplementedError", "DFM requires python-flint")),
    to_dfm_or_ddm = function() a,
    convert_to = function(K2) .ddm_convert_to(a, K2),
    copy = function() a,
    transpose = function() .ddm_new(.ddm_transpose(rows, n), c(n, m), K),
    add = function(b) { .ddm_check(a, "+", b, s, .subset2(b, "shape")); .ddm_new(.ddm_iadd(rows, .subset2(b, "rows")), s, K) },
    sub = function(b) { .ddm_check(a, "-", b, s, .subset2(b, "shape")); .ddm_new(.ddm_isub(rows, .subset2(b, "rows")), s, K) },
    neg = function() .ddm_new(.ddm_ineg(rows), s, K),
    mul = function(b) .ddm_new(.ddm_imul(rows, b), s, K),
    rmul = function(b) .ddm_new(.ddm_irmul(rows, b), s, K),
    matmul = function(b) {
      bs <- .subset2(b, "shape")
      .ddm_check(a, "*", b, n, bs[1L])
      .ddm_new(.ddm_imatmul(.ddm_zeros(c(m, bs[2L]), K)[["rows"]], rows, .subset2(b, "rows")), c(m, bs[2L]), K)
    },
    mul_elementwise = function(b) .ddm_new(Map(function(ai, bi) Map(`*`, ai, bi), rows, .subset2(b, "rows")), s, K),
    hstack = function(...) {
      new <- rows; cols <- n
      for (B in list(...)) { cols <- cols + .subset2(B, "shape")[2L]; new <- Map(c, new, .subset2(B, "rows")) }
      .ddm_new(new, c(m, cols), K)
    },
    vstack = function(...) {
      new <- rows; rr <- m
      for (B in list(...)) { rr <- rr + .subset2(B, "shape")[1L]; new <- c(new, .subset2(B, "rows")) }
      .ddm_new(new, c(rr, n), K)
    },
    applyfunc = function(func, domain) .ddm_new(lapply(rows, function(r) lapply(r, func)), s, domain),
    nnz = function() sum(vapply(unlist(rows, recursive = FALSE), .de_nz, NA)),
    scc = function() .sdm_method(.ddm_to_sdm(a), "scc")(),
    rref = function() {
      r <- .ddm_irref(rows, .de_kind(K) %in% c("RR", "CC"))
      .pytuple(list(.ddm_new(r$rows, s, K), as.list(r$pivots)))
    },
    rref_den = function() {
      r <- .ddm_irref_den(rows, K)
      .pytuple(list(.ddm_new(r$rows, s, K), r$den, as.list(r$pivots)))
    },
    nullspace = function() { r <- .ddm_method(a, "rref")(); .ddm_nullspace_from_rref(r[[1L]], unlist(r[[2L]])) },
    nullspace_from_rref = function(pivots = NULL) .ddm_nullspace_from_rref(a, if (is.null(pivots)) NULL else as.integer(unlist(pivots))),
    particular = function() .sdm_to_ddm(.sdm_method(.ddm_to_sdm(a), "particular")()),
    det = function() {
      if (m != n) stop(.dm_err("DMNonSquareMatrixError", "Determinant of non-square matrix"))
      .ddm_idet(rows, K)
    },
    inv = function() {
      if (m != n) stop(.dm_err("DMNonSquareMatrixError", "Determinant of non-square matrix"))
      .ddm_new(.ddm_iinv(rows, K), s, K)
    },
    lu = function() {
      r <- .ddm_ilu_split(.ddm_eye(c(m, m), K)[["rows"]], rows, K)
      .pytuple(list(.ddm_new(r$L, c(m, m), K), .ddm_new(r$U, s, K), lapply(r$swaps, function(z) .pytuple(as.list(z)))))
    },
    fflu = function() .ddm_fflu(a),
    qr = function() .ddm_qr(a),
    lu_solve = function(b) {
      bs <- .subset2(b, "shape")
      .ddm_check(a, "lu_solve", b, m, bs[1L])
      if (!.de_is_field(K)) stop(.dm_err("DMDomainError", "lu_solve requires a field"))
      r <- .ddm_ilu_split(.ddm_eye(c(m, m), K)[["rows"]], rows, K)
      x <- .ddm_zeros(c(n, bs[2L]), K)[["rows"]]
      .ddm_new(.ddm_ilu_solve(x, r$L, r$U, r$swaps, .subset2(b, "rows")), c(n, bs[2L]), K)
    },
    charpoly = function() {
      if (m != n) stop(.dm_err("DMNonSquareMatrixError", "Charpoly of non-square matrix"))
      vec <- .ddm_berk(rows, K)
      lapply(seq_len(n + 1L), function(i) vec[[i]][[1L]])
    },
    is_zero_matrix = function() all(vapply(unlist(rows, recursive = FALSE), .de_is_zero, NA)),
    is_upper = function() all(vapply(seq_len(m), function(i) all(vapply(rows[[i]][seq_len(min(i - 1L, n))], .de_is_zero, NA)), NA)),
    is_lower = function() all(vapply(seq_len(m), function(i) if (i >= n) TRUE else all(vapply(rows[[i]][(i + 1L):n], .de_is_zero, NA)), NA)),
    is_diagonal = function() .ddm_method(a, "is_upper")() && .ddm_method(a, "is_lower")(),
    diagonal = function() lapply(seq_len(min(m, n)), function(i) rows[[i]][[i]]),
    lll = function(delta = gmp::as.bigq(3, 4)) .ddm_lll(a, delta),
    lll_transform = function(delta = gmp::as.bigq(3, 4)) .ddm_lll_transform(a, delta),
    zeros = function(shape, domain) .ddm_zeros(.dm_shape(shape), domain),
    ones = function(shape, domain) { sh <- .dm_shape(shape); .ddm_new(rep(list(rep(list(.de_one(domain)), sh[2L])), sh[1L]), sh, domain) },
    eye = function(size, domain) .ddm_eye(if (is.list(size)) .dm_shape(size) else as.integer(size), domain),
    stop(.py_err("AttributeError", paste0("'DDM' object has no attribute '", name, "'"))))
}
# a 0-based Python index (negative from the end) to an R index
.dm_idx <- function(i, n) {
  i <- as.integer(as.character(i))
  if (i < 0L) i <- i + n
  if (i < 0L || i >= n) stop(.py_err("IndexError", "list index out of range"))
  i + 1L
}

#' @rdname DDM
#' @export
`$.DDM` <- function(x, name) .ddm_method(x, name)
#' @rdname DDM
#' @export
format.DDM <- function(x, ...) .dm_repr(x)
#' @rdname DDM
#' @export
print.DDM <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname DDM
#' @export
Ops.DDM <- function(e1, e2) .dm_ops(.Generic, e1, e2, "DDM")
