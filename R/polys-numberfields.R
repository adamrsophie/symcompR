# Algebraic number theory: modules in number fields, the Round Two
# algorithm for integral bases, prime decomposition and prime valuations,
# after sympy/polys/numberfields/{modules,utilities,basis,primes}.py.
# Integer matrices are list(nrow, cols) with cols a list of bigz column
# vectors; matrices over a field are lists of rows of field elements.

.bz <- function(x) gmp::as.bigz(x)
.bq <- function(x) gmp::as.bigq(x)
.zmat <- function(nrow, cols) list(nrow = as.integer(nrow), cols = cols)
.zmat_ncol <- function(M) length(M$cols)
.zmat_eye <- function(n) .zmat(n, lapply(seq_len(n), function(j) { v <- .bz(rep(0L, n)); v[j] <- 1; v }))
.zmat_get <- function(M, i, j) M$cols[[j]][i]
.zmat_hstack <- function(...) {
  ms <- list(...)
  .zmat(ms[[1L]]$nrow, do.call(c, lapply(ms, `[[`, "cols")))
}
.zmat_scale <- function(M, c) .zmat(M$nrow, lapply(M$cols, function(v) v * c))
.zmat_mul <- function(A, B) {
  # A (n x k) times B (k x m)
  .zmat(A$nrow, lapply(B$cols, function(b) {
    out <- .bz(rep(0L, A$nrow))
    for (k in seq_along(b)) if (b[k] != 0) out <- out + A$cols[[k]] * b[k]
    out
  }))
}
.zmat_to_rows <- function(M) lapply(seq_len(M$nrow), function(i) .bz(vapply(M$cols, function(v) as.character(v[i]), "")))
.zmat_from_rows <- function(rows, ncol) .zmat(length(rows), lapply(seq_len(ncol), function(j) .bz(vapply(rows, function(r) as.character(r[j]), ""))))
.zmat_hnf <- function(M, D = NULL) {
  rows <- .zmat_to_rows(M)
  n <- .zmat_ncol(M)
  if (!is.null(D)) { W <- .hnf_mod_D(rows, n, D); return(.zmat_from_rows(W, M$nrow)) }
  h <- .hnf(rows, n)
  .zmat_from_rows(h, if (length(h)) length(h[[1L]]) else 0L)
}
.zmat_flat <- function(M) do.call(c, lapply(seq_len(M$nrow), function(i) .bz(vapply(M$cols, function(v) as.character(v[i]), ""))))
.zmat_equal <- function(A, B) A$nrow == B$nrow && .zmat_ncol(A) == .zmat_ncol(B) && all(vapply(seq_along(A$cols), function(j) all(A$cols[[j]] == B$cols[[j]]), NA))
# determinant over ZZ by fraction-free (Bareiss) elimination
.zmat_det <- function(M) {
  n <- M$nrow
  if (n != .zmat_ncol(M)) stop("Determinant of non-square matrix", call. = FALSE)
  if (n == 0L) return(.bz(1))
  a <- .zmat_to_rows(M)
  sgn <- 1L; prev <- .bz(1)
  for (k in seq_len(n - 1L)) {
    if (a[[k]][k] == 0) {
      sw <- which(vapply((k + 1L):n, function(i) a[[i]][k] != 0, NA))
      if (!length(sw)) return(.bz(0))
      i <- k + sw[1L]; t <- a[[k]]; a[[k]] <- a[[i]]; a[[i]] <- t; sgn <- -sgn
    }
    for (i in (k + 1L):n) {
      for (j in (k + 1L):n) a[[i]][j] <- gmp::divq.bigz(a[[i]][j] * a[[k]][k] - a[[i]][k] * a[[k]][j], prev)
      a[[i]][k] <- 0
    }
    prev <- a[[k]][k]
  }
  a[[n]][n] * sgn
}

# fields for the dense elimination: GF(p) (values in [0, p)) and QQ
.nfK_ff <- function(p) {
  p <- .bz(p)
  list(kind = "FF", p = p, zero = .bz(0), one = .bz(1),
       add = function(a, b) gmp::mod.bigz(a + b, p), sub = function(a, b) gmp::mod.bigz(a - b, p),
       mul = function(a, b) gmp::mod.bigz(a * b, p), neg = function(a) gmp::mod.bigz(-a, p),
       quo = function(a, b) gmp::mod.bigz(a * gmp::inv.bigz(b, p), p), is_zero = function(a) a == 0,
       from_int = function(a) gmp::mod.bigz(.bz(a), p),
       to_int = function(a) { a <- gmp::mod.bigz(a, p); if (a <= gmp::divq.bigz(p, 2)) a else a - p })
}
.nfK_qq <- function() {
  list(kind = "QQ", zero = .bq(0), one = .bq(1),
       add = function(a, b) a + b, sub = function(a, b) a - b, mul = function(a, b) a * b, neg = function(a) -a,
       quo = function(a, b) a / b, is_zero = function(a) a == 0, from_int = function(a) .bq(a),
       to_int = function(a) { if (gmp::denominator(a) != 1) stop("CoercionFailed", call. = FALSE); gmp::numerator(a) })
}
# rows of field elements from an integer matrix
.fmat_from_zmat <- function(M, K) lapply(seq_len(M$nrow), function(i) lapply(M$cols, function(v) K$from_int(v[i])))

# sdm_rref_den: fraction-free Gauss-Jordan; rows of field elements
.fmat_rref_den <- function(rows, K, ncol) {
  # sparse rows: named lists keyed by column index
  A <- lapply(rows, function(r) { nz <- which(vapply(r, function(z) !K$is_zero(z), NA)); s <- r[nz]; names(s) <- as.character(nz); s })
  A <- Filter(length, A)
  if (!length(A)) return(list(rows = list(), den = K$one, pivots = integer()))
  minc <- function(r) min(as.integer(names(r)))
  if (length(A) == 1L) {
    Ai <- A[[1L]]; j <- minc(Ai)
    return(list(rows = list(Ai), den = Ai[[as.character(j)]], pivots = j))
  }
  A <- A[order(vapply(A, minc, 0L), method = "radix")]
  red <- list(); unred <- list()  # column (character) -> row index in out
  out <- list(); denom <- NULL; divisor <- NULL
  for (Ai in A) {
    Ai <- Ai[!(names(Ai) %in% names(red))]
    canc <- list()
    for (j in intersect(names(unred), names(Ai))) {
      Aij <- Ai[[j]]; Ai[[j]] <- NULL
      Aj <- out[[unred[[j]]]]
      for (k in names(Aj)) {
        v <- K$mul(Aij, Aj[[k]])
        if (is.null(canc[[k]])) canc[[k]] <- v
        else { s <- K$add(canc[[k]], v); if (K$is_zero(s)) canc[[k]] <- NULL else canc[[k]] <- s }
      }
    }
    d <- if (is.null(denom)) K$one else denom
    keys <- union(names(Ai), names(canc))
    newAi <- list()
    for (k in keys) {
      a <- Ai[[k]]; cc <- canc[[k]]
      v <- if (is.null(a)) K$neg(cc) else if (is.null(cc)) K$mul(a, d) else K$sub(K$mul(a, d), cc)
      if (!K$is_zero(v)) newAi[[k]] <- v
    }
    Ai <- newAi
    if (!length(Ai)) next
    j <- as.character(minc(Ai))
    Aij <- Ai[[j]]; Ai[[j]] <- NULL
    for (pk in names(unred)) {
      kidx <- unred[[pk]]
      Ak <- out[[kidx]]
      if (is.null(Ak[[j]])) {
        for (l in names(Ak)) { v <- K$mul(Ak[[l]], Aij); if (!is.null(divisor)) v <- K$quo(v, divisor); Ak[[l]] <- v }
        out[[kidx]] <- Ak
        next
      }
      Akj <- Ak[[j]]; Ak[[j]] <- NULL
      keys2 <- union(names(Ai), names(Ak))
      newAk <- list()
      for (l in keys2) {
        ai <- Ai[[l]]; ak <- Ak[[l]]
        v <- if (is.null(ak)) K$neg(K$mul(Akj, ai)) else if (is.null(ai)) K$mul(Aij, ak) else K$sub(K$mul(Aij, ak), K$mul(Akj, ai))
        if (!K$is_zero(v)) { if (!is.null(divisor)) v <- K$quo(v, divisor); newAk[[l]] <- v }
      }
      out[[kidx]] <- newAk
      if (!length(newAk)) { unred[[pk]] <- NULL; red[[pk]] <- kidx }
    }
    i <- length(out) + 1L
    out[[i]] <- if (length(Ai)) Ai else list()
    if (length(Ai)) unred[[j]] <- i else red[[j]] <- i
    if (!(Aij == K$one)) denom <- if (is.null(denom)) Aij else K$mul(denom, Aij)
    if (!is.null(divisor)) denom <- K$quo(denom, divisor)
    divisor <- denom
  }
  if (is.null(denom)) denom <- K$one
  c2r <- c(red, unred)
  pivots <- sort(as.integer(names(c2r)))
  rws <- lapply(pivots, function(pc) { r <- out[[c2r[[as.character(pc)]]]]; r[[as.character(pc)]] <- denom; r })
  list(rows = rws, den = denom, pivots = pivots)
}
# DomainMatrix.nullspace (divide_last = False): rows of the basis
.fmat_nullspace <- function(rows, K, ncol) {
  rr <- .fmat_rref_den(rows, K, ncol)
  pivots <- rr$pivots
  if (!length(pivots)) return(lapply(seq_len(ncol), function(j) { v <- rep(list(K$zero), ncol); v[[j]] <- K$one; v }))
  if (length(pivots) == ncol) return(list())
  pv <- rr$rows[[1L]][[as.character(pivots[1L])]]
  lapply(setdiff(seq_len(ncol), pivots), function(j) {
    v <- rep(list(K$zero), ncol); v[[j]] <- pv
    for (ii in seq_along(pivots)) { a <- rr$rows[[ii]][[as.character(j)]]; if (!is.null(a)) v[[pivots[ii]]] <- K$neg(a) }
    v
  })
}
# rref by Gauss-Jordan (unique): rows and pivots
.fmat_rref <- function(rows, K, ncol) {
  a <- rows; m <- length(a); pivots <- integer(); r <- 1L
  for (j in seq_len(ncol)) {
    if (r > m) break
    piv <- NULL
    for (i in r:m) if (!K$is_zero(a[[i]][[j]])) { piv <- i; break }
    if (is.null(piv)) next
    t <- a[[r]]; a[[r]] <- a[[piv]]; a[[piv]] <- t
    inv <- a[[r]][[j]]
    a[[r]] <- lapply(a[[r]], function(z) K$quo(z, inv))
    for (i in seq_len(m)) if (i != r && !K$is_zero(a[[i]][[j]])) {
      f <- a[[i]][[j]]
      a[[i]] <- lapply(seq_len(ncol), function(c) K$sub(a[[i]][[c]], K$mul(f, a[[r]][[c]])))
    }
    pivots <- c(pivots, j); r <- r + 1L
  }
  list(rows = a, pivots = pivots)
}
.fmat_cols_to_zmat <- function(vecs, K, nrow) .zmat(nrow, lapply(vecs, function(v) .bz(vapply(v, function(z) as.character(K$to_int(z)), ""))))

# ---------------------------------------------------------------------------
# utilities.py

# extract_fundamental_discriminant: list(D, F) of named integer vectors
# (names are the primes, in factorint order)
.extract_fundamental_discriminant <- function(a) {
  a <- .bz(a)
  if (!(gmp::mod.bigz(a, 4) %in% c(0, 1))) stop("To extract fundamental discriminant, number must be 0 or 1 mod 4.", call. = FALSE)
  if (a == 0) return(list(D = list(), F = list(`0` = 1L)))
  if (a == 1) return(list(D = list(), F = list()))
  fa <- factorint(as_sym(a))
  # SymPy's factorint inserts -1 after the prime factors
  nm <- names(fa)
  nm <- c(setdiff(nm, "-1"), intersect(nm, "-1"))
  D <- list(); F <- list(); n34 <- 0L
  for (pn in nm) {
    p <- .bz(pn); e <- as.integer(fa[[pn]])
    if (e %% 2L == 1L) {
      D[[pn]] <- 1L
      if (gmp::mod.bigz(p, 4) == 3) n34 <- n34 + 1L  # -1 counts: -1 % 4 == 3
      if (e >= 3L) F[[pn]] <- (e - 1L) %/% 2L
    } else F[[pn]] <- e %/% 2L
  }
  even <- !is.null(D[["2"]])
  if (even || n34 %% 2L == 1L) {
    e2 <- F[["2"]]
    if (e2 == 1L) F[["2"]] <- NULL else F[["2"]] <- e2 - 1L
    D[["2"]] <- if (even) 3L else 2L
  }
  list(D = D, F = F)
}

#' Fundamental discriminant
#'
#' `extract_fundamental_discriminant(a)` writes an integer `a` that is 0 or
#' 1 modulo 4 as `a = D*F^2` with `D` a fundamental discriminant, returning
#' the factorizations of `D` and `F` as named lists of exponents, as SymPy's
#' function of the same name.
#'
#' @param a An integer.
#' @return A list with elements `D` and `F`.
#' @examples
#' extract_fundamental_discriminant(-4 * 3^5 * 5^3)
#' @export
extract_fundamental_discriminant <- function(a) {
  r <- .extract_fundamental_discriminant(.bz(as.character(as_sym(a)$value)))
  lapply(r, function(d) lapply(d, as.integer))
}

# AlgIntPowers: the powers of theta in ZZ[theta] (T monic over ZZ), as
# bigz vectors of length n
.alg_int_powers <- function(Tc) {
  # Tc: bigz coefficients of T, leading first (monic)
  n <- length(Tc) - 1L
  env <- new.env()
  env$r <- list(-rev(Tc)[seq_len(n)])
  env$max <- n
  get <- function(e) {
    if (e < n) { v <- .bz(rep(0L, n)); v[e + 1L] <- 1; return(v) }
    if (e > env$max) {
      r <- env$r; c0 <- r[[1L]]
      for (k in (env$max + 1L):e) {
        prev <- r[[k - 1L - n + 1L]]
        b <- prev[n]
        nv <- .bz(rep(0L, n))
        nv[1L] <- c0[1L] * b
        if (n >= 2L) for (i in 2:n) nv[i] <- prev[i - 1L] + c0[i] * b
        r[[length(r) + 1L]] <- nv
      }
      env$r <- r; env$max <- e
    }
    env$r[[e - n + 1L]]
  }
  get
}

# coeff_search(m, R) as an iterator
.coeff_search <- function(m, R) {
  st <- new.env(); st$R0 <- R; st$R <- R; st$c <- rep(R, m); st$first <- TRUE
  function() {
    repeat {
      if (!st$first) {
        c <- st$c; j <- m
        while (c[j] == -st$R) j <- j - 1L
        c[j] <- c[j] - 1L
        if (j < m) c[(j + 1L):m] <- st$R
        if (all(c == 0L)) { st$R <- st$R + 1L; c <- rep(st$R, m) }
        st$c <- c
      }
      st$first <- FALSE
      c <- st$c
      if (st$R == st$R0 || (st$R %in% c) || ((-st$R) %in% c)) return(c)
    }
  }
}

# ---------------------------------------------------------------------------
# modules.py

.nf_power_basis <- function(Tc, gen) {
  e <- new.env()
  e$kind <- "PowerBasis"; e$Tc <- Tc; e$n <- length(Tc) - 1L; e$gen <- gen; e$parent <- NULL; e$mult_tab <- NULL
  e
}
.nf_submodule <- function(parent, M, denom = 1, mult_tab = NULL) {
  e <- new.env()
  e$kind <- "Submodule"; e$parent <- parent; e$M <- M; e$denom <- .bz(denom); e$n <- .zmat_ncol(M)
  e$mult_tab <- mult_tab; e$starts_with_unity <- NULL; e$is_sq_maxrank_HNF <- NULL
  e
}
.nf_elt <- function(module, col, denom = 1) list(module = module, col = .bz(col), denom = .bz(denom))
.nf_ancestors <- function(mod, include_self = FALSE) {
  a <- if (is.null(mod$parent)) list() else .nf_ancestors(mod$parent, TRUE)
  if (include_self) a <- c(a, list(mod))
  a
}
.nf_mod_equal <- function(A, B) {
  if (identical(A, B)) return(TRUE)
  if (A$kind == "PowerBasis" && B$kind == "PowerBasis") return(identical(as.character(A$Tc), as.character(B$Tc)))
  if (A$kind == "Submodule" && B$kind == "Submodule" && .nf_mod_equal(A$parent, B$parent))
    return(.nf_qq_matrix_equal(A, B))
  FALSE
}
.nf_qq_matrix_equal <- function(A, B) {
  if (A$n != B$n) return(FALSE)
  all(vapply(seq_len(A$n), function(j) all(.bq(A$M$cols[[j]]) / .bq(A$denom) == .bq(B$M$cols[[j]]) / .bq(B$denom)), NA))
}
.nf_nca <- function(A, B) {
  sa <- .nf_ancestors(A, TRUE); sb <- .nf_ancestors(B, TRUE)
  nca <- NULL
  for (i in seq_len(min(length(sa), length(sb)))) { if (.nf_mod_equal(sa[[i]], sb[[i]])) nca <- sa[[i]] else break }
  nca
}
.nf_basis_elt <- function(mod, j) { v <- .bz(rep(0L, mod$n)); v[j] <- 1; .nf_elt(mod, v) }
.nf_basis_elements <- function(mod) lapply(seq_len(mod$n), function(j) .nf_basis_elt(mod, j))
.nf_ilcm <- function(a, b) gmp::lcm.bigz(a, b)
.nf_reduced <- function(el) {
  if (el$denom == 1) return(el)
  g <- el$denom; for (c in as.list(el$col)) g <- gmp::gcd.bigz(g, c)
  if (g == 1) return(el)
  .nf_elt(el$module, gmp::divq.bigz(el$col, g), gmp::divq.bigz(el$denom, g))
}
.nf_to_parent <- function(el) {
  m <- el$module
  .nf_elt(m$parent, .zmat_mul(m$M, .zmat(m$n, list(el$col)))$cols[[1L]], m$denom * el$denom)
}
.nf_to_ancestor <- function(el, anc) if (.nf_mod_equal(el$module, anc)) el else .nf_to_ancestor(.nf_to_parent(el), anc)
.nf_over_power_basis <- function(el) { while (el$module$kind != "PowerBasis") el <- .nf_to_parent(el); el }
.nf_unify <- function(a, b) {
  if (.nf_mod_equal(a$module, b$module)) return(list(a, b))
  nca <- .nf_nca(a$module, b$module)
  if (is.null(nca)) stop("Cannot unify", call. = FALSE)
  list(.nf_to_ancestor(a, nca), .nf_to_ancestor(b, nca))
}
.nf_is_rat <- function(x) !is.list(x)
.nf_rat_num_den <- function(x) { q <- .bq(x); list(gmp::numerator(q), gmp::denominator(q)) }
.nf_add <- function(a, b) {
  if (.nf_is_rat(b)) return(.nf_add(a, .nf_element_from_rational(a$module, b)))
  if (.nf_is_rat(a)) return(.nf_add(b, a))
  if (!.nf_mod_equal(a$module, b$module)) { u <- .nf_unify(a, b); return(.nf_add(u[[1L]], u[[2L]])) }
  m <- .nf_ilcm(a$denom, b$denom)
  .nf_reduced(.nf_elt(a$module, a$col * gmp::divq.bigz(m, a$denom) + b$col * gmp::divq.bigz(m, b$denom), m))
}
.nf_scal <- function(a, q) {
  nd <- .nf_rat_num_den(q)
  if (nd[[1L]] == 1 && nd[[2L]] == 1) return(a)
  .nf_reduced(.nf_elt(a$module, a$col * nd[[1L]], a$denom * nd[[2L]]))
}
.nf_neg <- function(a) .nf_scal(a, -1)
.nf_sub <- function(a, b) .nf_add(a, if (.nf_is_rat(b)) -.bq(b) else .nf_neg(b))
.nf_mult_tab <- function(mod) {
  if (is.null(mod$mult_tab)) {
    if (mod$kind == "PowerBasis") {
      pw <- .alg_int_powers(mod$Tc)
      n <- mod$n
      mod$mult_tab <- lapply(seq_len(n), function(u) { l <- vector("list", n); for (v in u:n) l[[v]] <- pw(u + v - 2L); l })
    } else {
      gens <- lapply(.nf_basis_elements(mod), .nf_to_parent)
      n <- mod$n
      mod$mult_tab <- lapply(seq_len(n), function(u) { l <- vector("list", n); for (v in u:n) l[[v]] <- .nf_represent(mod, .nf_mul(gens[[u]], gens[[v]])); l })
    }
  }
  mod$mult_tab
}
.nf_mul <- function(a, b) {
  if (.nf_is_rat(b)) return(.nf_scal(a, b))
  if (.nf_is_rat(a)) return(.nf_scal(b, a))
  if (!.nf_mod_equal(a$module, b$module)) { u <- .nf_unify(a, b); return(.nf_mul(u[[1L]], u[[2L]])) }
  M <- .nf_mult_tab(a$module)
  A <- a$col; B <- b$col; n <- a$module$n
  C <- .bz(rep(0L, n))
  for (u in seq_len(n)) for (v in u:n) {
    c <- A[u] * B[v]
    if (v > u) c <- c + A[v] * B[u]
    if (c != 0) C <- C + M[[u]][[v]] * c
  }
  .nf_elt(a$module, C, a$denom * b$denom)
}
.nf_pow <- function(a, e) {
  # IntegerPowerable.__pow__: binary expansion, low bits first
  e <- as.integer(e)
  if (e == 0L) return(.nf_one(a$module))
  if (e == 1L) return(a)
  bits <- as.integer(intToBits(e))[seq_len(floor(log2(e)) + 1L)]
  p <- a; r <- NULL
  for (i in seq_along(bits)) {
    if (bits[i] == 1L) r <- if (is.null(r)) p else .nf_mul(r, p)
    if (i < length(bits)) p <- .nf_mul(p, p)
  }
  r
}
.nf_element_from_rational <- function(mod, a) {
  if (mod$kind == "PowerBasis") return(.nf_scal(.nf_basis_elt(mod, 1L), a))
  if (.nf_starts_with_unity(mod)) return(.nf_scal(.nf_basis_elt(mod, 1L), a))
  .nf_element_from_rational(mod$parent, a)
}
.nf_one <- function(mod) .nf_element_from_rational(mod, 1)
.nf_zero <- function(mod) .nf_scal(.nf_basis_elt(mod, 1L), 0)
.nf_qq_col <- function(el) .bq(el$col) / .bq(el$denom)
.nf_elt_equal <- function(a, b) .nf_mod_equal(a$module, b$module) && all(.nf_qq_col(a) == .nf_qq_col(b))
.nf_equiv <- function(a, b) {
  if (.nf_is_rat(b)) {
    if (a$module$kind == "PowerBasis") return(.nf_elt_equal(a, .nf_scal(.nf_basis_elt(a$module, 1L), b)))
    return(.nf_equiv(.nf_over_power_basis(a), b))
  }
  if (.nf_elt_equal(a, b)) return(TRUE)
  u <- .nf_unify(a, b)
  .nf_elt_equal(u[[1L]], u[[2L]])
}
.nf_starts_with_unity <- function(mod) {
  if (mod$kind == "PowerBasis") return(TRUE)
  if (is.null(mod$starts_with_unity)) mod$starts_with_unity <- .nf_equiv(.nf_basis_elt(mod, 1L), 1)
  mod$starts_with_unity
}
# Submodule.represent / PowerBasis.represent: integer coordinates or error
.nf_represent <- function(mod, el) {
  if (.nf_mod_equal(el$module, mod)) return(el$col)
  if (mod$kind == "PowerBasis") stop(structure(class = c("ClosureFailure", "error", "condition"), list(message = "Element not representable in ZZ[theta].", call = NULL)))
  if (.nf_mod_equal(el$module, mod$parent)) {
    rows <- lapply(seq_len(mod$M$nrow), function(i) c(lapply(mod$M$cols, function(v) .bq(v[i]) / .bq(mod$denom)), list(.nf_qq_col(el)[i])))
    K <- .nfK_qq()
    rr <- .fmat_rref(rows, K, mod$n + 1L)
    if ((mod$n + 1L) %in% rr$pivots) stop(structure(class = c("ClosureFailure", "error", "condition"), list(message = "Element outside QQ-span of this basis.", call = NULL)))
    x <- rep(list(.bq(0)), mod$n)
    for (i in seq_along(rr$pivots)) x[[rr$pivots[i]]] <- rr$rows[[i]][[mod$n + 1L]]
    if (any(vapply(x, function(z) gmp::denominator(z) != 1, NA))) stop(structure(class = c("ClosureFailure", "error", "condition"), list(message = "Element in QQ-span but not ZZ-span of this basis.", call = NULL)))
    return(.bz(vapply(x, function(z) as.character(gmp::numerator(z)), "")))
  }
  if (mod$parent$kind == "Submodule") return(.nf_represent(mod, .nf_elt(mod$parent, .nf_represent(mod$parent, el))))
  stop(structure(class = c("ClosureFailure", "error", "condition"), list(message = "Element outside ancestor chain of this module.", call = NULL)))
}
.nf_submodule_from_matrix <- function(mod, M, denom = 1) {
  if (M$nrow != mod$n) stop("Matrix row count must match base module.", call. = FALSE)
  .nf_submodule(mod, M, denom)
}
.nf_submodule_from_gens <- function(mod, gens, hnf = TRUE, hnf_modulus = NULL) {
  d <- if (length(gens) == 1L) gens[[1L]]$denom else Reduce(.nf_ilcm, lapply(gens, `[[`, "denom"))
  B <- .zmat(gens[[1L]]$module$n, lapply(gens, function(g) g$col * gmp::divq.bigz(d, g$denom)))
  if (hnf) B <- .zmat_hnf(B, hnf_modulus)
  .nf_submodule_from_matrix(mod, B, d)
}
# DomainMatrix.eye is sparse; SymPy keeps that format until an HNF replaces it
.nf_whole_submodule <- function(mod) { S <- .nf_submodule_from_matrix(mod, .zmat_eye(mod$n)); S$sparse <- TRUE; S }
.nf_sub_reduced <- function(S) {
  if (S$denom == 1) return(S)
  g <- S$denom; for (v in S$M$cols) for (c in as.list(v)) g <- gmp::gcd.bigz(g, c)
  if (g == 1) return(S)
  R <- .nf_submodule(S$parent, .zmat(S$M$nrow, lapply(S$M$cols, function(v) gmp::divq.bigz(v, g))), gmp::divq.bigz(S$denom, g), S$mult_tab)
  R$sparse <- S$sparse
  R
}
.nf_sub_add <- function(A, B, hnf = TRUE, hnf_modulus = NULL) {
  m <- .nf_ilcm(A$denom, B$denom)
  M <- .zmat_hstack(.zmat_scale(A$M, gmp::divq.bigz(m, A$denom)), .zmat_scale(B$M, gmp::divq.bigz(m, B$denom)))
  if (hnf) M <- .zmat_hnf(M, hnf_modulus)
  .nf_submodule_from_matrix(A$parent, M, m)
}
.nf_sub_mul <- function(S, other, hnf = TRUE, hnf_modulus = NULL) {
  if (.nf_is_rat(other)) {
    nd <- .nf_rat_num_den(other)
    if (nd[[1L]] == 1 && nd[[2L]] == 1) return(S)
    R <- .nf_submodule(S$parent, .zmat_scale(S$M, nd[[1L]]), S$denom * nd[[2L]])
    R$sparse <- S$sparse
    return(.nf_sub_reduced(R))
  }
  if (is.list(other) && !is.null(other$col) && .nf_mod_equal(other$module, S$parent)) {
    gens <- lapply(.nf_basis_elements(S), function(e) .nf_mul(other, .nf_to_parent(e)))
    return(.nf_submodule_from_gens(S$parent, gens, hnf, hnf_modulus))
  }
  alphas <- lapply(.nf_basis_elements(S), .nf_to_parent); betas <- lapply(.nf_basis_elements(other), .nf_to_parent)
  gens <- list()
  for (a in alphas) for (b in betas) gens[[length(gens) + 1L]] <- .nf_mul(a, b)
  .nf_submodule_from_gens(S$parent, gens, hnf, hnf_modulus)
}
.nf_is_sq_maxrank_HNF <- function(S) {
  if (!is.null(S$is_sq_maxrank_HNF)) return(S$is_sq_maxrank_HNF)
  M <- S$M; n <- M$nrow
  ok <- n == .zmat_ncol(M)
  if (ok) for (i in seq_len(n)) {
    if (i > 1L) for (j in seq_len(i - 1L)) if (.zmat_get(M, i, j) != 0) ok <- FALSE
    d <- .zmat_get(M, i, i)
    if (d <= 0) ok <- FALSE
    if (i < n) for (j in (i + 1L):n) { v <- .zmat_get(M, i, j); if (!(v >= 0 && v < d)) ok <- FALSE }
  }
  S$is_sq_maxrank_HNF <- ok
  ok
}
.nf_reduce_element <- function(S, el) {
  if (!.nf_is_sq_maxrank_HNF(S)) stop("Reduction not implemented unless matrix square max-rank HNF", call. = FALSE)
  B <- lapply(.nf_basis_elements(S), .nf_to_parent)
  a <- el
  for (i in rev(seq_len(S$n))) {
    b <- B[[i]]
    q <- gmp::divq.bigz(a$col[i] * b$denom, b$col[i] * a$denom)
    a <- .nf_sub(a, .nf_scal(b, q))
  }
  a
}
.nf_mod_rat <- function(el, m) .nf_reduce_element(.nf_sub_mul(.nf_whole_submodule(el$module), m), el)

# ModuleHomomorphism.matrix and kernel
.nf_hom_matrix <- function(domain, represent, mapping) {
  cols <- lapply(.nf_basis_elements(domain), function(e) represent(mapping(e)))
  .zmat(length(cols[[1L]]), cols)
}
.nf_hom_kernel <- function(domain, represent, mapping, modulus = NULL) {
  M <- .nf_hom_matrix(domain, represent, mapping)
  K <- if (is.null(modulus)) .nfK_qq() else .nfK_ff(modulus)
  rows <- .fmat_from_zmat(M, K)
  ns <- .fmat_nullspace(rows, K, .zmat_ncol(M))
  .nf_submodule_from_matrix(domain, .fmat_cols_to_zmat(ns, K, domain$n))
}
# EndomorphismRing.represent of x -> mult*x on C: the columns of its matrix stacked
.nf_end_represent <- function(C, mult) {
  M <- .nf_hom_matrix(C, function(e) .nf_represent(C, e), function(e) .nf_mul(mult, e))
  do.call(c, M$cols)
}

# PowerBasisElement helpers
.nf_poly_coeffs_from_poly <- function(mod, qcoeffs) {
  # element_from_poly: rational coefficients, leading first, reduced mod T
  n <- mod$n
  r <- up_trim(qcoeffs)
  if (length(r) - 1L >= n) r <- up_divmod(r, .bq(mod$Tc))$r
  if (!length(r) || all(r == 0)) return(.nf_zero(mod))
  d <- .lcm_den(r)
  c <- rev(gmp::numerator(r * .bq(d)))
  col <- .bz(rep(0L, n)); col[seq_along(c)] <- c
  .nf_elt(mod, col, d)
}
.nf_numerator_expr <- function(el, x) {
  cs <- el$col
  .add(lapply(seq_along(cs), function(i) .mul(list(.int(as.character(cs[i])), .pow(x, .int(i - 1L))))))
}
.nf_norm <- function(el, Tc) {
  # T.resultant(numerator) // denom^n; for monic T the resultant is the
  # determinant of multiplication by the numerator on the power basis
  pb <- el$module
  num <- .nf_elt(pb, el$col, 1)
  M <- .zmat(pb$n, lapply(.nf_basis_elements(pb), function(b) .nf_mul(num, b)$col))
  gmp::divq.bigz(.zmat_det(M), el$denom^length(el$col))
}
.nf_disc <- function(cs) {
  x <- Dummy("x")
  .bq(as.character(as_sym(discriminant(up_to_expr(.bq(cs), x), x))$value))
}

# ---------------------------------------------------------------------------
# basis.py

.gf_poly_from_bz <- function(c, p) .gf_trim(as.numeric(gmp::mod.bigz(c, p)))
.nf_dedekind <- function(Tc, p) {
  pn <- as.numeric(p)
  Tb <- .gf_poly_from_bz(Tc, pn)
  fl <- .gf_factor(Tb, pn)
  g <- 1
  for (fk in fl$factors) g <- .gf_mul(g, fk[[1L]], pn)
  h <- .gf_divmod(Tb, g, pn)$q
  gz <- .zp_from_gf_sym(g, pn); hz <- .zp_from_gf_sym(h, pn)
  f <- gmp::divq.bigz(gmp::numerator(up_sub(up_mul(.bq(gz), .bq(hz)), .bq(Tc))), p)
  fb <- .gf_poly_from_bz(f, pn)
  Z <- fb
  for (b in list(g, h)) Z <- .gf_gcd(Z, b, pn)
  U <- .gf_divmod(Tb, Z, pn)$q
  list(U = U, m = length(Z) - 1L)
}
.zp_from_gf_sym <- function(a, p) .bz(vapply(a, function(v) if (v > p %/% 2) as.character(v - p) else as.character(v), ""))

.nf_nilradical_mod_p <- function(H, p, q = NULL) {
  n <- H$n
  if (is.null(q)) { q <- p; while (q < n) q <- q * p }
  .nf_hom_kernel(H, function(e) .nf_represent(H, e), function(e) .nf_pow(e, as.integer(q)), modulus = p)
}
.nf_second_enlargement <- function(H, p, q) {
  Ip <- .nf_nilradical_mod_p(H, p, q)
  B <- .nf_submodule_from_matrix(H$parent, .zmat_mul(H$M, Ip$M), H$denom)
  C <- .nf_sub_add(B, .nf_sub_mul(H, p))
  gamma <- .nf_hom_kernel(H, function(v) v, function(e) .nf_end_represent(C, e), modulus = p)
  G <- .nf_submodule_from_matrix(H$parent, .zmat_mul(H$M, gamma$M), H$denom * p)
  list(.nf_sub_add(G, H), Ip)
}


# round_two(T): ZK as a Submodule of the power basis and the discriminant
.nf_round_two <- function(Tc, gen, radicals = NULL) {
  n <- length(Tc) - 1L
  D <- gmp::numerator(.nf_disc(Tc))
  ef <- .extract_fundamental_discriminant(D)
  F <- ef$F
  Ztheta <- .nf_power_basis(Tc, gen)
  H <- .nf_whole_submodule(Ztheta)
  nilrad <- NULL; plast <- NULL
  while (length(F)) {
    pn <- names(F)[length(F)]; e <- F[[pn]]; F[[pn]] <- NULL
    p <- .bz(pn); plast <- p
    dd <- .nf_dedekind(Tc, p)
    if (dd$m == 0L) next
    Uc <- .zp_from_gf_sym(dd$U, as.numeric(p))
    U <- .nf_poly_coeffs_from_poly(Ztheta, .bq(Uc))
    G <- .nf_sub_mul(H, .nf_scal(U, .bq(1) / .bq(p)))
    d <- .nf_ilcm(H$denom, G$denom)
    Dm <- gmp::divq.bigz(d, H$denom)^n
    for (i in seq_len(n)) Dm <- Dm * .zmat_get(H$M, i, i)
    H <- .nf_sub_add(H, G, hnf_modulus = Dm)
    if (e <= dd$m) next
    q <- p; while (q < n) q <- q * p
    r <- .nf_second_enlargement(H, p, q)
    H1 <- r[[1L]]; nilrad <- r[[2L]]
    while (!.nf_mod_equal(H1, H)) {
      H <- H1
      r <- .nf_second_enlargement(H, p, q); H1 <- r[[1L]]; nilrad <- r[[2L]]
    }
  }
  if (!is.null(nilrad) && is.environment(radicals)) radicals[[as.character(plast)]] <- nilrad
  ZK <- H
  ZK$starts_with_unity <- TRUE; ZK$is_sq_maxrank_HNF <- TRUE
  dK <- gmp::divq.bigz(D * .zmat_det(ZK$M)^2, ZK$denom^(2L * n))
  list(ZK = ZK, dK = dK, D = D)
}

# ---------------------------------------------------------------------------
# primes.py

.nf_check_maximal_order <- function(ZK) {
  pre <- "The submodule representing the maximal order should "
  cond <- NULL
  if (!(ZK$kind == "Submodule" && ZK$parent$kind == "PowerBasis")) cond <- "be a direct submodule of a power basis."
  else if (!.nf_starts_with_unity(ZK)) cond <- "have 1 as its first generator."
  else if (!.nf_is_sq_maxrank_HNF(ZK)) cond <- "have square matrix, of maximal rank, in Hermite Normal Form."
  if (!is.null(cond)) stop(structure(class = c("StructureError", "error", "condition"), list(message = paste0(pre, cond), call = NULL)))
}
.nf_prime_ideal <- function(ZK, p, alpha, f, e = NULL) {
  .nf_check_maximal_order(ZK)
  P <- new.env()
  P$ZK <- ZK; P$p <- .bz(p); P$alpha <- alpha; P$f <- as.integer(f); P$test_factor <- NULL
  class(P) <- "PrimeIdeal"
  P$e <- if (!is.null(e)) as.integer(e) else .nf_prime_valuation(.nf_sub_mul(ZK, P$p), P)
  P
}
.nf_prime_ideal_repr <- function(P, field_gen = NULL) {
  x <- field_gen %||% P$ZK$parent$gen
  a <- P$alpha
  ar <- sstr(.nf_numerator_expr(a, x))
  if (a$denom > 1) ar <- paste0("(", ar, ")/", as.character(a$denom))
  paste0("[ (", as.character(P$p), ", ", ar, ") e=", P$e, ", f=", P$f, " ]")
}
#' @export
format.PrimeIdeal <- function(x, ...) .nf_prime_ideal_repr(x)
#' @export
print.PrimeIdeal <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

.nf_compute_test_factor <- function(p, gens, ZK) {
  .nf_check_maximal_order(ZK)
  K <- .nfK_ff(p)
  rows <- list()
  for (g in gens) {
    M <- .nf_hom_matrix(ZK, function(e) .nf_represent(ZK, e), function(e) .nf_mul(g, e))
    rows <- c(rows, .fmat_from_zmat(M, K))
  }
  ns <- .fmat_nullspace(rows, K, ZK$n)
  x <- .bz(vapply(ns[[1L]], function(z) as.character(K$to_int(z)), ""))
  .nf_elt(ZK$parent, .zmat_mul(ZK$M, .zmat(ZK$n, list(x)))$cols[[1L]], ZK$denom)
}
.nf_qq_inv <- function(M) {
  # inverse of a square integer matrix over QQ, as columns of bigq
  n <- M$nrow; K <- .nfK_qq()
  rows <- lapply(seq_len(n), function(i) c(lapply(M$cols, function(v) .bq(v[i])), lapply(seq_len(n), function(j) if (i == j) .bq(1) else .bq(0))))
  rr <- .fmat_rref(rows, K, 2L * n)
  lapply(seq_len(n), function(j) .bq(vapply(seq_len(n), function(i) as.character(rr$rows[[i]][[n + j]]), "")))
}
# DomainMatrix.convert_to(ZZ) fails at the first non-integer entry: row by
# row for a dense matrix, in the dict order of rows and entries for a sparse
# one (ord = list(rows, keys))
.nf_coercion_failed <- function(cols, ord = NULL) {
  msg <- function(z) stop(.sympy_condition("CoercionFailed", sprintf("Cannot convert %s of type <class 'sympy.external.pythonmpq.PythonMPQ'> from QQ to ZZ", as.character(z))))
  if (is.null(ord)) {
    for (i in seq_along(cols[[1L]])) for (v in cols) if (gmp::denominator(v[i]) != 1) msg(v[i])
  } else {
    for (i in ord$rows) for (j in ord$keys[[i]]) if (gmp::denominator(cols[[j]][i]) != 1) msg(cols[[j]][i])
  }
}
# SDM.setitem over the entries of a new matrix, column by column as
# prime_valuation assigns them: existing keys keep their place, new keys go
# last, zeros are removed, and an emptied row leaves the row order
.nf_sdm_assign <- function(ord, A) {
  n <- A$nrow
  for (j in seq_len(.zmat_ncol(A))) for (i in seq_len(n)) {
    v <- A$cols[[j]][i]
    has <- j %in% ord$keys[[i]]
    if (v != 0) {
      if (!has) {
        if (!(i %in% ord$rows)) ord$rows <- c(ord$rows, i)
        ord$keys[[i]] <- c(ord$keys[[i]], j)
      }
    } else if (has) {
      ord$keys[[i]] <- setdiff(ord$keys[[i]], j)
      if (!length(ord$keys[[i]])) ord$rows <- setdiff(ord$rows, i)
    }
  }
  ord
}
.nf_assert <- function(cond) if (!isTRUE(cond)) stop(.sympy_condition("AssertionError", ""))
# prime_valuation(I, P)
.nf_prime_valuation <- function(I, P) {
  p <- P$p; ZK <- P$ZK
  n <- ZK$n; W <- ZK$M; d <- ZK$denom
  Wi <- .nf_qq_inv(W)
  # A = W^-1 * I.matrix * d / I.denom, over QQ then ZZ
  A <- lapply(I$M$cols, function(v) {
    out <- .bq(rep(0L, n))
    for (k in seq_along(v)) if (v[k] != 0) out <- out + Wi[[k]] * .bq(v[k])
    out * .bq(d) / .bq(I$denom)
  })
  # with W the sparse identity, A starts as the sparse diagonal of I
  sp <- isTRUE(ZK$sparse) && isTRUE(I$sparse)
  ord <- if (sp) list(rows = which(vapply(seq_len(n), function(i) A[[i]][i] != 0, NA)),
                      keys = lapply(seq_len(n), function(i) if (A[[i]][i] != 0) i else integer())) else NULL
  if (any(vapply(A, function(v) any(gmp::denominator(v) != 1), NA))) .nf_coercion_failed(A, ord)
  A <- .zmat(n, lapply(A, gmp::numerator))
  Dt <- .zmat_det(A)
  if (gmp::mod.bigz(Dt, p) != 0) return(0L)
  if (is.null(P$test_factor)) P$test_factor <- .nf_compute_test_factor(p, list(P$alpha), ZK)
  beta <- P$test_factor
  f <- gmp::divq.bigz(d^n, .zmat_det(W))
  need_complete <- gmp::mod.bigz(f, p) == 0
  v <- 0L
  repeat {
    A <- .zmat_mul(W, A)
    A <- .zmat(n, lapply(seq_len(n), function(j) {
      c <- .nf_mul(.nf_elt(ZK$parent, A$cols[[j]], d), beta)
      .nf_represent(ZK, c)
    }))
    if (sp) ord <- .nf_sdm_assign(ord, A)
    if (gmp::mod.bigz(.zmat_get(A, n, n), p) != 0) break
    Aq <- lapply(A$cols, function(v) .bq(v) / .bq(p))
    if (any(vapply(Aq, function(v) any(gmp::denominator(v) != 1), NA))) {
      if (need_complete) break
      .nf_coercion_failed(Aq, ord)
    }
    A <- .zmat(n, lapply(Aq, gmp::numerator))
    v <- v + 1L
  }
  v
}
# _two_elt_rep
.nf_two_elt_rep <- function(gens, ZK, p, f = NULL, Np = NULL) {
  .nf_check_maximal_order(ZK)
  pb <- ZK$parent
  if (all(vapply(gens, function(g) .nf_equiv(.nf_mod_rat(g, p), 0), NA))) return(.nf_zero(pb))
  if (is.null(Np)) Np <- if (!is.null(f)) p^f else abs(.zmat_det(.nf_submodule_from_gens(pb, gens)$M))
  omega <- lapply(.nf_basis_elements(ZK), .nf_to_parent)
  beta <- c(lapply(omega[-1L], function(om) .nf_scal(om, p)), gens)
  nxt <- .coeff_search(length(beta), 1L)
  repeat {
    c <- nxt()
    alpha <- NULL
    for (i in seq_along(beta)) { t <- .nf_scal(beta[[i]], c[i]); alpha <- if (is.null(alpha)) t else .nf_add(alpha, t) }
    nn <- gmp::divq.bigz(.nf_norm(alpha, pb$Tc), Np)
    if (gmp::mod.bigz(nn, p) != 0) return(.nf_mod_rat(alpha, p))
  }
}
.nf_gf_factor_list <- function(Tc, p) {
  pn <- as.numeric(p)
  fl <- .gf_factor(.gf_poly_from_bz(Tc, pn), pn)
  lapply(fl$factors, function(fk) list(fk[[1L]], fk[[2L]]))
}
.nf_prime_decomp_easy <- function(p, ZK) {
  Tc <- ZK$parent$Tc
  fl <- .nf_gf_factor_list(Tc, p)
  if (length(fl) == 1L && fl[[1L]][[2L]] == 1L) return(list(.nf_prime_ideal(ZK, p, .nf_zero(ZK$parent), ZK$n, 1L)))
  lapply(fl, function(t) .nf_prime_ideal(ZK, p, .nf_poly_coeffs_from_poly(ZK$parent, .bq(.zp_from_gf_sym(t[[1L]], as.numeric(p)))), length(t[[1L]]) - 1L, t[[2L]]))
}
# supplement_a_subspace over GF(p): an invertible matrix whose first columns are M
.nf_supplement <- function(M, p) {
  K <- .nfK_ff(p); n <- M$nrow; r <- .zmat_ncol(M)
  rows <- lapply(seq_len(n), function(i) c(lapply(M$cols, function(v) K$from_int(v[i])), lapply(seq_len(n), function(j) if (i == j) K$one else K$zero)))
  rr <- .fmat_rref(rows, K, r + n)
  if (!identical(rr$pivots[seq_len(r)], seq_len(r))) stop("M was not of maximal rank", call. = FALSE)
  Arows <- lapply(rr$rows, function(row) row[(r + 1L):(r + n)])
  # B = A^-1
  aug <- lapply(seq_len(n), function(i) c(Arows[[i]], lapply(seq_len(n), function(j) if (i == j) K$one else K$zero)))
  ri <- .fmat_rref(aug, K, 2L * n)
  .zmat(n, lapply(seq_len(n), function(j) .bz(vapply(seq_len(n), function(i) as.character(K$to_int(ri$rows[[i]][[n + j]])), ""))))
}
.nf_prime_decomp_compute_kernel <- function(I, p, ZK) {
  W <- I$M; n <- W$nrow; r <- .zmat_ncol(W)
  B <- if (r == 0L) .zmat_eye(n) else .zmat_hstack(W, .zmat(n, list(.zmat_eye(n)$cols[[1L]])))
  if (.zmat_ncol(B) < n) B <- .nf_supplement(B, p)
  G <- .nf_submodule_from_matrix(ZK, B)
  .nf_mult_tab(G)
  # discard_before(r)
  s <- G$n - r
  mt <- G$mult_tab
  Mnew <- lapply(seq_len(s), function(u) { l <- vector("list", s); for (v in u:s) l[[v]] <- mt[[r + u]][[r + v]][(r + 1L):G$n]; l })
  G <- .nf_submodule(ZK, .zmat(B$nrow, B$cols[(r + 1L):G$n]), 1, Mnew)
  N <- .nf_hom_kernel(G, function(e) .nf_represent(G, e), function(e) .nf_sub(.nf_pow(e, as.integer(p)), e), modulus = p)
  .nf_assert(.nf_starts_with_unity(N))
  list(N, G)
}
.nf_prime_decomp_maximal_ideal <- function(I, p, ZK) {
  m <- I$M$nrow; n <- .zmat_ncol(I$M)
  f <- m - n
  G <- .zmat_mul(ZK$M, I$M)
  gens <- lapply(G$cols, function(v) .nf_elt(ZK$parent, v, ZK$denom))
  alpha <- .nf_two_elt_rep(gens, ZK, p, f = f)
  .nf_prime_ideal(ZK, p, alpha, f)
}
# find_min_poly over GF(p): coefficients (leading first) and the powers of alpha
.nf_find_min_poly_ff <- function(alpha, p) {
  R <- alpha$module
  K <- .nfK_ff(p)
  one <- .nf_basis_elt(R, 1L)
  powers <- list(one)
  cols <- list(lapply(as.list(one$col), K$from_int))
  ak <- alpha
  for (k in seq_len(R$n)) {
    powers[[length(powers) + 1L]] <- ak
    akc <- lapply(as.list(ak$col), K$from_int)
    rows <- lapply(seq_len(R$n), function(i) c(lapply(cols, function(cc) cc[[i]]), list(akc[[i]])))
    rr <- .fmat_rref(rows, K, length(cols) + 1L)
    if ((length(cols) + 1L) %in% rr$pivots) { cols[[length(cols) + 1L]] <- akc; ak <- .nf_mul(ak, alpha); next }
    X <- rep(list(K$zero), length(cols))
    for (i in seq_along(rr$pivots)) X[[rr$pivots[i]]] <- rr$rows[[i]][[length(cols) + 1L]]
    coeffs <- c(list(K$one), lapply(rev(X), K$neg))
    return(list(coeffs = vapply(coeffs, function(z) as.numeric(z), 0), powers = powers))
  }
  NULL
}
.nf_prime_decomp_split_ideal <- function(I, p, N, G, ZK) {
  pn <- as.numeric(p)
  alpha <- .nf_to_parent(.nf_basis_elt(N, 2L))
  mp <- .nf_find_min_poly_ff(alpha, p)
  m <- mp$coeffs; alpha_powers <- mp$powers
  fl <- .gf_factor(m, pn)
  m1 <- fl$factors[[1L]][[1L]]
  m2 <- .gf_divmod(m, m1, pn)$q
  gx <- .gf_gcdex(m1, m2, pn)
  .nf_assert(length(gx$g) == 1L && gx$g[1L] == 1)
  U <- gx$s
  E <- rev(.zp_from_gf_sym(.gf_mul(U, m1, pn), pn))
  eps1 <- NULL
  for (i in seq_along(E)) { t <- .nf_scal(alpha_powers[[i]], E[i]); eps1 <- if (is.null(eps1)) t else .nf_add(eps1, t) }
  eps2 <- .nf_sub(.nf_element_from_rational(eps1$module, 1), eps1)
  lapply(list(eps1, eps2), function(eps) {
    e <- .nf_to_parent(eps)
    K <- .nfK_ff(p)
    cols <- c(lapply(I$M$cols, function(v) lapply(as.list(v), K$from_int)),
              lapply(.nf_basis_elements(ZK), function(om) lapply(as.list(.nf_mul(e, om)$col), K$from_int)))
    rows <- lapply(seq_len(ZK$n), function(i) lapply(cols, function(cc) cc[[i]]))
    rr <- .fmat_rref(rows, K, length(cols))
    W <- .zmat(ZK$n, lapply(rr$pivots, function(j) .bz(vapply(cols[[j]], function(z) as.character(K$to_int(z)), ""))))
    .nf_submodule_from_matrix(ZK, W)
  })
}

# ---------------------------------------------------------------------------
# public interface

.nf_poly_input <- function(T) {
  T <- as_sym(if (inherits(T, "Poly")) as_expr(T) else T)
  fs <- free_symbols(T)
  if (length(fs) != 1L) stop(.sympy_condition("ValueError", "Round 2 requires an irreducible univariate polynomial over ZZ or QQ."))
  x <- fs[[1L]]
  q <- expr_to_qpoly(expand(T), list(x), strict = TRUE)
  cs <- up_from_qpoly(q, 1L)
  list(x = x, cs = cs)
}
.nf_check_irreducible <- function(T, pi) {
  fl <- qp_factor_list(expr_to_qpoly(expand(as_sym(if (inherits(T, "Poly")) as_expr(T) else T)), list(pi$x), strict = TRUE))$factors
  if (length(fl) != 1L || fl[[1L]][[2L]] != 1L)
    stop(.sympy_condition("ValueError", "Round 2 requires an irreducible univariate polynomial over ZZ or QQ."))
}
# make_monic_over_integers_by_scaling_roots
.nf_monic_integer <- function(cs) {
  cs <- .bq(cs)
  lc <- cs[1L]
  cs <- cs / lc
  n <- length(cs) - 1L
  den <- .lcm_den(cs)
  # scale roots by den: coefficient of x^(n-i) times den^i
  if (den == 1) return(gmp::numerator(cs))
  gmp::numerator(.bq(vapply(seq_along(cs), function(i) as.character(cs[i] * .bq(den)^(i - 1L)), "")))
}

#' Integral bases, prime decomposition and valuations in number fields
#'
#' `round_two(T)` computes the maximal order of the number field defined by
#' the irreducible polynomial `T` with the Round Two algorithm, returning
#' the order as a submodule of the power basis and the field discriminant.
#' `prime_decomp(p, T)` factors the rational prime `p` into prime ideals and
#' `prime_valuation(I, P)` is the valuation of an ideal at a prime ideal, as
#' SymPy's functions of the same names.
#'
#' @param T An irreducible polynomial over ZZ or QQ.
#' @param p A rational prime.
#' @param I A submodule (ideal) of the maximal order.
#' @param P A prime ideal from `prime_decomp()`.
#' @return For `round_two()` a list `(ZK, dK)`; for `prime_decomp()` a
#'   list of prime ideals; for `prime_valuation()` an integer.
#' @examples
#' x <- Symbol("x")
#' round_two(x^3 + x^2 - 2*x + 8)
#' prime_decomp(2, x^3 + x^2 - 2*x + 8)
#' @export
round_two <- function(T) {
  pi <- .nf_poly_input(T)
  .nf_check_irreducible(T, pi)
  Tc <- .nf_monic_integer(pi$cs)
  r <- .nf_round_two(Tc, pi$x)
  structure(list(ZK = structure(list(S = r$ZK), class = "nf_Submodule"), dK = .int(as.character(r$dK))), class = "nf_round_two")
}
#' @export
format.nf_Submodule <- function(x, ...) {
  S <- x$S
  r <- paste0("Submodule[", paste(vapply(S$M$cols, function(v) paste0("[", paste(as.character(v), collapse = ", "), "]"), ""), collapse = ", "), "]")
  if (S$denom > 1) r <- paste0(r, "/", as.character(S$denom))
  r
}
#' @export
print.nf_Submodule <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
#' @export
format.nf_round_two <- function(x, ...) paste0("(", format(x$ZK), ", ", sstr(x$dK), ")")
#' @export
print.nf_round_two <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

#' @rdname round_two
#' @export
prime_decomp <- function(p, T) {
  pi <- .nf_poly_input(T)
  .nf_check_irreducible(T, pi)
  Tc <- .nf_monic_integer(pi$cs)
  p <- .bz(as.character(as_sym(p)$value))
  radicals <- new.env()
  r <- .nf_round_two(Tc, pi$x, radicals)
  ZK <- r$ZK; dK <- r$dK
  dT <- .nf_disc(pi$cs)
  fsq <- floor(dT / .bq(dK))
  if (gmp::numerator(fsq - .bq(p) * floor(fsq / .bq(p))) != 0) return(structure(.nf_prime_decomp_easy(p, ZK), class = "nf_prime_list"))
  radical <- radicals[[as.character(p)]] %||% .nf_nilradical_mod_p(ZK, p)
  stack <- list(radical); primes <- list()
  while (length(stack)) {
    I <- stack[[length(stack)]]; stack <- stack[-length(stack)]
    kg <- .nf_prime_decomp_compute_kernel(I, p, ZK)
    N <- kg[[1L]]; G <- kg[[2L]]
    if (N$n == 1L) primes[[length(primes) + 1L]] <- .nf_prime_decomp_maximal_ideal(I, p, ZK)
    else stack <- c(stack, .nf_prime_decomp_split_ideal(I, p, N, G, ZK))
  }
  structure(primes, class = "nf_prime_list")
}
#' @export
format.nf_prime_list <- function(x, ...) paste0("[", paste(vapply(unclass(x), format, ""), collapse = ", "), "]")
#' @export
print.nf_prime_list <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

#' @rdname round_two
#' @export
prime_valuation <- function(I, P) {
  S <- if (inherits(I, "nf_Submodule")) I$S else I
  .nf_prime_valuation(S, P)
}
