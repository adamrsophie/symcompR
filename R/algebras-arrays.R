# Quaternions, N-dimensional arrays and discrete transforms,
# after sympy/algebras, sympy/tensor/array and sympy/discrete.

#' Quaternions
#'
#' `Quaternion(a, b, c, d)` is `a + b*i + c*j + d*k`. `+`, `-`, `*` (Hamilton
#' product) and `/` by a scalar work; `quat_conjugate()`, `quat_norm()`,
#' `quat_inverse()` and `quat_to_rotation_matrix()` follow SymPy.
#'
#' @param a,b,c,d Components.
#' @param q A quaternion.
#' @return A `Quaternion`, expression or matrix.
#' @examples
#' Quaternion(1, 2, 3, 4) * Quaternion(5, 6, 7, 8)
#' @export
Quaternion <- function(a = 0, b = 0, c = 0, d = 0) {
  structure(list(c = lapply(list(a, b, c, d), as_sym)), class = "Quaternion")
}

#' @export
format.Quaternion <- function(x, ...) {
  paste0(.str(x$c[[1L]]), " + ", .str(x$c[[2L]]), "*i + ", .str(x$c[[3L]]), "*j + ", .str(x$c[[4L]]), "*k")
}

#' @export
print.Quaternion <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

.q_lift <- function(z) if (inherits(z, "Quaternion")) z else Quaternion(z, 0, 0, 0)

#' @export
Ops.Quaternion <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(structure(list(c = lapply(e1$c, .neg)), class = "Quaternion"))
    return(e1)
  }
  a <- .q_lift(e1)$c
  b <- .q_lift(e2)$c
  ex <- function(z) expand(z)
  switch(.Generic,
         "+" = do.call(Quaternion, Map(function(u, v) .add(list(u, v)), a, b)),
         "-" = do.call(Quaternion, Map(.sub, a, b)),
         "*" = {
           p <- function(i, j, sgn = 1L) {
             t <- .mul(list(a[[i]], b[[j]]))
             if (sgn < 0L) .neg(t) else t
           }
           Quaternion(ex(.add(list(p(1, 1), p(2, 2, -1L), p(3, 3, -1L), p(4, 4, -1L)))),
                      ex(.add(list(p(1, 2), p(2, 1), p(3, 4), p(4, 3, -1L)))),
                      ex(.add(list(p(1, 3), p(2, 4, -1L), p(3, 1), p(4, 2)))),
                      ex(.add(list(p(1, 4), p(2, 3), p(3, 2, -1L), p(4, 1)))))
         },
         "/" = {
           if (inherits(e2, "Quaternion")) return(e1 * quat_inverse(e2))
           do.call(Quaternion, lapply(a, function(u) .div(u, as_sym(e2))))
         },
         "==" = all(vapply(seq_len(4L), function(i) identical(expand(a[[i]])$key, expand(b[[i]])$key), TRUE)),
         "!=" = !all(vapply(seq_len(4L), function(i) identical(expand(a[[i]])$key, expand(b[[i]])$key), TRUE)),
         stop("unsupported quaternion operator"))
}

#' @rdname Quaternion
#' @export
quat_conjugate <- function(q) Quaternion(q$c[[1L]], .neg(q$c[[2L]]), .neg(q$c[[3L]]), .neg(q$c[[4L]]))

.q_norm2 <- function(q) .add(lapply(q$c, .sq))

#' @rdname Quaternion
#' @export
quat_norm <- function(q) .sqrt(.q_norm2(q))

#' @rdname Quaternion
#' @export
quat_inverse <- function(q) {
  n <- .q_norm2(q)
  if (is_zero_num(simplify(n))) stop("cannot invert a zero quaternion")
  do.call(Quaternion, lapply(quat_conjugate(q)$c, function(u) .div(u, n)))
}

#' @rdname Quaternion
#' @export
quat_to_rotation_matrix <- function(q) {
  n <- .q_norm2(q)
  s <- .inv(n)
  a <- q$c[[1]]; b <- q$c[[2]]; c <- q$c[[3]]; d <- q$c[[4]]
  m <- function(...) simplify(.mul(list(S_Two, s, .add(list(...)))))
  Matrix(list(
    list(simplify(.sub(S_One, m(.sq(c), .sq(d)))), m(.mul(list(b, c)), .neg(.mul(list(d, a)))), m(.mul(list(b, d)), .mul(list(c, a)))),
    list(m(.mul(list(b, c)), .mul(list(d, a))), simplify(.sub(S_One, m(.sq(b), .sq(d)))), m(.mul(list(c, d)), .neg(.mul(list(b, a))))),
    list(m(.mul(list(b, d)), .neg(.mul(list(c, a)))), m(.mul(list(c, d)), .mul(list(b, a))), simplify(.sub(S_One, m(.sq(b), .sq(c)))))))
}

# ---------------------------------------------------------------------------
# N-dimensional arrays

#' N-dimensional symbolic arrays
#'
#' `Array(data, shape)` stores expressions in row-major order.
#' `tensorproduct()`, `tensorcontraction()`, `permutedims()` and
#' `derive_by_array()` follow `sympy.tensor.array`. `arr_get()` reads one
#' entry using 1 based indices.
#'
#' @param data A list of entries (or a nested list with `shape` omitted).
#' @param shape An integer vector.
#' @param a,x Arrays or expressions.
#' @param axes A pair of 1 based axes to contract, or a permutation.
#' @param idx Indices.
#' @param ... Arrays; for `tensordiagonal()`, groups of 1 based axes to diagonalize.
#' @return A `SymArray`, or an expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' derive_by_array(x^2 * y, list(x, y))
#' tensorcontraction(Array(list(1, 2, 3, 4), c(2, 2)), c(1, 2))
#' @export
Array <- function(data, shape = NULL) {
  if (is.null(shape)) {
    sh <- integer()
    z <- data
    if (!is.list(data) && length(data) > 1L) sh <- length(data)
    while (is.list(z) && !inherits(z, "sym")) {
      sh <- c(sh, length(z))
      z <- z[[1L]]
    }
    shape <- sh
    flat <- function(z) {
      if (inherits(z, "sym") || !is.list(z)) return(if (inherits(z, "sym")) list(z) else as.list(z))
      do.call(c, lapply(z, flat))
    }
    data <- flat(data)
  }
  structure(list(data = lapply(as.list(data), as_sym), shape = as.integer(shape)), class = "SymArray")
}

#' @export
print.SymArray <- function(x, ...) {
  cat("Array of shape (", paste(x$shape, collapse = ", "), "):\n", sep = "")
  cat(vapply(x$data, .str, ""), sep = ", ")
  cat("\n")
  invisible(x)
}

#' @export
dim.SymArray <- function(x) x$shape

.arr_index <- function(shape, idx) {
  k <- 0L
  for (d in seq_along(shape)) k <- k * shape[d] + (idx[d] - 1L)
  k + 1L
}

.arr_indices <- function(shape) {
  g <- expand.grid(rev(lapply(shape, seq_len)))
  as.matrix(g[, rev(seq_along(shape)), drop = FALSE])
}

#' @rdname Array
#' @export
arr_get <- function(a, idx) a$data[[.arr_index(a$shape, idx)]]

#' @rdname Array
#' @export
tensorproduct <- function(...) {
  arrs <- list(...)
  Reduce(function(a, b) {
    if (!inherits(a, "SymArray")) a <- Array(list(a), integer())
    if (!inherits(b, "SymArray")) b <- Array(list(b), integer())
    d <- vector("list", length(a$data) * length(b$data))
    k <- 0L
    for (u in a$data) for (v in b$data) {
      k <- k + 1L
      d[[k]] <- .mul(list(u, v))
    }
    Array(d, c(a$shape, b$shape))
  }, arrs)
}

#' @rdname Array
#' @export
tensorcontraction <- function(a, axes) {
  i <- axes[1L]; j <- axes[2L]
  if (a$shape[i] != a$shape[j]) stop("contracted axes must have equal length")
  rest <- setdiff(seq_along(a$shape), c(i, j))
  sh <- a$shape[rest]
  if (!length(sh)) {
    return(.add(lapply(seq_len(a$shape[i]), function(t) {
      idx <- integer(length(a$shape))
      idx[c(i, j)] <- t
      arr_get(a, idx)
    })))
  }
  ind <- .arr_indices(sh)
  out <- vector("list", nrow(ind))
  for (r in seq_len(nrow(ind))) {
    out[[r]] <- .add(lapply(seq_len(a$shape[i]), function(t) {
      idx <- integer(length(a$shape))
      idx[rest] <- ind[r, ]
      idx[c(i, j)] <- t
      arr_get(a, idx)
    }))
  }
  Array(out, sh)
}

#' @rdname Array
#' @export
permutedims <- function(a, axes) {
  sh <- a$shape[axes]
  ind <- .arr_indices(sh)
  out <- vector("list", nrow(ind))
  for (r in seq_len(nrow(ind))) {
    idx <- integer(length(sh))
    idx[axes] <- ind[r, ]
    out[[r]] <- arr_get(a, idx)
  }
  Array(out, sh)
}

#' @rdname Array
#' @export
derive_by_array <- function(a, x) {
  if (!inherits(x, "SymArray") && !inherits(a, "SymArray") && (inherits(x, "sym") || !is.list(x)))
    return(.diff1(as_sym(a), as_sym(x)))
  vars <- if (inherits(x, "SymArray")) x$data else if (is.list(x) && !inherits(x, "sym")) lapply(x, as_sym) else list(as_sym(x))
  vshape <- if (inherits(x, "SymArray")) x$shape else if (inherits(x, "sym")) integer() else length(vars)
  items <- if (inherits(a, "SymArray")) a$data else list(as_sym(a))
  ashape <- if (inherits(a, "SymArray")) a$shape else integer()
  d <- list()
  for (v in vars) for (e in items) d[[length(d) + 1L]] <- .diff1(e, v)
  Array(d, c(vshape, ashape))
}

# ---------------------------------------------------------------------------
# Discrete transforms

#' Discrete transforms and convolutions
#'
#' Exact transforms over symbolic sequences: `sym_fft()` and `sym_ifft()`
#' (discrete Fourier transform with exact roots of unity),
#' `convolution()` (linear or `cycle`-length cyclic), `mobius_transform()`
#' and `inverse_mobius_transform()` over subset lattices, and
#' `walsh_hadamard_transform()`.
#'
#' @param seq A list or vector of values.
#' @param a,b Sequences.
#' @param cycle Length of a cyclic convolution (0 for linear).
#' @param inverse Compute the inverse transform.
#' @param subset Sum over subsets (`TRUE`) or supersets.
#' @return A list of expressions.
#' @examples
#' convolution(1:3, 4:5)
#' sym_fft(list(1, 2, 3, 4))
#' @export
convolution <- function(a, b, cycle = 0) {
  a <- lapply(as.list(a), as_sym)
  b <- lapply(as.list(b), as_sym)
  n <- length(a) + length(b) - 1L
  out <- rep(list(S_Zero), n)
  for (i in seq_along(a)) for (j in seq_along(b)) {
    k <- i + j - 1L
    out[[k]] <- .add(list(out[[k]], .mul(list(a[[i]], b[[j]]))))
  }
  if (cycle > 0L) {
    res <- rep(list(S_Zero), cycle)
    for (k in seq_along(out)) {
      m <- (k - 1L) %% cycle + 1L
      res[[m]] <- .add(list(res[[m]], out[[k]]))
    }
    out <- res
  }
  out
}

.pad2 <- function(seq) {
  seq <- lapply(as.list(seq), as_sym)
  n <- length(seq)
  m <- 2L^ceiling(log2(max(n, 1L)))
  c(seq, rep(list(S_Zero), m - n))
}

.dft <- function(seq, sgn) {
  n <- length(seq)
  w <- function(k) {
    ang <- .mul(list(.rat(2L * k * sgn, n), S_Pi))
    c1 <- fcall("cos", list(ang)); s1 <- fcall("sin", list(ang))
    if (c1$type == "Function" || s1$type == "Function") return(fcall("exp", list(.mul(list(S_I, ang)))))
    .add(list(c1, .mul(list(S_I, s1))))
  }
  lapply(0:(n - 1L), function(k) expand(.add(lapply(0:(n - 1L), function(j) .mul(list(seq[[j + 1L]], w((j * k) %% n)))))))
}

#' @rdname convolution
#' @export
sym_fft <- function(seq) .dft(.pad2(seq), 1L)

#' @rdname convolution
#' @export
sym_ifft <- function(seq) {
  s <- .pad2(seq)
  lapply(.dft(s, -1L), function(z) expand(.div(z, .int(length(s)))))
}

.lattice <- function(seq, subset, inverse) {
  a <- .pad2(seq)
  n <- length(a)
  sgn <- if (inverse) S_NegativeOne else S_One
  h <- 1L
  while (h < n) {
    for (i in 0:(n - 1L)) {
      if (bitwAnd(i, h) != 0L) {
        lo <- i - h + 1L; hi <- i + 1L
        if (subset) a[[hi]] <- .add(list(a[[hi]], .mul(list(sgn, a[[lo]]))))
        else a[[lo]] <- .add(list(a[[lo]], .mul(list(sgn, a[[hi]]))))
      }
    }
    h <- h * 2L
  }
  a
}

#' @rdname convolution
#' @export
mobius_transform <- function(seq, subset = TRUE) .lattice(seq, subset, FALSE)

#' @rdname convolution
#' @export
inverse_mobius_transform <- function(seq, subset = TRUE) .lattice(seq, subset, TRUE)

#' @rdname convolution
#' @export
walsh_hadamard_transform <- function(seq, inverse = FALSE) {
  a <- .pad2(seq)
  n <- length(a)
  h <- 1L
  while (h < n) {
    for (i in seq(0L, n - 1L, by = 2L * h)) for (j in i:(i + h - 1L)) {
      u <- a[[j + 1L]]; v <- a[[j + h + 1L]]
      a[[j + 1L]] <- .add(list(u, v))
      a[[j + h + 1L]] <- .sub(u, v)
    }
    h <- h * 2L
  }
  if (inverse) a <- lapply(a, function(z) .div(z, .int(n)))
  a
}

#' @rdname Array
#' @export
tensordiagonal <- function(a, ...) {
  groups <- lapply(list(...), as.integer)
  if (!length(groups)) return(a)
  if (any(vapply(groups, length, 0L) < 2L)) stop("need at least two axes to diagonalize", call. = FALSE)
  used <- unlist(groups)
  if (anyDuplicated(used) || any(used < 1L | used > length(a$shape))) stop("invalid axes", call. = FALSE)
  for (g in groups) if (length(unique(a$shape[g])) != 1L) stop("dimensions of diagonal axes must be equal", call. = FALSE)
  rest <- setdiff(seq_along(a$shape), used)
  # SymPy appends the diagonals, one axis per group, after the remaining axes
  sh <- c(a$shape[rest], vapply(groups, function(g) a$shape[g[1L]], 0L))
  ind <- .arr_indices(sh)
  out <- lapply(seq_len(nrow(ind)), function(r) {
    idx <- integer(length(a$shape))
    idx[rest] <- ind[r, seq_along(rest)]
    for (k in seq_along(groups)) idx[groups[[k]]] <- ind[r, length(rest) + k]
    arr_get(a, idx)
  })
  if (!length(sh)) return(out[[1L]])
  Array(out, sh)
}
