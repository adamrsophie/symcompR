# Simple Lie algebras, after sympy/liealgebras: Cartan types A to G, root
# systems and Weyl groups. Roots are lists of SymPy numbers (the
# coordinates SymPy uses, including its Floats in E's first simple root);
# positive roots are named lists keyed 1, 2, ... as SymPy's dicts.

#' Lie algebras
#'
#' `CartanType("A3")` (or `CartanType(list("A", 3))`) is a simple Lie algebra
#' of series A to G. `cartan_rank()`, `dimension()`, `basic_root()`,
#' `simple_root()`, `positive_roots()`, `highest_root()` (series A),
#' `roots()` (the number of roots), `cartan_matrix()`, `basis()`,
#' `lie_algebra()` and `dynkin_diagram()` are SymPy's methods.
#' `CartanMatrix(t)` and `DynkinDiagram(t)` take a type. `RootSystem(t)` has
#' `simple_roots()`, `all_roots()`, `root_space()`, `add_simple_roots()` and
#' `add_as_roots()`. `WeylGroup(t)` has `generators()`, `group_order()`,
#' `group_name()`, `element_order()`, `delete_doubles()`, `matrix_form()` and
#' `coxeter_diagram()`; Weyl group elements are strings such as `"r1*r2"`.
#'
#' @param ct,t A Cartan type, or a string such as `"B4"`.
#' @param x A Cartan type, root system or Weyl group.
#' @param i,j Indices from 1 (`simple_root`) or 0 (`basic_root`).
#' @param root1,root2 Roots or simple root indices.
#' @param weylelt A Weyl group element as a string.
#' @param reflections A vector of reflection numbers.
#' @return A Cartan type, root system or Weyl group, a list of roots, a
#'   matrix, a string or a number.
#' @examples
#' ct <- CartanType("B3")
#' simple_root(ct, 3)
#' cartan_matrix(ct)
#' cat(dynkin_diagram(CartanType("D5")))
#' element_order(WeylGroup("A3"), "r1*r2*r3")
#' @name liealgebras
NULL

#' @rdname liealgebras
#' @export
CartanType <- function(t) {
  if (inherits(t, "CartanType")) return(t)
  if (is.list(t)) { letter <- as.character(t[[1L]]); n <- as.integer(t[[2L]]) }
  else if (is.character(t)) { letter <- substr(t, 1L, 1L); n <- as.integer(substring(t, 2L)) }
  else stop("Argument must be a string (e.g. 'A3') or a list (e.g. ['A', 3])", call. = FALSE)
  if (n < 0) stop("Lie algebra rank cannot be negative", call. = FALSE)
  lim <- list(A = c(1, Inf), B = c(2, Inf), C = c(3, Inf), D = c(3, Inf))
  msg <- list(A = "n cannot be less than 1", B = "n cannot be less than 2", C = "n cannot be less than 3",
              D = "n cannot be less than 3")
  if (letter %in% names(lim)) {
    if (n < lim[[letter]][1L]) stop(msg[[letter]], call. = FALSE)
  } else if (letter == "E") {
    if (n < 6 || n > 8) return(NULL)
  } else if (letter == "F") {
    if (n != 4) return(NULL)
  } else if (letter == "G") {
    if (n != 2) return(NULL)
  } else return(NULL)
  structure(list(series = letter, n = n), class = "CartanType")
}

#' @export
print.CartanType <- function(x, ...) { cat(x$series, x$n, "\n", sep = ""); invisible(x) }

.zr <- function(n) rep(list(S_Zero), n)
.ri <- function(v) lapply(v, as_sym)

#' @rdname liealgebras
#' @export
cartan_rank <- function(x) CartanType(if (inherits(x, c("RootSystem", "WeylGroup"))) x$cartan_type else x)$n
#' @rdname liealgebras
#' @export
dimension <- function(ct) {
  ct <- CartanType(ct)
  switch(EXPR = ct$series, A = ct$n + 1, B = , C = , D = ct$n, E = 8, F = 4, G = 3)
}

#' @rdname liealgebras
#' @export
basic_root <- function(ct, i, j) {
  ct <- CartanType(ct); n <- ct$n
  if (ct$series == "E") { r <- .zr(8); r[[i + 1L]] <- S_NegativeOne; r[[j + 1L]] <- S_One; return(r) }
  if (ct$series == "G") stop("'TypeG' object has no attribute 'basic_root'", call. = FALSE)
  r <- .zr(if (ct$series == "A") n + 1L else n)
  r[[i + 1L]] <- S_One; r[[j + 1L]] <- S_NegativeOne
  r
}

#' @rdname liealgebras
#' @export
simple_root <- function(ct, i) {
  ct <- CartanType(ct); n <- ct$n
  switch(EXPR = ct$series,
    A = basic_root(ct, i - 1L, i),
    B = if (i < n) basic_root(ct, i - 1L, i) else { r <- .zr(n); r[[n]] <- S_One; r },
    C = if (i < n) basic_root(ct, i - 1L, i) else { r <- .zr(n); r[[n]] <- S_Two; r },
    D = if (i < n) basic_root(ct, i - 1L, i) else { r <- .zr(n); r[[n - 1L]] <- S_One; r[[n]] <- S_One; r },
    E = {
      # SymPy uses Python floats here, not Floats
      if (i == 1) { r <- rep(list(-0.5), 8); r[[1L]] <- 0.5; r[[8L]] <- 0.5; r }
      else if (i == 2) { r <- .zr(8); r[[2L]] <- S_One; r[[1L]] <- S_One; r }
      else {
        if (i %in% c(7, 8) && n == 6) stop("E6 only has six simple roots!", call. = FALSE)
        if (i == 8 && n == 7) stop("E7 only has seven simple roots!", call. = FALSE)
        basic_root(ct, i - 3L, i - 2L)
      }
    },
    F = if (i < 3) basic_root(ct, i - 1L, i) else if (i == 3) { r <- .zr(4); r[[4L]] <- S_One; r }
        else if (i == 4) rep(list(Rational(-1, 2)), 4),
    G = if (i == 1) .ri(c(0, 1, -1)) else .ri(c(1, -2, 1)))
}

#' @rdname liealgebras
#' @export
positive_roots <- function(ct) {
  ct <- CartanType(ct); n <- ct$n
  pr <- list()
  add <- function(r) pr[[as.character(length(pr) + 1L)]] <<- r
  pair_roots <- function(m, second) {
    for (i in 0:(m - 2L)) for (j in (i + 1L):(m - 1L)) {
      add(basic_root(ct, i, j))
      r <- basic_root(ct, i, j); r[[j + 1L]] <- second; add(r)
    }
  }
  switch(EXPR = ct$series,
    A = { for (i in 0:(n - 1L)) for (j in (i + 1L):n) add(basic_root(ct, i, j)) },
    B = , F = {
      pair_roots(n, S_One)
      for (i in seq_len(n)) { r <- .zr(n); r[[i]] <- S_One; add(r) }
      if (ct$series == "F") {
        h <- Rational(1, 2); mh <- Rational(-1, 2)
        add(rep(list(h), n))
        for (i in 2:4) { r <- rep(list(h), n); r[[i]] <- mh; add(r) }
        add(list(h, h, mh, mh)); add(list(h, mh, h, mh)); add(list(h, mh, mh, h)); add(list(h, mh, mh, mh))
      }
    },
    C = {
      pair_roots(n, S_One)
      for (i in seq_len(n)) { r <- .zr(n); r[[i]] <- S_Two; add(r) }
    },
    D = pair_roots(n, S_One),
    E = {
      h <- Rational(1, 2); mh <- Rational(-1, 2)
      m <- if (n == 8) 8L else n - 1L
      for (i in 0:(m - 1L)) if (i + 1L <= m - 1L) for (j in (i + 1L):(m - 1L)) {
        add(basic_root(ct, i, j))
        r <- basic_root(ct, i, j); r[[i + 1L]] <- S_One; add(r)
      }
      if (n == 7) add(.ri(c(0, 0, 0, 0, 0, 1, 1, 0)))
      root <- list(h, h, h, h, h, mh, mh, h)
      nb <- n - 1L
      grid <- as.matrix(expand.grid(rep(list(0:1), nb)))[, nb:1, drop = FALSE]
      for (r in seq_len(nrow(grid))) {
        bits <- grid[r, ]
        if (sum(bits) %% 2L) next
        # the root is updated in place, as SymPy's loop does
        for (b in seq_len(min(5L, nb))) if (bits[b] == 1L) root[[b]] <- mh
        if (nb >= 6L && bits[6L] == 1L) root[[6L]] <- h
        if (nb >= 7L && bits[7L] == 1L) root[[7L]] <- h
        add(root)
      }
    },
    G = { pr <- lapply(list(c(0, 1, -1), c(1, -2, 1), c(1, -1, 0), c(1, 0, 1), c(1, 1, -2), c(2, -1, -1)), .ri)
          names(pr) <- as.character(1:6) })
  pr
}

#' @rdname liealgebras
#' @export
highest_root <- function(ct) {
  ct <- CartanType(ct)
  if (ct$series != "A") stop("highest_root is defined for series A", call. = FALSE)
  basic_root(ct, 0L, ct$n)
}

.cartan_num_roots <- function(ct) {
  n <- ct$n
  switch(EXPR = ct$series, A = n * (n + 1), B = , C = 2 * n^2, D = 2 * n * (n - 1),
         E = c(`6` = 72, `7` = 126, `8` = 240)[[as.character(n)]], F = 48, G = 12)
}

#' @rdname liealgebras
#' @export
cartan_matrix <- function(ct) {
  if (inherits(ct, "RootSystem")) ct <- ct$cartan_type
  ct <- CartanType(ct); n <- ct$n
  if (ct$series == "F") return(Matrix(4, 4, f = c(2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2)))
  if (ct$series == "G") return(Matrix(2, 2, f = c(2, -1, -3, 2)))
  m <- diag(2, n)
  set <- function(i, j, v) m[i + 1L, j + 1L] <<- v
  switch(EXPR = ct$series,
    A = { if (n > 2) for (i in 1:(n - 2L)) { set(i, i + 1L, -1); set(i, i - 1L, -1) }
          if (n < 2) stop("index out of range", call. = FALSE)
          set(0, 1, -1); set(n - 1L, n - 2L, -1) },
    B = { if (n > 2) for (i in 1:(n - 2L)) { set(i, i + 1L, -1); set(i, i - 1L, -1) }
          set(0, 1, -1); set(n - 2L, n - 1L, -2); set(n - 1L, n - 2L, -1) },
    C = { if (n > 2) for (i in 1:(n - 2L)) { set(i, i + 1L, -1); set(i, i - 1L, -1) }
          set(0, 1, -1); set(n - 1L, n - 2L, -2) },
    D = { if (n > 3) for (i in 1:(n - 3L)) { set(i, i + 1L, -1); set(i, i - 1L, -1) }
          set(n - 2L, n - 3L, -1); set(n - 3L, n - 1L, -1); set(n - 1L, n - 3L, -1); set(0, 1, -1) },
    E = { if (n > 4) for (i in 3:(n - 2L)) { set(i, i + 1L, -1); set(i, i - 1L, -1) }
          set(0, 2, -1); set(2, 0, -1); set(1, 3, -1); set(3, 1, -1); set(2, 3, -1); set(n - 1L, n - 2L, -1) })
  Matrix(m)
}

#' @rdname liealgebras
#' @export
basis <- function(ct) {
  ct <- CartanType(ct); n <- ct$n
  # SymPy divides with Python's / for B and D, giving floats
  switch(EXPR = ct$series, A = n^2 - 1, B = (n^2 - n) / 2, C = n * (2 * n + 1), D = n * (n - 1) / 2,
         E = c(`6` = 78, `7` = 133, `8` = 248)[[as.character(n)]], F = 52, G = 14)
}

#' @rdname liealgebras
#' @export
lie_algebra <- function(ct) {
  ct <- CartanType(ct); n <- ct$n
  switch(EXPR = ct$series, A = paste0("su(", n + 1, ")"), B = , D = paste0("so(", 2 * n, ")"), C = paste0("sp(", 2 * n, ")"),
         stop("'Type", ct$series, "' object has no attribute 'lie_algebra'", call. = FALSE))
}

#' @rdname liealgebras
#' @export
dynkin_diagram <- function(ct) {
  if (inherits(ct, "RootSystem")) ct <- ct$cartan_type
  ct <- CartanType(ct); n <- ct$n
  nums <- function(r) paste(r, collapse = "   ")
  zeros <- function(k) paste(rep("0", k), collapse = "---")
  sp <- function(k) strrep(" ", k)
  switch(EXPR = ct$series,
    A = paste0(zeros(n), "\n", nums(seq_len(n))),
    B = paste0(zeros(n - 1L), "=>=0\n", nums(seq_len(n))),
    C = paste0(zeros(n - 1L), "=<=0\n", nums(seq_len(n))),
    D = paste0(sp(4 * (n - 3)), n - 1, "\n", sp(4 * (n - 3)), "0\n", sp(4 * (n - 3)), "|\n", sp(4 * (n - 3)), "|\n",
               zeros(n - 1L), "\n", nums(seq_len(n - 2L)), "   ", n),
    E = paste0(sp(8), 2, "\n", sp(8), "0\n", sp(8), "|\n", sp(8), "|\n", zeros(n - 1L), "\n1   ", nums(3:n)),
    F = paste0("0---0=>=0---0\n", nums(1:4)),
    G = "0\u2261<\u22610\n1   2")
}

#' @rdname liealgebras
#' @export
CartanMatrix <- function(t) cartan_matrix(CartanType(t))
#' @rdname liealgebras
#' @export
DynkinDiagram <- function(t) dynkin_diagram(CartanType(t))

# ---------------------------------------------------------------------------
# Root systems

#' @rdname liealgebras
#' @export
RootSystem <- function(t) structure(list(cartan_type = CartanType(t)), class = "RootSystem")

#' @rdname liealgebras
#' @export
simple_roots <- function(x) {
  ct <- x$cartan_type
  stats::setNames(lapply(seq_len(ct$n), function(i) simple_root(ct, i)), as.character(seq_len(ct$n)))
}

.rs_all_roots <- function(x) {
  alpha <- positive_roots(x$cartan_type)
  k <- max(as.integer(names(alpha)))
  for (val in names(alpha)[seq_along(alpha)]) {
    k <- k + 1L
    alpha[[as.character(k)]] <- lapply(alpha[[val]], function(z) if (is.numeric(z)) -z else .neg(z))
  }
  alpha
}

.lie_num <- function(z) if (is.numeric(z)) z else num_value_dbl(z)

#' @rdname liealgebras
#' @export
root_space <- function(x) paste(paste0("alpha[", seq_len(x$cartan_type$n), "]"), collapse = " + ")

#' @rdname liealgebras
#' @export
add_simple_roots <- function(x, root1, root2) {
  alpha <- simple_roots(x)
  if (root1 > length(alpha) || root2 > length(alpha)) stop("You've used a root that doesn't exist!", call. = FALSE)
  Map(function(a, b) if (is.numeric(a) || is.numeric(b)) .lie_num(a) + .lie_num(b) else .add(list(a, b)), alpha[[root1]], alpha[[root2]])
}

#' @rdname liealgebras
#' @export
add_as_roots <- function(x, root1, root2) {
  alpha <- .rs_all_roots(x)
  nr <- Map(function(a, b) .add(list(as_sym(a), as_sym(b))), root1, root2)
  key <- function(r) paste(vapply(r, function(z) as_sym(z)$key, ""), collapse = "|")
  if (key(nr) %in% vapply(alpha, key, "")) nr else "The sum of these two roots is not a root"
}

# ---------------------------------------------------------------------------
# Weyl groups

#' @rdname liealgebras
#' @export
WeylGroup <- function(t) structure(list(cartan_type = CartanType(t)), class = "WeylGroup")

#' @rdname liealgebras
#' @export
generators <- function(x) paste0("r", seq_len(x$cartan_type$n))

.weyl_group_order <- function(x) {
  ct <- x$cartan_type; n <- ct$n
  f <- function(k) as.numeric(gmp::factorialZ(k))
  switch(EXPR = ct$series, A = f(n + 1), B = , C = f(n) * 2^n, D = f(n) * 2^(n - 1),
         E = c(`6` = 51840, `7` = 2903040, `8` = 696729600)[[as.character(n)]], F = 1152, G = 12)
}

#' @rdname liealgebras
#' @export
group_name <- function(x) {
  ct <- x$cartan_type; n <- ct$n
  switch(EXPR = ct$series,
    A = paste0("S", n + 1, ": the symmetric group acting on ", n + 1, " elements."),
    B = , C = paste0("The hyperoctahedral group acting on ", 2 * n, " elements."),
    D = paste0("The symmetry group of the ", n, "-dimensional demihypercube."),
    E = paste0("The symmetry group of the ", n, "-polytope."),
    F = "The symmetry group of the 24-cell, or icositetrachoron.",
    G = "D6, the dihedral group of order 12, and symmetry group of the hexagon.")
}

# elts[1::3]: every third character from the second, the reflection numbers
.weyl_reflections <- function(weylelt) {
  ch <- strsplit(weylelt, "", fixed = TRUE)[[1L]]
  ch[seq(2L, length(ch), by = 3L)]
}

#' @rdname liealgebras
#' @export
delete_doubles <- function(x, reflections) {
  cp <- as.list(reflections)
  counter <- 1L
  # Python's enumerate over the list being shortened
  while (counter <= length(cp)) {
    if (counter < length(cp) && identical(cp[[counter + 1L]], cp[[counter]])) {
      cp <- cp[-c(counter, counter + 1L)]
    }
    counter <- counter + 1L
  }
  unlist(cp)
}

.mat_set <- function(M, i, j, v) { M$data[[(i - 1L) * M$ncol + j]] <- as_sym(v); M }

#' @rdname liealgebras
#' @export
matrix_form <- function(x, weylelt) {
  ct <- x$cartan_type; n <- ct$n
  refl <- as.integer(.weyl_reflections(weylelt))
  R <- function(...) Rational(...)
  E <- function(k) eye(k)
  mf <- switch(EXPR = ct$series, A = E(n + 1L), D = , B = , C = E(n), G = E(3L), F = E(4L), E = E(8L))
  for (a in refl) {
    mat <- switch(EXPR = ct$series,
      A = { m <- E(n + 1L); m <- .mat_set(m, a, a, 0); m <- .mat_set(m, a, a + 1L, 1); m <- .mat_set(m, a + 1L, a, 1); .mat_set(m, a + 1L, a + 1L, 0) },
      D = { m <- E(n)
            if (a < n) { m <- .mat_set(m, a, a, 0); m <- .mat_set(m, a, a + 1L, 1); m <- .mat_set(m, a + 1L, a, 1); .mat_set(m, a + 1L, a + 1L, 0) }
            else { m <- .mat_set(m, n - 1L, n, -1); m <- .mat_set(m, n - 1L, n - 1L, 0); m <- .mat_set(m, n, n - 1L, -1); .mat_set(m, n, n, 0) } },
      G = if (a == 1L) Matrix(list(list(1, 0, 0), list(0, 0, 1), list(0, 1, 0)))
          else Matrix(list(list(R(2, 3), R(2, 3), R(-1, 3)), list(R(2, 3), R(-1, 3), R(2, 3)), list(R(-1, 3), R(2, 3), R(2, 3)))),
      F = switch(as.character(a),
        "1" = Matrix(list(list(1, 0, 0, 0), list(0, 0, 1, 0), list(0, 1, 0, 0), list(0, 0, 0, 1))),
        "2" = Matrix(list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, 0, 1), list(0, 0, 1, 0))),
        "3" = Matrix(list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, -1))),
        { h <- R(1, 2); mh <- R(-1, 2)
          Matrix(list(list(h, h, h, h), list(h, h, mh, mh), list(h, mh, h, mh), list(h, mh, mh, h))) }),
      E = {
        if (a == 1L) {
          q <- R(1, 4); mq <- R(-1, 4); t <- R(3, 4)
          Matrix(list(list(t, q, q, q, q, q, q, mq), list(q, t, mq, mq, mq, mq, q, mq), list(q, mq, t, mq, mq, mq, mq, q),
                      list(q, mq, mq, t, mq, mq, mq, q), list(q, mq, mq, mq, t, mq, mq, q), list(q, mq, mq, mq, mq, t, mq, q),
                      list(q, mq, mq, mq, mq, mq, R(-3, 4), q), list(q, mq, mq, mq, mq, mq, mq, t)))
        } else if (a == 2L) {
          m <- E(8L); m <- .mat_set(m, 1L, 1L, 0); m <- .mat_set(m, 1L, 2L, -1); m <- .mat_set(m, 2L, 1L, -1); .mat_set(m, 2L, 2L, 0)
        } else {
          m <- E(8L); m <- .mat_set(m, a - 2L, a - 2L, 0); m <- .mat_set(m, a - 2L, a - 1L, 1); m <- .mat_set(m, a - 1L, a - 2L, 1); .mat_set(m, a - 1L, a - 1L, 0)
        }
      },
      B = , C = { m <- E(n)
        if (a == 1L) .mat_set(m, 1L, 1L, -1)
        else { m <- .mat_set(m, a - 1L, a - 1L, 0); m <- .mat_set(m, a - 1L, a, 1); m <- .mat_set(m, a, a - 1L, 1); .mat_set(m, a, a, 0) } })
    mf <- mf %*% mat
  }
  mf
}

.mat_is_identity <- function(M) {
  for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol))
    if (!identical(.mget(M, i, j)$key, (if (i == j) S_One else S_Zero)$key)) return(FALSE)
  TRUE
}

#' @rdname liealgebras
#' @export
element_order <- function(x, weylelt) {
  ct <- x$cartan_type
  if (ct$series == "G") {
    refl <- .weyl_reflections(weylelt)
    m <- delete_doubles(x, refl)
    while (!identical(delete_doubles(x, m), m)) { m <- delete_doubles(x, m); refl <- m }
    if (length(refl) %% 2L == 1L) return(2)
    if (!length(refl)) return(1)
    if (length(refl) == 1L) return(2)
    mm <- length(refl) %/% 2L
    lcm <- (6 * mm) / as.numeric(gmp::gcd.bigz(mm, 6))
    return(lcm / mm)
  }
  a <- matrix_form(x, weylelt); order <- 1
  mf <- matrix_form(x, weylelt)
  while (!.mat_is_identity(a)) { a <- a %*% mf; order <- order + 1 }
  order
}

#' @rdname liealgebras
#' @export
coxeter_diagram <- function(x) {
  ct <- x$cartan_type; n <- ct$n
  switch(EXPR = ct$series,
    A = , D = , E = dynkin_diagram(ct),
    B = , C = paste0(paste(rep("0", n - 1L), collapse = "---"), "===0\n", paste(seq_len(n), collapse = "   ")),
    F = paste0("0---0===0---0\n", paste(1:4, collapse = "   ")),
    G = "0\u2261\u2261\u22610\n1   2")
}
