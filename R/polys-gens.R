# Independent generators in composite domains. SymPy's PolyRing and
# FracField treat every generator as an indeterminate, so in ZZ[E, exp(1/2)]
# the element exp(1/2)*exp(1/2) is the monomial exp(1/2)**2, not E. The
# coefficients of a dpoly (R/polys-dpoly.R) are expressions, which would
# evaluate such products; a generator that is not a symbol is therefore
# stood in for by a Dummy inside the domain. One global table maps each
# generator to its Dummy, so domains built separately over the same
# generator share it. Elements enter the domain through .dn() (any mix of
# generators and their dummies is accepted) and leave it through .gd_out(),
# which is SymPy's to_sympy() and may let the products evaluate.

.gd_state <- new.env(parent = emptyenv())
.gd_state$fwd <- list()   # generator key -> Dummy
.gd_state$back <- list()  # Dummy key -> generator

.gen_dummy <- function(g) {
  if (is_symbol_like(g)) return(g)
  d <- .gd_state$fwd[[g$key]]
  if (is.null(d)) {
    d <- Dummy("gen")
    .gd_state$fwd[[g$key]] <- d
    .gd_state$back[[d$key]] <- g
  }
  d
}

# the generators of a composite domain as its elements see them
.dk_dgens <- function(K) lapply(.dk_gens(K), .gen_dummy)
.dk_has_dgens <- function(K) .dk_comp(K) && any(!vapply(.dk_gens(K), is_symbol_like, NA))

# TRUE when a contains no generator expression: numbers, I, symbols and
# integer powers of them only
.gd_clean <- function(a) {
  switch(a$type,
    Integer = , Rational = , Float = , Symbol = , Dummy = TRUE,
    Add = , Mul = { for (z in a$args) if (!.gd_clean(z)) return(FALSE); TRUE },
    Pow = a$args[[2L]]$type == "Integer" && .gd_clean(a$args[[1L]]),
    identical(a$key, S_I$key))
}

# an expression of the domain K (generators, their dummies, or both) with
# every generator that is not a symbol written as its Dummy
.gd_in <- function(K, a) {
  if (!.dk_has_dgens(K) || .gd_clean(a)) return(a)
  gens <- .dk_gens(K)
  if (.dk(K) == "FracField") {
    nd <- .as_numer_denom(a)
    return(.div(.gd_poly_in(K, gens, expand(nd[[1L]])), .gd_poly_in(K, gens, expand(nd[[2L]]))))
  }
  .gd_poly_in(K, gens, a)
}
.gd_poly_in <- function(K, gens, a) {
  ns <- which(!vapply(gens, is_symbol_like, NA))
  dg <- lapply(gens, .gen_dummy)
  allg <- c(gens, dg[ns])
  hasI <- .contains_key(a, S_I$key)
  if (hasI) allg <- c(allg, list(S_I))
  a <- .gd_rebuild(a, gens[ns], dg[ns])
  fl <- .has_type(a, "Float")
  p <- if (fl) NULL else tryCatch(expr_to_qpoly(a, allg, strict = TRUE), error = function(e) NULL)
  if (is.null(p)) {
    # not polynomial in the generators (or Float coefficients): replace the
    # generators that appear literally
    return(xreplace(a, lapply(ns, function(i) list(gens[[i]], dg[[i]]))))
  }
  n <- length(gens)
  e <- p$e
  e[, ns] <- e[, ns, drop = FALSE] + e[, n + seq_along(ns), drop = FALSE]
  keep <- c(seq_len(n), if (hasI) ncol(e))
  qpoly_to_expr(qp_norm(e[, keep, drop = FALSE], p$c), c(dg, if (hasI) list(S_I)))
}

# PolyRing._rebuild_expr: a generator becomes its Dummy, and a power b**e of
# the base of a generator bg**eg becomes Dummy**(e/eg) when eg divides e,
# so with the generator exp(1/2) the element E is its square
.gd_rebuild <- function(a, gens, dg) {
  pw <- Filter(Negate(is.null), lapply(seq_along(gens), function(i) {
    g <- gens[[i]]
    if (g$type == "Pow" || (g$type == "Function" && g$name == "exp") || identical(g$key, S_E$key)) c(.base_exp(g), list(dg[[i]]))
  }))
  keys <- vapply(gens, `[[`, "", "key")
  rb <- function(e) {
    k <- match(e$key, keys)
    if (!is.na(k)) return(dg[[k]])
    if (e$type == "Add") return(.add(lapply(e$args, rb)))
    if (e$type == "Mul") return(.mul(lapply(e$args, rb)))
    if (e$type == "Pow" || (e$type == "Function" && e$name == "exp") || identical(e$key, S_E$key)) {
      be <- .base_exp(e)
      for (g in pw) {
        if (identical(g[[1L]]$key, be[[1L]]$key)) {
          q <- .div(be[[2L]], g[[2L]])
          if (q$type == "Integer") return(.pow(g[[3L]], q))
        }
      }
      if (be[[2L]]$type == "Integer" && !is_one_num(be[[2L]])) return(.pow(rb(be[[1L]]), be[[2L]]))
    }
    e
  }
  if (!length(gens)) a else rb(a)
}

# to_sympy: the dummies back to their generators
.gd_out <- function(a) {
  if (!length(.gd_state$back) || !is_sym(a)) return(a)
  ks <- names(.gd_state$back)
  hit <- ks[vapply(ks, function(k) .contains_key(a, k), NA)]
  if (!length(hit)) return(a)
  xreplace(a, lapply(hit, function(k) list(.gd_key_sym(k), .gd_state$back[[k]])))
}
.gd_key_sym <- function(k) {
  for (d in .gd_state$fwd) if (identical(d$key, k)) return(d)
  NULL
}
# a dpoly as an expression with its coefficients kept in domain form
.dp_expr_raw <- function(p) {
  n <- length(p$c)
  if (!n) return(S_Zero)
  .add(lapply(seq_len(n), function(i) .mul(list(p$c[[i]], .pow(p$x, .int(n - i))))))
}
