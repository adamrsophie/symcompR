# Factoring with noncommutative symbols, after factor_nc and _mask_nc in
# sympy/core/exprtools.py.

# _mask_nc: one noncommutative symbol or object becomes a commutative dummy;
# several become noncommutative dummies. Returns the masked expression, the
# replacements back, and the noncommutative symbols.
.mask_nc <- function(e) {
  if (.is_comm(e)) return(list(expr = e, back = list(), nc = list()))
  nc_syms <- list(); nc_obj <- list()
  walk <- function(z) {
    if (.is_comm(z)) return(invisible())
    if (is_symbol_like(z)) { nc_syms[[z$key]] <<- z; return(invisible()) }
    if (!(z$type %in% c("Add", "Mul", "Pow"))) { nc_obj[[z$key]] <<- z; return(invisible()) }
    for (a in z$args) walk(a)
  }
  walk(e)
  rep <- list()
  if (length(nc_obj) == 1L && !length(nc_syms)) { rep[[1L]] <- list(nc_obj[[1L]], Dummy("mask0")); nc_obj <- list() }
  else if (length(nc_syms) == 1L && !length(nc_obj)) { rep[[1L]] <- list(nc_syms[[1L]], Dummy("mask0")); nc_syms <- list() }
  objs <- .key_sort(lapply(nc_obj, .sort_key), unname(nc_obj))
  for (n in objs) {
    d <- Dummy(paste0("mask", length(rep)), commutative = FALSE)
    rep[[length(rep) + 1L]] <- list(n, d)
    nc_syms[[d$key]] <- d
  }
  ex <- if (length(rep)) subs(e, rep) else e
  ns <- unname(nc_syms)
  list(expr = ex, back = lapply(rep, function(p) list(p[[2L]], p[[1L]])),
       nc = if (length(ns)) .key_sort(lapply(ns, .sort_key), ns) else list())
}

#' Factoring with noncommutative symbols
#'
#' `factor_nc(e)` factors an expression that contains noncommutative
#' symbols, as SymPy's `factor_nc()`: common commutative factors and common
#' noncommutative prefixes and suffixes come out, and the remaining sum is
#' factored with the noncommutative symbols kept in order.
#'
#' @param e An expression.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE)
#' factor_nc(x^2 + 2*A*x + A^2)
#' factor_nc(expand((x + A)*(x + B)))
#' factor_nc(A*B + A*Symbol("C", commutative = FALSE))
#' @export
factor_nc <- function(e) {
  e <- as_sym(e)
  if (!length(e$args) || is_relational(e) || e$type %in% c("And", "Or", "Not", "Xor")) return(e)
  if (e$type != "Add") return(sym_rebuild(e, lapply(e$args, factor_nc)))
  e <- .add(lapply(e$args, function(a) expand(a, deep = FALSE, mul = FALSE, multinomial = FALSE,
                                                power_base = FALSE, log = FALSE)))
  m <- .mask_nc(e)
  if (length(m$back)) return(subs(sym_factor(m$expr), m$back))
  # args_cnc of each term
  args <- lapply(if (e$type == "Add") e$args else list(e), function(t) {
    fs <- if (t$type == "Mul") t$args else list(t)
    list(c = Filter(.is_comm, fs), nc = Filter(Negate(.is_comm), fs))
  })
  c <- S_One; g <- S_One; l <- S_One; r <- S_One; hit <- FALSE
  for (i in seq_along(args)) {
    if (i == 1L) c <- .mul(args[[i]]$c)
    else if (length(args[[i]]$c)) c <- gcd(c, .mul(args[[i]]$c))
    else c <- S_One
  }
  if (!is_one_num(c)) {
    hit <- TRUE
    sp <- as_coeff_mul_split(c); c <- sp[[1L]]; g <- sp[[2L]]
    for (i in seq_along(args)) {
      cc <- args[[i]]$c
      if (!is_one_num(g)) { q <- .div(.mul(cc), g); cc <- if (q$type == "Mul") q$args else list(q) }
      cc <- if (length(cc)) c(list(.div(cc[[1L]], c)), cc[-1L]) else list(.inv(c))
      args[[i]]$c <- cc
    }
  }
  same <- function(a, b) identical(a$key, b$key)
  # a common noncommutative prefix, or a power of a common first base
  n <- NULL
  for (i in seq_along(args)) {
    nci <- args[[i]]$nc
    if (i == 1L) n <- nci
    else { k <- 0L; while (k < min(length(n), length(nci)) && same(n[[k + 1L]], nci[[k + 1L]])) k <- k + 1L; n <- n[seq_len(k)] }
    if (!length(n)) break
  }
  if (length(n)) {
    hit <- TRUE; l <- .mul(n)
    for (i in seq_along(args)) args[[i]]$nc <- args[[i]]$nc[-seq_along(n)]
  } else if (length(args[[1L]]$nc)) {
    be <- .base_exp(args[[1L]]$nc[[1L]]); b <- be[[1L]]; ex <- be[[2L]]
    if (ex$type == "Integer" && all(vapply(args, function(t) length(t$nc) > 0L, TRUE))) {
      bes <- lapply(args, function(t) .base_exp(t$nc[[1L]]))
      if (all(vapply(bes, function(z) z[[2L]]$type == "Integer" && same(z[[1L]], b), TRUE))) {
        ex <- .int(min(vapply(bes, function(z) as.numeric(z[[2L]]$value), 0)))
        hit <- TRUE; l <- .pow(b, ex)
        for (i in seq_along(args)) args[[i]]$nc[[1L]] <- .mul(list(.pow(b, .neg(ex)), args[[i]]$nc[[1L]]))
      }
    }
  }
  # a common noncommutative suffix, or a power of a common last base
  n <- NULL
  for (i in seq_along(args)) {
    nci <- rev(args[[i]]$nc)
    if (i == 1L) n <- nci
    else { k <- 0L; while (k < min(length(n), length(nci)) && same(n[[k + 1L]], nci[[k + 1L]])) k <- k + 1L; n <- n[seq_len(k)] }
    if (!length(n)) break
  }
  if (length(n)) {
    hit <- TRUE; r <- .mul(rev(n))
    for (i in seq_along(args)) { L <- length(args[[i]]$nc); args[[i]]$nc <- args[[i]]$nc[seq_len(L - length(n))] }
  } else if (all(vapply(args, function(t) length(t$nc) > 0L, TRUE))) {
    last <- function(t) t$nc[[length(t$nc)]]
    bes <- lapply(args, function(t) .base_exp(last(t)))
    b <- bes[[1L]][[1L]]
    if (all(vapply(bes, function(z) z[[2L]]$type == "Integer" && same(z[[1L]], b), TRUE))) {
      ex <- .int(min(vapply(bes, function(z) as.numeric(z[[2L]]$value), 0)))
      hit <- TRUE; r <- .pow(b, ex)
      for (i in seq_along(args)) { L <- length(args[[i]]$nc); args[[i]]$nc[[L]] <- .mul(list(args[[i]]$nc[[L]], .pow(b, .neg(ex)))) }
    }
  }
  mid <- if (hit) .add(lapply(args, function(t) .mul(c(t$c, t$nc)))) else e
  # factor the middle with the noncommutative symbols made commutative
  rep1 <- lapply(m$nc, function(s) list(s, Dummy(s$name)))
  unrep1 <- rev(lapply(rep1, function(p) list(p[[2L]], p[[1L]])))
  new_mid <- powsimp(sym_factor(if (length(rep1)) subs(mid, rep1) else mid))
  if (length(unrep1)) new_mid <- subs(new_mid, unrep1)
  if (new_mid$type == "Pow") return(.keep_coeff(c, .mul(list(g, l, new_mid, r))))
  if (new_mid$type == "Mul") {
    pem <- function(z) expand(z, deep = TRUE, mul = TRUE, power_exp = TRUE, power_base = FALSE, multinomial = TRUE, log = FALSE)
    cfac <- Filter(.is_comm, new_mid$args)
    ncfac <- list()
    for (f in Filter(Negate(.is_comm), new_mid$args)) {
      be <- .base_exp(f)
      if (be[[2L]]$type == "Integer" && be[[2L]]$value > 0) ncfac <- c(ncfac, rep(list(be[[1L]]), as.integer(be[[2L]]$value)))
      else ncfac[[length(ncfac) + 1L]] <- f
    }
    pre_mid <- .mul(c(list(g), cfac, list(l)))
    target <- pem(.div(e, c))
    # every ordering of the noncommutative factors, as itertools.variations
    perms <- function(v) if (length(v) <= 1L) list(v) else
      do.call(base::c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
    if (length(ncfac) <= 6L) for (s in perms(seq_along(ncfac))) {
      ok <- .mul(c(list(pre_mid), ncfac[s], list(r)))
      if (identical(pem(ok)$key, target$key)) return(.keep_coeff(c, ok))
    }
  }
  .keep_coeff(c, .mul(list(g, l, mid, r)))
}
