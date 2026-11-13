# radsimp, rad_rationalize, collect_sqrt and split_surds after
# sympy/simplify/radsimp.py.

# _unevaluated_Mul: numbers collected in front, products flattened, the
# other factors sorted, nothing combined
.umul <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) && !is_sym(args[[1L]])) args <- args[[1L]]
  flat <- list()
  for (a in args) flat <- c(flat, if (a$type == "Mul") a$args else list(a))
  co <- S_One; other <- list()
  for (a in flat) if (is_number_type(a)) co <- .mul(list(co, a)) else other[[length(other) + 1L]] <- a
  other <- .key_sort(lapply(other, .sort_key), other)
  if (!(is_one_num(co) && co$type == "Integer")) other <- c(list(co), other)
  if (!length(other)) return(S_One)
  if (length(other) == 1L) return(other[[1L]])
  new_sym("Mul", other)
}

.signsimp_ne <- function(e) {
  r <- signsimp(e)
  if (r$type == "Add" && could_extract_minus_sign(r)) return(new_sym("Mul", list(S_NegativeOne, .neg(r))))
  r
}

#' Rationalize denominators and denest square roots
#'
#' `radsimp()` removes square roots from denominators as SymPy's
#' `radsimp()`: at most `max_terms` radical terms are handled, and with
#' `symbolic = FALSE` denominators with symbols are left alone.
#' `rad_rationalize(num, den)` rationalizes `num/den` and returns the new
#' numerator and denominator. `collect_sqrt()` collects terms with common
#' square roots. `split_surds()` splits a sum of surds by their gcd.
#' `sqrtdenest()` rewrites `sqrt(a + b*sqrt(c))` as a sum of square roots
#' when possible.
#'
#' @param e,expr An expression.
#' @param symbolic Also handle denominators with symbols.
#' @param max_terms The largest number of radical terms handled.
#' @param num,den Numerator and denominator.
#' @param evaluate For `collect_sqrt()`: `FALSE` returns the terms and the
#'   number of radical terms instead of the sum.
#' @return An expression, or a list.
#' @examples
#' radsimp(1 / (1 + sqrt(S(2))))
#' sqrtdenest(sqrt(5 + 2 * sqrt(S(6))))
#' @export
radsimp <- function(e, symbolic = TRUE, max_terms = 4L) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) radsimp(z, symbolic, max_terms)))
  .radsimp(as_sym(e), symbolic, max_terms)
}

.ispow2 <- function(d, symbolic, log2 = FALSE) {
  if (d$type != "Pow") return(FALSE)
  e <- d$args[[2L]]
  if ((e$type == "Rational" && gmp::denominator(e$value) == 2) || (symbolic && identical(denom(e)$key, S_Two$key))) return(TRUE)
  if (log2) {
    q <- gmp::as.bigz(1)
    if (.is_rat(e)) q <- gmp::denominator(e$value)
    else if (symbolic) { dd <- denom(e); if (dd$type == "Integer") q <- gmp::as.bigz(dd$value) }
    if (q != 1 && q > 0) { lq <- log2(as.numeric(q)); if (abs(lq - round(lq)) < 1e-12 && gmp::as.bigz(2)^round(lq) == q) return(TRUE) }
  }
  FALSE
}

.rad_num <- function(rterms) {
  sq <- function(z) .sqrt(z)
  n <- length(rterms)
  A <- function(i) rterms[[i]][[1L]]; a <- function(i) rterms[[i]][[2L]]
  M <- function(...) .mul(list(...)); P <- function(b, k) .pow(b, .int(k))
  if (n == 2L) return(.sub(M(sq(A(1)), a(1)), M(sq(A(2)), a(2))))
  if (n == 3L) {
    return(M(.add(list(M(sq(A(1)), a(1)), M(sq(A(2)), a(2)), .neg(M(sq(A(3)), a(3))))),
             .add(list(M(S_Two, sq(A(1)), sq(A(2)), a(1), a(2)), .neg(M(A(1), P(a(1), 2))), .neg(M(A(2), P(a(2), 2))), M(A(3), P(a(3), 2))))))
  }
  if (n == 4L) {
    t1 <- .add(list(M(sq(A(1)), a(1)), M(sq(A(2)), a(2)), .neg(M(sq(A(3)), a(3))), .neg(M(sq(A(4)), a(4)))))
    t2 <- .add(list(M(S_Two, sq(A(1)), sq(A(2)), a(1), a(2)), .neg(M(A(1), P(a(1), 2))), .neg(M(A(2), P(a(2), 2))),
                    M(.int(-2), sq(A(3)), sq(A(4)), a(3), a(4)), M(A(3), P(a(3), 2)), M(A(4), P(a(4), 2))))
    t3 <- .add(list(M(.int(-8), sq(A(1)), sq(A(2)), sq(A(3)), sq(A(4)), a(1), a(2), a(3), a(4)),
                    M(P(A(1), 2), P(a(1), 4)), M(.int(-2), A(1), A(2), P(a(1), 2), P(a(2), 2)),
                    M(.int(-2), A(1), A(3), P(a(1), 2), P(a(3), 2)), M(.int(-2), A(1), A(4), P(a(1), 2), P(a(4), 2)),
                    M(P(A(2), 2), P(a(2), 4)), M(.int(-2), A(2), A(3), P(a(2), 2), P(a(3), 2)),
                    M(.int(-2), A(2), A(4), P(a(2), 2), P(a(4), 2)), M(P(A(3), 2), P(a(3), 4)),
                    M(.int(-2), A(3), A(4), P(a(3), 2), P(a(4), 2)), M(P(A(4), 2), P(a(4), 4))))
    return(M(t1, t2, t3))
  }
  if (n == 1L) return(sq(A(1)))
  stop("NotImplementedError")
}

# default_sort_key of a Python tuple or list of sort keys
.seq_sort_key <- function(keys, cls = "tuple") list(list(10, 0, cls), list(length(keys), keys), .ONE_SORT_KEY, 1)
.seq_nodes <- function(xs) 1 + sum(vapply(xs, function(z) if (is_sym(z)) .nodes(z) else .seq_nodes(z), 0))

.radsimp_handle <- function(expr, symbolic, max_terms) {
  handle <- function(z) .radsimp_handle(z, symbolic, max_terms)
  if (!length(expr$args)) return(expr)
  if (!.is_expr_node(expr)) return(sym_rebuild(expr, lapply(expr$args, handle)))
  fr <- fraction(expr); n <- fr[[1L]]; d <- fr[[2L]]
  if (!length(d$args) && !length(n$args)) return(expr)
  if (length(n$args)) {
    n <- sym_rebuild(n, lapply(n$args, handle))
    return(.umul(n, handle(.inv(d))))
  }
  if (!(is_one_num(n) && n$type == "Integer")) return(.umul(n, handle(.inv(d))))
  if (d$type == "Mul") return(.umul(lapply(d$args, function(di) handle(.inv(di)))))
  if (!symbolic && length(free_symbols(d))) return(expr)
  if (.ispow2(d, symbolic)) {
    d2 <- .pow(sqrtdenest(.sqrt(d$args[[1L]])), numer(d$args[[2L]]))
    if (!identical(d2$key, d$key)) return(handle(.inv(d2)))
  } else if (d$type == "Pow" && (isTRUE(ask_fact(d$args[[2L]], "integer")) || isTRUE(ask_fact(d$args[[1L]], "positive")))) {
    return(.pow(handle(.inv(d$args[[1L]])), d$args[[2L]]))
  }
  if (!(d$type == "Add" || .ispow2(d, symbolic))) return(.inv(sym_rebuild(d, lapply(d$args, handle))))
  keep <- TRUE
  d <- .mexpand(d)
  if (!length(d$args)) return(.inv(d))
  if (!length(free_symbols(d)) && .is_expr_node(d)) {
    d_ <- tryCatch(nsimplify(d), error = function(err) NULL)
    if (!is.null(d_) && is_number_type(d_) && isTRUE(.expr_equals(d_, d))) return(.inv(d_))
  }
  repeat {
    ck <- list(); cv <- list(); cp <- list()
    for (m in .aargs(d)) {
      p2 <- list(); other <- list()
      for (i in .margs(m)) {
        if (.ispow2(i, symbolic, log2 = TRUE)) {
          p2[[length(p2) + 1L]] <- if (identical(i$args[[2L]]$key, S_Half$key)) i$args[[1L]] else .pow(i$args[[1L]], .mul(list(S_Two, i$args[[2L]])))
        } else if (is_const(i, "I")) p2[[length(p2) + 1L]] <- S_NegativeOne
        else other[[length(other) + 1L]] <- i
      }
      p2 <- ordered(p2)
      k <- paste0("k:", paste(vapply(p2, .k, ""), collapse = "\u0001"))
      if (is.null(ck[[k]])) { ck[[k]] <- p2; cv[[k]] <- list() }
      cv[[k]][[length(cv[[k]]) + 1L]] <- .mul(other)
    }
    items <- lapply(names(ck), function(k) list(ck[[k]], cv[[k]]))
    keys1 <- vapply(items, function(it) 1 + .seq_nodes(it[[1L]]) + .seq_nodes(it[[2L]]), 0)
    keys2 <- lapply(items, function(it) .seq_sort_key(list(.seq_sort_key(lapply(it[[1L]], .sort_key)), .seq_sort_key(lapply(it[[2L]], .sort_key), "list"))))
    items <- .key_sort(lapply(seq_along(items), function(i) list(keys1[i], keys2[[i]])), items)
    rterms <- lapply(items, function(it) list(.mul(it[[1L]]), .add(it[[2L]])))
    nrad <- length(rterms) - (if (is_one_num(rterms[[1L]][[1L]]) && rterms[[1L]][[1L]]$type == "Integer") 1L else 0L)
    if (nrad < 1L) break
    if (nrad > max_terms) { keep <- FALSE; break }
    if (length(rterms) > 4L) {
      if (all(vapply(rterms, function(t) t[[1L]]$type == "Integer" && .is_rat(.pow(t[[2L]], S_Two)), TRUE))) {
        rr <- rad_rationalize(S_One, new_sym("Add", lapply(rterms, function(t) .mul(list(.sqrt(t[[1L]]), t[[2L]])))))
        n <- .mul(list(n, rr[[1L]])); d <- rr[[2L]]
      } else keep <- FALSE
      break
    }
    num <- powsimp(.rad_num(rterms))
    n <- .mul(list(n, num)); d <- .mul(list(d, num))
    d <- powdenest(.mexpand(d), force = symbolic)
    if (.contains_key(d, S_Zero$key) || .contains_key(d, S_NaN$key) || .contains_key(d, S_ComplexInfinity$key)) return(expr)
    if (!length(d$args)) break
  }
  if (!keep) return(expr)
  .umul(n, .inv(d))
}

.radsimp <- function(expr, symbolic, max_terms) {
  if (!.is_expr_node(expr)) return(sym_rebuild(expr, lapply(expr$args, function(a) .radsimp(a, symbolic, max_terms))))
  ca <- .coeff_add_split(expr); coeff <- ca[[1L]]; expr <- ca[[2L]]
  expr <- .expr_normal(expr)
  old <- fraction(expr)
  nd <- fraction(.radsimp_handle(expr, symbolic, max_terms)); n <- nd[[1L]]; d <- nd[[2L]]
  same <- function(a, b) identical(a[[1L]]$key, b[[1L]]$key) && identical(a[[2L]]$key, b[[2L]]$key)
  if (!same(old, list(n, d))) {
    if (length(d$args)) {
      was <- list(n, d)
      n <- .signsimp_ne(n); d <- .signsimp_ne(d)
      u <- .F_new(.umul(n, .inv(d)))
      u <- .umul(lapply(names(u$b), function(k) .pow(u$b[[k]], u$e[[k]])))
      nd <- fraction(u); n <- nd[[1L]]; d <- nd[[2L]]
      if (same(old, list(n, d))) { n <- was[[1L]]; d <- was[[2L]] }
    }
    n <- expand_mul(n)
    if (is_number_type(d) || d$type == "Add") {
      nd2 <- fraction(gcd_terms(.umul(n, .inv(d))))
      if (is_number_type(nd2[[2L]]) || count_ops(nd2[[2L]]) <= count_ops(d)) {
        n <- signsimp(nd2[[1L]]); d <- signsimp(nd2[[2L]])
        if (n$type == "Mul" && is_number_type(n$args[[1L]])) n <- .mul(n$args)
      }
    }
  }
  .add(list(coeff, .umul(n, .inv(d))))
}

#' @rdname radsimp
#' @export
rad_rationalize <- function(num, den) {
  num <- as_sym(num); den <- as_sym(den)
  if (den$type != "Add") return(list(num, den))
  sp <- split_surds(den); g <- sp[[1L]]; a <- sp[[2L]]; b <- sp[[3L]]
  a <- .mul(list(a, .sqrt(g)))
  num <- .mexpand(.mul(list(.sub(a, b), num)))
  den <- .mexpand(.sub(.pow(a, S_Two), .pow(b, S_Two)))
  rad_rationalize(num, den)
}

#' @rdname radsimp
#' @export
split_surds <- function(expr) {
  expr <- as_sym(expr)
  args <- .key_sort(lapply(expr$args, .sort_key), expr$args)
  cm <- lapply(args, as_coeff_mul_split)
  surds <- lapply(Filter(function(x) x[[2L]]$type == "Pow", cm), function(x) .pow(x[[2L]], S_Two))
  surds <- .key_sort(lapply(surds, .sort_key), surds)
  sg <- .split_gcd(surds); g <- sg[[1L]]; b1 <- sg[[2L]]; b2 <- sg[[3L]]
  g2 <- g
  if (!length(b2) && length(b1) >= 2L) {
    b1n <- lapply(b1, function(x) .div(x, g))
    b1n <- Filter(function(x) !is_one_num(x), b1n)
    sg <- .split_gcd(b1n); g1 <- sg[[1L]]; b2 <- sg[[3L]]
    g2 <- .mul(list(g, g1))
  }
  b1k <- vapply(b1, .k, "")
  a1v <- list(); a2v <- list()
  for (cs in cm) {
    c <- cs[[1L]]; s <- cs[[2L]]
    if (s$type == "Pow" && identical(s$args[[2L]]$key, S_Half$key) && s$args[[1L]]$key %in% b1k)
      a1v[[length(a1v) + 1L]] <- .mul(list(c, .sqrt(.div(s$args[[1L]], g2))))
    else a2v[[length(a2v) + 1L]] <- .mul(list(c, s))
  }
  list(g2, .add(a1v), .add(a2v))
}

.split_gcd <- function(a) {
  g <- a[[1L]]; b1 <- list(g); b2 <- list()
  for (x in a[-1L]) {
    g1 <- gcd(g, x)
    if (is_one_num(g1)) b2[[length(b2) + 1L]] <- x
    else { g <- g1; b1[[length(b1) + 1L]] <- x }
  }
  list(g, b1, b2)
}

#' @rdname radsimp
#' @export
collect_sqrt <- function(expr, evaluate = TRUE) {
  expr <- as_sym(expr)
  cp <- .acp(expr); coeff <- cp[[1L]]; expr <- cp[[2L]]
  vars <- list()
  for (a in .aargs(expr)) for (m in Filter(.is_comm, .margs(a))) {
    if (!length(free_symbols(m)) && .is_expr_node(m) &&
        ((m$type == "Pow" && m$args[[2L]]$type == "Rational" && gmp::denominator(m$args[[2L]]$value) == 2) || is_const(m, "I")))
      vars[[m$key]] <- m
  }
  d <- do.call(collect_const, c(list(expr), unname(vars), list(Numbers = FALSE)))
  hit <- !identical(expr$key, d$key)
  if (!evaluate) {
    nrad <- 0L
    args <- ordered(.aargs(d))
    for (i in seq_along(args)) {
      for (ci in Filter(.is_comm, .margs(args[[i]]))) {
        if ((ci$type == "Pow" && ci$args[[2L]]$type == "Rational" && gmp::denominator(ci$args[[2L]]$value) == 2) || is_const(ci, "I")) {
          nrad <- nrad + 1L; break
        }
      }
      args[[i]] <- .mul(list(args[[i]], coeff))
    }
    if (!(hit || nrad)) args <- list(.add(args))
    return(list(args, nrad))
  }
  .mul(list(coeff, d))
}
