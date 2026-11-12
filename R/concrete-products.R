# Product.doit and Product._eval_product of sympy/concrete/products.py, and
# deltaproduct of sympy/concrete/delta.py.

.prod_node <- function(f, lims) .sum_new(f, lims, "Product")

# Product.doit
.prod_doit <- function(e, deep = TRUE) {
  lims <- lapply(e$args[-1L], function(t) t$args)
  reps <- list()
  for (l in lims) {
    d <- .dummy_inherited(l[[1L]], l[[2L]], l[[3L]])
    if (!is.null(d)) reps[[length(reps) + 1L]] <- list(l[[1L]], d)
  }
  if (length(reps)) {
    undo <- lapply(reps, function(p) list(p[[2L]], p[[1L]]))
    return(xreplace(.prod_doit(xreplace(e, reps), deep), undo))
  }
  f <- e$args[[1L]]
  for (n in seq_along(lims)) {
    i <- lims[[n]][[1L]]; a <- lims[[n]][[2L]]; b <- lims[[n]][[3L]]
    dif <- .sub(b, a)
    if (isTRUE(ask_fact(dif, "integer")) && isTRUE(ask_fact(dif, "negative"))) {
      a0 <- a; a <- .add(list(b, S_One)); b <- .sub(a0, S_One); f <- .inv(f)
    }
    g <- .eval_product(f, i, a, b)
    if (is.null(g) || is_nan(g)) return(.prod_node(powsimp(f), lims[n:length(lims)]))
    f <- g
  }
  if (deep) doit(f) else powsimp(f)
}

.prod_direct <- function(term, k, a, n) {
  m <- as.integer(as.character(.sub(n, a)$value))
  .mul(lapply(seq_len(m + 1L) - 1L, function(j) subs(term, k, .add(list(a, .int(j))))))
}

# Product._eval_product
.eval_product <- function(term, k, a, n) {
  len <- .add(list(.sub(n, a), S_One))
  if (!.has1(term, k)) {
    if (isTRUE(ask_fact(.sub(term, S_One), "zero"))) return(S_One)
    return(.pow(term, len))
  }
  if (identical(a$key, n$key)) return(subs(term, k, a))
  if (.has_function(term, "KroneckerDelta") && .has_simple_delta(term, k)) return(deltaproduct(term, list(k, a, n)))
  dif <- .sub(n, a)
  definite <- dif$type == "Integer"
  if (definite && dif$value < 100) return(.prod_direct(term, k, a, n))
  if (isTRUE(is_polynomial(term, k))) {
    A <- S_One; B <- S_One; Q <- S_One
    rs <- roots(term, k)
    mult <- attr(rs, "multiplicity")
    M <- 0L
    for (j in seq_along(rs)) {
      r <- rs[[j]]; m <- .int(mult[j])
      M <- M + as.integer(mult[j])
      A <- .mul(list(A, .pow(fcall("RisingFactorial", list(.sub(a, r), len)), m)))
      Q <- .mul(list(Q, .pow(.sub(n, r), m)))
    }
    if (M < .deg_in(term, k)) {
      # quo(poly, Q.as_poly(k)): Q has no k, so this divides by it
      arg <- expand(.div(term, Q))
      B <- .prod_doit(.prod_node(arg, list(list(k, a, n))))
    }
    return(.mul(list(.pow(.lc_in(term, k), len), A, B)))
  }
  if (term$type == "Add") {
    factored <- factor_terms(term, fraction = TRUE)
    if (factored$type == "Mul") return(.eval_product(factored, k, a, n))
  } else if (term$type == "Mul") {
    cm <- .as_coeff_mul_dep(term, k); without <- cm[[1L]]; with <- cm[[2L]]
    if (length(with) >= 2L) {
      exclude <- list(); include <- list()
      for (t in with) {
        p <- .eval_product(t, k, a, n)
        if (!is.null(p)) exclude[[length(exclude) + 1L]] <- p else include[[length(include) + 1L]] <- t
      }
      if (!length(exclude)) return(NULL)
      arg <- if (!length(include)) S_One else if (length(include) == 1L) include[[1L]] else new_sym("Mul", include)
      B <- .prod_doit(.prod_node(arg, list(list(k, a, n))))
      return(.mul(list(.pow(without, len), .mul(exclude), B)))
    }
    p <- .eval_product(with[[1L]], k, a, n)
    if (is.null(p)) p <- .prod_doit(.prod_node(with[[1L]], list(list(k, a, n))))
    return(.mul(list(.pow(without, len), p)))
  } else if (term$type == "Pow" || (term$type == "Function" && term$name == "exp")) {
    be <- .abe(term); base <- be[[1L]]; ex <- be[[2L]]
    if (!.has1(base, k)) return(.pw(base, summation(ex, list(k, a, n))))
    if (!.has1(ex, k)) {
      p <- .eval_product(base, k, a, n)
      if (!is.null(p)) return(.pow(p, ex))
    }
  } else if (term$type == "Product") {
    evaluated <- doit(term)
    f <- .eval_product(evaluated, k, a, n)
    return(if (is.null(f)) .prod_node(evaluated, list(list(k, a, n))) else f)
  }
  if (definite) return(.prod_direct(term, k, a, n))
  NULL
}

# ---------------------------------------------------------------------------
# delta.py: deltaproduct and its helpers

.delta_solve <- function(eqs) {
  syms <- unique(unlist(lapply(eqs, free_symbols), recursive = FALSE))
  syms <- syms[!duplicated(vapply(syms, `[[`, "", "key"))]
  if (!length(syms)) return(list())
  sl <- tryCatch(solve(eqs, syms, dict = TRUE), error = function(err) NULL)
  if (is.null(sl)) return(NULL)
  bykey <- stats::setNames(syms, vapply(syms, `[[`, "", "name"))
  lapply(sl, function(d) lapply(names(d), function(nm) list(bykey[[nm]], d[[nm]])))
}

.remove_multiple_delta <- function(expr) {
  if (expr$type == "Add") return(.add(lapply(expr$args, .remove_multiple_delta)))
  if (expr$type != "Mul") return(expr)
  eqs <- list(); newargs <- list()
  for (a in expr$args) {
    if (a$type == "Function" && a$name == "KroneckerDelta") eqs[[length(eqs) + 1L]] <- .sub(a$args[[1L]], a$args[[2L]])
    else newargs[[length(newargs) + 1L]] <- a
  }
  if (!length(eqs)) return(expr)
  solns <- .delta_solve(eqs)
  if (is.null(solns)) return(expr)
  if (!length(solns)) return(S_Zero)
  if (length(solns) == 1L) {
    newargs <- c(newargs, lapply(solns[[1L]], function(p) KroneckerDelta(p[[1L]], p[[2L]])))
    expr2 <- .mul(newargs)
    if (!identical(expr2$key, expr$key)) return(.remove_multiple_delta(expr2))
  }
  expr
}

.simplify_delta <- function(expr) {
  if (expr$type == "Function" && expr$name == "KroneckerDelta") {
    slns <- .delta_solve(list(.sub(expr$args[[1L]], expr$args[[2L]])))
    if (!is.null(slns) && length(slns) == 1L)
      return(.mul(lapply(slns[[1L]], function(p) KroneckerDelta(p[[1L]], p[[2L]]))))
  }
  expr
}

#' @rdname summation
#' @export
deltaproduct <- function(f, limit) {
  f <- as_sym(f); limit <- lapply(limit, as_sym)
  x <- limit[[1L]]; a <- limit[[2L]]; b <- limit[[3L]]
  if (isTRUE(ask_fact(.sub(b, a), "negative"))) return(S_One)
  if (!.has_function(f, "KroneckerDelta")) return(product(f, limit))
  if (f$type == "Add") {
    delta <- NULL; terms <- list()
    for (arg in .key_sort(lapply(f$args, .sort_key), f$args)) {
      if (is.null(delta) && .has_simple_delta(arg, x)) delta <- arg else terms[[length(terms) + 1L]] <- arg
    }
    newexpr <- .add(terms)
    k <- Dummy("kprime", integer = TRUE)
    inner <- .mul(list(deltaproduct(newexpr, list(x, a, .sub(k, S_One))), subs(delta, x, k),
                       deltaproduct(newexpr, list(x, .add(list(k, S_One)), b))))
    result <- .add(list(deltaproduct(newexpr, limit),
                        deltasummation(inner, list(k, a, b), no_piecewise = .has_simple_delta(newexpr, x))))
    return(.remove_multiple_delta(result))
  }
  delta <- .extract_delta(f, x)[[1L]]
  if (is.null(delta)) {
    g <- .expand_delta(f, x)
    if (!identical(g$key, f$key)) return(sym_factor(deltaproduct(g, limit)))
    return(product(f, limit))
  }
  .add(list(.remove_multiple_delta(.mul(list(subs(f, x, a), KroneckerDelta(b, a)))),
            .simplify_delta(KroneckerDelta(b, .sub(a, S_One)))))
}
