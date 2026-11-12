# Sum.doit and the summation algorithms of sympy/concrete/summations.py,
# gosper.py and delta.py: direct sums, eval_sum_symbolic (telescoping,
# Faulhaber, geometric, Gosper), hypergeometric sums through hyperexpand,
# sums of rational functions by residues, zeta sums and Kronecker deltas.

# Python truthiness of a result: None and numeric zero are false
.truthy <- function(s) !is.null(s) && !(is_number_type(s) && num_is_zero(s))
.sum_node <- function(f, lims) .sum_new(f, lims)

# as_two_terms: the first argument in SymPy's order and the rest
.as_two_terms <- function(e) {
  args <- .sympy_args(e)
  if (length(args) == 1L) return(list(if (e$type == "Add") S_Zero else S_One, e))
  op <- if (e$type == "Add") .add else .mul
  list(args[[1L]], if (length(args) == 2L) args[[2L]] else op(args[-1L]))
}

# Expr.as_independent(i) of a product: noncommutative factors after the first
# dependent one stay dependent
.as_indep_mul <- function(e, i) {
  if (e$type != "Mul") return(if (.has1(e, i)) list(S_One, e) else list(e, S_One))
  comm <- vapply(e$args, .is_comm, TRUE)
  ca <- e$args[comm]; nc <- e$args[!comm]
  dep <- vapply(ca, function(t) .has1(t, i), TRUE)
  indep <- ca[!dep]; depend <- ca[dep]
  for (k in seq_along(nc)) {
    if (.has1(nc[[k]], i)) { depend <- c(depend, nc[k:length(nc)]); break }
    indep <- c(indep, nc[k])
  }
  list(.mul(indep), if (length(depend) > 1L && length(nc)) new_sym("Mul", depend) else .mul(depend))
}

# _dummy_with_inherited_properties_concrete
.SUM_ASSUMPTIONS <- c("extended_nonnegative", "nonnegative", "extended_nonpositive", "nonpositive",
                      "extended_positive", "positive", "extended_negative", "negative", "integer",
                      "rational", "finite", "zero", "real", "extended_real")
.dummy_inherited <- function(x, a, b) {
  if (!(x$type %in% c("Symbol", "Dummy"))) return(NULL)
  keep <- list(); add <- FALSE
  for (s in .SUM_ASSUMPTIONS) {
    if (isTRUE(ask_fact(x, s))) keep[[s]] <- TRUE
    else if (isTRUE(ask_fact(a, s)) && isTRUE(ask_fact(b, s))) { keep[[s]] <- TRUE; add <- TRUE }
  }
  if (!add) return(NULL)
  do.call(Dummy, c(list("d"), keep))
}

# Sum.doit
.sum_doit <- function(e, deep = TRUE) {
  f0 <- e$args[[1L]]
  f <- if (deep) doit(f0) else f0
  lims <- lapply(e$args[-1L], function(t) t$args)
  reps <- list()
  for (l in lims) {
    d <- .dummy_inherited(l[[1L]], l[[2L]], l[[3L]])
    if (!is.null(d)) reps[[length(reps) + 1L]] <- list(l[[1L]], d)
  }
  if (length(reps)) {
    undo <- lapply(reps, function(p) list(p[[2L]], p[[1L]]))
    did <- .sum_doit(xreplace(e, reps), deep)
    return(xreplace(did, undo))
  }
  for (n in seq_along(lims)) {
    i <- lims[[n]][[1L]]; a <- lims[[n]][[2L]]; b <- lims[[n]][[3L]]
    dif <- .sub(b, a)
    if (is_negone_num(dif)) return(S_Zero)
    if (isTRUE(ask_fact(dif, "integer")) && isTRUE(ask_fact(dif, "negative"))) {
      a0 <- a; a <- .add(list(b, S_One)); b <- .sub(a0, S_One); f <- .neg(f)
    }
    newf <- .eval_sum(f, i, a, b)
    if (is.null(newf)) {
      if (identical(f$key, f0$key)) {
        z <- .eval_zeta_function(e, f, i, a, b)
        if (!is.null(z)) return(z)
        return(e)
      }
      return(.sum_node(f, lims[n:length(lims)]))
    }
    f <- newf
  }
  if (deep && f$type != "Piecewise") return(doit(f))
  f
}

.eval_sum <- function(f, i, a, b) {
  if (isTRUE(ask_fact(f, "zero"))) return(S_Zero)
  if (!.has1(f, i)) return(.mul(list(f, .add(list(.sub(b, a), S_One)))))
  if (identical(a$key, b$key)) return(subs(f, i, a))
  if (.is_comparable(a) && .is_comparable(b) && isTRUE(ask_fact(.sub(a, b), "extended_positive")))
    return(.eval_sum(f, i, .add(list(b, S_One)), .sub(a, S_One)))
  if (f$type == "Piecewise") {
    pairs <- .pw_pairs(f)
    if (!any(vapply(pairs, function(p) .has1(p[[2L]], i), TRUE))) {
      newargs <- list()
      for (p in pairs) {
        ne <- .eval_sum(p[[1L]], i, a, b)
        if (is.null(ne)) return(NULL)
        newargs[[length(newargs) + 1L]] <- list(ne, p[[2L]])
      }
      return(.piecewise(newargs))
    }
  }
  if (.has_function(f, "KroneckerDelta")) {
    f <- .map_bottom_up(f, function(z) if (z$type == "Sum") sym_factor(z) else z)
    if (.has_simple_delta(f, i)) return(deltasummation(f, list(i, a, b)))
  }
  dif <- .sub(b, a)
  definite <- dif$type == "Integer"
  if (definite && dif$value < 100) return(.eval_sum_direct(f, i, a, b))
  if (f$type == "Piecewise") return(NULL)
  value <- .eval_sum_symbolic(expand(f), i, a, b)
  if (!is.null(value)) return(value)
  if (definite) return(.eval_sum_direct(f, i, a, b))
  NULL
}

# apart(f, i), or f when it is not a rational function of i (PolynomialError)
.sum_apart <- function(f, i) {
  if (!length(f$args)) return(f)
  nd <- .as_numer_denom(f)
  if (!isTRUE(is_polynomial(nd[[1L]], i)) || !isTRUE(is_polynomial(nd[[2L]], i))) return(f)
  tryCatch(apart(f, i), error = function(err) f)
}

.eval_sum_direct <- function(expr, i, a, b) {
  dif <- .sub(b, a)
  if (expr$type == "Mul") {
    sp <- .as_indep_mul(expr, i); without <- sp[[1L]]; with <- sp[[2L]]
    if (!is_one_num(without)) {
      s <- .eval_sum_direct(with, i, a, b)
      if (.truthy(s)) {
        r <- .mul(list(without, s))
        if (!is_nan(r)) return(r)
      }
    } else {
      lr <- .as_two_terms(expr); L <- lr[[1L]]; R <- lr[[2L]]
      if (!.has1(L, i)) {
        sR <- .eval_sum_direct(R, i, a, b)
        if (.truthy(sR)) return(.mul(list(L, sR)))
      }
      if (!.has1(R, i)) {
        sL <- .eval_sum_direct(L, i, a, b)
        if (.truthy(sL)) return(.mul(list(sL, R)))
      }
    }
  }
  expr <- .sum_apart(expr, i)
  if (expr$type == "Add") {
    sp <- .as_indep(expr, list(i), as_add = TRUE); without <- sp[[1L]]; with <- sp[[2L]]
    if (!is_zero_num(without)) {
      s <- .eval_sum_direct(with, i, a, b)
      if (.truthy(s)) {
        r <- .add(list(.mul(list(without, .add(list(dif, S_One)))), s))
        if (!is_nan(r)) return(r)
      }
    } else {
      lr <- .as_two_terms(expr)
      lsum <- .eval_sum_direct(lr[[1L]], i, a, b)
      rsum <- .eval_sum_direct(lr[[2L]], i, a, b)
      if (!is.null(lsum) && !is.null(rsum)) {
        r <- .add(list(lsum, rsum))
        if (!is_nan(r)) return(r)
      }
    }
  }
  n <- as.integer(as.character(dif$value))
  .add(lapply(seq_len(n + 1L) - 1L, function(j) subs(expr, i, .add(list(a, .int(j))))))
}

.eval_sum_symbolic <- function(f, i, a, b) {
  f_orig <- f
  if (!.has1(f, i)) return(.mul(list(f, .add(list(.sub(b, a), S_One)))))
  if (f$type == "Mul") {
    sp <- .as_indep_mul(f, i); without <- sp[[1L]]; with <- sp[[2L]]
    if (!is_one_num(without)) {
      s <- .eval_sum_symbolic(with, i, a, b)
      if (.truthy(s)) {
        r <- .mul(list(without, s))
        if (!is_nan(r)) return(r)
      }
    } else {
      lr <- .as_two_terms(f); L <- lr[[1L]]; R <- lr[[2L]]
      if (!.has1(L, i)) {
        sR <- .eval_sum_symbolic(R, i, a, b)
        if (.truthy(sR)) return(.mul(list(L, sR)))
      }
      if (!.has1(R, i)) {
        sL <- .eval_sum_symbolic(L, i, a, b)
        if (.truthy(sL)) return(.mul(list(sL, R)))
      }
    }
  }
  f <- .sum_apart(f, i)
  if (f$type == "Add") {
    lr <- .as_two_terms(f); L <- lr[[1L]]; R <- lr[[2L]]
    lrsum <- .telescopic(L, R, i, a, b)
    if (.truthy(lrsum)) return(lrsum)
    sp <- .as_indep(f, list(i), as_add = TRUE); without <- sp[[1L]]; with <- sp[[2L]]
    if (!is_zero_num(without)) {
      s <- .eval_sum_symbolic(with, i, a, b)
      if (.truthy(s)) {
        r <- .add(list(.mul(list(without, .add(list(.sub(b, a), S_One)))), s))
        if (!is_nan(r)) return(r)
      }
    } else {
      lsum <- .eval_sum_symbolic(L, i, a, b)
      rsum <- .eval_sum_symbolic(R, i, a, b)
      if (!is.null(lsum) && !is.null(rsum)) {
        r <- .add(list(lsum, rsum))
        if (!is_nan(r)) return(r)
      }
    }
  }
  m <- sym_match(f, .pow(i, Wild("n")))
  if (!is.null(m)) {
    n <- m[["n"]]
    if (n$type == "Integer") {
      if (n$value >= 0) {
        if ((is_const(b, "oo") && !is_const(a, "-oo")) || (is_const(a, "-oo") && !is_const(b, "oo"))) return(S_Infinity)
        n1 <- .add(list(n, S_One))
        return(expand(.div(.sub(fcall("bernoulli", list(n1, .add(list(b, S_One)))), fcall("bernoulli", list(n1, a))), n1)))
      } else if (a$type == "Integer" && a$value >= 1) {
        if (is_negone_num(n)) return(.sub(fcall("harmonic", list(b)), fcall("harmonic", list(.sub(a, S_One)))))
        an <- .neg(n)
        return(.sub(fcall("harmonic", list(b, an)), fcall("harmonic", list(.sub(a, S_One), an))))
      }
    }
  }
  infs <- function(z) .has1(z, S_Infinity) || .has1(z, S_NegativeInfinity)
  if (!(infs(a) || infs(b))) {
    c1 <- Wild("c1", exclude = list(i)); c2 <- Wild("c2", exclude = list(i)); c3 <- Wild("c3", exclude = list(i))
    m <- sym_match(powsimp(f), .pow(c1, Wild("wexp")))
    if (!is.null(m)) {
      me <- sym_match(expand(m[["wexp"]]), .add(list(.mul(list(c2, i)), c3)))
      if (!is.null(me)) {
        cv <- m[["c1"]]
        p <- .pow(cv, me[["c3"]]); q <- .pow(cv, me[["c2"]])
        r <- .div(.mul(list(p, .sub(.pow(q, a), .pow(q, .add(list(b, S_One)))))), .sub(S_One, q))
        l <- .mul(list(p, .add(list(.sub(b, a), S_One))))
        return(.piecewise(list(list(l, Eq(q, S_One)), list(r, S_true))))
      }
    }
    r <- .gosper_sum_f(f, i, a, b)
    if (!is.null(r) && r$type %in% c("Mul", "Add")) {
      lfs <- vapply(c(free_symbols(a), free_symbols(b)), `[[`, "", "key")
      non_limit <- Filter(function(s) !(s$key %in% lfs), free_symbols(r))
      den <- denom(together(r))
      dk <- vapply(free_symbols(den), `[[`, "", "key")
      den_sym <- Filter(function(s) s$key %in% dk, non_limit)
      args <- list()
      for (v in ordered(den_sym)) {
        s <- tryCatch(solve(den, v), error = function(err) NULL)
        if (is.null(s)) next
        if (length(s)) {
          mm <- Eq(v, s[[1L]])
          if (!.is_false(mm))
            args[[length(args) + 1L]] <- list(.sum_doit(.sum_node(subs(f_orig, v, s[[1L]]), list(list(i, a, b)))), mm)
        }
        break
      }
      args[[length(args) + 1L]] <- list(r, S_true)
      return(.piecewise(args))
    }
    if (!is.null(r) && !is_nan(r)) return(r)
  }
  h <- .eval_sum_hyper(f_orig, i, a, b)
  if (!is.null(h)) return(h)
  r <- tryCatch(.eval_sum_residue(f_orig, i, a, b), error = function(err) NULL)
  if (!is.null(r)) return(r)
  factored <- sym_factor(f_orig)
  if (!identical(factored$key, f_orig$key)) return(.eval_sum_symbolic(factored, i, a, b))
  NULL
}

.telescopic_direct <- function(L, R, n, i, a, b)
  .add(lapply(seq_len(n) - 1L, function(m) .add(list(subs(L, i, .add(list(a, .int(m)))), subs(R, i, .sub(b, .int(m)))))))

.telescopic <- function(L, R, i, a, b) {
  if (L$type == "Add" || R$type == "Add") return(NULL)
  k <- Wild("k")
  sol <- sym_match(.neg(R), subs(L, i, .add(list(i, k))))
  s <- NULL
  if (!is.null(sol) && !is.null(sol[["k"]])) {
    s <- sol[["k"]]
    if (!(s$type == "Integer" && is_zero_num(.add(list(subs(L, i, .add(list(i, s))), R))))) s <- NULL
  }
  if (is.null(s)) {
    m <- Dummy("m")
    sl <- tryCatch(solve(.add(list(subs(L, i, .add(list(i, m))), R)), m), error = function(err) NULL)
    if (is.null(sl)) return(NULL)
    sl <- Filter(function(si) si$type == "Integer" && isTRUE(ask_fact(expand(.add(list(subs(L, i, .add(list(i, si))), R))), "zero")), sl)
    if (length(sl) != 1L) return(NULL)
    s <- sl[[1L]]
  }
  if (s$value < 0) return(.telescopic_direct(R, L, as.integer(as.character(abs(s$value))), i, a, b))
  if (s$value > 0) return(.telescopic_direct(L, R, as.integer(as.character(s$value)), i, a, b))
  NULL
}

# ---------------------------------------------------------------------------
# hypersimp and Gosper's algorithm

# rewrite(gamma) with SymPy's Piecewise forms for the Pochhammer symbols
.rewrite_gamma_sympy <- function(e) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  g <- function(v) fcall("gamma", list(v))
  a <- z$args
  switch(z$name,
    factorial = g(.add(list(a[[1L]], S_One))),
    binomial = .div(g(.add(list(a[[1L]], S_One))), .mul(list(g(.add(list(a[[2L]], S_One))), g(.add(list(.sub(a[[1L]], a[[2L]]), S_One)))))),
    RisingFactorial = .piecewise(list(
      list(.div(g(.add(list(a[[1L]], a[[2L]]))), g(a[[1L]])), Gt(a[[1L]], 0)),
      list(.div(.mul(list(.pow(S_NegativeOne, a[[2L]]), g(.sub(S_One, a[[1L]])))), g(.add(list(.neg(a[[2L]]), .neg(a[[1L]]), S_One)))), S_true))),
    FallingFactorial = .piecewise(list(
      list(.div(g(.add(list(a[[1L]], S_One))), g(.add(list(.neg(a[[2L]]), a[[1L]], S_One)))), Ge(a[[1L]], 0)),
      list(.div(.mul(list(.pow(S_NegativeOne, a[[2L]]), g(.sub(a[[2L]], a[[1L]])))), g(.neg(a[[1L]]))), S_true))),
    beta = .div(.mul(list(g(a[[1L]]), g(a[[2L]]))), g(.add(list(a[[1L]], a[[2L]])))),
    z)
})

.hypersimp <- function(f, k) {
  g <- .div(subs(f, k, .add(list(k, S_One))), f)
  g <- .rewrite_gamma_sympy(g)
  if (.has_type(g, "Piecewise")) {
    g <- piecewise_fold(g)
    if (g$type == "Piecewise") g <- g$args[[length(g$args)]]$args[[1L]]
  }
  g <- expand_func(g)
  g <- powsimp(g, deep = TRUE, combine = "exp")
  if (.is_ratfun(g, k)) simplify(g, ratio = Inf) else NULL
}

.monic_in <- function(p, n) {
  cl <- .coeff_list(p, n)
  if (!length(cl)) return(S_Zero)
  lc <- cl[[1L]]
  expand(.add(lapply(seq_along(cl), function(j) .mul(list(cancel(.div(cl[[j]], lc)), .pow(n, .int(length(cl) - j)))))))
}
.lc_in <- function(p, n) { cl <- .coeff_list(p, n); if (length(cl)) cl[[1L]] else S_Zero }
.deg_in <- function(p, n) length(.coeff_list(p, n)) - 1L
.quo_in <- function(p, d, n) {
  r <- .expr_poly_divmod(.coeff_list(p, n), .coeff_list(d, n))
  .coeffs_to_expr(r$q, n)
}

# gosper_normal over the field of the coefficients
.gosper_normal <- function(p, q, n) {
  a <- .lc_in(p, n); A <- .monic_in(p, n)
  b <- .lc_in(q, n); B <- .monic_in(q, n)
  C <- S_One; Z <- cancel(.div(a, b))
  h <- Dummy("h")
  R <- resultant(A, subs(B, n, .add(list(n, h))), n)
  rs <- tryCatch(roots(R, h), error = function(err) list())
  rts <- Filter(function(r) r$type == "Integer" && r$value >= 0, as.list(rs))
  if (length(rts)) rts <- rts[order(vapply(rts, function(r) as.numeric(r$value), 0))]
  for (iv in rts) {
    d <- .monic_in(gcd(A, subs(B, n, .add(list(n, iv))), n), n)
    A <- .quo_in(A, d, n)
    B <- .quo_in(B, subs(d, n, .sub(n, iv)), n)
    for (j in seq_len(as.integer(as.character(iv$value)))) C <- expand(.mul(list(C, subs(d, n, .sub(n, .int(j))))))
  }
  list(expand(.mul(list(A, Z))), B, C)
}

.gosper_term_f <- function(f, n) {
  r <- .hypersimp(f, n)
  if (is.null(r)) return(NULL)
  nd <- .as_numer_denom(r)
  abc <- .gosper_normal(expand(nd[[1L]]), expand(nd[[2L]]), n)
  A <- abc[[1L]]; B <- expand(subs(abc[[2L]], n, .sub(n, S_One))); C <- abc[[3L]]
  N <- .deg_in(A, n); M <- .deg_in(B, n); K <- .deg_in(C, n)
  lcA <- .lc_in(A, n); lcB <- .lc_in(B, n)
  if (N != M || !isTRUE(ask_fact(cancel(.sub(lcA, lcB)), "zero"))) {
    d <- K - max(N, M)
  } else if (N == 0L) {
    d <- K + 1L
  } else {
    d1 <- K - N + 1L
    nth <- function(p, k) { cl <- .coeff_list(p, n); cl[[length(cl) - k]] }
    d2 <- cancel(.div(.sub(nth(B, N - 1L), nth(A, N - 1L)), lcA))
    d2 <- if (d2$type == "Integer") as.integer(as.character(d2$value)) else -1L
    d <- max(d1, d2)
  }
  if (d < 0L) return(NULL)
  cs <- lapply(0:d, function(k) Dummy(paste0("c", k)))
  x <- .add(lapply(0:d, function(k) .mul(list(cs[[k + 1L]], .pow(n, .int(d - k))))))
  H <- expand(.sub(.sub(.mul(list(A, subs(x, n, .add(list(n, S_One))))), .mul(list(B, x))), C))
  eqs <- Filter(function(c) !is_zero_num(c), .coeff_list(H, n))
  sol <- .linsolve_exprs(eqs, cs)
  if (is.null(sol)) return(NULL)
  x <- xreplace(x, sol)
  solved <- vapply(sol, function(p) p[[1L]]$key, "")
  for (c in cs) if (!(c$key %in% solved)) x <- subs(x, c, S_Zero)
  if (isTRUE(ask_fact(x, "zero")) || is_zero_num(x)) return(NULL)
  .div(.mul(list(B, x)), C)
}

.gosper_sum_f <- function(f, k, a = NULL, b = NULL) {
  g <- tryCatch(.gosper_term_f(f, k), error = function(err) NULL)
  if (is.null(g)) return(NULL)
  if (is.null(a)) return(sym_factor(.mul(list(f, g))))
  F1 <- .mul(list(f, .add(list(g, S_One)))); F0 <- .mul(list(f, g))
  result <- .sub(subs(F1, k, b), subs(F0, k, a))
  if (is_nan(result)) {
    result <- tryCatch(.sub(limit(F1, k, b), limit(F0, k, a)), error = function(err) NULL)
    if (is.null(result)) return(NULL)
  }
  sym_factor(result)
}

# ---------------------------------------------------------------------------
# hypergeometric sums

.as_coeff_mul_dep <- function(e, i) {
  if (e$type == "Mul") {
    args <- .sympy_args(e)
    dep <- vapply(args, function(t) .has1(t, i), TRUE)
    return(list(.mul(args[!dep]), args[dep]))
  }
  if (.has1(e, i)) list(S_One, list(e)) else list(e, list())
}

.eval_sum_hyper1 <- function(f, i, a) {
  if (!is_zero_num(a)) return(.eval_sum_hyper1(subs(f, i, .add(list(i, a))), i, S_Zero))
  if (is_zero_num(subs(f, i, S_Zero))) {
    if (is_zero_num(simplify(subs(f, i, Dummy("i", integer = TRUE, positive = TRUE))))) return(list(S_Zero, S_true))
    return(.eval_sum_hyper1(subs(f, i, .add(list(i, S_One))), i, S_Zero))
  }
  hs <- .hypersimp(f, i)
  if (is.null(hs)) return(NULL)
  if (hs$type == "Float") hs <- nsimplify(hs)
  nd <- fraction(sym_factor(hs))
  ab <- list(); fl <- list()
  for (k in 1:2) { cm <- .as_coeff_mul_dep(nd[[k]], i); ab[[k]] <- cm[[1L]]; fl[[k]] <- cm[[2L]] }
  params <- list(list(), list())
  for (k in 1:2) for (fac in fl[[k]]) {
    mul <- S_One
    if (fac$type == "Pow") {
      mul <- fac$args[[2L]]; fac <- fac$args[[1L]]
      if (mul$type != "Integer") return(NULL)
    }
    cl <- tryCatch(.coeff_list(fac, i), error = function(err) NULL)
    if (is.null(cl) || length(cl) != 2L) return(NULL)
    m <- cl[[1L]]; nn <- cl[[2L]]
    ab[[k]] <- .mul(list(ab[[k]], .pow(m, mul)))
    params[[k]] <- c(params[[k]], rep(list(.div(nn, m)), as.integer(as.character(mul$value))))
  }
  ap <- c(params[[1L]], list(S_One)); bq <- params[[2L]]
  x <- .div(ab[[1L]], ab[[2L]])
  h <- hyper(ap, bq, x)
  f <- combsimp(f)
  conv <- if (h$type == "Function" && h$name == "hyper") .hyper_convergence(h) else S_true
  list(.mul(list(subs(f, i, S_Zero), hyperexpand(h))), conv)
}

.ILLEGAL <- c("nan", "oo", "-oo", "zoo")
.n_illegal <- function(e) {
  n <- 0L
  walk <- function(z) {
    if (z$type == "Constant" && z$name %in% .ILLEGAL) n <<- n + 1L
    for (a in z$args) if (is_sym(a)) walk(a)
  }
  walk(e)
  n
}
.n_illegal_pair <- function(p) sum(vapply(p, function(z) z$type == "Constant" && z$name %in% .ILLEGAL, TRUE))

.eval_sum_hyper <- function(f, i, a, b) {
  if (is.null(.hypersimp(f, i))) return(NULL)
  if (.sub(b, a)$type == "Integer") return(NULL)
  old_sum <- .sum_node(f, list(list(i, a, b)))
  if (!is_const(b, "oo")) {
    if (is_const(a, "-oo")) {
      res <- .eval_sum_hyper1(subs(f, i, .neg(i)), i, .neg(b))
      if (!is.null(res)) return(.piecewise(list(res, list(old_sum, S_true))))
      return(NULL)
    }
    had <- .n_illegal(f)
    res1 <- .eval_sum_hyper1(f, i, a)
    if (is.null(res1) || .n_illegal_pair(res1) > had) return(NULL)
    res2 <- .eval_sum_hyper1(f, i, .add(list(b, S_One)))
    if (is.null(res2) || .n_illegal_pair(res2) > had) return(NULL)
    cond <- And(res1[[2L]], res2[[2L]])
    if (.is_false(cond)) return(NULL)
    return(.piecewise(list(list(.sub(res1[[1L]], res2[[1L]]), cond), list(old_sum, S_true))))
  }
  if (is_const(a, "-oo")) {
    res1 <- .eval_sum_hyper1(subs(f, i, .neg(i)), i, S_One)
    res2 <- .eval_sum_hyper1(f, i, S_Zero)
    if (is.null(res1) || is.null(res2)) return(NULL)
    cond <- And(res1[[2L]], res2[[2L]])
    if (.is_false(cond)) return(NULL)
    fsc <- free_symbols(cond)
    empty <- length(fsc) == 1L && tryCatch(identical(.cond_set(cond, fsc[[1L]])$key, S_EmptySet$key), error = function(err) FALSE)
    if (empty) return(NULL)
    return(.piecewise(list(list(.add(list(res1[[1L]], res2[[1L]])), cond), list(old_sum, S_true))))
  }
  res <- .eval_sum_hyper1(f, i, a)
  if (is.null(res)) return(NULL)
  if (.is_false(res[[2L]])) {
    if (!length(free_symbols(res[[1L]]))) {
      fd <- subs(f, i, .add(list(Dummy("i", integer = TRUE, positive = TRUE), a)))
      if (isTRUE(ask_fact(fd, "positive")) || isTRUE(ask_fact(fd, "zero"))) return(S_Infinity)
      if (isTRUE(ask_fact(fd, "negative"))) return(S_NegativeInfinity)
    }
    return(NULL)
  }
  .piecewise(list(res, list(old_sum, S_true)))
}

# ---------------------------------------------------------------------------
# sums of rational functions by residues

.eval_sum_residue <- function(f, i, a, b) {
  fs <- free_symbols(f)
  if (any(vapply(fs, function(s) !identical(s$key, i$key), TRUE))) return(NULL)
  okb <- function(z) z$type == "Integer" || is_const(z, "oo") || is_const(z, "-oo")
  if (!okb(a) || !okb(b)) return(NULL)
  if (!(is_inf_const(a) || is_inf_const(b))) return(NULL)
  if (isTRUE(ask_fact(.sub(a, b), "extended_positive")))
    return(.neg(.eval_sum_residue(f, i, .add(list(b, S_One)), .sub(a, S_One))))
  monoms <- function(p) { cl <- .coeff_list(p, i); d <- length(cl) - 1L; d - which(!vapply(cl, is_zero_num, TRUE)) + 1L }
  is_even_function <- function(nm, dn) {
    a1 <- monoms(nm); a2 <- monoms(dn)
    (all(a1 %% 2 == 0) && all(a2 %% 2 == 0)) || (all(a1 %% 2 == 1) && all(a2 %% 2 == 1))
  }
  match_rational <- function(g) {
    nd <- .as_numer_denom(g)
    if (!isTRUE(is_polynomial(nd[[1L]], i)) || !isTRUE(is_polynomial(nd[[2L]], i))) return(NULL)
    list(expand(nd[[1L]]), expand(nd[[2L]]))
  }
  m <- match_rational(f); alternating <- FALSE
  if (is.null(m)) {
    m <- match_rational(.div(f, .pow(S_NegativeOne, i)))
    if (is.null(m)) return(NULL)
    alternating <- TRUE
  }
  numer <- m[[1L]]; denom <- m[[2L]]
  if (.deg_in(denom, i) - .deg_in(numer, i) < 2L) return(NULL)
  both_inf <- is_const(a, "-oo") && is_const(b, "oo")
  if (!is_even_function(numer, denom)) {
    n <- .deg_in(denom, i); cl <- .coeff_list(denom, i)
    shift <- .div(.div(.neg(cl[[2L]]), cl[[1L]]), .int(n))
    if (is_zero_num(shift)) return(NULL)
    if (shift$type == "Integer") {
      numer <- expand(subs(numer, i, .add(list(i, shift)))); denom <- expand(subs(denom, i, .add(list(i, shift))))
      if (!is_even_function(numer, denom)) return(NULL)
      f <- if (alternating) .mul(list(.pow(S_NegativeOne, i), .div(.mul(list(.pow(S_NegativeOne, shift), numer)), denom)))
           else .div(numer, denom)
      return(.eval_sum_residue(f, i, .sub(a, shift), .sub(b, shift)))
    } else if (!both_inf) return(NULL)
  }
  get_poles <- function() {
    rts <- all_roots(sqf_part(denom, i), i)
    ints <- list(); nonints <- list()
    for (r in rts) {
      t <- ask_fact(r, "integer")
      if (is.na(t)) return(NULL)
      if (t) ints[[length(ints) + 1L]] <- r else nonints[[length(nonints) + 1L]] <- r
    }
    list(ints, nonints)
  }
  residues <- function(poles) lapply(poles, function(p) .sum_get_residue(numer, denom, p, alternating, i))
  sumlist <- function(l) Reduce(function(x, y) .add(list(x, y)), l, S_Zero)
  if (both_inf) {
    poles <- get_poles()
    if (is.null(poles) || length(poles[[1L]])) return(NULL)
    return(.mul(list(.neg(S_Pi), sumlist(residues(poles[[2L]])))))
  }
  if (is_const(a, "-oo")) { a0 <- a; a <- .neg(b); b <- .neg(a0) }
  poles <- get_poles()
  if (is.null(poles)) return(NULL)
  ints <- poles[[1L]]; nonints <- poles[[2L]]
  av <- as.integer(as.character(a$value))
  if (length(ints)) {
    iv <- vapply(ints, function(r) as.integer(as.character(r$value)), 0L)
    if (length(iv) != max(iv) - min(iv) + 1L) return(NULL)
    if (av <= max(iv)) return(NULL)
  }
  full_sum <- .mul(list(.neg(S_Pi), sumlist(residues(c(ints, nonints)))))
  fat <- function(k) xreplace(f, list(list(i, .int(k))))
  if (!length(ints)) {
    half_sum <- .div(.add(list(full_sum, fat(0L))), S_Two)
    neg <- if (av < 0) lapply(av:-1L, fat) else list()
    pos <- if (av > 0) lapply(0:(av - 1L), fat) else list()
    return(.sub(.add(list(half_sum, sumlist(neg))), sumlist(pos)))
  }
  half_sum <- .div(full_sum, S_Two)
  lo <- max(iv) + 1L
  extr <- if (lo <= av - 1L) lapply(lo:(av - 1L), fat) else list()
  .sub(half_sum, sumlist(extr))
}

.sum_get_residue <- function(numer, denom, a, alternating, x) {
  z <- Dummy("z")
  nz <- subs(numer, x, z); dz <- subs(denom, x, z)
  if (.mul(list(S_Two, a))$type == "Integer") {
    trig <- if (!alternating) fcall("cot", list(.mul(list(S_Pi, z)))) else fcall("csc", list(.mul(list(S_Pi, z))))
    return(residue(.mul(list(trig, .div(nz, dz))), z, a))
  }
  p <- 0L; d <- dz
  repeat {
    r <- .expr_poly_divmod(.coeff_list(d, z), .coeff_list(.sub(z, a), z))
    if (!all(vapply(r$r, function(c) isTRUE(ask_fact(c, "zero")) || is_zero_num(c), TRUE)) || is_zero_num(d)) break
    d <- .coeffs_to_expr(r$q, z); p <- p + 1L
  }
  trig <- if (!alternating) .inv(fcall("tan", list(.mul(list(S_Pi, z))))) else .inv(fcall("sin", list(.mul(list(S_Pi, z)))))
  expr <- .mul(list(trig, .div(nz, d)))
  if (p > 1L) expr <- .div(diff(expr, z, p - 1L), fcall("factorial", list(.int(p - 1L))))
  xreplace(expr, list(list(.mul(list(S_Pi, z)), expand_mul(.mul(list(S_Pi, a)))), list(z, a)))
}

# Sum.eval_zeta_function
.eval_zeta_function <- function(self, f, i, a, b) {
  if (.is_comparable(a) && .is_comparable(b) && isTRUE(ask_fact(.sub(a, b), "extended_positive")))
    return(.eval_zeta_function(self, f, i, .add(list(b, S_One)), .sub(a, S_One)))
  if (!is_const(b, "oo")) return(NULL)
  w <- Wild("w", exclude = list(i)); y <- Wild("y", exclude = list(i)); z <- Wild("z", exclude = list(i))
  r <- sym_match(f, .pow(.add(list(.mul(list(w, i)), y)), .neg(z)))
  if (is.null(r) || !length(r)) return(NULL)
  coeff <- .div(S_One, .pow(r[["w"]], r[["z"]]))
  s <- r[["z"]]
  q <- .add(list(.div(r[["y"]], r[["w"]]), a))
  cond <- And(Not(Contains(.neg(q), Naturals0)), Gt(fcall("re", list(s)), S_One))
  .piecewise(list(list(.mul(list(coeff, fcall("zeta", list(s, q)))), cond), list(self, S_true)))
}

# ---------------------------------------------------------------------------
# Kronecker deltas (sympy/concrete/delta.py)

.is_simple_delta <- function(delta, index) {
  if (!(delta$type == "Function" && delta$name == "KroneckerDelta" && .has1(delta, index))) return(FALSE)
  d <- .sub(delta$args[[1L]], delta$args[[2L]])
  isTRUE(is_polynomial(d, index)) && .deg_in(d, index) == 1L
}
.has_simple_delta <- function(expr, index) {
  if (!.has_function(expr, "KroneckerDelta")) return(FALSE)
  if (.is_simple_delta(expr, index)) return(TRUE)
  if (expr$type %in% c("Add", "Mul")) return(any(vapply(expr$args, .has_simple_delta, TRUE, index = index)))
  FALSE
}
.extract_delta <- function(expr, index) {
  if (!.has_simple_delta(expr, index)) return(list(NULL, expr))
  if (expr$type == "Function" && expr$name == "KroneckerDelta") return(list(expr, S_One))
  if (expr$type != "Mul") stop("Incorrect expr")
  delta <- NULL; terms <- list()
  for (a in expr$args) {
    if (is.null(delta) && .is_simple_delta(a, index)) delta <- a else terms[[length(terms) + 1L]] <- a
  }
  list(delta, .mul(terms))
}
.expand_delta <- function(expr, index) {
  if (expr$type != "Mul") return(expr)
  delta <- FALSE; terms <- list(S_One); add <- TRUE
  for (h in expr$args) {
    if (!delta && h$type == "Add" && .has_simple_delta(h, index)) {
      delta <- TRUE
      terms <- lapply(h$args, function(t) .mul(list(terms[[1L]], t)))
    } else terms <- lapply(terms, function(t) .mul(list(t, h)))
  }
  .add(terms)
}

#' @rdname summation
#' @param limit A limit `list(i, a, b)`.
#' @param no_piecewise Do not wrap the result in a Piecewise.
#' @export
deltasummation <- function(f, limit, no_piecewise = FALSE) {
  f <- as_sym(f); limit <- lapply(limit, as_sym)
  x <- limit[[1L]]; a <- limit[[2L]]; b <- limit[[3L]]
  if (isTRUE(ask_fact(.sub(b, a), "negative"))) return(S_Zero)
  if (!.has_function(f, "KroneckerDelta")) return(summation(f, limit))
  g <- .expand_delta(f, x)
  if (g$type == "Add") return(piecewise_fold(.add(lapply(g$args, deltasummation, limit = limit, no_piecewise = no_piecewise))))
  de <- .extract_delta(g, x); delta <- de[[1L]]; expr <- de[[2L]]
  if (is.null(delta)) return(summation(f, limit))
  solns <- solve(.sub(delta$args[[1L]], delta$args[[2L]]), x)
  if (!length(solns)) return(S_Zero)
  if (length(solns) != 1L) return(.sum_node(f, list(limit)))
  value <- solns[[1L]]
  if (no_piecewise) return(subs(expr, x, value))
  .piecewise(list(list(subs(expr, x, value), as_relational(Interval(a, b), value)), list(S_Zero, S_true)))
}
