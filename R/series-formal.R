# Formal power series, after sympy/series/formal.py: the rational
# algorithm, the hypergeometric algorithm (simple differential equations,
# their recurrence equations and rsolve_hypergeometric), compute_fps and
# the FormalPowerSeries object with its truncation, arithmetic, derivative,
# integral, product, composition and inverse.

.fps_add_args <- function(e) if (e$type == "Add") .sympy_args(e) else list(e)
.fps_coeff_add <- function(e, x) {
  # as_coeff_add(x): the part free of x and the terms that depend on x
  ts <- .fps_add_args(e)
  dep <- vapply(ts, function(t) .has1(t, x), NA)
  list(.add(ts[!dep]), ts[dep])
}
.fps_is_ratfun <- function(e, x) isTRUE(.is_rational_function_in(e, list(x)))
.fps_fact <- function(z) fcall("factorial", list(z))

# rational_algorithm
.fps_rational_algorithm <- function(f, x, k, order = 4L, full = FALSE) {
  d <- f
  ds <- list()
  for (i in 0:order) {
    if (i) d <- diff(d, x)
    if (.fps_is_ratfun(d, x)) {
      coeff <- S_Zero; sep <- S_Zero
      terms <- apart(d, x, full = full)
      if (.has_type(terms, "RootSum")) terms <- doit(terms)
      for (t in .fps_add_args(terms)) {
        nd <- .as_numer_denom(t)
        num <- nd[[1L]]; den <- nd[[2L]]
        if (!.has1(den, x)) { sep <- .add(list(sep, t)); next }
        if (den$type == "Mul") {
          ind <- .as_independent(den, x)
          den <- ind[[2L]]; num <- .div(num, ind[[1L]])
        }
        be <- .as_base_exp(den)
        den <- be[[1L]]; j <- be[[2L]]
        ca <- .fps_coeff_add(den, x)
        a <- ca[[1L]]
        if (is_zero_num(a)) { sep <- .add(list(sep, t)); next }
        xc <- coeff(ca[[2L]][[1L]], x)
        a <- .div(a, .neg(xc))
        num <- .div(num, .pow(xc, j))
        bin <- .div(.fps_fact(.add(list(j, k, S_NegativeOne))), .mul(list(.fps_fact(k), .fps_fact(.sub(j, S_One)))))
        ak <- .div(.mul(list(.pow(S_NegativeOne, j), num, bin)), .pow(a, .add(list(j, k))))
        coeff <- .add(list(coeff, ak))
      }
      if (isTRUE(ask_fact(coeff, "zero")) || is_zero_num(coeff)) return(NULL)
      if (.has1(coeff, x) || .contains_key(coeff, S_ComplexInfinity$key) || .contains_key(coeff, S_Infinity$key) || is_nan(coeff) || .contains_key(coeff, S_NaN$key)) return(NULL)
      for (j in seq_len(i) - 1L) {
        coeff <- .div(coeff, .add(list(k, .int(j + 1L))))
        sep <- integrate(sep, x)
        last <- ds[[length(ds)]]; ds <- ds[-length(ds)]
        sep <- .add(list(sep, limit(.sub(last, sep), x, S_Zero)))
      }
      return(list(subs(coeff, k, .sub(k, .int(i))), sep, .int(i)))
    }
    ds[[length(ds) + 1L]] <- d
  }
  NULL
}

# rational_independent
.fps_rational_independent <- function(terms, x) {
  if (!length(terms)) return(list())
  ind <- terms[1L]
  for (t in terms[-1L]) {
    n <- .as_independent(t, x)[[2L]]
    placed <- FALSE
    for (i in seq_along(ind)) {
      dd <- .as_independent(ind[[i]], x)[[2L]]
      q <- cancel(.div(n, dd))
      if (.fps_is_ratfun(q, x)) { ind[[i]] <- .add(list(ind[[i]], t)); placed <- TRUE; break }
    }
    if (!placed) ind[[length(ind) + 1L]] <- t
  }
  ind
}

# exp_re and hyper_re: the recurrence equation of a differential equation
.fps_g_parts <- function(DE, gname) {
  lapply(.fps_add_args(DE), function(t) {
    fs <- if (t$type == "Mul") t$args else list(t)
    isg <- vapply(fs, function(z) .has_fun(z, gname), NA)
    d <- .mul(fs[isg]); c <- .mul(fs[!isg])
    j <- if (d$type == "Derivative") sum(vapply(d$args[-1L], function(v) if (v$type == "Tuple") as.integer(num_value_dbl(v$args[[2L]])) else 1L, 0L)) else 0L
    list(c, j)
  })
}
.fps_exp_re <- function(DE, gname, k) {
  RE <- S_Zero; mini <- NULL
  for (p in .fps_g_parts(DE, gname)) {
    j <- p[[2L]]
    if (is.null(mini) || j < mini) mini <- j
    RE <- .add(list(RE, .mul(list(p[[1L]], fcall(gname, list(.add(list(k, .int(j)))))))))
  }
  if (!is.null(mini) && mini != 0L) RE <- subs(RE, k, .sub(k, .int(mini)))
  RE
}
.fps_hyper_re <- function(DE, gname, k, x) {
  RE <- S_Zero; mini <- NULL
  for (p in .fps_g_parts(expand(DE), gname)) {
    ci <- .as_independent(p[[1L]], x)
    c <- ci[[1L]]; v <- ci[[2L]]
    l <- .so_as_coeff_exponent(v, x)[[2L]]
    j <- p[[2L]]
    RE <- .add(list(RE, .mul(list(c, fcall("RisingFactorial", list(.add(list(k, S_One, .neg(l))), .int(j))),
                                  fcall(gname, list(.add(list(k, .int(j), .neg(l)))))))))
    jl <- .sub(.int(j), l)
    if (is.null(mini) || .is_true(Lt(jl, mini))) mini <- jl
  }
  RE <- subs(RE, k, .sub(k, mini))
  apps <- .fps_g_apps(RE, gname)
  if (length(apps)) RE <- collect(RE, apps)
  RE
}
.fps_g_apps <- function(e, gname) {
  out <- list()
  walk <- function(z) { if (z$type == "Function" && z$name == gname) out[[z$key]] <<- z; for (a in z$args) walk(a) }
  walk(e)
  unname(out)
}

# the hypergeometric recurrence solver
.fps_transformation_a <- function(f, x, P, Q, k, m, shift)
  list(.mul(list(f, .pow(x, .neg(shift)))), subs(P, k, .add(list(k, shift))), subs(Q, k, .add(list(k, shift))), m)
.fps_transformation_c <- function(f, x, P, Q, k, m, scale)
  list(subs(f, x, .pow(x, scale)), subs(P, k, .div(k, scale)), subs(Q, k, .div(k, scale)), .mul(list(m, scale)))
.fps_transformation_e <- function(f, x, P, Q, k, m)
  list(diff(f, x), .mul(list(subs(P, k, .add(list(k, S_One))), .add(list(k, m, S_One)))), .mul(list(subs(Q, k, .add(list(k, S_One))), .add(list(k, S_One)))), m)
.fps_roots <- function(p, k) {
  r <- roots(p, k)
  if (is.null(r) || !length(r)) return(list())
  r
}
.fps_compute_formula <- function(f, x, P, Q, k, m, k_max) {
  sol <- list()
  mi <- as.integer(num_value_dbl(m)); km <- as.integer(num_value_dbl(k_max))
  for (i in seq.int(km + 1L, km + mi, length.out = max(mi, 0L))) {
    if (i < 0L) next
    r <- .div(limit(diff(f, x, i), x, S_Zero), .fps_fact(.int(i)))
    if (isTRUE(ask_fact(r, "zero")) || is_zero_num(r)) next
    kterm <- .add(list(.mul(list(m, k)), .int(i)))
    res <- r
    p <- subs(P, k, kterm); q <- subs(Q, k, kterm)
    c1 <- .so_leadterm(subs(p, k, .inv(k)), k)[[1L]]
    c2 <- .so_leadterm(subs(q, k, .inv(k)), k)[[1L]]
    res <- .mul(list(res, .pow(.div(.neg(c1), c2), k)))
    for (rt in .fps_root_items(p, k)) res <- .mul(list(res, .pow(fcall("RisingFactorial", list(.neg(rt[[1L]]), k)), rt[[2L]])))
    for (rt in .fps_root_items(q, k)) res <- .div(res, .pow(fcall("RisingFactorial", list(.neg(rt[[1L]]), k)), rt[[2L]]))
    sol[[length(sol) + 1L]] <- list(res, kterm)
  }
  sol
}
.fps_root_items <- function(p, k) {
  r <- roots(p, k)
  if (is.null(r) || !length(r)) return(list())
  mult <- attr(r, "multiplicity")
  lapply(seq_along(r), function(i) list(r[[i]], .int(mult[i])))
}
# linsolve as a substitution: SymPy's dict(zip(syms, solution))
.fps_linsolve_subs <- function(eqs, syms) {
  sol <- tryCatch(linsolve(eqs, syms), error = function(err) NULL)
  if (is.null(sol) || !length(sol)) return(list())
  nm <- vapply(syms, function(s) s$name, "")
  out <- list()
  for (i in seq_along(syms)) if (!is.null(sol[[nm[i]]])) out[[length(out) + 1L]] <- list(syms[[i]], sol[[nm[i]]])
  out
}
.fps_rsolve_hypergeometric_ <- function(f, x, P, Q, k, m) {
  proots <- .fps_root_items(P, k); qroots <- .fps_root_items(Q, k)
  allr <- c(proots, qroots)
  dens <- lapply(Filter(function(r) isTRUE(ask_fact(r[[1L]], "rational")), allr), function(r) .as_numer_denom(r[[1L]])[[2L]])
  scale <- if (length(dens)) Reduce(function(a, b) .int_lcm(a, b), dens) else S_One
  tr <- .fps_transformation_c(f, x, P, Q, k, m, scale); f <- tr[[1L]]; P <- tr[[2L]]; Q <- tr[[3L]]; m <- tr[[4L]]
  qroots <- .fps_root_items(Q, k)
  k_min <- if (length(qroots)) do.call(Min, lapply(qroots, `[[`, 1L)) else S_Zero
  shift <- .add(list(k_min, m))
  tr <- .fps_transformation_a(f, x, P, Q, k, m, shift); f <- tr[[1L]]; P <- tr[[2L]]; Q <- tr[[3L]]; m <- tr[[4L]]
  l <- limit(.mul(list(x, f)), x, S_Zero)
  if (!(l$type == "Limit") && !is_zero_num(l)) return(NULL)
  qroots <- .fps_root_items(Q, k)
  k_max <- if (length(qroots)) do.call(Max, lapply(qroots, `[[`, 1L)) else S_Zero
  ind <- S_Zero; mp <- S_NegativeInfinity
  top <- as.integer(num_value_dbl(.add(list(k_max, m))))
  for (i in seq_len(top + 1L) - 1L) {
    r <- .div(limit(diff(f, x, i), x, S_Zero), .fps_fact(.int(i)))
    if (identical(ask_fact(r, "finite"), FALSE)) {
      old_f <- f
      tr <- .fps_transformation_a(f, x, P, Q, k, m, .int(i)); f <- tr[[1L]]; P <- tr[[2L]]; Q <- tr[[3L]]; m <- tr[[4L]]
      tr <- .fps_transformation_e(f, x, P, Q, k, m); f <- tr[[1L]]; P <- tr[[2L]]; Q <- tr[[3L]]; m <- tr[[4L]]
      rr <- .fps_rsolve_hypergeometric_(f, x, P, Q, k, m)
      sol <- rr[[1L]]; ind <- rr[[2L]]; mp <- rr[[3L]]
      sol <- lapply(sol, function(sc) list(.div(sc[[1L]], .mul(list(.add(list(sc[[2L]], S_One)), coeff(as_coeff_add_split(sc[[2L]])[[2L]], k)))), .add(list(sc[[2L]], S_One))))
      sol <- lapply(sol, function(sc) list(sc[[1L]], .add(list(sc[[2L]], .int(i)))))
      ind <- integrate(ind, x)
      ind <- .add(list(ind, limit(.sub(old_f, ind), x, S_Zero)))
      mp <- .add(list(mp, S_One))
      return(list(sol, ind, mp))
    }
    if (!is_zero_num(r)) {
      ind <- .add(list(ind, .mul(list(r, .pow(x, .add(list(.int(i), shift)))))))
      pow_x <- .div(.add(list(.int(i), shift)), scale)
      if (.is_true(Gt(pow_x, mp))) mp <- pow_x
    }
  }
  ind <- subs(ind, x, .pow(x, .inv(scale)))
  sol <- .fps_compute_formula(f, x, P, Q, k, m, k_max)
  sol <- lapply(sol, function(sc) list(sc[[1L]], .add(list(sc[[2L]], shift))))
  sol <- lapply(sol, function(sc) list(sc[[1L]], .div(sc[[2L]], scale)))
  list(sol, ind, mp)
}
.int_lcm <- function(a, b) as_sym(gmp::lcm.bigz(gmp::as.bigz(as.character(num_value_dbl(a))), gmp::as.bigz(as.character(num_value_dbl(b)))))
as_coeff_add_split <- function(e) {
  # Expr.as_coeff_Add: the numeric term and the rest
  if (e$type == "Add" && is_number_type(e$args[[1L]])) return(list(e$args[[1L]], .add(e$args[-1L])))
  if (is_number_type(e)) return(list(e, S_Zero))
  list(S_Zero, e)
}

# rsolve_hypergeometric
.fps_rsolve_hypergeometric <- function(f, x, P, Q, k, m) {
  result <- .fps_rsolve_hypergeometric_(f, x, P, Q, k, m)
  if (is.null(result)) return(NULL)
  sol_list <- result[[1L]]; ind <- result[[2L]]; mp <- result[[3L]]
  conds <- list(); ress <- list()
  for (sc in sol_list) {
    res <- sc[[1L]]
    ca <- as_coeff_add_split(sc[[2L]])
    j <- ca[[1L]]; mk <- ca[[2L]]
    c <- coeff(mk, k)
    if (identical(ask_fact(j, "integer"), FALSE)) {
      res <- .mul(list(res, .pow(x, fcall("frac", list(j)))))
      j <- fcall("floor", list(j))
    }
    res <- subs(res, k, .div(.sub(k, j), c))
    cond <- Eq(fcall("Mod", list(k, c)), fcall("Mod", list(j, c)))
    key <- cond$key
    if (is.null(conds[[key]])) { conds[[key]] <- cond; ress[[key]] <- res } else ress[[key]] <- .add(list(ress[[key]], res))
  }
  pairs <- lapply(names(conds), function(kk) list(ress[[kk]], conds[[kk]]))
  pairs[[length(pairs) + 1L]] <- list(S_Zero, S_true)
  sol <- .piecewise(pairs)
  s <- if (is_const(mp, "-oo")) S_Zero else if (identical(ask_fact(mp, "integer"), FALSE)) fcall("ceiling", list(mp)) else .add(list(mp, S_One))
  if (.is_true(Lt(s, 0))) {
    si <- as.integer(num_value_dbl(s))
    for (kk in seq.int(si, -1L)) ind <- .add(list(ind, subs(.mul(list(sol, .pow(x, k))), k, .int(kk))))
    s <- S_Zero
  }
  list(sol, ind, s)
}

.fps_solve_hyper_RE <- function(f, x, RE, gname, k) {
  terms <- .fps_add_args(RE)
  if (length(terms) == 2L) {
    gs <- .fps_g_apps(RE, gname)
    if (length(gs) != 2L) return(NULL)
    sh <- vapply(gs, function(g) num_value_dbl(.sub(g$args[[1L]], k)), 0)
    gs <- gs[order(sh)]
    P <- coeff(RE, gs[[1L]]); Q <- coeff(RE, gs[[2L]])
    m <- .sub(gs[[2L]]$args[[1L]], gs[[1L]]$args[[1L]])
    return(.fps_rsolve_hypergeometric(f, x, P, Q, k, m))
  }
  NULL
}
.fps_solve_explike_DE <- function(f, x, DE, gname, k) {
  for (p in .fps_g_parts(DE, gname)) if (length(free_symbols(p[[1L]]))) return(NULL)
  RE <- .fps_exp_re(DE, gname, k)
  init <- list()
  nt <- length(.fps_add_args(RE))
  for (i in seq_len(nt) - 1L) {
    if (i) f <- diff(f, x)
    init[[length(init) + 1L]] <- list(.int(i), limit(f, x, S_Zero))
  }
  sol <- tryCatch(rsolve(RE, fcall(gname, list(k)), k, init), error = function(err) NULL)
  if (!is.null(sol) && !is_zero_num(sol)) return(list(.div(sol, .fps_fact(k)), S_Zero, S_Zero))
  NULL
}
.fps_solve_simple <- function(f, x, DE, gname, k) {
  RE <- .fps_hyper_re(DE, gname, k, x)
  init <- list()
  nt <- length(.fps_add_args(RE))
  for (i in seq_len(nt) - 1L) {
    if (i) f <- diff(f, x)
    init[[length(init) + 1L]] <- list(.int(i), .div(limit(f, x, S_Zero), .fps_fact(.int(i))))
  }
  sol <- tryCatch(rsolve(RE, fcall(gname, list(k)), k, init), error = function(err) NULL)
  if (!is.null(sol) && !is_zero_num(sol)) return(list(sol, S_Zero, S_Zero))
  NULL
}

# simpleDE and hyper_algorithm
.fps_hyper_algorithm <- function(f, x, k, order = 4L) {
  gname <- "g"
  g <- Function(gname)
  a <- lapply(seq_len(order) - 1L, function(i) Symbol(paste0("a", i)))
  dn <- function(e, n) if (n == 0L) e else diff(e, x, n)
  des <- list()
  found <- FALSE
  for (kk in seq_len(order)) {
    eq <- .add(c(list(dn(f, kk)), lapply(seq_len(kk) - 1L, function(i) .mul(list(a[[i + 1L]], dn(f, i))))))
    gx <- g(x)
    DE <- .add(c(list(.fps_dg(gx, x, kk)), lapply(seq_len(kk) - 1L, function(i) .mul(list(a[[i + 1L]], .fps_dg(gx, x, i))))))
    eq <- expand(eq)
    terms <- .add_print_order(.fps_add_args(eq))
    ind <- .fps_rational_independent(terms, x)
    if (!(found || length(ind) == kk)) next
    sub <- .fps_linsolve_subs(ind, a[seq_len(kk)])
    if (length(sub)) {
      found <- TRUE
      DE <- xreplace(DE, sub)
    }
    DE <- .as_numer_denom(DE)[[1L]]
    DE <- .fps_first_derivative_factor(sym_factor(DE))
    DE <- collect(DE, .fps_dg(gx, x, 1L))
    sol <- .fps_solve_de(f, x, DE, kk, gname, k)
    if (!is.null(sol)) return(sol)
    if (!length(setdiff(vapply(free_symbols(DE), .k, ""), x$key))) des[[length(des) + 1L]] <- DE
  }
  for (DE in des) {
    sol <- .fps_solve_simple(f, x, DE, gname, k)
    if (!is.null(sol)) return(sol)
  }
  NULL
}
.fps_dg <- function(gx, x, n) if (n == 0L) gx else Derivative(gx, list(x, n))
.fps_first_derivative_factor <- function(e) {
  # as_coeff_mul(Derivative)[1][0]
  fs <- if (e$type == "Mul") .sympy_args(e) else list(e)
  dep <- Filter(function(z) .has_type(z, "Derivative"), fs)
  if (length(dep)) dep[[1L]] else e
}
.fps_solve_de <- function(f, x, DE, order, gname, k) {
  sol <- NULL
  syms <- Filter(function(s) !identical(s$key, x$key), free_symbols(DE))
  RE <- if (length(syms)) .fps_transform_DE_RE(DE, gname, k, order, syms, x) else .fps_hyper_re(DE, gname, k, x)
  if (!length(Filter(function(s) !identical(s$key, k$key), free_symbols(RE)))) sol <- .fps_solve_hyper_RE(f, x, RE, gname, k)
  if (!is.null(sol)) return(sol)
  if (length(syms)) DE <- .fps_transform_explike_DE(DE, gname, x, order, syms)
  if (!length(Filter(function(s) !identical(s$key, x$key), free_symbols(DE)))) sol <- .fps_solve_explike_DE(f, x, DE, gname, k)
  sol
}
.fps_transform_explike_DE <- function(DE, gname, x, order, syms) {
  g <- Function(gname); gx <- g(x)
  eq <- list()
  hc <- coeff(DE, .fps_dg(gx, x, order))
  for (i in seq_len(order) - 1L) {
    c <- coeff(DE, .fps_dg(gx, x, i))
    c <- collect(expand(.div(c, hc)), x)
    eq <- c(eq, .fps_add_args(c))
  }
  temp <- list(); brk <- FALSE
  for (e in eq) {
    if (.has1(e, x)) { brk <- TRUE; break }
    if (length(free_symbols(e))) temp[[length(temp) + 1L]] <- e
  }
  if (!brk) eq <- temp
  if (length(eq)) {
    sub <- .fps_linsolve_subs(eq, syms)
    if (length(sub)) {
      DE <- xreplace(DE, sub)
      DE <- .fps_first_derivative_factor(sym_factor(DE))
      DE <- collect(DE, .fps_dg(gx, x, 1L))
    }
  }
  DE
}
.fps_transform_DE_RE <- function(DE, gname, k, order, syms, x) {
  RE <- .fps_hyper_re(DE, gname, k, x)
  eq <- lapply(seq_len(max(order - 1L, 0L)), function(i) coeff(RE, fcall(gname, list(.add(list(k, .int(i)))))))
  sub <- .fps_linsolve_subs(eq, syms)
  if (length(sub)) {
    RE <- xreplace(RE, sub)
    RE <- .as_numer_denom(sym_factor(RE))[[1L]]
    apps <- .fps_g_apps(RE, gname)
    if (length(apps)) RE <- collect(RE, apps)
    fs <- if (RE$type == "Mul") .sympy_args(RE) else list(RE)
    dep <- Filter(function(z) .has_fun(z, gname), fs)
    if (length(dep)) RE <- dep[[1L]]
    for (i in seq_len(order) - 1L) {
      if (!is_zero_num(coeff(RE, fcall(gname, list(.add(list(k, .int(i))))))) && i) { RE <- subs(RE, k, .sub(k, .int(i))); break }
    }
  }
  RE
}

# Coeff(p, x, n): the n-th coefficient of the polynomial p in x
def_function("Coeff", nargs = 3L, eval = function(args) {
  p <- args[[1L]]; x <- args[[2L]]; n <- args[[3L]]
  if (isTRUE(is_polynomial(p, x)) && isTRUE(ask_fact(n, "integer"))) return(coeff(p, x, as.integer(num_value_dbl(n))))
  NULL
})

# _compute_fps
.fps_compute <- function(f, x, x0, dir, hyper, order, rational, full) {
  if (is_const(x0, "oo") || is_const(x0, "-oo")) {
    dir <- if (is_const(x0, "oo")) S_One else S_NegativeOne
    r <- .fps_compute(subs(f, x, .inv(x)), x, S_Zero, dir, hyper, order, rational, full)
    if (is.null(r)) return(NULL)
    return(list(r[[1L]], subs(r[[2L]], x, .inv(x)), subs(r[[3L]], x, .inv(x))))
  }
  if (!is_zero_num(x0) || identical(dir$key, S_NegativeOne$key)) {
    if (identical(dir$key, S_NegativeOne$key)) { rep <- .add(list(.neg(x), x0)); rep2 <- .neg(x); rep2b <- x0 }
    else { rep <- .add(list(x, x0)); rep2 <- x; rep2b <- .neg(x0) }
    r <- .fps_compute(subs(f, x, rep), x, S_Zero, S_One, hyper, order, rational, full)
    if (is.null(r)) return(NULL)
    back <- .add(list(rep2, rep2b))
    return(list(r[[1L]], subs(r[[2L]], x, back), subs(r[[3L]], x, back)))
  }
  if (isTRUE(is_polynomial(f, x))) {
    k <- Dummy("k")
    ak <- SeqFormula(fcall("Coeff", list(f, x, k)), list(k, 1, oo))
    xk <- SeqFormula(.pow(x, k), list(k, 0, oo))
    return(list(ak, xk, coeff(f, x, 0L)))
  }
  if (f$type == "Add") {
    result <- FALSE
    ak <- SeqFormula(S_Zero, list(Dummy("k"), 0, oo))
    ind <- S_Zero; xk <- NULL
    for (t in .fps_add_args(f)) {
      res <- .fps_compute(t, x, S_Zero, S_One, hyper, order, rational, full)
      if (!is.null(res)) {
        if (!result) { result <- TRUE; xk <- res[[2L]] }
        s1 <- .seq_bounds(ak)[[1L]]; s2 <- .seq_bounds(res[[1L]])[[1L]]
        if (.is_true(Gt(s2, s1))) { seqq <- ak; s <- s1; ff <- s2 } else { seqq <- res[[1L]]; s <- s2; ff <- s1 }
        si <- as.integer(num_value_dbl(s)); fi <- as.integer(num_value_dbl(ff))
        save <- if (fi > si) .add(lapply(seq_len(fi - si) - 1L, function(i) .mul(list(.seq_coeff(seqq, .int(si + i)), .seq_coeff(xk, .int(si + i)))))) else S_Zero
        ak <- SeqAdd(ak, res[[1L]])
        ind <- .add(list(ind, res[[3L]], save))
      } else ind <- .add(list(ind, t))
    }
    if (result) return(list(ak, xk, ind))
    return(NULL)
  }
  syms <- Filter(function(s) !identical(s$key, x$key), free_symbols(f))
  ef <- expand(f)
  si <- if (length(syms)) .fps_as_independent_many(ef, syms) else list(ef, S_One)
  f <- si[[1L]]; symb <- si[[2L]]
  k <- Dummy("k")
  result <- NULL
  if (rational) result <- .fps_rational_algorithm(f, x, k, order, full)
  if (is.null(result) && hyper) result <- .fps_hyper_algorithm(f, x, k, order)
  if (is.null(result)) return(NULL)
  symb <- if (isTRUE(ask_fact(symb, "zero"))) S_One else powsimp(symb)
  ak <- SeqFormula(result[[1L]], list(k, result[[3L]], oo))
  xk <- SeqFormula(powsimp(.mul(list(.pow(x, k), symb))), list(k, 0, oo))
  ind <- powsimp(.mul(list(result[[2L]], symb)))
  list(ak, xk, ind)
}
# Expr.as_independent(*syms): the part free of all syms and the rest
.fps_as_independent_many <- function(e, syms) {
  dep <- function(t) any(vapply(syms, function(s) .has1(t, s), NA))
  if (e$type == "Add") { ts <- e$args; d <- vapply(ts, dep, NA); return(list(.add(ts[!d]), .add(ts[d]))) }
  fs <- if (e$type == "Mul") e$args else list(e)
  d <- vapply(fs, dep, NA)
  list(.mul(fs[!d]), .mul(fs[d]))
}

#' Formal power series
#'
#' `fps(f, x, x0, dir)` is SymPy's `fps()`: the formal power series of `f`
#' about `x0`, as a `FormalPowerSeries` (`f` itself when none is found).
#' `compute_fps()` returns the coefficient sequence, the power sequence and
#' the independent terms. `fps_truncate()` (also `truncate()`) and
#' `fps_polynomial()` give the truncated series; `fps_term()` the term of
#' a given order.
#'
#' @param f An expression.
#' @param x The variable.
#' @param x0 The expansion point.
#' @param dir `1` or `"+"`, `-1` or `"-"`.
#' @param hyper,order,rational,full Options of SymPy's algorithms.
#' @param s A `FormalPowerSeries`.
#' @param n Number of terms.
#' @return A `FormalPowerSeries` or an expression.
#' @examples
#' x <- Symbol("x")
#' fps_truncate(fps(exp(x), x))
#' fps_truncate(fps(log(1 + x), x), 4)
#' @export
fps <- function(f, x = NULL, x0 = 0, dir = 1, hyper = TRUE, order = 4, rational = TRUE, full = FALSE) {
  f <- as_sym(f)
  if (is.null(x)) {
    fr <- free_symbols(f)
    if (length(fr) == 1L) x <- fr[[1L]]
    else if (!length(fr)) return(f)
    else .notimpl_error("multivariate formal power series")
  }
  x <- as_sym(x)
  r <- compute_fps(f, x, x0, dir, hyper, order, rational, full)
  if (is.null(r)) return(f)
  dir <- if (identical(dir, "-") || identical(dir, -1)) S_NegativeOne else S_One
  new_sym("FormalPowerSeries", list(f, x, as_sym(x0), dir, Tuple(r[[1L]], r[[2L]], r[[3L]])))
}

#' @rdname fps
#' @export
compute_fps <- function(f, x, x0 = 0, dir = 1, hyper = TRUE, order = 4, rational = TRUE, full = FALSE) {
  f <- as_sym(f); x <- as_sym(x)
  if (!.has1(f, x)) return(NULL)
  dir <- if (identical(dir, "+") || identical(dir, 1)) S_One else if (identical(dir, "-") || identical(dir, -1)) S_NegativeOne else stop("Dir must be '+' or '-'", call. = FALSE)
  .fps_compute(f, x, as_sym(x0), dir, hyper, as.integer(order), rational, full)
}

.fps_parts <- function(s) list(f = s$args[[1L]], x = s$args[[2L]], x0 = s$args[[3L]], dir = s$args[[4L]],
                               ak = s$args[[5L]]$args[[1L]], xk = s$args[[5L]]$args[[2L]], ind = s$args[[5L]]$args[[3L]])
.fps_get_pow_x <- function(term, x) {
  be <- .as_base_exp(.as_independent(term, x)[[2L]])
  if (!.has1(be[[1L]], x)) return(S_Zero)
  be[[2L]]
}
# FormalPowerSeries._eval_term
.fps_eval_term <- function(s, pt) {
  p <- .fps_parts(s)
  term <- tryCatch(.mul(list(simplify(.seq_coeff(p$ak, pt)), .seq_coeff(p$xk, pt))), error = function(err) S_Zero)
  if (!is_zero_num(p$ind)) {
    ind <- S_Zero
    sym <- free_symbols(s)
    for (t in .fps_add_args(p$ind)) {
      px <- .fps_get_pow_x(t, p$x)
      if (length(sym) && any(vapply(sym, function(v) .has1(px, v), NA))) px <- as_coeff_add_split(px)[[1L]]
      if (is_zero_num(as_sym(pt)) && .is_true(Lt(px, 1))) ind <- .add(list(ind, t))
      else if (.is_true(Ge(px, pt)) && .is_true(Lt(px, .add(list(as_sym(pt), S_One))))) ind <- .add(list(ind, t))
    }
    term <- .add(list(term, ind))
  }
  collect(term, p$x)
}

#' @rdname fps
#' @param pt The index of the term.
#' @export
fps_term <- function(s, pt) .fps_eval_term(s, .int(pt))

#' @rdname fps
#' @export
fps_polynomial <- function(s, n = 6) {
  p <- .fps_parts(s)
  terms <- list()
  sym <- free_symbols(s)
  i <- 0L
  repeat {
    t <- .fps_eval_term(s, .int(i))
    xp <- .fps_get_pow_x(t, p$x)
    if (length(sym) && any(vapply(sym, function(v) .has1(xp, v), NA))) xp <- as_coeff_add_split(xp)[[1L]]
    if (.is_true(Ge(xp, n))) break
    if (isTRUE(ask_fact(xp, "integer")) && i == n + 1L) break
    if (!is_zero_num(t)) terms[[length(terms) + 1L]] <- t
    i <- i + 1L
    if (i > 10L * n + 50L) break
  }
  .add(terms)
}

#' @rdname fps
#' @export
fps_truncate <- function(s, n = 6) {
  if (.ffps_is(s)) return(.ffps_truncate(s, n))
  if (!(is_sym(s) && s$type == "FormalPowerSeries"))
    stop(.sympy_condition("AttributeError", sprintf("'%s' object has no attribute 'truncate'", as_sym(s)$type)))
  p <- .fps_parts(s)
  pt_xk <- .seq_coeff(p$xk, .int(n))
  x0 <- if (is_const(p$x0, "-oo")) S_Infinity else p$x0
  .add(list(fps_polynomial(s, n), .so_new(pt_xk, list(p$x), list(x0))))
}


register_str_printer("FormalPowerSeries", function(e)
  paste0("FormalPowerSeries(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))
