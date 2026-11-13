# The constant coefficient recurrence solver used by fps() and its helpers.
# The ODE solvers are in the files solvers-ode-*.R; rsolve() is in
# solvers-recurr.R.

# Solve linear equations for the constants that they determine; constants
# left free stay symbolic.
.solve_constants <- function(eqs, cs) {
  eqs <- lapply(eqs, expand)
  lin <- all(vapply(eqs, function(e) .is_linear_in(e, cs), TRUE))
  if (!lin) {
    # nonlinear in the constants: general solver, first verified branch
    sols <- tryCatch(solve(eqs, cs, dict = TRUE), error = function(e) list())
    for (sl in sols) {
      if (!is.list(sl) || is.null(names(sl))) next
      res <- vapply(eqs, function(e) is_zero_num(simplify(subs(e, sl))), TRUE)
      if (all(res)) return(lapply(Filter(function(c) !is.null(sl[[c$name]]), cs), function(c) list(c, sl[[c$name]])))
    }
    stop("initial conditions are inconsistent")
  }
  s <- .solve_linear_system(eqs, cs)
  if (is.null(s)) stop("initial conditions are inconsistent")
  lapply(Filter(function(c) !is.null(s[[c$name]]), cs), function(c) list(c, s[[c$name]]))
}


# The constant coefficient recurrence solver used by fps(); rsolve() is in
# solvers-recurr.R.
.rsolve_cc <- function(eq, yn, n, init = NULL) {
  e <- expand(.eq_to_expr(eq))
  yn <- as_sym(yn)
  n <- as_sym(n)
  name <- yn$name
  shifts <- list()
  walk <- function(z) {
    if (z$type == "Function" && z$name == name) {
      d <- expand(.sub(z$args[[1L]], n))
      if (d$type != "Integer") stop("rsolve() needs integer shifts of the index")
      shifts[[z$key]] <<- list(node = z, k = .int_val(d))
    }
    for (a in z$args) walk(a)
  }
  walk(e)
  shifts <- unname(shifts)
  if (!length(shifts)) stop("the unknown does not appear in the recurrence")
  ks <- vapply(shifts, `[[`, 0, "k")
  lo <- min(ks); hi <- max(ks)
  ds <- lapply(shifts, function(s) Dummy("s"))
  er <- xreplace(e, Map(function(s, d) list(s$node, d), shifts, ds))
  coefs <- rep(list(S_Zero), hi - lo + 1)
  for (i in seq_along(shifts)) coefs[[ks[i] - lo + 1]] <- .linear_coeff(er, ds[[i]])
  if (any(vapply(coefs, function(c) .has1(c, n), TRUE))) {
    rest <- expand(.sub(er, .add(Map(function(c, d) .mul(list(c, d)), lapply(seq_along(shifts), function(i) coefs[[ks[i] - lo + 1]]), ds))))
    return(.rsolve_general(coefs, lo, rest, n, init))
  }
  rest <- expand(.sub(er, .add(Map(function(c, d) .mul(list(c, d)), lapply(seq_along(shifts), function(i) coefs[[ks[i] - lo + 1]]), ds))))
  m <- Dummy("m")
  poly <- .add(lapply(seq_along(coefs), function(k) .mul(list(coefs[[k]], .pow(m, .int(k - 1L))))))
  rs <- .roots_legacy(poly, m)
  mult <- attr(rs, "multiplicity")
  if (sum(mult) != length(coefs) - 1L) stop("could not factor the characteristic polynomial")
  terms <- list()
  idx <- 0L
  cs <- list()
  for (i in seq_along(rs)) for (j in 0:(mult[i] - 1L)) {
    c <- Symbol(paste0("C", idx))
    cs[[length(cs) + 1L]] <- c
    terms[[length(terms) + 1L]] <- .mul(list(c, .pow(n, .int(j)), .pow(rs[[i]], n)))
    idx <- idx + 1L
  }
  yh <- .add(terms)
  yp <- S_Zero
  if (!is_zero_num(rest)) yp <- .rsolve_particular(coefs, lo, .neg(rest), n, rs, mult)
  sol <- .add(list(yh, yp))
  if (!is.null(init)) {
    eqs <- lapply(init, function(p) .sub(subs(sol, n, as_sym(p[[1L]])), as_sym(p[[2L]])))
    # inconsistent initial values: no solution, as SymPy's None
    k <- tryCatch(.solve_constants(eqs, cs), error = function(e) NULL)
    if (is.null(k)) return(NULL)
    sol <- simplify_if_cheap(xreplace(sol, k))
  }
  sol
}

# SymPy's rsolve() driver for polynomial or rational coefficients: clear
# denominators, shift the lowest term to y(n) and call rsolve_hyper().
.rsolve_general <- function(coefs, lo, rest, n, init) {
  coefs <- lapply(coefs, simplify)
  common <- S_One
  for (c in coefs) if (!is_polynomial(c, n)) common <- lcm(common, fraction(together(c))[[2L]], n)
  ind <- fraction(together(rest))
  if (is_polynomial(ind[[2L]], n)) common <- lcm(common, ind[[2L]], n)
  if (!identical(common$key, S_One$key)) {
    coefs <- lapply(coefs, function(c) {
      nd <- fraction(together(c))
      .mul(list(nd[[1L]], quo(common, nd[[2L]], n)))
    })
    rest <- .mul(list(ind[[1L]], quo(common, ind[[2L]], n)))
  }
  sh <- .add(list(n, .int(-lo)))
  coefs <- lapply(coefs, function(c) expand(subs(c, n, sh)))
  rest <- expand(subs(rest, n, sh))
  r <- .rsolve_hyper(coefs, .neg(rest), n)
  if (is.null(r)) return(NULL)
  sol <- r$result
  if (length(r$C) && length(init)) {
    eqs <- lapply(init, function(p) .sub(subs(sol, n, as_sym(p[[1L]])), as_sym(p[[2L]])))
    cs <- tryCatch(.solve_constants(eqs, r$C), error = function(e) NULL)
    if (is.null(cs)) return(NULL)
    sol <- simplify_if_cheap(xreplace(sol, cs))
  }
  sol
}

.rsolve_particular <- function(coefs, lo, r, n, rs, mult) {
  terms <- if (r$type == "Add") r$args else list(r)
  groups <- list()
  for (t in terms) {
    base <- S_One; deg <- 0L
    facs <- if (t$type == "Mul") t$args else list(t)
    for (f in facs) {
      if (!.has1(f, n)) next
      d <- .degree_in(f, n)
      if (!is.na(d)) {
        deg <- deg + d
      } else if (f$type == "Pow" && !.has1(f$args[[1L]], n) && identical(f$args[[2L]]$key, n$key)) {
        base <- f$args[[1L]]
      } else stop("rsolve() cannot handle this inhomogeneous term")
    }
    k <- base$key
    groups[[k]] <- list(base = base, deg = max(deg, if (is.null(groups[[k]])) 0L else groups[[k]]$deg))
  }
  unknowns <- list()
  trial <- list()
  for (g in groups) {
    s <- 0L
    for (i in seq_along(rs)) if (identical(rs[[i]]$key, g$base$key)) s <- mult[i]
    for (d in 0:g$deg) {
      u <- Dummy(paste0("A", length(unknowns)))
      unknowns[[length(unknowns) + 1L]] <- u
      trial[[length(trial) + 1L]] <- .mul(list(u, .pow(n, .int(d + s)), .pow(g$base, n)))
    }
  }
  Y <- .add(trial)
  L <- .add(lapply(seq_along(coefs), function(k) .mul(list(coefs[[k]], subs(Y, n, .add(list(n, .int(lo + k - 1L))))))))
  res <- expand(.sub(L, r))
  res <- expand_power_exp(res)
  gens <- .find_gens(list(res))
  pw <- Filter(function(g) g$type == "Pow" && identical(g$args[[2L]]$key, n$key), gens)
  bases <- unique(c(list(n), pw))
  dm <- lapply(pw, function(g) Dummy("g"))
  rs2 <- xreplace(res, Map(list, pw, dm))
  eqs <- list()
  for (t in if (rs2$type == "Add") rs2$args else list(rs2)) {
    mono <- S_One
    for (f in if (t$type == "Mul") t$args else list(t)) {
      b <- .base_exp(f)[[1L]]
      if (identical(b$key, n$key) || any(vapply(dm, function(d) identical(d$key, b$key), TRUE))) mono <- .mul(list(mono, f))
    }
    eqs[[mono$key]] <- c(eqs[[mono$key]], list(.div(t, mono)))
  }
  sol <- .linsolve_exprs(unname(lapply(eqs, .add)), unknowns)
  if (is.null(sol)) stop("could not find a particular solution")
  xreplace(Y, sol)
}
