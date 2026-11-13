# The nth order linear ODE helpers, after sympy/solvers/ode/nonhomogeneous.py:
# the characteristic equations of constant coefficient and Euler equations,
# variation of parameters and undetermined coefficients.

# nonhomogeneous._test_term
.ode_test_term <- function(coeff, func, order) {
  x <- func$args[[1L]]
  if (order < 0L) .ode_stop("ValueError", "order should be greater than 0")
  if (is_zero_num(coeff)) return(TRUE)
  if (order == 0L) return(!.sp_in(x, free_symbols(coeff)))
  if (coeff$type == "Mul") {
    if (.has1(coeff, func)) return(FALSE)
    xo <- .pow(x, .int(order))
    return(any(vapply(coeff$args, function(a) identical(a$key, xo$key), TRUE)))
  }
  if (coeff$type == "Pow") {
    be <- .base_exp(coeff)
    return(identical(be[[1L]]$key, x$key) && identical(be[[2L]]$key, .int(order)$key))
  }
  if (order == 1L) return(identical(x$key, coeff$key))
  FALSE
}

# group equal roots, keeping their first order: a list of list(root, count)
.ode_counter <- function(xs) {
  keys <- vapply(xs, .k, "")
  u <- unique(keys)
  lapply(u, function(k) list(xs[[match(k, keys)]], sum(keys == k)))
}

# rootof(chareq, k) for k in range(degree): explicit roots ordered as
# CRootOf orders them (real roots increasing, then complex roots by real and
# imaginary part), or CRootOf objects when roots() cannot find them all
.ode_rootofs <- function(expr, symbol, deg) {
  rs <- tryCatch(roots(expr, symbol, multiple = TRUE, cubics = FALSE, quartics = FALSE), error = function(e) list())
  if (length(rs) != deg) {
    # rootof() of a reducible polynomial: the roots of each irreducible
    # factor, a rational root for a linear one and CRootOf objects otherwise
    rs <- list()
    for (fm in factor_list(expr, symbol)$factors) {
      g <- fm[[1L]]
      if (!.has1(g, symbol)) next
      gd <- degree(Poly(g, symbol))
      gr <- tryCatch(roots(g, symbol, multiple = TRUE, cubics = FALSE, quartics = FALSE), error = function(e) list())
      if (length(gr) != gd) gr <- lapply(seq_len(gd) - 1L, function(k) CRootOf(g, symbol, k))
      for (j in seq_len(as.integer(as.character(fm[[2L]])))) rs <- c(rs, gr)
    }
  }
  vals <- vapply(rs, function(r) as.complex(.neval(r, list())), 0i)
  isreal <- vapply(rs, function(r) isTRUE(ask_fact(r, "real")), TRUE) | abs(Im(vals)) < 1e-12
  o <- order(!isreal, Re(vals), ifelse(isreal, 0, Im(vals)))
  rs[o]
}

# nonhomogeneous._get_euler_characteristic_eq_sols
.ode_get_euler_characteristic_eq_sols <- function(eq, func, match_obj) {
  x <- func$args[[1L]]
  chareq <- S_Zero
  symbol <- Dummy("x")
  for (i in .rk_ge0(match_obj)) {
    ii <- as.integer(i)
    chareq <- .add(list(chareq, expand(.mul(list(match_obj[[i]], .ode_d(.pow(x, symbol), x, ii), .pow(x, .neg(symbol)))))))
  }
  deg <- degree(Poly(chareq, symbol))
  chareqroots <- .ode_rootofs(chareq, symbol, deg)
  collectterms <- list()
  constants <- rev(.ode_numbered_constants(list(eq), deg * 2L))
  pop <- function() {
    v <- constants[[length(constants)]]
    constants <<- constants[-length(constants)]
    v
  }
  charroots <- .ode_counter(chareqroots)
  gsol <- S_Zero
  lnx <- fcall("log", list(x))
  for (rm in charroots) {
    root <- rm[[1L]]
    multiplicity <- rm[[2L]]
    for (i in seq_len(multiplicity) - 1L) {
      if (root$type == "CRootOf") {
        gsol <- .add(list(gsol, .mul(list(.pow(x, root), pop()))))
        if (multiplicity != 1L) .ode_stop("ValueError", "Value should be 1")
        collectterms <- c(list(list(S_Zero, root, S_Zero)), collectterms)
      } else if (isTRUE(ask_fact(root, "real"))) {
        gsol <- .add(list(gsol, .mul(list(.pow(lnx, .int(i)), .pow(x, root), pop()))))
        collectterms <- c(list(list(.int(i), root, S_Zero)), collectterms)
      } else {
        reroot <- fcall("re", list(root))
        imroot <- fcall("im", list(root))
        c1 <- pop()
        c2 <- pop()
        gsol <- .add(list(gsol, .mul(list(.pow(lnx, .int(i)), .pow(x, reroot),
          .add(list(.mul(list(c1, fcall("sin", list(.mul(list(Abs(imroot), lnx)))))), .mul(list(c2, fcall("cos", list(.mul(list(imroot, lnx))))))))))))
        collectterms <- c(list(list(.int(i), reroot, imroot)), collectterms)
      }
    }
  }
  gsol <- Eq(func, gsol)
  gensols <- list()
  for (ct in collectterms) {
    i <- ct[[1L]]; reroot <- ct[[2L]]; imroot <- ct[[3L]]
    if (is_zero_num(imroot)) gensols[[length(gensols) + 1L]] <- .mul(list(.pow(lnx, i), .pow(x, reroot)))
    else {
      sin_form <- .mul(list(.pow(lnx, i), .pow(x, reroot), fcall("sin", list(.mul(list(Abs(imroot), lnx))))))
      if (.sp_in(sin_form, gensols)) gensols[[length(gensols) + 1L]] <- .mul(list(.pow(lnx, i), .pow(x, reroot), fcall("cos", list(.mul(list(imroot, lnx))))))
      else gensols[[length(gensols) + 1L]] <- sin_form
    }
  }
  list(gsol, gensols)
}

# nonhomogeneous._solve_variation_of_parameters
.ode_solve_variation_of_parameters <- function(eq, func, roots, homogen_sol, order, match_obj, simplify_flag = TRUE) {
  x <- func$args[[1L]]
  r <- match_obj
  psol <- S_Zero
  wr <- wronskian(roots, x)
  if (simplify_flag) {
    wr <- simplify(wr)
    wr <- trigsimp(wr, deep = TRUE, recursive = TRUE)
  }
  if (!.ode_true(wr))
    .ode_nie(paste0("Cannot find ", order, " solutions to the homogeneous equation necessary to apply variation of parameters to ",
                    .sp_py_str(eq), " (Wronskian == 0)"))
  if (length(roots) != order)
    .ode_nie(paste0("Cannot find ", order, " solutions to the homogeneous equation necessary to apply variation of parameters to ",
                    .sp_py_str(eq), " (number of terms != order)"))
  negoneterm <- .pow(S_NegativeOne, .int(order))
  for (i in roots) {
    others <- Filter(function(s) !identical(s$key, i$key), roots)
    integrand <- .div(.mul(list(wronskian(others, x), .rk(r, -1L))), wr)
    psol <- .add(list(psol, .div(.mul(list(negoneterm, Integral(integrand, x), i)), .rk(r, order))))
    negoneterm <- .neg(negoneterm)
  }
  if (simplify_flag) {
    psol <- simplify(psol)
    psol <- trigsimp(psol, deep = TRUE)
  }
  Eq(func, .add(list(homogen_sol$args[[2L]], psol)))
}

# nonhomogeneous._get_const_characteristic_eq_sols
.ode_get_const_characteristic_eq_sols <- function(r, func, order) {
  x <- func$args[[1L]]
  chareq <- S_Zero
  symbol <- Dummy("x")
  for (i in .rk_ge0(r)) chareq <- .add(list(chareq, .mul(list(r[[i]], .pow(symbol, .int(as.integer(i)))))))
  P <- Poly(chareq, symbol)
  chareqroots <- roots(chareq, symbol, multiple = TRUE)
  if (length(chareqroots) != order) chareqroots <- .ode_rootofs(chareq, symbol, degree(P))
  chareq_is_complex <- !all(vapply(all_coeffs(P), function(c) isTRUE(ask_fact(c, "real")), TRUE))
  counts <- .ode_counter(chareqroots)
  ckeys <- vapply(counts, function(p) p[[1L]]$key, "")
  used <- character()
  collectterms <- list()
  gensols <- list()
  conjugate_roots <- list()
  for (root in chareqroots) {
    if (root$key %in% used) next
    used <- c(used, root$key)
    multiplicity <- counts[[match(root$key, ckeys)]][[2L]]
    for (i in seq_len(multiplicity) - 1L) {
      xi <- .pow(x, .int(i))
      if (chareq_is_complex) {
        gensols[[length(gensols) + 1L]] <- .mul(list(xi, fcall("exp", list(.mul(list(root, x))))))
        collectterms <- c(list(list(.int(i), root, S_Zero)), collectterms)
        next
      }
      reroot <- fcall("re", list(root))
      imroot <- fcall("im", list(root))
      if (.has_function(imroot, "atan2") && .has_function(reroot, "atan2")) {
        gensols[[length(gensols) + 1L]] <- .mul(list(xi, fcall("exp", list(.mul(list(root, x))))))
        collectterms <- c(list(list(.int(i), root, S_Zero)), collectterms)
      } else {
        if (.sp_in(root, conjugate_roots)) {
          collectterms <- c(list(list(.int(i), reroot, imroot)), collectterms)
          next
        }
        if (is_zero_num(imroot)) {
          gensols[[length(gensols) + 1L]] <- .mul(list(xi, fcall("exp", list(.mul(list(reroot, x))))))
          collectterms <- c(list(list(.int(i), reroot, S_Zero)), collectterms)
          next
        }
        conjugate_roots[[length(conjugate_roots) + 1L]] <- .ode_conjugate(root)
        ex <- fcall("exp", list(.mul(list(reroot, x))))
        gensols[[length(gensols) + 1L]] <- .mul(list(xi, ex, fcall("sin", list(.mul(list(.ode_abs_im(imroot), x))))))
        gensols[[length(gensols) + 1L]] <- .mul(list(xi, ex, fcall("cos", list(.mul(list(imroot, x))))))
        collectterms <- c(list(list(.int(i), reroot, imroot)), collectterms)
      }
    }
  }
  list(gensols, collectterms)
}

# conjugate(root) and Abs(im(root)): a CRootOf evaluates them from its
# numerical value, the conjugate being the root of the same polynomial
.ode_crootof_values <- function(root) vapply(.numeric_roots(root$args[[1L]], Symbol("x")), as.complex, 0i)

.ode_conjugate <- function(root) {
  if (root$type != "CRootOf") return(fcall("conjugate", list(root)))
  v <- .ode_crootof_values(root)
  k <- as.integer(as.character(root$args[[2L]]$value)) + 1L
  new_sym("CRootOf", list(root$args[[1L]], .int(which.min(Mod(v - Conj(v[k]))) - 1L)))
}

.ode_abs_im <- function(imroot) {
  if (imroot$type == "Function" && imroot$name == "im" && imroot$args[[1L]]$type == "CRootOf") {
    v <- Im(as.complex(.neval(imroot$args[[1L]], list())))
    if (v > 0) return(imroot)
    if (v < 0) return(.neg(imroot))
  }
  Abs(imroot)
}

# nonhomogeneous._get_simplified_sol
.ode_get_simplified_sol <- function(sol, func, collectterms) {
  x <- func$args[[1L]]
  collectterms <- rev(.key_sort(lapply(collectterms, function(t) lapply(t, .sort_key)), collectterms))
  sol <- sol[[1L]]$args[[2L]]
  sol <- expand_mul(sol)
  for (ct in collectterms) {
    i <- ct[[1L]]; reroot <- ct[[2L]]; imroot <- ct[[3L]]
    base <- .mul(list(.pow(x, i), fcall("exp", list(.mul(list(reroot, x))))))
    sol <- collect(sol, .mul(list(base, fcall("sin", list(.mul(list(.ode_abs_im(imroot), x)))))))
    sol <- collect(sol, .mul(list(base, fcall("cos", list(.mul(list(imroot, x)))))))
  }
  for (ct in collectterms) sol <- collect(sol, .mul(list(.pow(x, ct[[1L]]), fcall("exp", list(.mul(list(ct[[2L]], x)))))))
  sol <- powsimp(sol)
  Eq(func, sol)
}

# nonhomogeneous._undetermined_coefficients_match. The trial set function
# keeps SymPy's mutable default set: terms it adds to it stay there for the
# rest of the call of _undetermined_coefficients_match, which defines it.
.ode_state$trial_default <- list()

.ode_ucm_test_term <- function(expr, x, a, b) {
  if (!.has1(expr, x)) return(TRUE)
  if (expr$type == "Add") return(all(vapply(.sympy_args(expr), .ode_ucm_test_term, TRUE, x = x, a = a, b = b)))
  if (expr$type == "Mul") {
    if (.has_function(expr, c("sin", "cos"))) {
      foundtrig <- FALSE
      for (i in .sympy_args(expr)) {
        if (.has_function(i, c("sin", "cos"))) {
          if (foundtrig) return(FALSE)
          foundtrig <- TRUE
        }
      }
    }
    return(all(vapply(.sympy_args(expr), .ode_ucm_test_term, TRUE, x = x, a = a, b = b)))
  }
  if (expr$type == "Function") return(expr$name %in% c("sin", "cos", "exp", "sinh", "cosh") && !is.null(sym_match(expr$args[[1L]], a * x + b)))
  if (expr$type == "Pow" && expr$args[[1L]]$type %in% c("Symbol", "Dummy") && expr$args[[2L]]$type == "Integer" && expr$args[[2L]]$value >= 0)
    return(TRUE)
  if (expr$type == "Pow" && .expr_is_number(expr$args[[1L]])) return(!is.null(sym_match(expr$args[[2L]], a * x + b)))
  expr$type %in% c("Symbol", "Dummy") || isTRUE(.expr_is_number(expr))
}

.ode_remove_coefficient <- function(expr, x) {
  term <- S_One
  if (expr$type == "Mul") {
    for (i in .sympy_args(expr)) if (.has1(i, x)) term <- .mul(list(term, i))
  } else if (.has1(expr, x)) term <- expr
  term
}

.ode_set_add <- function(s, e) {
  s[[e$key]] <- e
  s
}
.ode_set_union <- function(a, b) {
  for (k in names(b)) a[[k]] <- b[[k]]
  a
}
.ode_set_eq <- function(a, b) length(a) == length(b) && all(names(a) %in% names(b))

.ode_get_trial_set <- function(expr, x, exprs = NULL) {
  is_default <- is.null(exprs)
  if (is_default) exprs <- .ode_state$trial_default
  expr <- expand_mul(expr)
  if (expr$type == "Add") {
    for (term in .sympy_args(expr)) {
      rc <- .ode_remove_coefficient(term, x)
      if (rc$key %in% names(exprs)) next
      exprs <- .ode_set_add(exprs, rc)
      if (is_default) {
        .ode_state$trial_default <- exprs
        is_default <- FALSE
      }
      exprs <- .ode_set_union(exprs, .ode_get_trial_set(term, x, exprs))
    }
  } else {
    term <- .ode_remove_coefficient(expr, x)
    tmpset <- .ode_set_add(exprs, term)
    oldset <- NULL
    while (is.null(oldset) || !.ode_set_eq(tmpset, oldset)) {
      oldset <- tmpset
      expr <- .diff1(expr, x)
      term <- .ode_remove_coefficient(expr, x)
      if (term$type == "Add") tmpset <- .ode_set_union(tmpset, .ode_get_trial_set(term, x, tmpset))
      else tmpset <- .ode_set_add(tmpset, term)
    }
    exprs <- tmpset
  }
  exprs
}

.ode_undetermined_coefficients_match <- function(expr, x, func = NULL, eq_homogeneous = S_Zero) {
  .ode_state$trial_default <- list()
  a <- Wild("a", exclude = list(x))
  b <- Wild("b", exclude = list(x))
  expr <- powsimp(expr, combine = "exp")
  test <- .ode_ucm_test_term(expr, x, a, b)
  if (!test) return(.pydict(list("test"), list(FALSE)))
  temp_set <- list()
  for (i in if (expr$type == "Add") .sympy_args(expr) else list(expr)) {
    act <- .ode_get_trial_set(i, x)
    if (!is_zero_num(eq_homogeneous) || !identical(eq_homogeneous$key, S_Zero$key)) {
      is_hom <- function(term) isTRUE(.sp_fact(expand(.ode_sub_func_doit(eq_homogeneous, func, term)), "zero"))
      while (any(vapply(act, is_hom, TRUE))) {
        act2 <- list()
        for (ts in act) act2 <- .ode_set_add(act2, .mul(list(x, ts)))
        act <- act2
      }
    }
    temp_set <- .ode_set_union(temp_set, act)
  }
  .pydict(list("test", "trialset"), list(TRUE, structure(unname(temp_set), pyset = TRUE)))
}

# nonhomogeneous._solve_undetermined_coefficients
.ode_solve_undetermined_coefficients <- function(eq, func, order, match, trialset) {
  r <- match
  gensols <- r$list
  gsol <- r$sol
  x <- func$args[[1L]]
  if (length(gensols) != order)
    .ode_nie(paste0("Cannot find ", order, " solutions to the homogeneous equation necessary to apply undetermined coefficients to ",
                    .sp_py_str(eq), " (number of terms != order)"))
  trialfunc <- S_Zero
  coefflist <- list()
  k <- 0L
  for (i in trialset) {
    c <- Dummy(paste0("a", k))
    k <- k + 1L
    coefflist[[length(coefflist) + 1L]] <- c
    trialfunc <- .add(list(trialfunc, .mul(list(c, i))))
  }
  eqs <- .ode_sub_func_doit(eq, func, trialfunc)
  ck <- vapply(trialset, .k, "")
  cv <- rep(list(S_Zero), length(trialset))
  ckeys <- as.list(trialset)
  eqs <- .mexpand(eqs)
  for (i in if (eqs$type == "Add") .sympy_args(eqs) else list(eqs)) {
    s <- .ode_separatevars(i, dict = TRUE, symbols = list(x))
    sx <- .sv(s, x)
    j <- match(sx$key, ck)
    if (is.na(j)) {
      ck <- c(ck, sx$key)
      ckeys[[length(ckeys) + 1L]] <- sx
      cv[[length(cv) + 1L]] <- s$coeff
    } else cv[[j]] <- .add(list(cv[[j]], s$coeff))
  }
  coeffvals <- .sp_solve(cv, list(coefflist), list(dict = TRUE, `_raw` = TRUE))
  if (!length(coeffvals))
    .ode_nie(paste0("Could not solve `", .sp_py_str(eq), "` using the method of undetermined coefficients (unable to solve for coefficients)."))
  cvd <- coeffvals[[1L]]
  psol <- subs(trialfunc, lapply(seq_along(cvd$k), function(i) list(cvd$k[[i]], cvd$v[[i]])))
  Eq(func, .add(list(gsol$args[[2L]], psol)))
}
