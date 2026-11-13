# Lie group symmetries of first order ODEs, after
# sympy/solvers/ode/lie_group.py: the heuristics that find infinitesimals
# xi(x, f(x)) and eta(x, f(x)), and the lie_group solver built on them.

.ODE_LIE_HEURISTICS <- c("abaco1_simple", "abaco1_product", "abaco2_similar", "abaco2_unique_unknown",
                         "abaco2_unique_general", "linear", "function_sum", "bivariate", "chi")

# an infinitesimal dictionary {xi: a, eta: b} in the given key order
.ode_inf <- function(k1, v1, k2, v2) .pydict(list(k1, k2), list(as_sym(v1), as_sym(v2)))

.ode_inf_same <- function(a, b) {
  length(a$keys) == length(b$keys) &&
    all(vapply(seq_along(a$keys), function(i) {
      j <- match(a$keys[[i]]$key, vapply(b$keys, .k, ""))
      !is.na(j) && identical(a$values[[i]]$key, b$values[[j]]$key)
    }, TRUE))
}
.ode_inf_in <- function(inf, lst) any(vapply(lst, .ode_inf_same, TRUE, b = inf))

.ode_xi_eta <- function(x, func) list(fcall("xi", list(x, func)), fcall("eta", list(x, func)))

# separatevars(e, dict=True, symbols=[x, y]) with its truth value: NULL when
# the expression does not separate or the coefficient is zero
.ode_sep_ok <- function(e, syms) {
  s <- .ode_separatevars(e, dict = TRUE, symbols = syms)
  if (is.null(s) || !.ode_true(s$coeff)) return(NULL)
  s
}

.ode_lie_abaco1_simple <- function(match, comp = FALSE) {
  xieta <- list()
  y <- match$y; h <- match$h; func <- match$func; x <- func$args[[1L]]
  hx <- match$hx; hy <- match$hy
  xe <- .ode_xi_eta(x, func)
  add <- function(inf) {
    if (!comp) return(list(inf))
    if (!.ode_inf_in(inf, xieta)) xieta[[length(xieta) + 1L]] <<- inf
    NULL
  }
  if (!.sp_in(y, free_symbols(hy))) {
    fx <- fcall("exp", list(integrate(hy, x)))
    r <- add(.ode_inf(xe[[1L]], S_Zero, xe[[2L]], fx))
    if (!is.null(r)) return(r)
  }
  factor <- .div(hy, h)
  if (!.sp_in(x, free_symbols(factor))) {
    fy <- fcall("exp", list(integrate(factor, y)))
    r <- add(.ode_inf(xe[[1L]], S_Zero, xe[[2L]], subs(fy, y, func)))
    if (!is.null(r)) return(r)
  }
  factor <- .div(.neg(hx), h)
  if (!.sp_in(y, free_symbols(factor))) {
    fx <- fcall("exp", list(integrate(factor, x)))
    r <- add(.ode_inf(xe[[1L]], fx, xe[[2L]], S_Zero))
    if (!is.null(r)) return(r)
  }
  factor <- .div(.neg(hx), .pow(h, S_Two))
  if (!.sp_in(x, free_symbols(factor))) {
    fy <- fcall("exp", list(integrate(factor, y)))
    r <- add(.ode_inf(xe[[1L]], subs(fy, y, func), xe[[2L]], S_Zero))
    if (!is.null(r)) return(r)
  }
  if (length(xieta)) xieta else NULL
}

.ode_lie_abaco1_product <- function(match, comp = FALSE) {
  xieta <- list()
  y <- match$y; h <- match$h; hinv <- match$hinv; func <- match$func; x <- func$args[[1L]]
  xe <- .ode_xi_eta(x, func)
  inf <- .ode_sep_ok(.div(.diff1(.diff1(fcall("log", list(h)), y), x), .pow(h, S_Two)), list(x, y))
  if (!is.null(inf)) {
    fx <- .sv(inf, x)
    gy <- simplify(.mul(list(fx, .diff1(.inv(.mul(list(fx, h))), x))))
    if (!.sp_in(x, free_symbols(gy))) {
      gy <- fcall("exp", list(integrate(gy, y)))
      inf <- .ode_inf(xe[[2L]], S_Zero, xe[[1L]], subs(.mul(list(fx, gy)), y, func))
      if (!comp) return(list(inf))
      if (!.ode_inf_in(inf, xieta)) xieta[[length(xieta) + 1L]] <- inf
    }
  }
  u1 <- Dummy("u1")
  inf <- .ode_sep_ok(.div(.diff1(.diff1(fcall("log", list(hinv)), y), x), .pow(hinv, S_Two)), list(x, y))
  if (!is.null(inf)) {
    fx <- .sv(inf, x)
    gy <- simplify(.mul(list(fx, .diff1(.inv(.mul(list(fx, hinv))), x))))
    if (!.sp_in(x, free_symbols(gy))) {
      gy <- fcall("exp", list(integrate(gy, y)))
      etaval <- .mul(list(fx, gy))
      etaval <- subs(subs(etaval, list(list(x, u1), list(y, x))), u1, y)
      inf <- .ode_inf(xe[[2L]], subs(etaval, y, func), xe[[1L]], S_Zero)
      if (!comp) return(list(inf))
      if (!.ode_inf_in(inf, xieta)) xieta[[length(xieta) + 1L]] <- inf
    }
  }
  if (length(xieta)) xieta else NULL
}

# the coefficients of a polynomial in x and y, as Poly(e, x, y).as_dict()
.ode_xy_coeffs <- function(e, x, y) {
  p <- Poly(e, x, y)
  cs <- coeffs(p)
  lapply(cs, as_sym)
}

# solve(eqs, *syms) for the dictionary of a linear system, as a list of
# (symbol, value) pairs; NULL when there is no solution
.ode_solve_dict <- function(eqs, syms) {
  sl <- tryCatch(.sp_solve(eqs, syms, list(dict = TRUE, `_raw` = TRUE)), error = function(e) list())
  if (!length(sl)) return(NULL)
  d <- sl[[1L]]
  lapply(seq_along(d$k), function(i) list(d$k[[i]], d$v[[i]]))
}

.ode_lie_bivariate <- function(match, comp = FALSE) {
  h <- match$h; hx <- match$hx; hy <- match$hy; func <- match$func; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  if (!.is_rational_function_in(h, free_symbols(h))) return(NULL)
  sy <- lapply(c("etax", "etay", "etad", "xix", "xiy", "xid"), Symbol)
  ipde <- .add(list(sy[[1L]], .mul(list(.sub(sy[[2L]], sy[[4L]]), h)), .neg(.mul(list(sy[[5L]], .pow(h, S_Two)))),
                    .neg(.mul(list(sy[[6L]], hx))), .neg(.mul(list(sy[[3L]], hy)))))
  num <- .as_numer_denom(cancel(ipde))[[1L]]
  deg <- total_degree(Poly(num, x, y))
  deta <- fcall("deta", list(x, y))
  dxi <- fcall("dxi", list(x, y))
  ipde <- .add(list(.diff1(deta, x), .mul(list(.sub(.diff1(deta, y), .diff1(dxi, x)), h)),
                    .neg(.mul(list(.diff1(dxi, y), .pow(h, S_Two)))), .neg(.mul(list(dxi, hx))), .neg(.mul(list(deta, hy)))))
  xieq <- Symbol("xi0")
  etaeq <- Symbol("eta0")
  for (i in 0:deg) {
    if (i) {
      xieq <- .add(list(xieq, .add(lapply(0:i, function(p) .mul(list(Symbol(paste0("xi_", p, "_", i - p)), .pow(x, .int(p)), .pow(y, .int(i - p))))))))
      etaeq <- .add(list(etaeq, .add(lapply(0:i, function(p) .mul(list(Symbol(paste0("eta_", p, "_", i - p)), .pow(x, .int(p)), .pow(y, .int(i - p))))))))
    }
    pden <- expand(.as_numer_denom(doit(subs(ipde, list(list(dxi, xieq), list(deta, etaeq)), simultaneous = TRUE)))[[1L]])
    if (isTRUE(is_polynomial(pden, x, y)) && pden$type == "Add") {
      polyy <- .ode_xy_coeffs(pden, x, y)
      if (length(polyy)) {
        xk <- c(x$key, y$key)
        symset <- Filter(function(s) !(s$key %in% xk), .sp_sorted(c(free_symbols(xieq), free_symbols(etaeq))))
        soldict <- .ode_solve_dict(polyy, symset)
        if (is.null(soldict)) .ode_stop("IndexError", "list index out of range")
        if (any(vapply(soldict, function(p) .ode_true(p[[2L]]), TRUE))) {
          xired <- subs(xieq, soldict)
          etared <- subs(etaeq, soldict)
          ones <- lapply(symset, function(s) list(s, S_One))
          return(list(.ode_inf(xe[[2L]], subs(subs(etared, ones), y, func), xe[[1L]], subs(subs(xired, ones), y, func))))
        }
      }
    }
  }
  NULL
}

.ode_lie_chi <- function(match, comp = FALSE) {
  h <- match$h; hy <- match$hy; func <- match$func; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  if (!.is_rational_function_in(h, free_symbols(h))) return(NULL)
  schi <- Symbol("schi"); schix <- Symbol("schix"); schiy <- Symbol("schiy")
  cpde <- .add(list(schix, .mul(list(h, schiy)), .neg(.mul(list(hy, schi)))))
  num <- .as_numer_denom(cancel(cpde))[[1L]]
  deg <- total_degree(Poly(num, x, y))
  chi <- fcall("chi", list(x, y))
  cpde <- .add(list(.diff1(chi, x), .mul(list(h, .diff1(chi, y))), .neg(.mul(list(hy, chi)))))
  chieq <- Symbol("chi")
  if (deg >= 1) for (i in 1:deg) {
    chieq <- .add(list(chieq, .add(lapply(0:i, function(p) .mul(list(Symbol(paste0("chi_", p, "_", i - p)), .pow(x, .int(p)), .pow(y, .int(i - p))))))))
    cnum <- expand(.as_numer_denom(cancel(doit(subs(cpde, chi, chieq))))[[1L]])
    if (isTRUE(is_polynomial(cnum, x, y)) && cnum$type == "Add") {
      cpoly <- .ode_xy_coeffs(cnum, x, y)
      if (length(cpoly)) {
        xk <- c(x$key, y$key)
        solsyms <- Filter(function(s) !(s$key %in% xk), free_symbols(chieq))
        soldict <- .ode_solve_dict(cpoly, solsyms)
        if (is.null(soldict)) .ode_stop("IndexError", "list index out of range")
        if (any(vapply(soldict, function(p) .ode_true(p[[2L]]), TRUE))) {
          chieq <- subs(chieq, soldict)
          chieq <- subs(chieq, lapply(solsyms, function(s) list(s, S_One)))
          qr <- div(chieq, h)
          return(list(.ode_inf(xe[[2L]], subs(qr$r, y, func), xe[[1L]], .neg(subs(qr$q, y, func)))))
        }
      }
    }
  }
  NULL
}

.ode_lie_function_sum <- function(match, comp = FALSE) {
  xieta <- list()
  h <- match$h; func <- match$func; hinv <- match$hinv; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  for (odefac in list(h, hinv)) {
    factor <- .mul(list(odefac, .diff_spec_eval(.inv(odefac), list(list(x, 2L)))))
    sep <- .ode_sep_ok(.diff1(.inv(factor), y), list(x, y))
    if (!is.null(sep) && .has1(.sv(sep, x), x) && .has1(.sv(sep, y), y)) {
      k <- Dummy("k")
      gy <- .mul(list(k, integrate(.sv(sep, y), y)))
      fdd <- .inv(.mul(list(k, .sv(sep, x), sep$coeff)))
      fx <- simplify(.sub(.div(fdd, factor), gy))
      check <- simplify(.sub(.diff_spec_eval(fx, list(list(x, 2L))), fdd))
      if (.ode_true(fx)) {
        if (!.ode_true(check)) {
          fx <- subs(fx, k, S_One)
          gy <- .div(gy, k)
        } else {
          sol <- .sp_solve(check, list(k), list())
          if (length(sol)) {
            sol <- sol[[1L]]
            fx <- subs(fx, k, sol)
            gy <- .mul(list(.div(gy, k), sol))
          } else next
        }
        isinv <- identical(odefac$key, hinv$key)
        if (isinv) {
          fx <- subs(fx, x, y)
          gy <- subs(gy, y, x)
        }
        etaval <- factor_terms(.add(list(fx, gy)))
        if (etaval$type == "Mul") etaval <- .mul(Filter(function(a) .has1(a, x) || .has1(a, y), etaval$args))
        inf <- if (isinv) .ode_inf(xe[[2L]], subs(etaval, y, func), xe[[1L]], S_Zero)
               else .ode_inf(xe[[1L]], subs(etaval, y, func), xe[[2L]], S_Zero)
        if (!comp) return(list(inf))
        xieta[[length(xieta) + 1L]] <- inf
      }
    }
    if (length(xieta)) return(xieta)
  }
  NULL
}

.ode_lie_abaco2_similar <- function(match, comp = FALSE) {
  h <- match$h; hx <- match$hx; hy <- match$hy; func <- match$func; hinv <- match$hinv
  x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  for (inverse in c(FALSE, TRUE)) {
    g <- if (inverse) hinv else h
    gx <- if (inverse) .diff1(hinv, x) else hx
    gyd <- if (inverse) .diff1(hinv, y) else hy
    factor <- cancel(.div(.diff1(g, y), .diff_spec_eval(g, list(list(y, 2L)))))
    factorx <- .diff1(factor, x)
    factory <- .diff1(factor, y)
    out <- function(tau, gxv) {
      if (inverse) list(.ode_inf(xe[[2L]], subs(tau, x, func), xe[[1L]], subs(gxv, x, func)))
      else list(.ode_inf(xe[[1L]], tau, xe[[2L]], gxv))
    }
    if (!.has1(factor, x) && !.has1(factor, y)) {
      A <- Wild("A", exclude = list(y))
      B <- Wild("B", exclude = list(y))
      C <- Wild("C", exclude = list(x, y))
      m <- sym_match(g, A + B * fcall("exp", list(.div(y, C))))
      if (!is.null(m) && .ode_true(m$C %||% S_Zero) && .ode_true(m$B %||% S_Zero)) {
        A_ <- m$A %||% S_Zero
        tau <- .div(fcall("exp", list(.neg(integrate(.div(A_, m$C), x)))), m$B)
        return(out(tau, .mul(list(A_, tau))))
      }
    } else {
      gamma <- cancel(.div(factorx, factory))
      if (!.has1(gamma, y)) {
        tauint <- cancel(.div(.add(list(.mul(list(gamma, gyd)), .neg(.diff1(gamma, x)), .neg(gx))), .add(list(g, gamma))))
        if (!.has1(tauint, y)) {
          tau <- fcall("exp", list(integrate(tauint, x)))
          return(out(tau, .mul(list(.neg(tau), gamma))))
        }
      }
    }
  }
  NULL
}

.ode_lie_abaco2_unique_unknown <- function(match, comp = FALSE) {
  h <- match$h; hx <- match$hx; hy <- match$hy; func <- match$func; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  funclist <- list()
  for (atom in .sp_sorted(atoms(h, "Pow"))) {
    be <- .base_exp(atom)
    if (.has1(be[[1L]], x) && .has1(be[[1L]], y) && be[[2L]]$type != "Integer") funclist[[length(funclist) + 1L]] <- atom
  }
  for (fn in .sp_sorted(Filter(.is_undefined_application, atoms(h, "Function")))) {
    fs <- free_symbols(fn)
    if (.sp_in(x, fs) && .sp_in(y, fs)) funclist[[length(funclist) + 1L]] <- fn
  }
  for (f in funclist) {
    frac <- cancel(.div(.diff1(f, y), .diff1(f, x)))
    sep <- .ode_sep_ok(frac, list(x, y))
    if (!is.null(sep)) {
      xitry1 <- .sv(sep, x)
      etatry1 <- .div(S_NegativeOne, .mul(list(.sv(sep, y), sep$coeff)))
      pde1 <- .add(list(.mul(list(.diff1(etatry1, y), h)), .neg(.mul(list(.diff1(xitry1, x), h))),
                        .neg(.mul(list(xitry1, hx))), .neg(.mul(list(etatry1, hy)))))
      if (!.ode_true(simplify(pde1))) return(list(.ode_inf(xe[[1L]], xitry1, xe[[2L]], subs(etatry1, y, func))))
      xitry2 <- .inv(etatry1)
      etatry2 <- .inv(xitry1)
      pde2 <- .add(list(.diff1(etatry2, x), .neg(.mul(list(.diff1(xitry2, y), .pow(h, S_Two)))),
                        .neg(.mul(list(xitry2, hx))), .neg(.mul(list(etatry2, hy)))))
      if (!.ode_true(simplify(expand(pde2)))) return(list(.ode_inf(xe[[1L]], subs(xitry2, y, func), xe[[2L]], etatry2)))
    } else {
      etatry <- .div(S_NegativeOne, frac)
      pde <- .add(list(.diff1(etatry, x), .mul(list(.diff1(etatry, y), h)), .neg(hx), .neg(.mul(list(etatry, hy)))))
      if (!.ode_true(simplify(pde))) return(list(.ode_inf(xe[[1L]], S_One, xe[[2L]], subs(etatry, y, func))))
      xitry <- .neg(frac)
      pde <- .add(list(.mul(list(.neg(.diff1(xitry, x)), h)), .neg(.mul(list(.diff1(xitry, y), .pow(h, S_Two)))),
                       .neg(.mul(list(xitry, hx))), .neg(hy)))
      if (!.ode_true(simplify(expand(pde)))) return(list(.ode_inf(xe[[1L]], subs(xitry, y, func), xe[[2L]], S_One)))
    }
  }
  NULL
}

.ode_lie_abaco2_unique_general <- function(match, comp = FALSE) {
  hx <- match$hx; hy <- match$hy; func <- match$func; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  d <- function(e, v) .diff1(e, v)
  A <- d(hx, y)
  B <- .add(list(d(hy, y), .pow(hy, S_Two)))
  C <- .sub(d(hx, x), .pow(hx, S_Two))
  if (!(.ode_true(A) && .ode_true(B) && .ode_true(C))) return(NULL)
  Ax <- d(A, x); Ay <- d(A, y); Axy <- d(Ax, y); Axx <- d(Ax, x); Ayy <- d(Ay, y)
  D <- .sub(.mul(list(simplify(.add(list(.mul(list(S_Two, Axy)), .mul(list(hx, Ay)), .neg(.mul(list(Ax, hy))),
                                          .mul(list(.add(list(.mul(list(hx, hy)), .mul(list(S_Two, A)))), A))))), A)),
            .mul(list(.int(3), Ax, Ay)))
  if (!.ode_true(D)) {
    E1 <- simplify(.add(list(.mul(list(.int(3), .pow(Ax, S_Two))),
                             .mul(list(.sub(.mul(list(.add(list(.pow(hx, S_Two), .mul(list(S_Two, C)))), A)), .mul(list(S_Two, Axx))), A)))))
    if (.ode_true(E1)) {
      E2 <- simplify(.sub(.mul(list(.add(list(.mul(list(S_Two, Ayy)), .mul(list(.sub(.mul(list(S_Two, B)), .pow(hy, S_Two)), A)))), A)),
                          .mul(list(.int(3), .pow(Ay, S_Two)))))
      if (!.ode_true(E2)) {
        E3 <- simplify(.sub(.mul(list(E1, .sub(.mul(list(.add(list(.mul(list(.int(28), Ax)), .mul(list(.int(4), hx, A)))), .pow(A, .int(3)))),
                                               .mul(list(E1, .add(list(.mul(list(hy, A)), Ay))))))),
                            .mul(list(d(E1, x), .int(8), .pow(A, .int(4))))))
        if (!.ode_true(E3)) {
          etaval <- cancel(.div(.add(list(.mul(list(.int(4), .pow(A, .int(3)), .sub(Ax, .mul(list(hx, A))))),
                                          .mul(list(E1, .sub(.mul(list(hy, A)), Ay))))), .mul(list(S_Two, A, E1))))
          if (!.has_free(etaval, list(x))) {
            etaval <- fcall("exp", list(integrate(etaval, y)))
            xival <- .div(.mul(list(.int(-4), .pow(A, .int(3)), etaval)), E1)
            if (!.has_free(xival, list(y))) return(list(.ode_inf(xe[[1L]], xival, xe[[2L]], subs(etaval, y, func))))
          }
        }
      }
    }
  } else {
    E1 <- simplify(.sub(.mul(list(.add(list(.mul(list(S_Two, Ayy)), .mul(list(.sub(.mul(list(S_Two, B)), .pow(hy, S_Two)), A)))), A)),
                        .mul(list(.int(3), .pow(Ay, S_Two)))))
    if (.ode_true(E1)) {
      E2 <- simplify(.add(list(.mul(list(.int(4), .pow(A, .int(3)), D)), .neg(.pow(D, S_Two)),
                               .mul(list(E1, .sub(.mul(list(.sub(.mul(list(S_Two, Axx)), .mul(list(.add(list(.pow(hx, S_Two), .mul(list(S_Two, C)))), A))), A)),
                                                  .mul(list(.int(3), .pow(Ax, S_Two)))))))))
      if (!.ode_true(E2)) {
        E3 <- simplify(.add(list(.mul(list(.neg(.mul(list(A, D))), d(E1, y))),
                                 .mul(list(.add(list(.mul(list(.sub(d(E1, x), .mul(list(hy, D))), A)), .mul(list(.int(3), Ay, D)),
                                                     .mul(list(.sub(.mul(list(A, hx)), .mul(list(.int(3), Ax))), E1)))), E1)))))
        if (!.ode_true(E3)) {
          etaval <- cancel(.div(.sub(.mul(list(.sub(.mul(list(A, hx)), Ax), E1)), .mul(list(.add(list(Ay, .mul(list(A, hy)))), D))),
                                .mul(list(S_Two, A, D))))
          if (!.has_free(etaval, list(x))) {
            etaval <- fcall("exp", list(integrate(etaval, y)))
            xival <- .div(.mul(list(.neg(E1), etaval)), D)
            if (!.has_free(xival, list(y))) return(list(.ode_inf(xe[[1L]], xival, xe[[2L]], subs(etaval, y, func))))
          }
        }
      }
    }
  }
  NULL
}

.ode_lie_linear <- function(match, comp = FALSE) {
  h <- match$h; hx <- match$hx; hy <- match$hy; func <- match$func; x <- func$args[[1L]]; y <- match$y
  xe <- .ode_xi_eta(x, func)
  # numbered_symbols("c") read through islice: every second symbol
  symlist <- lapply(c(1L, 3L, 5L, 7L, 9L, 11L), function(i) Dummy(paste0("c", i)))
  C0 <- symlist[[1L]]; C1 <- symlist[[2L]]; C2 <- symlist[[3L]]; C3 <- symlist[[4L]]; C4 <- symlist[[5L]]; C5 <- symlist[[6L]]
  pde <- .add(list(C3, .mul(list(.sub(C4, C0), h)), .neg(.mul(list(.add(list(.mul(list(C0, x)), .mul(list(C1, y)), C2)), hx))),
                   .neg(.mul(list(.add(list(.mul(list(C3, x)), .mul(list(C4, y)), C5)), hy))), .neg(.mul(list(C1, .pow(h, S_Two))))))
  pde <- .as_numer_denom(pde)[[1L]]
  pde <- powsimp(expand(pde))
  ck <- character()
  cv <- list()
  if (pde$type == "Add") {
    for (term in .sympy_args(pde)) {
      if (term$type == "Mul") {
        rem <- .mul(Filter(function(m) !(.has1(m, x) || .has1(m, y)), term$args))
        xypart <- .div(term, rem)
      } else {
        rem <- S_One
        xypart <- term
      }
      j <- match(xypart$key, ck)
      if (is.na(j)) {
        ck <- c(ck, xypart$key)
        cv[[length(cv) + 1L]] <- rem
      } else cv[[j]] <- .add(list(cv[[j]], rem))
    }
  }
  soldict <- .ode_solve_dict(cv, symlist)
  if (length(soldict) && any(vapply(soldict, function(p) .ode_true(p[[2L]]), TRUE))) {
    ones <- lapply(symlist, function(s) list(s, S_One))
    xival <- subs(subs(.add(list(.mul(list(C0, x)), .mul(list(C1, func)), C2)), soldict), ones)
    etaval <- subs(subs(.add(list(.mul(list(C3, x)), .mul(list(C4, func)), C5)), soldict), ones)
    return(list(.ode_inf(xe[[1L]], xival, xe[[2L]], etaval)))
  }
  NULL
}

#' Infinitesimals of a first order ODE
#'
#' `infinitesimals(eq, func)` finds infinitesimals `xi(x, f(x))` and
#' `eta(x, f(x))` of the Lie point symmetries of the first order ODE `eq`,
#' as SymPy's `infinitesimals()`: the heuristics `abaco1_simple`,
#' `abaco1_product`, `abaco2_similar`, `abaco2_unique_unknown`,
#' `abaco2_unique_general`, `linear`, `function_sum`, `bivariate` and `chi`
#' are tried in this order.
#'
#' @param eq The ODE, an expression or an equation.
#' @param func The unknown function; detected when `NULL`.
#' @param order The order of the ODE; computed when `NULL`.
#' @param hint `"default"` for the first heuristic that succeeds, `"all"` for
#'   the infinitesimals of every heuristic, or the name of one heuristic.
#' @param match Internal: the right hand side `h` and the variable `y` used
#'   by the `lie_group` solver.
#' @return A list of dictionaries keyed by `xi(x, f(x))` and `eta(x, f(x))`.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' infinitesimals(Eq(Derivative(f(x), x), x^2 * f(x)), f(x))
#' @export
infinitesimals <- function(eq, func = NULL, order = NULL, hint = "default", match = NULL) {
  eq <- as_sym(eq)
  if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  if (is.null(func)) {
    pp <- .ode_preprocess(eq)
    eq <- pp[[1L]]
    func <- pp[[2L]]
  }
  func <- as_sym(func)
  variables <- func$args
  if (length(variables) != 1L) .ode_stop("ValueError", "ODE's have only one independent variable")
  x <- variables[[1L]]
  if (is.null(order)) order <- .ode_order(eq, func)
  if (order != 1L) .ode_nie("Infinitesimals for only first order ODE's have been implemented")
  df <- .ode_d(func, x)
  if (!is.null(match)) {
    h <- match$h
    y <- match$y
  } else {
    a <- Wild("a", exclude = list(df))
    b <- Wild("b", exclude = list(df))
    m <- sym_match(collect(expand(eq), df), a * df + b)
    if (!is.null(m)) h <- .neg(simplify(.div(m$b %||% S_Zero, m$a %||% S_Zero)))
    else {
      sol <- tryCatch(.sp_solve(eq, list(df), list()), NotImplementedError = function(e)
        .ode_nie("Infinitesimals for the first order ODE could not be found"))
      h <- sol[[1L]]
    }
    y <- Dummy("y")
    h <- subs(h, func, y)
  }
  u <- Dummy("u")
  hx <- .diff1(h, x)
  hy <- .diff1(h, y)
  hinv <- subs(subs(.inv(h), list(list(x, u), list(y, x))), u, y)
  match <- list(h = h, func = func, hx = hx, hy = hy, y = y, hinv = hinv)
  heur <- function(nm) get(paste0(".ode_lie_", nm), mode = "function")
  if (hint == "all") {
    xieta <- list()
    for (nm in .ODE_LIE_HEURISTICS) {
      inflist <- heur(nm)(match, comp = TRUE)
      if (length(inflist)) for (inf in inflist) if (!.ode_inf_in(inf, xieta)) xieta[[length(xieta) + 1L]] <- inf
    }
    if (length(xieta)) return(xieta)
    .ode_nie("Infinitesimals could not be found for the given ODE")
  } else if (hint == "default") {
    for (nm in .ODE_LIE_HEURISTICS) {
      xieta <- heur(nm)(match, comp = FALSE)
      if (length(xieta)) return(xieta)
    }
    .ode_nie("Infinitesimals could not be found for the given ODE")
  } else if (!(hint %in% .ODE_LIE_HEURISTICS)) {
    .ode_stop("ValueError", paste0("Heuristic not recognized: ", hint))
  } else {
    xieta <- heur(hint)(match, comp = TRUE)
    if (length(xieta)) return(xieta)
    .ode_stop("ValueError", "Infinitesimals could not be found using the given heuristic")
  }
}

# lie_group._lie_group_remove
.ode_lie_group_remove <- function(coords) {
  if (.is_undefined_application(coords)) return(coords$args[[1L]])
  if (coords$type == "Add") {
    for (fn in Filter(.is_undefined_application, atoms(coords, "Function"))) coords <- subs(coords, fn, S_Zero)
    return(coords)
  }
  if (coords$type == "Pow") {
    be <- .base_exp(coords)
    return(.pow(.ode_lie_group_remove(be[[1L]]), .ode_lie_group_remove(be[[2L]])))
  }
  if (coords$type == "Mul") return(.mul(lapply(coords$args, .ode_lie_group_remove)))
  coords
}

# lie_group._ode_lie_group_try_heuristic
.ode_lie_group_try_heuristic <- function(eq, heuristic, func, match, inf) {
  xi <- function(...) fcall("xi", list(...))
  eta <- function(...) fcall("eta", list(...))
  fname <- func$name
  F <- function(...) fcall(fname, list(...))
  x <- func$args[[1L]]
  y <- match$y
  h <- match$h
  tempsol <- list()
  if (!length(inf)) {
    inf <- tryCatch(infinitesimals(eq, hint = heuristic, func = func, order = 1L, match = match),
                    ValueError = function(e) NULL)
    if (is.null(inf)) return(NULL)
  }
  C1 <- Symbol("C1")
  for (infsim in inf) {
    xiinf <- subs(.ode_inf_get(infsim, xi(x, func), "xi"), func, y)
    etainf <- subs(.ode_inf_get(infsim, eta(x, func), "eta"), func, y)
    if (identical(simplify(.div(etainf, xiinf))$key, h$key)) next
    fxy <- F(x, y)
    rpde <- .add(list(.mul(list(.diff1(fxy, x), xiinf)), .mul(list(.diff1(fxy, y), etainf))))
    r <- pdsolve(rpde, func = fxy)$args[[2L]]
    s <- pdsolve(.sub(rpde, S_One), func = fxy)$args[[2L]]
    newcoord <- lapply(list(r, s), .ode_lie_group_remove)
    r <- Dummy("r")
    s <- Dummy("s")
    rcoord <- newcoord[[1L]]
    scoord <- newcoord[[2L]]
    sol <- tryCatch(.sp_solve(list(.sub(r, rcoord), .sub(s, scoord)), list(x, y), list(dict = TRUE, `_raw` = TRUE)),
                    NotImplementedError = function(e) NULL)
    if (!length(sol)) next
    sd <- sol[[1L]]
    getv <- function(v) sd$v[[match(v$key, vapply(sd$k, .k, ""))]]
    xsub <- getv(x)
    ysub <- getv(y)
    num <- simplify(.add(list(.diff1(scoord, x), .mul(list(.diff1(scoord, y), h)))))
    denom <- simplify(.add(list(.diff1(rcoord, x), .mul(list(.diff1(rcoord, y), h)))))
    if (.ode_true(num) && .ode_true(denom)) {
      diffeq <- simplify(subs(.div(num, denom), list(list(x, xsub), list(y, ysub))))
      sep <- .ode_separatevars(diffeq, dict = TRUE, symbols = list(r, s))
      if (!is.null(sep)) {
        deq <- .add(list(integrate(.inv(.sv(sep, s)), s), C1, .neg(integrate(.mul(list(sep$coeff, .sv(sep, r))), r))))
        deq <- subs(deq, list(list(r, rcoord), list(s, scoord)))
        sdeq <- tryCatch(.sp_solve(deq, list(y), list()), NotImplementedError = function(e) NULL)
        if (is.null(sdeq)) tempsol[[length(tempsol) + 1L]] <- deq
        else return(lapply(sdeq, function(v) Eq(F(x), v)))
      }
    } else if (.ode_true(denom)) {
      return(list(Eq(F(x), .sp_solve(.sub(scoord, C1), list(y), list())[[1L]])))
    } else if (.ode_true(num)) {
      return(list(Eq(F(x), .sp_solve(.sub(rcoord, C1), list(y), list())[[1L]])))
    }
  }
  if (length(tempsol)) return(lapply(tempsol, function(v) Eq(subs(v, y, F(x)), S_Zero)))
  NULL
}

# lie_group._ode_lie_group
.ode_lie_group <- function(s, func, order, match) {
  heuristics <- .ODE_LIE_HEURISTICS
  inf <- list()
  x <- func$args[[1L]]
  df <- .ode_d(func, x)
  xis <- match$xi
  etas <- match$eta
  y <- match$y
  if (!is.null(y)) h <- .neg(simplify(.div(match$d, match$e)))
  else {
    y <- Dummy("y")
    h <- subs(s, func, y)
  }
  if (!is.null(xis) && !is.null(etas)) {
    xe <- .ode_xi_eta(x, func)
    inf <- list(.ode_inf(xe[[1L]], xis, xe[[2L]], etas))
    if (isTRUE(checkinfsol(Eq(df, s), inf, func = func, order = 1L)[[1L]][[1L]])) heuristics <- c("user_defined", heuristics)
  }
  m <- list(h = h, y = y)
  sol <- NULL
  for (heuristic in heuristics) {
    sol <- .ode_lie_group_try_heuristic(Eq(df, s), heuristic, func, m, inf)
    if (length(sol)) return(sol)
  }
  sol
}
