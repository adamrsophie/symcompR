# The solvers of single ODEs, after sympy/solvers/ode/single.py: the
# SingleODEProblem and each SingleODESolver subclass, in SymPy's order.
#
# A problem is an environment with cached properties. A solver class is a
# list of methods, a subclass copying its parent's methods before setting its
# own; a solver is an environment of class "ode_solver" holding its class,
# its problem and the attributes its methods set, as the Python objects do.

.ode_classes <- new.env()

.ode_defclass <- function(name, parent = NULL, ...) {
  def <- if (is.null(parent)) list() else get(parent, envir = .ode_classes)
  new <- list(...)
  for (k in names(new)) def[[k]] <- new[[k]]
  def$name <- name
  assign(name, def, envir = .ode_classes)
  invisible(def)
}

.ode_cls <- function(name) get(name, envir = .ode_classes)

.ode_solver_new <- function(cls, prob) {
  self <- new.env()
  self$def <- .ode_cls(cls)
  self$name <- cls
  self$prob <- prob
  self$matched <- NULL
  class(self) <- "ode_solver"
  self
}

#' @export
format.ode_solver <- function(x, ...) paste0("<", x$name, ">")
#' @export
print.ode_solver <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# SingleODESolver.matches
.ode_matches <- function(self) {
  ord <- self$def$order
  if (!is.null(ord) && !(.op_order(self$prob) %in% ord)) {
    self$matched <- FALSE
    return(FALSE)
  }
  if (is.null(self$matched)) self$matched <- isTRUE(self$def$matches(self))
  self$matched
}

# SingleODESolver.get_general_solution
.ode_get_general_solution <- function(self, simplify = TRUE) {
  if (!.ode_matches(self))
    .ode_stop("ODEMatchError", sprintf("%s solver cannot solve:\n%s", self$def$hint, .str(self$prob$eq)), "NotImplementedError")
  self$def$gen(self, simplify)
}

# ---------------------------------------------------------------------------
# SingleODEProblem

.ode_problem <- function(eq, func, sym, prep = TRUE, params = list()) {
  p <- new.env()
  p$eq <- eq
  p$func <- func
  p$sym <- sym
  p$prep <- prep
  p$params <- params
  p
}

.op_order <- function(p) {
  if (is.null(p$order_)) p$order_ <- .ode_order(p$eq, p$func)
  p$order_
}

.op_eq_preprocessed <- function(p) {
  if (is.null(p$eq_pre_)) {
    if (p$prep) {
      pp <- .ode_preprocess(p$eq, p$func)
      if (!identical(pp[[2L]]$key, p$func$key)) .ode_stop("ValueError", "")
      p$eq_pre_ <- pp[[1L]]
    } else p$eq_pre_ <- p$eq
  }
  p$eq_pre_
}

.op_eq_high_order_free <- function(p) {
  if (is.null(p$eq_hof_)) {
    a <- Wild("a", exclude = list(p$func))
    c1 <- Wild("c1", exclude = list(p$sym))
    reduced_eq <- NULL
    if (p$eq$type == "Add") {
      deriv_coef <- coeff(p$eq, .ode_d(p$func, p$sym, .op_order(p)))
      if (!is_one_num(deriv_coef) && !is_zero_num(deriv_coef)) {
        r <- sym_match(deriv_coef, a * .pow(p$func, c1))
        if (!is.null(r) && .ode_true(r$c1)) {
          den <- .pow(p$func, r$c1)
          reduced_eq <- .add(lapply(p$eq$args, function(arg) .div(arg, den)))
        }
      }
    }
    if (is.null(reduced_eq)) reduced_eq <- expand(p$eq)
    p$eq_hof_ <- reduced_eq
  }
  p$eq_hof_
}

.op_eq_expanded <- function(p) {
  if (is.null(p$eq_exp_)) p$eq_exp_ <- expand(.op_eq_preprocessed(p))
  p$eq_exp_
}

.op_constants <- function(p, num = 1L, start = 1L, prefix = "C") .ode_numbered_constants(list(p$eq), num, start, prefix)

# Basic.free_symbols with SymPy's Derivative.free_symbols: the symbols of
# the differentiated expression and of the counts, not the variables
.ode_free_symbols <- function(e) {
  # a Tuple stands for the Derivative so that no term cancels against it
  e <- .map_bottom_up(e, function(z) if (z$type == "Derivative") new_sym("Tuple", c(list(z$args[[1L]]), lapply(z$args[-1L], function(t) t$args[[2L]]))) else z)
  free_symbols(e)
}

.op_is_autonomous <- function(p) {
  if (is.null(p$auto_)) {
    u <- Dummy("u")
    syms <- .ode_free_symbols(subs(p$eq, p$func, u))
    p$auto_ <- !.sp_in(p$sym, syms)
  }
  p$auto_
}

# SingleODEProblem.get_linear_coefficients: the coefficients of f(x) and its
# derivatives keyed "0", ..., "order", and the right hand side keyed "-1"
.op_linear_coefficients <- function(eq, func, order) {
  x <- func$args[[1L]]
  ders <- lapply(0:order, function(i) .ode_d(func, x, i))
  r <- tryCatch(.pv_lin_eq2dict(eq, ders), NonlinearError = function(e) NULL)
  if (is.null(r)) return(NULL)
  rhs <- r[[1L]]
  lhs <- r[[2L]]
  if (.has1(rhs, func) || any(vapply(lhs$values, function(c) .has1(c, func), TRUE))) return(NULL)
  keys <- vapply(lhs$keys, .k, "")
  terms <- list()
  for (i in 0:order) {
    j <- match(ders[[i + 1L]]$key, keys)
    terms[[as.character(i)]] <- if (is.na(j)) S_Zero else lhs$values[[j]]
  }
  terms[["-1"]] <- rhs
  terms
}

.rk <- function(r, i) r[[as.character(i)]]
.rk_ge0 <- function(r) names(r)[as.integer(names(r)) >= 0L]

# ---------------------------------------------------------------------------
# SinglePatternODESolver

.ode_wilds <- function(self) {
  p <- self$prob
  self$def$wilds_(self, p$func, p$sym, .op_order(p))
}

.ode_wilds_match <- function(self) {
  lapply(.ode_wilds(self), function(w) {
    v <- self$wm[[w$name]]
    if (is.null(v)) S_Zero else v
  })
}

.ode_pattern_matches <- function(self) {
  p <- self$prob
  eq <- .op_eq_expanded(p)
  fx <- p$func
  x <- p$sym
  order <- .op_order(p)
  df <- .ode_d(fx, x, order)
  if (!(order %in% c(1L, 2L))) return(FALSE)
  pattern <- self$def$equation(self, fx, x, order)
  if (!.has_type(coeff(pattern, df), "Wild")) eq <- expand(.div(eq, coeff(eq, df)))
  eq <- collect(eq, list(.ode_d(fx, x), fx), func = cancel)
  m <- sym_match(eq, pattern)
  self$wm <- m
  if (!is.null(m)) return(isTRUE(self$def$verify(self, fx)))
  FALSE
}

.ode_defclass("SingleODESolver",
  matches = function(self) .ode_nie("Subclasses of SingleODESolver should implement matches."),
  gen = function(self, simplify_flag = TRUE) .ode_nie("Subclasses of SingleODESolver should implement get_general_solution."))

.ode_defclass("SinglePatternODESolver", "SingleODESolver",
  matches = .ode_pattern_matches,
  verify = function(self, fx) TRUE)

# ---------------------------------------------------------------------------
# NthAlgebraic

.ode_defclass("NthAlgebraic", "SingleODESolver",
  hint = "nth_algebraic", has_integral = TRUE,
  matches = function(self) {
    eq <- self$prob$eq
    func <- self$prob$func
    var <- self$prob$sym
    diffx <- function(e) .ode_diffx(e, var)
    expand_diffx <- function(e) {
      if (e$type != "Derivative") return(e)
      toreplace <- e$args[[1L]]
      for (vc in .ode_variable_count(e)) {
        for (k in seq_len(vc[[2L]])) {
          if (identical(vc[[1L]]$key, var$key)) toreplace <- diffx(toreplace)
          else toreplace <- Derivative(toreplace, vc[[1L]])
        }
      }
      toreplace
    }
    unreplace <- function(e) .map_bottom_up(e, function(z) if (z$type == "Function" && identical(z$name, "diffx")) Derivative(z$args[[1L]], var) else z)
    subs_eqn <- .map_bottom_up(eq, expand_diffx)
    old <- .ode_state$diffx_var
    .ode_state$diffx_var <- var
    on.exit(.ode_state$diffx_var <- old)
    solns <- tryCatch(.sp_solve(subs_eqn, list(func), list(simplify = FALSE)), NotImplementedError = function(e) list())
    solns <- lapply(solns, function(s) simplify(unreplace(s)))
    self$solutions <- lapply(solns, function(s) Eq(func, s))
    length(self$solutions) != 0L
  },
  gen = function(self, simplify_flag = TRUE) self$solutions)

# NthAlgebraic._get_diffx: an undefined function that is "invertible" with
# the inverse expr -> Integral(expr, var) + Dummy('C')
.ode_diffx <- function(e, var) fcall("diffx", list(e))
.ode_diffx_inverse <- function(rhs) {
  var <- .ode_state$diffx_var
  .add(list(Integral(rhs, var), Dummy("C")))
}

# ---------------------------------------------------------------------------
# FirstExact

.ode_defclass("FirstExact", "SinglePatternODESolver",
  hint = "1st_exact", has_integral = TRUE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    df <- .ode_d(fx, x)
    list(Wild("P", exclude = list(df)), Wild("Q", exclude = list(df)))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    w[[1L]] + w[[2L]] * .ode_d(fx, x)
  },
  verify = function(self, fx) {
    x <- self$prob$sym
    y <- Dummy("y")
    mn <- .ode_wilds_match(self)
    m <- subs(mn[[1L]], fx, y)
    n <- subs(mn[[2L]], fx, y)
    numerator <- cancel(.sub(.diff1(m, y), .diff1(n, x)))
    if (isTRUE(.sp_fact(numerator, "zero"))) return(TRUE)
    factor_n <- cancel(.div(numerator, n))
    factor_m <- cancel(.div(.neg(numerator), m))
    if (!.sp_in(y, free_symbols(factor_n))) {
      fac <- factor_n
      iv <- x
    } else if (!.sp_in(x, free_symbols(factor_m))) {
      fac <- factor_m
      iv <- y
    } else return(FALSE)
    fac <- fcall("exp", list(Integral(fac, iv)))
    m <- .mul(list(m, fac))
    n <- .mul(list(n, fac))
    self$wm$P <- subs(m, y, fx)
    self$wm$Q <- subs(n, y, fx)
    TRUE
  },
  gen = function(self, simplify_flag = TRUE) {
    mn <- .ode_wilds_match(self)
    fx <- self$prob$func
    x <- self$prob$sym
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    y <- Dummy("y")
    m <- subs(mn[[1L]], fx, y)
    n <- subs(mn[[2L]], fx, y)
    im <- Integral(m, x)
    inner <- .add(list(im, Integral(.sub(n, .diff1(im, y)), y)))
    list(Eq(Subs(inner, y, fx), C1))
  })

# ---------------------------------------------------------------------------
# FirstLinear

.ode_defclass("FirstLinear", "SinglePatternODESolver",
  hint = "1st_linear", has_integral = TRUE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    list(Wild("P", exclude = list(fx)), Wild("Q", exclude = list(fx, .ode_d(fx, x))))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    .ode_d(fx, x) + w[[1L]] * fx - w[[2L]]
  },
  gen = function(self, simplify_flag = TRUE) {
    pq <- .ode_wilds_match(self)
    P <- pq[[1L]]; Q <- pq[[2L]]
    fx <- self$prob$func
    x <- self$prob$sym
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    ip <- Integral(P, x)
    list(Eq(fx, .mul(list(.add(list(C1, Integral(.mul(list(Q, fcall("exp", list(ip)))), x))),
                          fcall("exp", list(.neg(ip)))))))
  })

# ---------------------------------------------------------------------------
# AlmostLinear

.ode_defclass("AlmostLinear", "SinglePatternODESolver",
  hint = "almost_linear", has_integral = TRUE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    df <- .ode_d(fx, x)
    list(Wild("P", exclude = list(df)), Wild("Q", exclude = list(df)))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    w[[1L]] * .ode_d(fx, x) + w[[2L]]
  },
  verify = function(self, fx) {
    ab <- .ode_wilds_match(self)
    a <- ab[[1L]]; b <- ab[[2L]]
    if (b$type == "Add") {
      cb <- .sp_indep(b, list(fx))
      c <- cb[[1L]]; b <- cb[[2L]]
    } else c <- S_Zero
    db <- .diff1(b, fx)
    if (!is_zero_num(db) && !.has1(simplify(.div(db, a)), fx)) {
      fb <- factor_terms(b)
      self$ly <- .sp_indep(fb, list(fx), as_add = FALSE)[[2L]]
      self$ax <- .div(a, .diff1(self$ly, fx))
      self$cx <- .neg(c)
      self$bx <- .div(fb, self$ly)
      return(TRUE)
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    x <- self$prob$sym
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    ib <- Integral(.div(self$bx, self$ax), x)
    list(Eq(self$ly, .mul(list(.add(list(C1, Integral(.mul(list(.div(self$cx, self$ax), fcall("exp", list(ib)))), x))),
                               fcall("exp", list(.neg(ib)))))))
  })

# ---------------------------------------------------------------------------
# Bernoulli

.ode_defclass("Bernoulli", "SinglePatternODESolver",
  hint = "Bernoulli", has_integral = TRUE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    list(Wild("P", exclude = list(fx)), Wild("Q", exclude = list(fx)), Wild("n", exclude = list(x, fx, .ode_d(fx, x))))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    .ode_d(fx, x) + w[[1L]] * fx - w[[2L]] * .pow(fx, w[[3L]])
  },
  gen = function(self, simplify_flag = TRUE) {
    w <- .ode_wilds_match(self)
    P <- w[[1L]]; Q <- w[[2L]]; n <- w[[3L]]
    fx <- self$prob$func
    x <- self$prob$sym
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    if (is_one_num(n)) {
      gensol <- Eq(fcall("log", list(fx)), .add(list(C1, Integral(.add(list(.neg(P), Q)), x))))
    } else {
      ip <- Integral(P, x)
      inner <- Integral(.mul(list(Q, fcall("exp", list(.mul(list(.neg(n), ip)))), fcall("exp", list(ip)))), x)
      gensol <- Eq(.pow(fx, .sub(S_One, n)),
                   .mul(list(.sub(C1, .mul(list(.sub(n, S_One), inner))),
                             fcall("exp", list(.mul(list(.neg(.sub(S_One, n)), ip)))))))
    }
    list(gensol)
  })

# ---------------------------------------------------------------------------
# Factorable

.ode_defclass("Factorable", "SingleODESolver",
  hint = "factorable", has_integral = FALSE,
  matches = function(self) {
    eq_orig <- self$prob$eq
    fx <- self$prob$func
    nd <- .as_numer_denom(eq_orig)
    fl <- factor_list(nd[[1L]])
    ms <- vapply(fl$factors, function(p) as.numeric(p[[2L]]), 0)
    is_reduced <- sum(ms) > 1 || .has1(nd[[2L]], fx)
    self$eqs <- lapply(Filter(function(p) .has1(p[[1L]], fx), fl$factors), `[[`, 1L)
    is_reduced
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    sols <- list()
    eq <- NULL
    for (eq in self$eqs) {
      sol <- tryCatch(dsolve(eq, fx), NotImplementedError = function(e) NULL)
      if (is.null(sol)) next
      sols <- c(sols, .ode_sol_list(sol))
    }
    if (!length(sols)) .ode_nie(paste0("The given ODE ", .str(eq), " cannot be solved by the factorable method"))
    sols
  })

# ---------------------------------------------------------------------------
# RiccatiSpecial

.ode_defclass("RiccatiSpecial", "SinglePatternODESolver",
  hint = "Riccati_special_minus2", has_integral = FALSE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    df <- .ode_d(fx, x)
    list(Wild("a", exclude = list(x, fx, df, S_Zero)), Wild("b", exclude = list(x, fx, df, S_Zero)),
         Wild("c", exclude = list(x, fx, df)), Wild("d", exclude = list(x, fx, df)))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    w[[1L]] * .ode_d(fx, x) + w[[2L]] * .pow(fx, S_Two) + w[[3L]] * fx / x + w[[4L]] / .pow(x, S_Two)
  },
  gen = function(self, simplify_flag = TRUE) {
    w <- .ode_wilds_match(self)
    a <- w[[1L]]; b <- w[[2L]]; c <- w[[3L]]; d <- w[[4L]]
    fx <- self$prob$func
    x <- self$prob$sym
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    mu <- .sqrt(.sub(.mul(list(S_Two, S_Two, d, b)), .pow(.sub(a, c), S_Two)))
    tn <- fcall("tan", list(.add(list(.mul(list(.div(mu, .mul(list(S_Two, a))), fcall("log", list(x)))), C1))))
    list(Eq(fx, .div(.sub(.sub(a, c), .mul(list(mu, tn))), .mul(list(S_Two, b, x)))))
  })

# ---------------------------------------------------------------------------
# RationalRiccati

.ode_defclass("RationalRiccati", "SinglePatternODESolver",
  hint = "1st_rational_riccati", has_integral = FALSE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    ex <- list(fx, .ode_d(fx, x))
    list(Wild("b0", exclude = ex), Wild("b1", exclude = ex), Wild("b2", exclude = ex))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    .ode_d(fx, x) - w[[1L]] - w[[2L]] * fx - w[[3L]] * .pow(fx, S_Two)
  },
  matches = function(self) {
    eq <- .op_eq_expanded(self$prob)
    fx <- self$prob$func
    x <- self$prob$sym
    if (.op_order(self$prob) != 1L) return(FALSE)
    mf <- .ode_match_riccati(eq, fx, x)
    if (!isTRUE(mf[[1L]])) return(FALSE)
    funcs <- mf[[2L]]
    self$wm <- list(b0 = funcs[[1L]], b1 = funcs[[2L]], b2 = funcs[[3L]])
    TRUE
  },
  gen = function(self, simplify_flag = TRUE) {
    b <- .ode_wilds_match(self)
    .ode_solve_riccati(self$prob$func, self$prob$sym, b[[1L]], b[[2L]], b[[3L]], gensol = TRUE)
  })

# ---------------------------------------------------------------------------
# SecondNonlinearAutonomousConserved

.ode_defclass("SecondNonlinearAutonomousConserved", "SinglePatternODESolver",
  hint = "2nd_nonlinear_autonomous_conserved", has_integral = TRUE, order = 2L,
  wilds_ = function(self, fx, x, order) list(Wild("fy", exclude = list(S_Zero, .ode_d(fx, x), .ode_d(fx, x, 2L)))),
  equation = function(self, fx, x, order) .ode_d(fx, x, 2L) + .ode_wilds(self)[[1L]],
  verify = function(self, fx) .op_is_autonomous(self$prob),
  gen = function(self, simplify_flag = TRUE) {
    g <- .ode_wilds_match(self)[[1L]]
    fx <- self$prob$func
    x <- self$prob$sym
    u <- Dummy("u")
    g <- subs(g, fx, u)
    cs <- .op_constants(self$prob, 2L)
    inside <- .add(list(.mul(list(.int(-2L), Integral(g, u))), cs[[1L]]))
    lhs <- Integral(.inv(.sqrt(inside)), list(u, fx))
    list(Eq(lhs, .add(list(cs[[2L]], x))), Eq(lhs, .sub(cs[[2L]], x)))
  })

# ---------------------------------------------------------------------------
# Liouville

.ode_defclass("Liouville", "SinglePatternODESolver",
  hint = "Liouville", has_integral = TRUE, order = 2L,
  wilds_ = function(self, fx, x, order) {
    df <- .ode_d(fx, x)
    list(Wild("d", exclude = list(df, .ode_d(fx, x, 2L))), Wild("e", exclude = list(df)), Wild("k", exclude = list(df)))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    w[[1L]] * .ode_d(fx, x, 2L) + w[[2L]] * .pow(.ode_d(fx, x), S_Two) + w[[3L]] * .ode_d(fx, x)
  },
  verify = function(self, fx) {
    w <- .ode_wilds_match(self)
    self$y <- Dummy("y")
    x <- self$prob$sym
    self$g <- subs(simplify(.div(w[[2L]], w[[1L]])), fx, self$y)
    self$h <- subs(simplify(.div(w[[3L]], w[[1L]])), fx, self$y)
    if (.sp_in(self$y, free_symbols(self$h)) || .sp_in(x, free_symbols(self$g))) return(FALSE)
    TRUE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    x <- self$prob$sym
    cs <- .op_constants(self$prob, 2L)
    int <- Integral(fcall("exp", list(Integral(self$g, self$y))), list(self$y, fx))
    list(Eq(.add(list(int, .mul(list(cs[[1L]], Integral(fcall("exp", list(.neg(Integral(self$h, x)))), x))), cs[[2L]])), S_Zero))
  })

# ---------------------------------------------------------------------------
# Separable and SeparableReduced

.ode_defclass("Separable", "SinglePatternODESolver",
  hint = "separable", has_integral = TRUE, order = 1L,
  wilds_ = function(self, fx, x, order) {
    df <- .ode_d(fx, x)
    list(Wild("d", exclude = list(df, .ode_d(fx, x, 2L))), Wild("e", exclude = list(df)))
  },
  equation = function(self, fx, x, order) {
    w <- .ode_wilds(self)
    w[[1L]] + w[[2L]] * .ode_d(fx, x)
  },
  verify = function(self, fx) {
    w <- .ode_wilds_match(self)
    self$y <- Dummy("y")
    x <- self$prob$sym
    d <- .ode_separatevars(subs(w[[1L]], fx, self$y))
    e <- .ode_separatevars(subs(w[[2L]], fx, self$y))
    self$m1 <- .ode_separatevars(d, dict = TRUE, symbols = list(x, self$y))
    self$m2 <- .ode_separatevars(e, dict = TRUE, symbols = list(x, self$y))
    !is.null(self$m1) && !is.null(self$m2)
  },
  match_object = function(self) list(self$m1, self$m2, self$prob$sym, self$prob$func),
  gen = function(self, simplify_flag = TRUE) {
    mo <- self$def$match_object(self)
    m1 <- mo[[1L]]; m2 <- mo[[2L]]; x <- mo[[3L]]; fx <- mo[[4L]]
    C1 <- .op_constants(self$prob, 1L)[[1L]]
    y <- self$y
    int <- Integral(.div(.mul(list(.sv(m2, "coeff"), .sv(m2, y))), .sv(m1, y)), list(y, fx))
    gen_sol <- Eq(int, .add(list(Integral(.div(.mul(list(.neg(.sv(m1, "coeff")), .sv(m1, x))), .sv(m2, x)), x), C1)))
    list(gen_sol)
  })

.ode_sr_degree <- function(expr, x) {
  for (val in expr) {
    if (.has1(val, x)) {
      if (val$type == "Pow" && identical(.base_exp(val)[[1L]]$key, x$key)) return(.base_exp(val)[[2L]])
      if (identical(val$key, x$key)) return(S_One)
      return(.ode_sr_degree(.sympy_args(val), x))
    }
  }
  S_Zero
}

.ode_defclass("SeparableReduced", "Separable",
  hint = "separable_reduced", has_integral = TRUE, order = 1L,
  powers = function(self, expr) {
    fx <- self$prob$func
    x <- self$prob$sym
    self$y <- Dummy("y")
    exprs <- if (expr$type %in% c("Add", "Mul", "Pow")) .ode_atoms_type(expr, expr$type) else list(expr)
    pows <- list()
    for (arg in exprs) {
      if (.has1(arg, x)) {
        u <- .sp_indep(arg, list(x, fx))[[2L]]
        uy <- list(subs(u, fx, self$y))
        pw <- .div(.ode_sr_degree(uy, x), .ode_sr_degree(uy, self$y))
        pows[[pw$key]] <- pw
      }
    }
    pows
  },
  verify = function(self, fx) {
    w <- .ode_wilds_match(self)
    x <- self$prob$sym
    fac <- simplify(.mul(list(.div(x, fx), .div(w[[1L]], w[[2L]]))))
    nd <- .as_numer_denom(fac)
    num <- expand(nd[[1L]])
    dem <- expand(nd[[2L]])
    pows <- self$def$powers(self, num)
    for (p in self$def$powers(self, dem)) pows[[p$key]] <- p
    pows <- unname(pows)
    if (length(pows) == 1L && !identical(pows[[1L]]$key, S_ComplexInfinity$key)) {
      self$t <- Dummy("t")
      self$r2 <- list(t = self$t)
      pat <- .mul(list(.pow(x, pows[[1L]]), fx))
      num <- subs(num, pat, self$t)
      dem <- subs(dem, pat, self$t)
      test <- .div(num, dem)
      free <- free_symbols(test)
      if (length(free) == 1L && identical(free[[1L]]$key, self$t$key)) {
        self$r2$power <- pows[[1L]]
        self$r2$u <- test
        return(TRUE)
      }
      return(FALSE)
    }
    FALSE
  },
  match_object = function(self) {
    fx <- self$prob$func
    x <- self$prob$sym
    u <- subs(self$r2$u, self$r2$t, self$y)
    ycoeff <- .inv(.mul(list(self$y, .sub(self$r2$power, u))))
    m1 <- .sv_new(list(self$y, x), list(S_One, .neg(.inv(x))), S_One)
    m2 <- .sv_new(list(self$y, x), list(ycoeff, S_One), S_One)
    list(m1, m2, x, .mul(list(.pow(x, self$r2$power), fx)))
  })

.ode_atoms_type <- function(e, type) {
  out <- list()
  walk <- function(z) {
    if (z$type == type) out[[z$key]] <<- z
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

# ---------------------------------------------------------------------------
# Homogeneous coefficients

.ode_hc_wilds <- function(self, fx, x, order) {
  df <- .ode_d(fx, x)
  list(Wild("d", exclude = list(df, .ode_d(fx, x, 2L))), Wild("e", exclude = list(df)))
}
.ode_hc_equation <- function(self, fx, x, order) {
  w <- .ode_wilds(self)
  w[[1L]] + w[[2L]] * .ode_d(fx, x)
}

.ode_hc_match_object <- function(self) {
  self$u1 <- Dummy("u1")
  list(d = self$d, e = self$e, fx = self$prob$func, x = self$prob$sym, u = self$u, u1 = self$u1, y = self$y,
       xarg = S_Zero, yarg = S_Zero)
}

.ode_hc_finish <- function(sol, mo) {
  fx <- mo$fx; u <- mo$u; x <- mo$x
  sol <- subs(sol, fx, u)
  sol <- subs(sol, list(list(u, .sub(u, mo$yarg)), list(x, .sub(x, mo$xarg)), list(u, fx)))
  sol
}

.ode_hc_verify_dep <- function(self, fx) {
  w <- .ode_wilds_match(self)
  self$y <- Dummy("y")
  x <- self$prob$sym
  self$d <- .ode_separatevars(subs(w[[1L]], fx, self$y))
  self$e <- .ode_separatevars(subs(w[[2L]], fx, self$y))
  ordera <- homogeneous_order(self$d, x, self$y)
  orderb <- homogeneous_order(self$e, x, self$y)
  if (!is.null(ordera) && !is.null(orderb) && identical(ordera$key, orderb$key)) {
    self$u <- Dummy("u")
    v <- simplify(subs(.add(list(self$d, .mul(list(self$u, self$e)))), list(list(x, S_One), list(self$y, self$u)), simultaneous = TRUE))
    return(!is_zero_num(v))
  }
  FALSE
}

.ode_hc_verify_indep <- function(self, fx) {
  w <- .ode_wilds_match(self)
  self$y <- Dummy("y")
  x <- self$prob$sym
  self$d <- .ode_separatevars(subs(w[[1L]], fx, self$y))
  self$e <- .ode_separatevars(subs(w[[2L]], fx, self$y))
  ordera <- homogeneous_order(self$d, x, self$y)
  orderb <- homogeneous_order(self$e, x, self$y)
  if (!is.null(ordera) && !is.null(orderb) && identical(ordera$key, orderb$key)) {
    self$u <- Dummy("u")
    v <- simplify(subs(.add(list(self$e, .mul(list(self$u, self$d)))), list(list(x, self$u), list(self$y, S_One)), simultaneous = TRUE))
    return(!is_zero_num(v))
  }
  FALSE
}

.ode_hc_gen_dep <- function(self, simplify_flag = TRUE) {
  mo <- self$def$match_object(self)
  C1 <- .op_constants(self$prob, 1L)[[1L]]
  d <- mo$d; e <- mo$e; u1 <- mo$u1; x <- mo$x; y <- mo$y
  integrand <- subs(.neg(.div(e, .add(list(d, .mul(list(u1, e)))))), list(list(x, S_One), list(y, u1)), simultaneous = TRUE)
  int <- Integral(integrand, list(u1, .div(mo$fx, x)))
  sol <- .ode_eq_map(Eq(fcall("log", list(x)), .add(list(int, fcall("log", list(C1))))), function(z) logcombine(z, force = TRUE))
  list(.ode_hc_finish(sol, mo))
}

.ode_hc_gen_indep <- function(self, simplify_flag = TRUE) {
  mo <- self$def$match_object(self)
  C1 <- .op_constants(self$prob, 1L)[[1L]]
  d <- mo$d; e <- mo$e; u1 <- mo$u1; x <- mo$x; y <- mo$y
  integrand <- simplify(subs(.neg(.div(d, .add(list(e, .mul(list(u1, d)))))), list(list(x, u1), list(y, S_One)), simultaneous = TRUE))
  int <- Integral(integrand, list(u1, .div(x, mo$fx)))
  sol <- .ode_eq_map(Eq(fcall("log", list(mo$fx)), .add(list(int, fcall("log", list(C1))))), function(z) logcombine(z, force = TRUE))
  list(.ode_hc_finish(sol, mo))
}

.ode_defclass("HomogeneousCoeffSubsDepDivIndep", "SinglePatternODESolver",
  hint = "1st_homogeneous_coeff_subs_dep_div_indep", has_integral = TRUE, order = 1L,
  wilds_ = .ode_hc_wilds, equation = .ode_hc_equation, verify = .ode_hc_verify_dep,
  match_object = .ode_hc_match_object, gen = .ode_hc_gen_dep)

.ode_defclass("HomogeneousCoeffSubsIndepDivDep", "SinglePatternODESolver",
  hint = "1st_homogeneous_coeff_subs_indep_div_dep", has_integral = TRUE, order = 1L,
  wilds_ = .ode_hc_wilds, equation = .ode_hc_equation, verify = .ode_hc_verify_indep,
  match_object = .ode_hc_match_object, gen = .ode_hc_gen_indep)

.ode_hc_best_gen <- function(self, simplify_flag = TRUE) {
  sol1 <- .ode_hc_gen_indep(self)[[1L]]
  sol2 <- .ode_hc_gen_dep(self)[[1L]]
  fx <- self$prob$func
  if (simplify_flag) {
    sol1 <- .ode_odesimp(self$prob$eq, sol1, fx, "1st_homogeneous_coeff_subs_indep_div_dep")
    sol2 <- .ode_odesimp(self$prob$eq, sol2, fx, "1st_homogeneous_coeff_subs_dep_div_indep")
  }
  sol1 <- .ode_sol_list(sol1)
  sol2 <- .ode_sol_list(sol2)
  c1 <- max(vapply(sol1, .ode_sol_simplicity, 0, func = fx, trysolving = !simplify_flag))
  c2 <- max(vapply(sol2, .ode_sol_simplicity, 0, func = fx, trysolving = !simplify_flag))
  if (c2 < c1) sol2 else sol1
}

# HomogeneousCoeffBest(HomogeneousCoeffSubsIndepDivDep, HomogeneousCoeffSubsDepDivIndep)
.ode_defclass("HomogeneousCoeffBest", "HomogeneousCoeffSubsIndepDivDep",
  hint = "1st_homogeneous_coeff_best", has_integral = FALSE, order = 1L,
  verify = function(self, fx) .ode_hc_verify_indep(self, fx) && .ode_hc_verify_dep(self, fx),
  gen = .ode_hc_best_gen)

# ---------------------------------------------------------------------------
# LinearCoefficients

.ode_lc_match <- function(expr, func) {
  x <- func$args[[1L]]
  abc <- function(eq) {
    eq <- .mexpand(eq)
    c <- .sp_indep(eq, list(x, func), as_add = TRUE)[[1L]]
    if (!(c$type %in% c("Integer", "Rational"))) return(NULL)
    a <- coeff(eq, x)
    if (!(a$type %in% c("Integer", "Rational"))) return(NULL)
    b <- coeff(eq, func)
    if (!(b$type %in% c("Integer", "Rational"))) return(NULL)
    if (identical(eq$key, .add(list(.mul(list(a, x)), .mul(list(b, func)), c))$key)) return(list(a, b, c))
    NULL
  }
  matchf <- function(arg) {
    nd <- .as_numer_denom(together(arg))
    m <- abc(nd[[1L]])
    if (!is.null(m)) {
      a1 <- m[[1L]]; b1 <- m[[2L]]; c1 <- m[[3L]]
      m <- abc(nd[[2L]])
      if (!is.null(m)) {
        a2 <- m[[1L]]; b2 <- m[[2L]]; c2 <- m[[3L]]
        d <- .sub(.mul(list(a2, b1)), .mul(list(a1, b2)))
        if ((.sp_bool(c1) || .sp_bool(c2)) && .sp_bool(d)) return(list(a1, b1, c1, a2, b2, c2, d))
      }
    }
    NULL
  }
  fs <- Filter(function(fi) fi$name != func$name && length(fi$args) == 1L && fi$args[[1L]]$type != "Function", atoms(expr, "Function"))
  m <- if (length(fs)) lapply(fs, function(fi) fi$args[[1L]]) else list(expr)
  m1 <- matchf(m[[length(m)]])
  m <- m[-length(m)]
  if (!is.null(m1) && all(vapply(m, function(mi) .ode_same(matchf(mi), m1), TRUE))) {
    a1 <- m1[[1L]]; b1 <- m1[[2L]]; c1 <- m1[[3L]]; a2 <- m1[[4L]]; b2 <- m1[[5L]]; c2 <- m1[[6L]]; den <- m1[[7L]]
    return(list(.div(.sub(.mul(list(b2, c1)), .mul(list(b1, c2))), den), .div(.sub(.mul(list(a1, c2)), .mul(list(a2, c1))), den)))
  }
  NULL
}

.ode_defclass("LinearCoefficients", "HomogeneousCoeffBest",
  hint = "linear_coefficients", has_integral = TRUE, order = 1L,
  wilds_ = .ode_hc_wilds, equation = .ode_hc_equation,
  verify = function(self, fx) {
    w <- .ode_wilds_match(self)
    self$d <- w[[1L]]; self$e <- w[[2L]]
    wl <- .ode_wilds(self)
    a <- wl[[1L]]; b <- wl[[2L]]
    F <- .div(self$d, self$e)
    x <- self$prob$sym
    params <- .ode_lc_match(F, fx)
    if (!is.null(params)) {
      self$xarg <- params[[1L]]; self$yarg <- params[[2L]]
      u <- Dummy("u"); t <- Dummy("t")
      self$y <- Dummy("y")
      df <- .ode_d(fx, x)
      dummy_eq <- subs(self$prob$eq, list(list(df, t), list(fx, u)))
      reps <- list(list(x, .add(list(x, self$xarg))), list(u, .add(list(u, self$yarg))), list(t, df), list(u, fx))
      dummy_eq <- simplify(subs(dummy_eq, reps))
      r2 <- sym_match(collect(expand(dummy_eq), list(df, fx)), a * df + b)
      if (length(r2)) {
        self$d <- r2$d %||% S_Zero
        self$e <- r2$e %||% S_Zero
        orderd <- homogeneous_order(self$d, x, fx)
        ordere <- homogeneous_order(self$e, x, fx)
        if (!is.null(orderd) && !is.null(ordere) && identical(orderd$key, ordere$key)) {
          self$d <- subs(self$d, fx, self$y)
          self$e <- subs(self$e, fx, self$y)
          return(TRUE)
        }
        return(FALSE)
      }
      return(FALSE)
    }
    FALSE
  },
  match_object = function(self) {
    self$u1 <- Dummy("u1")
    u <- Dummy("u")
    list(d = self$d, e = self$e, fx = self$prob$func, x = self$prob$sym, u = u, u1 = self$u1, y = self$y,
         xarg = self$xarg, yarg = self$yarg)
  })

# ---------------------------------------------------------------------------
# NthOrderReducible

.ode_defclass("NthOrderReducible", "SingleODESolver",
  hint = "nth_order_reducible", has_integral = FALSE,
  matches = function(self) {
    eq <- .op_eq_preprocessed(self$prob)
    func <- self$prob$func
    x <- self$prob$sym
    if (length(func$args) != 1L) .ode_stop("ValueError", "Function must have exactly one argument")
    vc <- list()
    for (d in atoms(eq, "Derivative")) {
      v <- .ode_variable_count(d)
      if (identical(d$args[[1L]]$key, func$key) && length(v) == 1L) vc[[length(vc) + 1L]] <- v[[1L]]
    }
    ords <- unlist(lapply(vc, function(p) if (identical(p[[1L]]$key, x$key)) p[[2L]]))
    if (length(ords) < 2L) return(FALSE)
    self$smallest <- min(ords)
    D <- Dummy()
    if (.has1(subs(eq, .ode_d(func, x, self$smallest), D), func)) return(FALSE)
    TRUE
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- self$prob$eq
    fx <- self$prob$func
    x <- self$prob$sym
    n <- self$smallest
    names <- vapply(Filter(.is_undefined_application, atoms(eq, "Function")), function(a) a$name, "")
    repeat {
      nm <- Dummy()$name
      if (!(nm %in% names)) break
    }
    gx <- fcall(nm, list(x))
    w <- .ode_d(fx, x, n)
    geq <- subs(eq, w, gx)
    gsol <- .ode_sol_list(dsolve(geq, gx))
    lapply(gsol, function(gsoli) dsolve(subs(gsoli, gx, w), fx))
  })

# ---------------------------------------------------------------------------
# SecondHypergeometric

.ode_defclass("SecondHypergeometric", "SingleODESolver",
  hint = "2nd_hypergeometric", has_integral = TRUE,
  matches = function(self) {
    eq <- .op_eq_preprocessed(self$prob)
    func <- self$prob$func
    r <- .ode_match_2nd_hypergeometric(eq, func)
    self$match_object <- NULL
    if (!is.null(r)) {
      d <- .ode_equivalence_hypergeometric(r[[1L]], r[[2L]], func)
      if (!is.null(d) && identical(d$type, "2F1")) {
        mo <- .ode_match_2nd_2F1_hypergeometric(d$I0, d$k, d$sing_point, func)
        if (!is.null(mo)) {
          mo$A <- r[[1L]]
          mo$B <- r[[2L]]
          self$match_object <- mo
        }
      }
    }
    !is.null(self$match_object)
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- self$prob$eq
    func <- self$prob$func
    sol <- NULL
    if (identical(self$match_object$type, "2F1")) {
      sol <- .ode_get_sol_2F1_hypergeometric(eq, func, self$match_object)
      if (is.null(sol)) .ode_nie(paste0("The given ODE ", .str(eq), " cannot be solved by the hypergeometric method"))
    }
    list(sol)
  })

# ---------------------------------------------------------------------------
# The nth order linear solvers

.ode_cc_ok <- function(self, x) {
  order <- .op_order(self$prob)
  order > 0L && !is.null(self$r) && !any(vapply(.rk_ge0(self$r), function(i) .has1(self$r[[i]], x), TRUE))
}

.ode_defclass("NthLinearConstantCoeffHomogeneous", "SingleODESolver",
  hint = "nth_linear_constant_coeff_homogeneous", has_integral = FALSE,
  matches = function(self) {
    p <- self$prob
    self$r <- .op_linear_coefficients(.op_eq_high_order_free(p), p$func, .op_order(p))
    if (.ode_cc_ok(self, p$sym)) return(!.ode_true(.rk(self$r, -1L)))
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    order <- .op_order(self$prob)
    rc <- .ode_get_const_characteristic_eq_sols(self$r, fx, order)
    constants <- .op_constants(self$prob, length(rc[[1L]]))
    gsol <- Eq(fx, .add(lapply(seq_along(rc[[1L]]), function(i) .mul(list(constants[[i]], rc[[1L]][[i]])))))
    if (simplify_flag) gsol <- .ode_get_simplified_sol(list(gsol), fx, rc[[2L]])
    list(gsol)
  })

.ode_defclass("NthLinearConstantCoeffVariationOfParameters", "SingleODESolver",
  hint = "nth_linear_constant_coeff_variation_of_parameters", has_integral = TRUE,
  matches = function(self) {
    p <- self$prob
    self$r <- .op_linear_coefficients(.op_eq_high_order_free(p), p$func, .op_order(p))
    if (.ode_cc_ok(self, p$sym)) return(.ode_true(.rk(self$r, -1L)))
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- .op_eq_high_order_free(self$prob)
    fx <- self$prob$func
    order <- .op_order(self$prob)
    rc <- .ode_get_const_characteristic_eq_sols(self$r, fx, order)
    constants <- .op_constants(self$prob, length(rc[[1L]]))
    homogen_sol <- Eq(fx, .add(lapply(seq_along(rc[[1L]]), function(i) .mul(list(constants[[i]], rc[[1L]][[i]])))))
    homogen_sol <- .ode_solve_variation_of_parameters(eq, fx, rc[[1L]], homogen_sol, order, self$r, simplify_flag)
    if (simplify_flag) homogen_sol <- .ode_get_simplified_sol(list(homogen_sol), fx, rc[[2L]])
    list(homogen_sol)
  })

.ode_defclass("NthLinearConstantCoeffUndeterminedCoefficients", "SingleODESolver",
  hint = "nth_linear_constant_coeff_undetermined_coefficients", has_integral = FALSE,
  matches = function(self) {
    p <- self$prob
    eq <- .op_eq_high_order_free(p)
    x <- p$sym
    self$r <- .op_linear_coefficients(eq, p$func, .op_order(p))
    if (.ode_cc_ok(self, x) && .ode_true(.rk(self$r, -1L))) {
      eq_homogeneous <- .add(list(eq, .neg(.rk(self$r, -1L))))
      u <- .ode_undetermined_coefficients_match(.rk(self$r, -1L), x, p$func, eq_homogeneous)
      if (isTRUE(u$values[[1L]])) {
        self$trialset <- u$values[[2L]]
        return(TRUE)
      }
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- self$prob$eq
    fx <- self$prob$func
    order <- .op_order(self$prob)
    rc <- .ode_get_const_characteristic_eq_sols(self$r, fx, order)
    constants <- .op_constants(self$prob, length(rc[[1L]]))
    homogen_sol <- Eq(fx, .add(lapply(seq_along(rc[[1L]]), function(i) .mul(list(constants[[i]], rc[[1L]][[i]])))))
    self$r$list <- rc[[1L]]
    self$r$sol <- homogen_sol
    self$r$simpliy_flag <- simplify_flag
    gsol <- .ode_solve_undetermined_coefficients(eq, fx, order, self$r, self$trialset)
    if (simplify_flag) gsol <- .ode_get_simplified_sol(list(gsol), fx, rc[[2L]])
    list(gsol)
  })

# x**order/coeff times the linear coefficients, for the Euler solvers
.ode_euler_r <- function(self, eq) {
  p <- self$prob
  fx <- p$func
  x <- p$sym
  order <- .op_order(p)
  m <- .op_linear_coefficients(eq, fx, order)
  self$r <- NULL
  if (order > 0L && !is.null(m)) {
    factor <- .div(.pow(x, .int(order)), .rk(m, order))
    self$r <- lapply(m, function(v) .mul(list(factor, v)))
  }
  !is.null(self$r) && all(vapply(.rk_ge0(self$r), function(i) .ode_test_term(self$r[[i]], fx, as.integer(i)), TRUE))
}

.ode_defclass("NthLinearEulerEqHomogeneous", "SingleODESolver",
  hint = "nth_linear_euler_eq_homogeneous", has_integral = FALSE,
  matches = function(self) {
    if (.ode_euler_r(self, .op_eq_preprocessed(self$prob))) return(!.ode_true(.rk(self$r, -1L)))
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    list(.ode_get_euler_characteristic_eq_sols(self$prob$eq, self$prob$func, self$r)[[1L]])
  })

.ode_defclass("NthLinearEulerEqNonhomogeneousVariationOfParameters", "SingleODESolver",
  hint = "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters", has_integral = TRUE,
  matches = function(self) {
    if (.ode_euler_r(self, .op_eq_preprocessed(self$prob))) return(.ode_true(.rk(self$r, -1L)))
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- self$prob$eq
    fx <- self$prob$func
    order <- .op_order(self$prob)
    hr <- .ode_get_euler_characteristic_eq_sols(eq, fx, self$r)
    homogen_sol <- hr[[1L]]
    self$r[["-1"]] <- .div(.rk(self$r, -1L), .rk(self$r, order))
    sol <- .ode_solve_variation_of_parameters(eq, fx, hr[[2L]], homogen_sol, order, self$r, simplify_flag)
    list(Eq(fx, .add(list(homogen_sol$args[[2L]], .mul(list(.sub(sol$args[[2L]], homogen_sol$args[[2L]]), .rk(self$r, order)))))))
  })

.ode_defclass("NthLinearEulerEqNonhomogeneousUndeterminedCoefficients", "SingleODESolver",
  hint = "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients", has_integral = FALSE,
  matches = function(self) {
    x <- self$prob$sym
    if (.ode_euler_r(self, .op_eq_high_order_free(self$prob)) && .ode_true(.rk(self$r, -1L))) {
      pe <- posify(subs(.rk(self$r, -1L), x, fcall("exp", list(x))))
      u <- .ode_undetermined_coefficients_match(subs(pe$expr, pe$reps), x)
      return(isTRUE(u$values[[1L]]))
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    x <- self$prob$sym
    chareq <- S_Zero
    eq <- S_Zero
    symbol <- Dummy("x")
    for (i in .rk_ge0(self$r)) {
      ii <- as.integer(i)
      term <- .mul(list(self$r[[i]], .ode_d(.pow(x, symbol), x, ii), .pow(x, .neg(symbol))))
      chareq <- .add(list(chareq, expand(term)))
    }
    for (i in seq_len(degree(Poly(chareq, symbol)))) eq <- .add(list(eq, .mul(list(coeff(chareq, .pow(symbol, .int(i))), .ode_d(fx, x, i)))))
    c0 <- .sp_indep(chareq, list(symbol), as_add = TRUE)[[1L]]
    if (.ode_true(c0)) eq <- .add(list(eq, .mul(list(c0, fx))))
    pe <- posify(subs(.rk(self$r, -1L), x, fcall("exp", list(x))))
    eq <- .add(list(eq, subs(pe$expr, pe$reps)))
    inst <- .ode_solver_new("NthLinearConstantCoeffUndeterminedCoefficients", .ode_problem(eq, fx, x))
    self$const_undet_instance <- inst
    sol <- .ode_get_general_solution(inst, simplify = simplify_flag)[[1L]]
    sol <- subs(sol, x, fcall("log", list(x)))
    sol <- expand(subs(sol, fcall(fx$name, list(fcall("log", list(x)))), fx))
    list(sol)
  })

# ---------------------------------------------------------------------------
# Bessel and Airy equations

.ode_defclass("SecondLinearBessel", "SingleODESolver",
  hint = "2nd_linear_bessel", has_integral = FALSE,
  matches = function(self) {
    eq <- .op_eq_high_order_free(self$prob)
    f <- self$prob$func
    order <- .op_order(self$prob)
    x <- self$prob$sym
    df <- .ode_d(f, x)
    d2 <- .ode_d(f, x, 2L)
    a <- Wild("a", exclude = list(f, df))
    b <- Wild("b", exclude = list(x, f, df))
    a4 <- Wild("a4", exclude = list(x, f, df))
    b4 <- Wild("b4", exclude = list(x, f, df))
    c4 <- Wild("c4", exclude = list(x, f, df))
    d4 <- Wild("d4", exclude = list(x, f, df))
    a3 <- Wild("a3", exclude = list(f, df, d2))
    b3 <- Wild("b3", exclude = list(f, df, d2))
    c3 <- Wild("c3", exclude = list(f, df, d2))
    deq <- a3 * d2 + b3 * df + c3 * f
    r <- sym_match(collect(eq, list(d2, df, f)), deq)
    if (order == 2L && length(r)) {
      if (!all(vapply(r, function(v) is_polynomial(v), TRUE))) {
        eq <- expand(.as_numer_denom(eq)[[1L]])
        r <- sym_match(collect(eq, list(d2, df, f)), deq)
      }
    }
    get <- function(m, k) if (is.null(m[[k]])) S_Zero else m[[k]]
    if (length(r) && !is_zero_num(get(r, "a3"))) {
      cf <- sym_match(factor_sym(get(r, "a3")), a4 * .pow(.sub(x, b), b4))
      if (!is.null(cf)) {
        if (is_zero_num(get(cf, "b4"))) return(FALSE)
        point <- get(cf, "b")
      } else return(FALSE)
      ra3 <- get(r, "a3"); rb3 <- get(r, "b3"); rc3 <- get(r, "c3")
      if (.ode_true(point)) {
        sh <- .add(list(x, point))
        ra3 <- simplify(subs(ra3, x, sh))
        rb3 <- simplify(subs(rb3, x, sh))
        rc3 <- simplify(subs(rc3, x, sh))
      }
      sc <- .mul(list(get(cf, "a4"), .pow(x, .add(list(.int(-2L), get(cf, "b4"))))))
      ra3 <- cancel(.div(ra3, sc))
      rb3 <- cancel(.div(rb3, sc))
      rc3 <- cancel(.div(rc3, sc))
      coeff1 <- sym_match(factor_sym(rb3), a4 * x)
      if (is.null(coeff1)) return(FALSE)
      c2 <- sym_match(expand(rc3), a - b)
      if (is.null(c2)) return(FALSE)
      coeff2 <- sym_match(factor_sym(get(c2, "a")), .pow(c4, S_Two) * .pow(x, S_Two * a4))
      if (is.null(coeff2)) return(FALSE)
      if (is_zero_num(get(c2, "b"))) coeff2$d4 <- S_Zero
      else coeff2$d4 <- get(sym_match(factor_sym(get(c2, "b")), .pow(d4, S_Two)), "d4")
      self$rn <- list(n = coeff2$d4, a4 = get(coeff2, "c4"), d4 = get(coeff2, "a4"), c4 = get(coeff1, "a4"), b4 = point)
      return(TRUE)
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    x <- self$prob$sym
    rn <- self$rn
    n <- .sqrt(.add(list(.pow(rn$n, S_Two), .mul(list(.rat(1L, 4L), .pow(.sub(rn$c4, S_One), S_Two))))))
    cs <- .op_constants(self$prob, 2L)
    a4 <- rn$a4; d4 <- rn$d4
    arg <- .div(.mul(list(a4, .pow(x, d4))), d4)
    nu <- .div(n, d4)
    body <- .mul(list(.pow(x, .div(.sub(S_One, rn$c4), S_Two)),
                      .add(list(.mul(list(cs[[1L]], fcall("besselj", list(nu, arg)))), .mul(list(cs[[2L]], fcall("bessely", list(nu, arg))))))))
    list(Eq(fx, subs(body, x, .sub(x, rn$b4))))
  })

.ode_defclass("SecondLinearBesselTransform", "SingleODESolver",
  hint = "2nd_linear_bessel_transform", has_integral = FALSE, order = 2L,
  matches = function(self) {
    eq <- .op_eq_high_order_free(self$prob)
    f <- self$prob$func
    x <- self$prob$sym
    order <- .op_order(self$prob)
    if (order != 2L) return(FALSE)
    m <- .op_linear_coefficients(eq, f, order)
    if (is.null(m)) return(FALSE)
    a3 <- .rk(m, 2L); b3 <- .rk(m, 1L); c3 <- .rk(m, 0L); rhs <- .rk(m, -1L)
    if (!is_zero_num(b3) || is_zero_num(a3) || !is_zero_num(rhs)) return(FALSE)
    c3 <- .div(c3, a3)
    A_sym <- Wild("A_sym", exclude = list(x))
    k_sym <- Wild("k_sym", exclude = list(x))
    pm <- sym_match(c3, A_sym * .pow(x, k_sym))
    if (!length(pm) || is.null(pm$k_sym)) return(FALSE)
    if (identical(pm$k_sym$key, .int(-2L)$key)) return(FALSE)
    self$match_dict <- list(A = pm$A_sym %||% S_Zero, k = pm$k_sym)
    TRUE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    x <- self$prob$sym
    cs <- .op_constants(self$prob, 2L)
    A <- self$match_dict$A
    k <- self$match_dict$k
    nu <- .inv(.add(list(k, S_Two)))
    z <- .mul(list(.div(.mul(list(S_Two, .sqrt(Abs(A)))), .add(list(k, S_Two))), .pow(x, .div(.add(list(k, S_Two)), S_Two))))
    list(Eq(fx, .mul(list(.sqrt(x), .add(list(.mul(list(cs[[1L]], fcall("besselj", list(nu, z)))),
                                                .mul(list(cs[[2L]], fcall("bessely", list(nu, z))))))))))
  })

.ode_defclass("SecondLinearAiry", "SingleODESolver",
  hint = "2nd_linear_airy", has_integral = FALSE,
  matches = function(self) {
    eq <- .op_eq_high_order_free(self$prob)
    f <- self$prob$func
    order <- .op_order(self$prob)
    x <- self$prob$sym
    df <- .ode_d(f, x)
    a4 <- Wild("a4", exclude = list(x, f, df))
    b4 <- Wild("b4", exclude = list(x, f, df))
    m <- .op_linear_coefficients(eq, f, order)
    if (order == 2L && !is.null(m) && !is_zero_num(.rk(m, 2L))) {
      if (isTRUE(.sp_fact(.rk(m, 1L), "zero"))) {
        rn <- sym_match(cancel(.div(.rk(m, 0L), .rk(m, 2L))), a4 + b4 * x)
        if (length(rn) && !is.null(rn$b4) && !is_zero_num(rn$b4)) {
          self$rn <- list(b = rn$a4 %||% S_Zero, m = rn$b4)
          return(TRUE)
        }
      }
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    fx <- self$prob$func
    x <- self$prob$sym
    cs <- .op_constants(self$prob, 2L)
    b <- self$rn$b
    m <- self$rn$m
    if (isTRUE(.sp_fact(m, "positive"))) arg <- .sub(.neg(.div(b, .pow(cbrt(m), S_Two))), .mul(list(cbrt(m), x)))
    else arg <- .add(list(.neg(.div(b, .pow(cbrt(.neg(m)), S_Two))), .mul(list(cbrt(.neg(m)), x))))
    list(Eq(fx, .add(list(.mul(list(cs[[1L]], fcall("airyai", list(arg)))), .mul(list(cs[[2L]], fcall("airybi", list(arg))))))))
  })

# ---------------------------------------------------------------------------
# LieGroup

.ode_defclass("LieGroup", "SingleODESolver",
  hint = "lie_group", has_integral = FALSE,
  matches = function(self) {
    eq <- self$prob$eq
    fx <- self$prob$func
    order <- .op_order(self$prob)
    x <- self$prob$sym
    df <- .ode_d(fx, x)
    y <- Dummy("y")
    d <- Wild("d", exclude = list(df, .ode_d(fx, x, 2L)))
    e <- Wild("e", exclude = list(df))
    if (all(c("xi", "eta") %in% names(self$prob$params)) && order == 1L) {
      self$r3 <- list(xi = self$prob$params$xi, eta = self$prob$params$eta)
      r <- sym_match(collect(eq, df, exact = TRUE), d + e * df)
      if (length(r)) {
        self$r3$d <- subs(r$d %||% S_Zero, fx, y)
        self$r3$e <- subs(r$e %||% S_Zero, fx, y)
        self$r3$y <- y
      }
      return(TRUE)
    }
    FALSE
  },
  gen = function(self, simplify_flag = TRUE) {
    eq <- self$prob$eq
    x <- self$prob$sym
    func <- self$prob$func
    order <- .op_order(self$prob)
    df <- .ode_d(func, x)
    eqsol <- tryCatch(.sp_solve(eq, list(df), list()), NotImplementedError = function(e) list())
    desols <- list()
    for (s in eqsol) {
      sol <- .ode_lie_group(s, func, order, self$r3)
      if (length(sol)) desols <- c(desols, sol)
    }
    if (!length(desols)) .ode_nie(paste0("The given ODE ", .str(eq), " cannot be solved by the lie group method"))
    desols
  })

# single.solver_map, in SymPy's order
.ODE_SOLVER_MAP <- list(
  factorable = "Factorable",
  nth_linear_constant_coeff_homogeneous = "NthLinearConstantCoeffHomogeneous",
  nth_linear_euler_eq_homogeneous = "NthLinearEulerEqHomogeneous",
  nth_linear_constant_coeff_undetermined_coefficients = "NthLinearConstantCoeffUndeterminedCoefficients",
  nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients = "NthLinearEulerEqNonhomogeneousUndeterminedCoefficients",
  separable = "Separable",
  `1st_exact` = "FirstExact",
  `1st_linear` = "FirstLinear",
  Bernoulli = "Bernoulli",
  Riccati_special_minus2 = "RiccatiSpecial",
  `1st_rational_riccati` = "RationalRiccati",
  `1st_homogeneous_coeff_best` = "HomogeneousCoeffBest",
  `1st_homogeneous_coeff_subs_indep_div_dep` = "HomogeneousCoeffSubsIndepDivDep",
  `1st_homogeneous_coeff_subs_dep_div_indep` = "HomogeneousCoeffSubsDepDivIndep",
  almost_linear = "AlmostLinear",
  linear_coefficients = "LinearCoefficients",
  separable_reduced = "SeparableReduced",
  nth_linear_constant_coeff_variation_of_parameters = "NthLinearConstantCoeffVariationOfParameters",
  nth_linear_euler_eq_nonhomogeneous_variation_of_parameters = "NthLinearEulerEqNonhomogeneousVariationOfParameters",
  Liouville = "Liouville",
  `2nd_linear_airy` = "SecondLinearAiry",
  `2nd_linear_bessel` = "SecondLinearBessel",
  `2nd_linear_bessel_transform` = "SecondLinearBesselTransform",
  `2nd_hypergeometric` = "SecondHypergeometric",
  nth_order_reducible = "NthOrderReducible",
  `2nd_nonlinear_autonomous_conserved` = "SecondNonlinearAutonomousConserved",
  nth_algebraic = "NthAlgebraic",
  lie_group = "LieGroup")

# ---------------------------------------------------------------------------
# sympy.simplify.simplify.separatevars, kept here in SymPy's form

.ode_args_c <- function(e) if (e$type == "Mul") e$args else list(e)

.ode_separatevars <- function(expr, dict = FALSE, symbols = list(), force = FALSE) {
  expr <- as_sym(expr)
  if (dict) return(.ode_separatevars_dict(.ode_separatevars_(expr, force), symbols))
  .ode_separatevars_(expr, force)
}

.ode_separatevars_ <- function(expr, force) {
  if (expr$type == "Function" && expr$name == "Abs") {
    arg <- expr$args[[1L]]
    if (arg$type == "Mul" && !.expr_is_number(arg)) {
      s <- .ode_separatevars(arg, dict = TRUE, force = force)
      if (!is.null(s)) return(.mul(lapply(.sv_values(s), function(v) Abs(v))))
      return(expr)
    }
  }
  if (length(free_symbols(expr)) < 2L) return(expr)
  if (expr$type == "Mul") {
    args <- .sympy_args(expr)
    changed <- FALSE
    for (i in seq_along(args)) {
      a <- args[[i]]
      args[[i]] <- .ode_separatevars(a, force = force)
      changed <- changed || !identical(args[[i]]$key, a$key)
    }
    if (changed) expr <- .mul(args)
    return(expr)
  }
  if (expr$type == "Pow" && !is_const(expr$args[[1L]], "E")) expr <- .pow(.ode_separatevars(expr$args[[1L]], force = force), expr$args[[2L]])
  expr <- expand(expr, mul = FALSE, multinomial = FALSE, force = force)
  if (force) {
    pe <- posify(expr)
    expr <- subs(factor_sym(pe$expr), pe$reps)
  } else expr <- factor_sym(expr)
  if (expr$type != "Add") return(expr)
  args <- .sympy_args(expr)
  ckeys <- function(a) vapply(Filter(.is_comm, .ode_args_c(a)), .k, "")
  commonc_keys <- ckeys(args[[1L]])
  for (a in args[-1L]) commonc_keys <- intersect(commonc_keys, ckeys(a))
  allc <- Filter(.is_comm, .ode_args_c(args[[1L]]))
  commonc <- .mul(Filter(function(z) z$key %in% commonc_keys, allc))
  commonc <- .ode_as_coeff_Mul(commonc)[[2L]]
  commonc_set <- vapply(Filter(.is_comm, .ode_args_c(commonc)), .k, "")
  args <- lapply(args, function(a) {
    cs <- .ode_args_c(a)
    c <- Filter(function(z) .is_comm(z) && !(z$key %in% commonc_set), cs)
    nc <- Filter(Negate(.is_comm), cs)
    .mul(c(list(.mul(c)), list(.mul(nc))))
  })
  nonsepar <- .add(args)
  if (length(free_symbols(nonsepar)) > 1L) {
    if (force) {
      pe <- posify(nonsepar)
      e2 <- subs(factor_sym(pe$expr), pe$reps)
    } else e2 <- factor_sym(nonsepar)
    if (e2$type != "Add") nonsepar <- e2
  }
  .mul(list(commonc, nonsepar))
}

.ode_as_coeff_Mul <- function(e) {
  if (e$type %in% .NUMBER_TYPES) return(list(e, S_One))
  if (e$type == "Mul" && e$args[[1L]]$type %in% .NUMBER_TYPES) return(list(e$args[[1L]], .mul(e$args[-1L])))
  list(S_One, e)
}

# the separated factors: list(coeff, keys, values) with one value per symbol
.sv_new <- function(syms, vals, coeff) list(coeff = coeff, keys = vapply(syms, .k, ""), vals = vals)
.sv <- function(m, s) if (is.character(s)) m$coeff else m$vals[[match(s$key, m$keys)]]
.sv_values <- function(m) c(m$vals, list(m$coeff))

.ode_separatevars_dict <- function(expr, symbols) {
  if (is.null(symbols)) return(list(coeff = expr, keys = character(), vals = list()))
  if (!length(symbols)) {
    symbols <- free_symbols(expr)
    if (!length(symbols)) return(NULL)
  }
  keys <- vapply(symbols, .k, "")
  parts <- rep(list(list()), length(symbols))
  cf <- list()
  for (i in .ode_args_c(expr)) {
    fk <- vapply(free_symbols(i), .k, "")
    inter <- which(keys %in% fk)
    if (length(inter) > 1L) return(NULL)
    if (!length(inter)) cf[[length(cf) + 1L]] <- i
    else parts[[inter]][[length(parts[[inter]]) + 1L]] <- i
  }
  list(coeff = .mul(cf), keys = keys, vals = lapply(parts, .mul))
}
