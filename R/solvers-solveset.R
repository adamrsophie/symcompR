# Port of sympy/solvers/solveset.py: the univariate solver solveset(), its
# inverters and the equation class solvers. The systems part (linsolve,
# nonlinsolve, substitution) is in solvers-solveset-systems.R.

# ---------------------------------------------------------------------------
# Small helpers

.sset_err <- function(cls, msg) stop(.sympy_condition(cls, msg))
# a Python 2-tuple, as the inverters return
.sset_pair <- function(a, b) .pytuple(list(a, b))
.sset_is_fn <- function(e, nm) is_sym(e) && e$type == "Function" && e$name %in% nm
.SSET_TRIG <- c("sin", "cos", "tan", "cot", "sec", "csc")
.SSET_HYP <- c("sinh", "cosh", "tanh", "coth", "sech", "csch")
.sset_is_trig <- function(e) .sset_is_fn(e, .SSET_TRIG)
.sset_is_hyp <- function(e) .sset_is_fn(e, .SSET_HYP)
.sset_is_class <- function(e, cls) {
  switch(cls,
         TrigonometricFunction = .sset_is_trig(e),
         HyperbolicFunction = .sset_is_hyp(e),
         .sset_is_fn(e, cls))
}
# isinstance(f, exp) or (f.is_Pow and f.base == S.Exp1)
.sset_is_exp <- function(e) .sset_is_fn(e, "exp") || (e$type == "Pow" && identical(e$args[[1L]]$key, S_E$key))
.sset_exp_arg <- function(e) if (e$type == "Pow") e$args[[2L]] else e$args[[1L]]
.sset_set_is <- function(s, nm) is_sym(s) && s$type == "SpecialSet" && s$name == nm
.sset_is_empty <- function(s) is_sym(s) && s$type == "EmptySet"
.sset_is_fs <- function(s) is_sym(s) && s$type == "FiniteSet"
.sset_eq <- function(a, b) identical(a$key, b$key)
# Set.is_subset as a fuzzy value
.sset_subset <- function(a, b) {
  r <- tryCatch(.is_subset(a, b), error = function(e) NA)
  if (length(r) != 1L) NA else r
}
.sset_sub_reals <- function(d) isTRUE(.sset_subset(d, Reals))
# fuzzy assumption query
.sset_fact <- function(e, f) {
  v <- tryCatch(ask_fact(e, f), error = function(err) NA)
  if (length(v) != 1L) NA else v
}
.sset_lam <- function(v, e) Lambda(v, e)
.sset_img <- function(v, e, s) .imageset(Lambda(v, e), s)
# Add.make_args and Mul.make_args in SymPy order
.sset_add_args <- function(e) if (e$type == "Add") .sympy_args(e) else list(e)
.sset_mul_args <- function(e) if (e$type == "Mul") .sympy_args(e) else list(e)
.sset_union <- function(sets) if (!length(sets)) S_EmptySet else do.call(Union, sets)
.sset_has_type <- function(e, nm) .has_function(e, nm)
.sset_is_expr <- function(e) is_sym(e) && .is_expr_s(e)
.sset_free_has <- function(e, s) .has1(e, s)

# _masked(f, *atoms): replace the atoms by real dummies
.sset_masked <- function(f, fname) {
  ats <- Filter(function(a) .sset_is_fn(a, fname), atoms(f, "Function"))
  ats <- .sp_ordered(ats)
  mask <- list()
  k <- 0L
  for (a in ats) {
    for (i in mask) a <- xreplace(a, list(list(i[[1L]], i[[2L]])))
    mask[[length(mask) + 1L]] <- list(a, Dummy(paste0("a", k), real = TRUE))
    k <- k + 1L
  }
  for (i in seq_along(mask)) {
    f <- xreplace(f, list(list(mask[[i]][[1L]], mask[[i]][[2L]])))
    mask[[i]] <- list(mask[[i]][[2L]], mask[[i]][[1L]])
  }
  list(f, rev(mask))
}

# the iteration order of a FiniteSet: FiniteSet.__new__ orders its args by
# Set._infimum_key, then node count and default_sort_key
.sset_elems <- function(s) if (.sset_is_fs(s)) .set_print_order(s$args) else list()

# ---------------------------------------------------------------------------
# Inverters

# Function.inverse() for the classes SymPy defines it on
.SSET_INVERSE <- c(log = "exp", asin = "sin", acos = "cos", atan = "tan", acot = "cot", asec = "sec",
                   acsc = "csc", asinh = "sinh", acosh = "cosh", atanh = "tanh", acoth = "coth",
                   asech = "sech", acsch = "csch", erf = "erfinv", erfc = "erfcinv", erfinv = "erf",
                   erfcinv = "erfc", conjugate = "conjugate", tan = "atan", cot = "acot",
                   sinh = "asinh", cosh = "acosh", tanh = "atanh", coth = "acoth", sech = "asech",
                   csch = "acsch", exp = "log")
.sset_inverse <- function(f) {
  if (!is_sym(f) || f$type != "Function") return(NULL)
  r <- .SSET_INVERSE[f$name]
  if (is.na(r)) NULL else r[[1L]]
}

.sset_invert <- function(f_x, y, x, domain = Complexes) {
  x <- as_sym(x)
  if (!(x$type %in% c("Symbol", "Dummy"))) .sset_err("ValueError", "x must be a symbol")
  f_x <- as_sym(f_x)
  if (!.has1(f_x, x)) .sset_err("ValueError", "Inverse of constant function doesn't exist")
  y <- as_sym(y)
  if (.has1(y, x)) .sset_err("ValueError", "y should be independent of x ")
  r <- if (.sset_sub_reals(domain)) .sset_invert_real(f_x, FiniteSet(y), x)
       else .sset_invert_complex(f_x, FiniteSet(y), x)
  x1 <- r[[1L]]; s <- r[[2L]]
  if (!.sset_eq(x1, x)) return(r)
  if (.sset_set_is(domain, "Complexes")) return(r)
  if (.sset_is_fs(s)) return(.sset_pair(x1, Intersection(s, domain)))
  if (.sset_set_is(domain, "Reals")) return(r)
  .sset_pair(x1, Intersection(s, domain))
}

#' Invert real and complex valued functions
#'
#' `invert_real(f_x, y, x)` and `invert_complex(f_x, y, x)` reduce the
#' equation `f_x = y` to `g(x) = h` for the simplest `g(x)` they can reach,
#' and return the pair `(g(x), set of h)`, as SymPy's functions of the same
#' names. When `g(x)` is `x` the set holds the solutions.
#'
#' @param f_x An expression in `x`.
#' @param y The right hand side, free of `x`.
#' @param x The symbol to invert for.
#' @param domain The domain, `Complexes` by default for `invert_complex()`.
#' @return A tuple `(g(x), set)`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' invert_real(exp(x), y, x)
#' invert_complex(exp(x), y, x)
#' @export
invert_real <- function(f_x, y, x) .sset_invert(f_x, y, x, Reals)

#' @rdname invert_real
#' @export
invert_complex <- function(f_x, y, x, domain = Complexes) .sset_invert(f_x, y, x, as_sym(domain))

.sset_invert_real <- function(f, g_ys, symbol) {
  if (.sset_eq(f, symbol) || .sset_is_empty(g_ys)) return(.sset_pair(symbol, g_ys))
  n <- Dummy("n", real = TRUE)
  if (.sset_is_exp(f))
    return(.sset_invert_real(.sset_exp_arg(f), .sset_img(n, log(n), g_ys), symbol))
  inv <- .sset_inverse(f)
  if (!is.null(inv) && !.sset_is_trig(f) && !.sset_is_hyp(f)) {
    if (length(f$args) > 1L) .sset_err("ValueError", "Only functions with one argument are supported.")
    return(.sset_invert_real(f$args[[1L]], .sset_img(n, fcall(inv, list(n)), g_ys), symbol))
  }
  if (.sset_is_fn(f, "Abs")) return(.sset_invert_abs(f$args[[1L]], g_ys, symbol))
  if (f$type == "Add") {
    gh <- .sp_indep(f, list(symbol))
    g <- gh[[1L]]; h <- gh[[2L]]
    if (!is_zero_num(g)) return(.sset_invert_real(h, .sset_img(n, .sub(n, g), g_ys), symbol))
  }
  if (f$type == "Mul") {
    gh <- .sp_indep(f, list(symbol))
    g <- gh[[1L]]; h <- gh[[2L]]
    if (!is_one_num(g)) return(.sset_invert_real(h, .sset_img(n, .div(n, g), g_ys), symbol))
  }
  if (f$type == "Pow") {
    base <- f$args[[1L]]; expo <- f$args[[2L]]
    base_has_sym <- .has1(base, symbol)
    expo_has_sym <- .has1(expo, symbol)
    if (!expo_has_sym) {
      if (isTRUE(.sset_fact(expo, "rational"))) {
        nd <- .as_numer_denom(expo)
        num <- nd[[1L]]; den <- nd[[2L]]
        par <- function(e, r) {
          m <- Mod(e, .int(2))
          if (is_exact_number(m)) m$value == r else FALSE
        }
        if (par(den, 0) && par(num, 1) && identical(.sset_fact(den, "zero"), FALSE)) {
          root <- Lambda(n, .sset_real_root(n, expo))
          g_ys_pos <- Intersection(g_ys, Interval(0, S_Infinity))
          res <- .imageset(root, g_ys_pos)
          return(.sset_invert_real(base, res, symbol))
        }
        if (par(den, 1)) {
          root <- Lambda(n, .sset_real_root(n, expo))
          res <- .imageset(root, g_ys)
          if (par(num, 0)) {
            neg_res <- .sset_img(n, .neg(n), res)
            return(.sset_invert_real(base, Union(res, neg_res), symbol))
          }
          if (par(num, 1)) return(.sset_invert_real(base, res, symbol))
        }
      } else if (isTRUE(.sset_fact(expo, "irrational"))) {
        root <- Lambda(n, .sset_real_root(n, expo))
        g_ys_pos <- Intersection(g_ys, Interval(0, S_Infinity))
        res <- .imageset(root, g_ys_pos)
        return(.sset_invert_real(base, res, symbol))
      }
    }
    if (!base_has_sym) {
      rhs <- g_ys$args[[1L]]
      if (isTRUE(.sset_fact(base, "positive"))) {
        # the unevaluated log(n, base) evaluates at each element of a finite g_ys
        img <- if (g_ys$type == "FiniteSet") FiniteSet(lapply(g_ys$args, function(el) fcall("log", list(el, base))))
               else .sset_img(n, .sset_log_base(n, base), g_ys)
        return(.sset_invert_real(expo, img, symbol))
      } else if (isTRUE(.sset_fact(base, "negative"))) {
        il <- .sset_integer_log(rhs, base)
        if (il[[2L]]) return(.sset_invert_real(expo, FiniteSet(il[[1L]]), symbol))
        return(.sset_pair(expo, S_EmptySet))
      } else if (isTRUE(.sset_fact(base, "zero"))) {
        one <- Eq(rhs, 1)
        if (.is_true(one)) return(.sset_invert_real(expo, FiniteSet(S_Zero), symbol))
        if (.is_false(one)) return(.sset_pair(expo, S_EmptySet))
      }
    }
  }
  if (.sset_is_trig(f) || .sset_is_hyp(f)) return(.sset_invert_trig_hyp_real(f, g_ys, symbol))
  .sset_pair(f, g_ys)
}

# sympy.functions.elementary.miscellaneous.real_root(arg, n)
.sset_real_root <- function(arg, n) {
  Piecewise(list(.pow(arg, .inv(n)), Or(Eq(n, S_One), Eq(n, S_NegativeOne))),
            list(.mul(list(fcall("sign", list(arg)), .pow(fcall("Abs", list(arg)), .inv(n)))),
                 And(Eq(im(arg), S_Zero), Eq(Mod(n, .int(2)), S_One))),
            list(.pow(arg, .inv(n)), S_true))
}

# log(n, base, evaluate=False) for a symbol n: log(n)/log(base)
.sset_log_base <- function(n, base) {
  if (identical(base$key, S_E$key)) return(log(n))
  .div(log(n), log(base))
}

# integer_log(rhs, base) for a negative integer base: list(value, exact)
.sset_integer_log <- function(y, x) {
  r <- tryCatch(integer_log(y, x), error = function(err) NULL)
  if (is.null(r)) return(list(S_Zero, FALSE))
  list(as_sym(r[[1L]]), isTRUE(as.logical(if (is_sym(r[[2L]])) .is_true(r[[2L]]) else r[[2L]])))
}

.sset_trig_inverses <- function(fn) {
  pi_ <- S_Pi
  switch(fn,
    sin = list(list(function(y) fcall("asin", list(y)), function(y) .sub(pi_, fcall("asin", list(y)))),
               .mul(list(S_Two, pi_)), Interval(-1, 1)),
    cos = list(list(function(y) fcall("acos", list(y)), function(y) .neg(fcall("acos", list(y)))),
               .mul(list(S_Two, pi_)), Interval(-1, 1)),
    tan = list(list(function(y) fcall("atan", list(y))), pi_, Reals),
    cot = list(list(function(y) fcall("acot", list(y))), pi_, Reals),
    sec = list(list(function(y) fcall("asec", list(y)), function(y) .neg(fcall("asec", list(y)))),
               .mul(list(S_Two, pi_)), Union(Interval(S_NegativeInfinity, -1), Interval(1, S_Infinity))),
    csc = list(list(function(y) fcall("acsc", list(y)), function(y) .sub(pi_, fcall("acsc", list(y)))),
               .mul(list(S_Two, pi_)), Union(Interval(S_NegativeInfinity, -1), Interval(1, S_Infinity))))
}

.sset_invert_trig_hyp_real <- function(f, g_ys, symbol) {
  if (.sset_is_hyp(f)) {
    n <- Dummy("n", real = TRUE)
    a0 <- f$args[[1L]]
    dom_or <- function(g_ys_dom, isect = TRUE) {
      bad <- if (isect) g_ys_dom$type == "Intersection" else g_ys_dom$type == "Complement"
      if (bad) {
        if (.sset_is_fs(g_ys)) return(g_ys)
        return(NULL)
      }
      g_ys_dom
    }
    if (f$name == "sinh") return(.sset_invert_real(a0, .sset_img(n, fcall("asinh", list(n)), g_ys), symbol))
    if (f$name == "cosh") {
      d <- dom_or(Intersection(g_ys, Interval(1, S_Infinity)))
      if (is.null(d)) return(.sset_pair(f, g_ys))
      return(.sset_invert_real(a0, Union(.sset_img(n, fcall("acosh", list(n)), d),
                                         .sset_img(n, .neg(fcall("acosh", list(n))), d)), symbol))
    }
    if (f$name == "sech") {
      d <- dom_or(Intersection(g_ys, Interval(0, 1, TRUE, FALSE)))
      if (is.null(d)) return(.sset_pair(f, g_ys))
      return(.sset_invert_real(a0, Union(.sset_img(n, fcall("asech", list(n)), d),
                                         .sset_img(n, .neg(fcall("asech", list(n))), d)), symbol))
    }
    if (f$name == "tanh") {
      d <- dom_or(Intersection(g_ys, Interval(-1, 1, TRUE, TRUE)))
      if (is.null(d)) return(.sset_pair(f, g_ys))
      return(.sset_invert_real(a0, .sset_img(n, fcall("atanh", list(n)), d), symbol))
    }
    if (f$name == "coth") {
      d <- dom_or(Complement(g_ys, Interval(-1, 1)), FALSE)
      if (is.null(d)) return(.sset_pair(f, g_ys))
      return(.sset_invert_real(a0, .sset_img(n, fcall("acoth", list(n)), d), symbol))
    }
    if (f$name == "csch") {
      d <- dom_or(Complement(g_ys, FiniteSet(S_Zero)), FALSE)
      if (is.null(d)) return(.sset_pair(f, g_ys))
      return(.sset_invert_real(a0, .sset_img(n, fcall("acsch", list(n)), d), symbol))
    }
  } else if (.sset_is_trig(f) && .sset_is_fs(g_ys)) {
    ti <- .sset_trig_inverses(f$name)
    invs <- ti[[1L]]; period <- ti[[2L]]; rng <- ti[[3L]]
    n <- Dummy("n", integer = TRUE)
    create_return_set <- function(g) {
      invsimg <- .sset_union(lapply(invs, function(inv) .sset_img(n, .add(list(.mul(list(period, n)), inv(g))), Integers)))
      r <- .sset_invert_real(f$args[[1L]], invsimg, symbol)
      if (.sset_eq(r[[1L]], symbol)) {
        conds <- .contains(rng, g)
        return(ConditionSet(symbol, conds, r[[2L]]))
      }
      ConditionSet(symbol, Eq(f, g), Reals)
    }
    retset <- .sset_union(lapply(.sset_elems(g_ys), create_return_set))
    return(.sset_pair(symbol, retset))
  }
  .sset_pair(f, g_ys)
}

.sset_hyp_inverses <- function(fn) {
  ipi <- .mul(list(S_I, S_Pi)); i2pi <- .mul(list(S_Two, S_I, S_Pi))
  switch(fn,
    sinh = list(list(function(y) fcall("asinh", list(y)), function(y) .sub(ipi, fcall("asinh", list(y)))), i2pi, list()),
    cosh = list(list(function(y) fcall("acosh", list(y)), function(y) .neg(fcall("acosh", list(y)))), i2pi, list()),
    tanh = list(list(function(y) fcall("atanh", list(y))), ipi, list(S_NegativeOne, S_One)),
    coth = list(list(function(y) fcall("acoth", list(y))), ipi, list(S_NegativeOne, S_One)),
    sech = list(list(function(y) fcall("asech", list(y)), function(y) .neg(fcall("asech", list(y)))), i2pi, list(S_Zero)),
    csch = list(list(function(y) fcall("acsch", list(y)), function(y) .sub(ipi, fcall("acsch", list(y)))), i2pi, list(S_Zero)))
}

.sset_invert_trig_hyp_complex <- function(f, g_ys, symbol) {
  if (.sset_is_trig(f) && .sset_is_fs(g_ys)) {
    n <- Dummy("n", integer = TRUE)
    tn <- .mul(list(S_Two, n, S_Pi))
    Ls <- if (f$name %in% c("sin", "csc")) {
      F <- if (f$name == "sin") "asin" else "acsc"
      list(function(a) .add(list(tn, fcall(F, list(a)))),
           function(a) .sub(.add(list(tn, S_Pi)), fcall(F, list(a))))
    } else if (f$name %in% c("cos", "sec")) {
      F <- if (f$name == "cos") "acos" else "asec"
      list(function(a) .add(list(tn, fcall(F, list(a)))),
           function(a) .sub(tn, fcall(F, list(a))))
    } else {
      F <- .sset_inverse(f)
      list(function(a) .add(list(.mul(list(n, S_Pi)), fcall(F, list(a)))))
    }
    invs <- S_EmptySet
    for (L in Ls)
      invs <- Union(invs, .sset_union(lapply(.sset_elems(g_ys), function(g) .sset_img(n, L(g), Integers))))
    return(.sset_invert_complex(f$args[[1L]], invs, symbol))
  }
  if (.sset_is_hyp(f) && .sset_is_fs(g_ys)) {
    hi <- .sset_hyp_inverses(f$name)
    invs <- hi[[1L]]; period <- hi[[2L]]; excl <- hi[[3L]]
    n <- Dummy("n", integer = TRUE)
    create_return_set <- function(g) {
      invsimg <- .sset_union(lapply(invs, function(inv) .sset_img(n, .add(list(.mul(list(period, n)), inv(g))), Integers)))
      r <- .sset_invert_complex(f$args[[1L]], invsimg, symbol)
      if (.sset_eq(r[[1L]], symbol)) {
        conds <- if (length(excl)) do.call(And, lapply(excl, function(e) Ne(g, e))) else S_true
        return(ConditionSet(symbol, conds, r[[2L]]))
      }
      ConditionSet(symbol, Eq(f, g), Complexes)
    }
    retset <- .sset_union(lapply(.sset_elems(g_ys), create_return_set))
    return(.sset_pair(symbol, retset))
  }
  .sset_pair(f, g_ys)
}

.sset_invert_complex <- function(f, g_ys, symbol) {
  if (.sset_eq(f, symbol) || .sset_is_empty(g_ys)) return(.sset_pair(symbol, g_ys))
  n <- Dummy("n")
  if (f$type == "Add") {
    gh <- .sp_indep(f, list(symbol))
    g <- gh[[1L]]; h <- gh[[2L]]
    if (!is_zero_num(g)) return(.sset_invert_complex(h, .sset_img(n, .sub(n, g), g_ys), symbol))
  }
  if (f$type == "Mul") {
    gh <- .sp_indep(f, list(symbol))
    g <- gh[[1L]]; h <- gh[[2L]]
    if (!is_one_num(g)) {
      if (is_const(g, "-oo") || is_const(g, "zoo") || is_const(g, "oo")) return(.sset_pair(h, S_EmptySet))
      return(.sset_invert_complex(h, .sset_img(n, .div(n, g), g_ys), symbol))
    }
  }
  if (f$type == "Pow") {
    expo <- f$args[[2L]]
    if (is_exact_number(expo) && .sset_eq(g_ys, FiniteSet(S_Zero))) {
      if (isTRUE(.sset_fact(expo, "positive"))) return(.sset_invert_complex(f$args[[1L]], g_ys, symbol))
    }
  }
  inv <- .sset_inverse(f)
  if (!is.null(inv) && !.sset_is_trig(f) && !.sset_is_hyp(f) && !.sset_is_exp(f)) {
    if (length(f$args) > 1L) .sset_err("ValueError", "Only functions with one argument are supported.")
    return(.sset_invert_complex(f$args[[1L]], .sset_img(n, fcall(inv, list(n)), g_ys), symbol))
  }
  if (.sset_is_exp(f)) {
    if (is_sym(g_ys) && g_ys$type == "ImageSet") {
      lam <- g_ys$args[[1L]]
      g_ys_expr <- lam$args[[2L]]
      g_ys_vars <- .lambda_vars(lam$args[[1L]]$args)
      k <- Dummy(paste0("k", length(g_ys_vars)))
      vars1 <- c(list(k), g_ys_vars)
      body <- .add(list(.mul(list(S_I, .add(list(.mul(list(S_Two, k, S_Pi)), arg(g_ys_expr))))),
                        log(Abs(g_ys_expr))))
      lam1 <- new_sym("Lambda", list(Tuple(do.call(Tuple, vars1)), body))
      base <- do.call(ProductSet, rep(list(Integers), length(vars1)))
      exp_invs <- .imageset(lam1, base)
      return(.sset_invert_complex(.sset_exp_arg(f), exp_invs, symbol))
    } else if (.sset_is_fs(g_ys)) {
      sets <- lapply(Filter(function(g_y) !is_zero_num(g_y), .sset_elems(g_ys)), function(g_y)
        .sset_img(n, .add(list(.mul(list(S_I, .add(list(.mul(list(S_Two, n, S_Pi)), arg(g_y))))), log(Abs(g_y)))), Integers))
      return(.sset_invert_complex(.sset_exp_arg(f), .sset_union(sets), symbol))
    }
  }
  if (.sset_is_trig(f) || .sset_is_hyp(f)) return(.sset_invert_trig_hyp_complex(f, g_ys, symbol))
  .sset_pair(f, g_ys)
}

.sset_invert_abs <- function(f, g_ys, symbol) {
  if (!.sset_is_fs(g_ys)) {
    pos <- Intersection(g_ys, Interval(0, S_Infinity))
    parg <- .sset_invert_real(f, pos, symbol)
    narg <- .sset_invert_real(.neg(f), pos, symbol)
    if (!.sset_eq(parg[[1L]], narg[[1L]])) .sset_err("NotImplementedError", "")
    return(.sset_pair(parg[[1L]], Union(narg[[2L]], parg[[2L]])))
  }
  unknown <- list()
  for (a in g_ys$args) {
    ok <- if (is_number_type(a) || is_inf_const(a)) .sset_fact(a, "nonnegative") else .sset_fact(a, "positive")
    if (is.na(ok)) unknown[[length(unknown) + 1L]] <- a
    else if (!ok) return(.sset_pair(symbol, S_EmptySet))
  }
  conditions <- if (length(unknown)) do.call(And, lapply(unknown, function(i) Contains(i, Interval(0, S_Infinity)))) else S_true
  n <- Dummy("n", real = TRUE)
  r <- .sset_invert_real(f, Union(.sset_img(n, n, g_ys), .sset_img(n, .neg(n), g_ys)), symbol)
  .sset_pair(r[[1L]], ConditionSet(r[[1L]], conditions, r[[2L]]))
}

# ---------------------------------------------------------------------------
# Domain checks

#' Check that a point is in the domain of a function
#'
#' `domain_check(f, symbol, p)` returns `FALSE` when `f` or one of its
#' subexpressions becomes infinite at `symbol = p`, as SymPy's
#' `domain_check()`.
#'
#' @param f An expression.
#' @param symbol The variable.
#' @param p The point.
#' @return `TRUE` or `FALSE`.
#' @examples
#' x <- Symbol("x")
#' domain_check(1/x, x, 0)
#' domain_check(x^2, x, 3)
#' @export
domain_check <- function(f, symbol, p) {
  f <- as_sym(f); p <- as_sym(p); symbol <- as_sym(symbol)
  if (isTRUE(.sset_fact(p, "infinite"))) return(FALSE)
  .sset_domain_check(f, symbol, p)
}

.sset_domain_check <- function(f, symbol, p) {
  if (!length(f$args) && isTRUE(.sset_fact(f, "finite"))) return(TRUE)
  fs <- tryCatch(subs(f, symbol, p), error = function(e) NULL)
  if (!is.null(fs) && isTRUE(.sset_fact(fs, "infinite"))) return(FALSE)
  if (f$type == "Piecewise") {
    for (pc in .sps_pw_pairs(f)) {
      condsubs <- subs(pc[[2L]], symbol, p)
      if (.is_false(condsubs)) next
      if (.is_true(condsubs)) return(.sset_domain_check(pc[[1L]], symbol, p))
      return(TRUE)
    }
    # SymPy falls off the loop and returns None
    return(FALSE)
  }
  for (g in f$args) if (!.sset_domain_check(g, symbol, p)) return(FALSE)
  TRUE
}

.sset_is_finite_with_finite_vars <- function(f, domain = Complexes) {
  f <- as_sym(f)
  real <- .sset_sub_reals(domain)
  reps <- lapply(free_symbols(f), function(s) {
    A <- .sset_assumptions0(s)
    if (is.null(A[["finite"]])) A[["finite"]] <- TRUE
    if (real) { if (is.null(A[["real"]])) A[["real"]] <- TRUE }
    else if (is.null(A[["complex"]])) A[["complex"]] <- TRUE
    list(s, do.call(Dummy, c(list(NULL), as.list(unlist(A)))))
  })
  g <- if (length(reps)) xreplace(f, reps) else f
  .sset_fact(g, "finite")
}

# Symbol.assumptions0: the full assumptions of a symbol as a named list
.sset_assumptions0 <- function(s) {
  fx <- s$extra$facts
  if (is.null(fx)) return(list())
  fx <- fx[!is.na(fx)]
  as.list(fx)
}

# symbol._assumptions_orig == {name: TRUE}
.sset_orig_is <- function(s, nm) {
  if (s$type == "Symbol") {
    g <- s$extra$given
    return(length(g) == 1L && identical(names(g), nm) && isTRUE(g[[1L]]))
  }
  ref <- do.call(Symbol, stats::setNames(list("t", TRUE), c("name", nm)))$extra$facts
  fx <- s$extra$facts
  !is.null(fx) && identical(fx[order(names(fx))], ref[order(names(ref))])
}

.sset_is_function_class_equation <- function(func_class, f, symbol) {
  f <- as_sym(f); symbol <- as_sym(symbol)
  if (f$type %in% c("Mul", "Add"))
    return(all(vapply(.sympy_args(f), function(a) .sset_is_function_class_equation(func_class, a, symbol), NA)))
  if (f$type == "Pow") {
    if (!.has1(f$args[[2L]], symbol)) return(.sset_is_function_class_equation(func_class, f$args[[1L]], symbol))
    return(FALSE)
  }
  if (!.has1(f, symbol)) return(TRUE)
  if (.sset_is_class(f, func_class)) {
    a <- f$args[[1L]]
    if (!.sset_poly_ok(a, symbol)) return(FALSE)
    return(.sset_degree(a, symbol) <= 1)
  }
  FALSE
}

# Poly(e, symbol) succeeds: e is a polynomial in symbol
.sset_poly_ok <- function(e, symbol) isTRUE(tryCatch(is_polynomial(e, symbol), error = function(err) FALSE))
.sset_degree <- function(e, symbol) {
  d <- degree(e, symbol)
  if (is_sym(d)) { if (is_const(d, "-oo")) -Inf else as.numeric(as.character(d$value)) } else as.numeric(d)
}
# Poly(e, symbol).all_coeffs()[0]
.sset_lead_coeff <- function(e, symbol) {
  p <- .sp_as_poly(e, list(symbol))
  if (is.null(p) || !length(p$terms)) return(S_Zero)
  degs <- vapply(p$terms, function(t) t$monom[[1L]], 0L)
  .add(lapply(p$terms[degs == max(degs)], `[[`, "coeff"))
}

# ---------------------------------------------------------------------------
# Equation classes

.sset_solve_as_rational <- function(f, symbol, domain) {
  f <- together(.sp_mexpand_rec(f), deep = TRUE)
  gh <- fraction(f)
  g <- gh[[1L]]; h <- gh[[2L]]
  if (!.has1(h, symbol)) {
    r <- tryCatch(.sset_solve_as_poly(g, symbol, domain),
                  NotImplementedError = function(e) ConditionSet(symbol, Eq(f, 0), domain),
                  CoercionFailed = function(e) S_EmptySet)
    return(r)
  }
  valid <- .sset_solveset(g, symbol, domain)
  invalid <- .sset_solveset(h, symbol, domain)
  Complement(valid, invalid)
}

.sset_solve_trig <- function(f, symbol, domain) {
  trig_expr <- NULL; count <- 0L
  for (ex in preorder_traversal(f)) {
    if ((.sset_is_trig(ex) || .sset_is_hyp(ex)) && .has1(ex, symbol)) {
      if (is.null(trig_expr)) { trig_expr <- ex; count <- 1L }
      else if (isFALSE(trig_expr)) next
      else if (.sset_eq(ex, trig_expr)) count <- count + 1L
      else { trig_expr <- FALSE; count <- 0L; break }
    }
  }
  if (count == 1L) {
    r <- .sset_invert(f, S_Zero, symbol, domain)
    if (.sset_eq(r[[1L]], symbol)) {
      cond <- S_true
      if (length(setdiff(vapply(free_symbols(trig_expr), .k, ""), symbol$key))) {
        ah <- .sp_indep(trig_expr$args[[1L]], list(symbol), as_add = TRUE)
        mh <- .sp_indep(ah[[2L]], list(symbol), as_add = FALSE)
        nd <- .as_numer_denom(mh[[1L]])
        cond <- And(Ne(nd[[1L]], 0), Ne(nd[[2L]], 0))
      }
      return(ConditionSet(symbol, cond, r[[2L]]))
    }
    return(ConditionSet(symbol, Eq(f, 0), domain))
  } else if (count > 0L) {
    y <- Dummy("y")
    f_cov <- subs(f, trig_expr, y)
    sol_cov <- solveset(f_cov, y, domain)
    if (.sset_is_fs(sol_cov))
      return(.sset_union(lapply(.sset_elems(sol_cov), function(s) .sset_solve_trig(.sub(trig_expr, s), symbol, domain))))
  }
  sol <- tryCatch(.sset_solve_trig1(f, symbol, domain), SolveTrig1Error = function(e) NULL)
  if (is.null(sol)) {
    sol <- tryCatch(.sset_solve_trig2(f, symbol, domain), ValueError = function(e) NULL)
    if (is.null(sol))
      .sset_err("NotImplementedError", .filldedent("Solution to this kind of trigonometric equations is yet to be implemented"))
  }
  sol
}

.sset_trig1_fail <- function(msg) .sset_err("SolveTrig1Error", msg)

.sset_lcm_list <- function(xs) {
  if (!length(xs)) return(S_One)
  Reduce(function(a, b) lcm(a, b), xs)
}
.sset_gcd_list <- function(xs) {
  if (!length(xs)) return(S_Zero)
  Reduce(function(a, b) gcd(a, b), xs)
}

.sset_solve_trig1 <- function(f, symbol, domain) {
  x <- Dummy("x")
  if (.sset_is_function_class_equation("HyperbolicFunction", f, symbol)) {
    cov <- exp(x)
    inverter <- if (.sset_sub_reals(domain)) function(a, b, c) .sset_invert(a, b, c, Reals)
                else function(a, b, c) .sset_invert(a, b, c, Complexes)
  } else {
    cov <- exp(.mul(list(S_I, x)))
    inverter <- function(a, b, c) .sset_invert(a, b, c, Complexes)
  }
  f <- trigsimp(f)
  f_original <- f
  trig_functions <- Filter(function(a) .sset_is_trig(a) || .sset_is_hyp(a), atoms(f, "Function"))
  trig_arguments <- lapply(trig_functions, function(e) e$args[[1L]])
  if (!any(vapply(trig_arguments, function(a) .has1(a, symbol), NA)))
    return(solveset(f_original, symbol, domain))
  denominators <- list(); numerators <- list()
  for (ar in trig_arguments) {
    if (!.sset_poly_ok(ar, symbol)) .sset_trig1_fail("trig argument is not a polynomial")
    dg <- .sset_degree(ar, symbol)
    if (dg > 1) .sset_trig1_fail("degree of variable must not exceed one")
    if (dg == 0) next
    c <- .sset_lead_coeff(ar, symbol)
    fr <- fraction(c)
    numerators[[length(numerators) + 1L]] <- fr[[1L]]
    denominators[[length(denominators) + 1L]] <- fr[[2L]]
  }
  mu <- .div(.sset_lcm_list(denominators), .sset_gcd_list(numerators))
  f <- subs(f, symbol, .mul(list(mu, x)))
  f <- rewrite(f, "exp")
  f <- together(f)
  gh <- fraction(f)
  y <- Dummy("y")
  g <- expand(gh[[1L]]); h <- expand(gh[[2L]])
  g <- subs(g, cov, y); h <- subs(h, cov, y)
  if (.has1(g, x) || .has1(h, x)) .sset_trig1_fail("change of variable not possible")
  solns <- Complement(solveset(g, y, Complexes), solveset(h, y, Complexes))
  if (solns$type == "ConditionSet") .sset_trig1_fail("polynomial has ConditionSet solution")
  if (.sset_is_fs(solns)) {
    if (any(vapply(solns$args, function(s) s$type == "CRootOf", NA)))
      .sset_trig1_fail("polynomial results in RootOf object")
    cov <- subs(cov, x, .div(symbol, mu))
    result <- .sset_union(lapply(.sset_elems(solns), function(s) inverter(cov, s, symbol)[[2L]]))
    if (length(Filter(function(a) a$type %in% c("Symbol", "Dummy"), atoms(mu)))) {
      syms <- Filter(function(a) a$type %in% c("Symbol", "Dummy"), atoms(mu))
      fr <- fraction(mu)
      condnum <- .sp_indep(fr[[1L]], syms, as_add = FALSE)[[2L]]
      condden <- .sp_indep(fr[[2L]], syms, as_add = FALSE)[[2L]]
      cond <- And(Ne(condnum, 0), Ne(condden, 0))
    } else cond <- S_true
    if (.sset_set_is(domain, "Complexes")) return(ConditionSet(symbol, cond, result))
    return(ConditionSet(symbol, cond, Intersection(result, domain)))
  }
  if (.sset_is_empty(solns)) return(S_EmptySet)
  .sset_trig1_fail("polynomial solutions must form FiniteSet")
}

.sset_solve_trig2 <- function(f, symbol, domain) {
  f <- trigsimp(f)
  f_original <- f
  trig_functions <- Filter(function(a) .sset_is_trig(a), atoms(f, "Function"))
  trig_arguments <- lapply(trig_functions, function(e) e$args[[1L]])
  denominators <- list(); numerators <- list()
  if (!length(trig_functions)) return(ConditionSet(symbol, Eq(f_original, 0), domain))
  for (ar in trig_arguments) {
    if (!.sset_poly_ok(ar, symbol)) .sset_err("ValueError", "give up, we cannot solve if this is not a polynomial in x")
    dg <- .sset_degree(ar, symbol)
    if (dg > 1) .sset_err("ValueError", "degree of variable inside polynomial should not exceed one")
    if (dg == 0) next
    c <- .sset_lead_coeff(ar, symbol)
    if (!is_exact_number(c)) return(ConditionSet(symbol, Eq(f_original, 0), domain))
    q <- gmp::as.bigq(c$value)
    numerators[[length(numerators) + 1L]] <- gmp::numerator(q)
    denominators[[length(denominators) + 1L]] <- gmp::denominator(q)
  }
  x <- Dummy("x")
  L <- Reduce(function(a, b) a * b / gmp::gcd(a, b), denominators, gmp::as.bigz(1))
  G <- Reduce(function(a, b) gmp::gcd(a, b), numerators, gmp::as.bigz(0))
  mu <- make_rational(gmp::as.bigq(2L * L, G))
  f <- subs(f, symbol, .mul(list(mu, x)))
  f <- rewrite(f, "tan")
  f <- expand_trig(f)
  f <- together(f)
  gh <- fraction(f)
  y <- Dummy("y")
  g <- expand(gh[[1L]]); h <- expand(gh[[2L]])
  tx <- fcall("tan", list(x))
  g <- subs(g, tx, y); h <- subs(h, tx, y)
  if (.has1(g, x) || .has1(h, x)) return(ConditionSet(symbol, Eq(f_original, 0), domain))
  solns <- Complement(solveset(g, y, Reals), solveset(h, y, Reals))
  if (.sset_is_fs(solns)) {
    tsm <- fcall("tan", list(.div(symbol, mu)))
    result <- .sset_union(lapply(.sset_elems(solns), function(s) .sset_invert(tsm, s, symbol, Reals)[[2L]]))
    dsol <- .sset_invert(tsm, S_Infinity, symbol, Reals)[[2L]]
    if (.sset_degree(h, y) > .sset_degree(g, y)) result <- Union(result, dsol)
    return(Intersection(result, domain))
  }
  if (.sset_is_empty(solns)) return(S_EmptySet)
  ConditionSet(symbol, Eq(f_original, 0), Reals)
}

# roots(f, symbol, cubics=True, quartics=True, quintics=True, domain='EX')
# as a list of (root, multiplicity) in dict order
.sset_roots <- function(f, symbol) {
  r <- roots(f, symbol, cubics = TRUE, quartics = TRUE, quintics = TRUE)
  if (is.null(r) || !length(r)) return(list(vals = list(), mult = integer()))
  m <- attr(r, "multiplicity")
  if (is.null(m)) m <- rep(1L, length(r))
  list(vals = unclass(r)[seq_along(r)], mult = m)
}

.sset_expand_complex_set <- function(result) {
  if (!.sset_is_fs(result)) return(result)
  ok <- all(vapply(result$args, function(s)
    !length(free_symbols(s)) && !.has_undef(s) && s$type != "CRootOf", NA))
  if (ok) {
    s <- Dummy("s")
    result <- .imageset(Lambda(s, expand_complex(s)), result)
  }
  result
}

.sset_solve_as_poly <- function(f, symbol, domain = Complexes) {
  result <- NULL
  if (.sset_poly_ok(f, symbol)) {
    rs <- .sset_roots(f, symbol)
    num_roots <- sum(rs$mult)
    if (.sset_degree(f, symbol) <= num_roots) {
      result <- FiniteSet(rs$vals)
    } else {
      solns <- tryCatch(all_roots(f, symbol), error = function(e) NULL)
      if (!is.null(solns) && .sset_degree(f, symbol) <= length(solns)) result <- FiniteSet(solns)
      else result <- ConditionSet(symbol, Eq(f, 0), domain)
    }
  } else {
    P <- tryCatch(Poly(f), error = function(e) NULL)
    if (is.null(P)) return(ConditionSet(symbol, Eq(f, 0), domain))
    gens <- Filter(function(g) .has1(g, symbol), P$gens[seq_len(.poly_ngens(P))])
    if (length(gens) == 1L) {
      gen <- gens[[1L]]
      pe <- .sset_gen_poly(f, gen)
      if (is.null(pe)) return(ConditionSet(symbol, Eq(f, 0), domain))
      deg <- .sset_degree(pe$expr, pe$dummy)
      rs <- .sset_roots(pe$expr, pe$dummy)
      poly_solns <- FiniteSet(rs$vals)
      if (length(.sset_elems(poly_solns)) < deg) result <- ConditionSet(symbol, Eq(f, 0), domain)
      gen <- pe$gen
      if (!.sset_eq(gen, symbol)) {
        y <- Dummy("y")
        inv <- if (.sset_sub_reals(domain)) .sset_invert(gen, y, symbol, Reals) else .sset_invert(gen, y, symbol, Complexes)
        if (.sset_eq(inv[[1L]], symbol)) {
          result <- .sset_union(lapply(.sset_elems(poly_solns), function(s) subs(inv[[2L]], y, s)))
          if (.sset_is_fs(result) && gen$type == "Pow" && is_exact_number(gen$args[[1L]]))
            result <- FiniteSet(lapply(result$args, expand_log))
        } else result <- ConditionSet(symbol, Eq(f, 0), domain)
      }
    } else result <- ConditionSet(symbol, Eq(f, 0), domain)
  }
  if (!is.null(result)) {
    result <- .sset_expand_complex_set(result)
    if (.sset_is_fs(result) && !.sset_set_is(domain, "Complexes")) result <- Intersection(result, domain)
    return(result)
  }
  ConditionSet(symbol, Eq(f, 0), domain)
}

# f as a polynomial in the generator gen: the expression with gen replaced
# by a dummy, or NULL
.sset_gen_poly <- function(f, gen) {
  d <- Dummy("g")
  e <- .sset_gen_subs(expand(f), gen, d)
  if (.has1(e, gen) || length(Filter(function(s) !.sset_eq(s, d) && .has1(gen, s), free_symbols(e)))) return(NULL)
  if (!.sset_poly_ok(e, d)) return(NULL)
  list(expr = e, dummy = d, gen = gen)
}

# replace every power of gen by a power of d
.sset_gen_subs <- function(e, gen, d) {
  be <- .sset_base_exp(gen)
  walk <- function(z) {
    if (.sset_eq(z, gen)) return(d)
    zb <- .sset_base_exp(z)
    if (.sset_eq(zb[[1L]], be[[1L]])) {
      r <- .div(zb[[2L]], be[[2L]])
      if (r$type == "Integer" && r$value > 0) return(.pow(d, r))
    }
    if (!length(z$args)) return(z)
    sym_rebuild(z, lapply(z$args, walk))
  }
  walk(e)
}
.sset_base_exp <- function(e) {
  if (e$type == "Pow") return(e$args)
  if (.sset_is_fn(e, "exp")) {
    a <- e$args[[1L]]
    cm <- .sset_coeff_mul(a)
    return(list(fcall("exp", list(cm[[2L]])), cm[[1L]]))
  }
  list(e, S_One)
}
# as_coeff_Mul(rational=True)
.sset_coeff_mul <- function(e) {
  if (is_exact_number(e)) return(list(e, S_One))
  if (e$type == "Mul") {
    a <- .sympy_args(e)
    if (is_exact_number(a[[1L]])) return(list(a[[1L]], .mul(a[-1L])))
  }
  list(S_One, e)
}

.sset_solve_radical <- function(f, unradf, symbol, solveset_solver) {
  if (is.null(unradf)) { eq <- f; cov <- list() } else { eq <- unradf[[1L]]; cov <- unradf[[2L]] }
  if (!length(cov)) {
    result <- Complement(solveset_solver(eq, symbol),
                         .sset_union(lapply(denoms(f, symbol), function(g) solveset_solver(g, symbol))))
  } else {
    y <- cov[[1L]]; yeq <- cov[[2L]]
    if (.sset_is_empty(solveset_solver(.sub(y, S_I), y))) {
      yreal <- Dummy("yreal", real = TRUE)
      yeq <- xreplace(yeq, list(list(y, yreal)))
      eq <- xreplace(eq, list(list(y, yreal)))
      y <- yreal
    }
    g_y_s <- solveset_solver(yeq, symbol)
    f_y_sols <- solveset_solver(eq, y)
    result <- .sset_union(lapply(.sset_elems(g_y_s), function(g_y) .imageset(Lambda(y, g_y), f_y_sols)))
  }
  check_finiteset <- function(solutions) {
    f_set <- list(); c_set <- list()
    for (s in .sset_elems(solutions)) {
      ok <- .sp_checksol(f, .sd1(symbol, s), list())
      if (isTRUE(ok)) f_set[[length(f_set) + 1L]] <- s else c_set[[length(c_set) + 1L]] <- s
    }
    Union(FiniteSet(f_set), ConditionSet(symbol, Eq(f, 0), FiniteSet(c_set)))
  }
  check_set <- function(solutions) {
    if (.sset_is_empty(solutions)) return(solutions)
    if (solutions$type == "ConditionSet") return(solutions)
    if (.sset_is_fs(solutions)) return(check_finiteset(solutions))
    if (solutions$type == "Complement") return(Complement(check_set(solutions$args[[1L]]), solutions$args[[2L]]))
    if (solutions$type == "Union") return(.sset_union(lapply(solutions$args, check_set)))
    solutions
  }
  check_set(result)
}

.sset_solve_abs <- function(f, symbol, domain) {
  if (!.sset_sub_reals(domain))
    .sset_err("ValueError", .filldedent("Absolute values cannot be inverted in the complex domain."))
  pm <- .sset_match_abs(f, symbol)
  f_p <- pm[[1L]]; f_q <- pm[[2L]]; f_r <- pm[[3L]]
  if (!(isTRUE(.sset_fact(f_p, "zero")) || isTRUE(.sset_fact(f_q, "zero")))) {
    domain <- continuous_domain(f_q, symbol, domain)
    q_pos_cond <- solve_univariate_inequality(Ge(f_q, 0), symbol, relational = FALSE, domain = domain)
    q_neg_cond <- Complement(domain, q_pos_cond)
    sols_q_pos <- Intersection(solveset(.add(list(.mul(list(f_p, f_q)), f_r)), symbol, Reals), q_pos_cond)
    sols_q_neg <- Intersection(solveset(.add(list(.mul(list(f_p, .neg(f_q))), f_r)), symbol, Reals), q_neg_cond)
    return(Union(sols_q_pos, sols_q_neg))
  }
  ConditionSet(symbol, Eq(f, 0), domain)
}

# f.match(p*Abs(q) + r): list(p, q, r), zeros when there is no match. As
# AssocOp._matches_commutative, the terms are tried last first in ordered()
# order, stably sorted by a leading numeric coefficient.
.sset_match_abs <- function(f, symbol) {
  terms <- list(f)
  if (f$type == "Add") {
    terms <- .sp_ordered(.sympy_args(f))
    key <- vapply(terms, function(t)
      if (t$type == "Mul" && is_number_type(t$args[[1L]])) as.numeric(N(t$args[[1L]])) else 0, 0)
    terms <- terms[order(key, method = "radix")]
  }
  for (i in rev(seq_along(terms))) {
    pq <- .sset_match_pabs(terms[[i]])
    if (!is.null(pq)) return(list(pq[[1L]], pq[[2L]], .add(terms[-i])))
  }
  list(S_Zero, S_Zero, S_Zero)
}

# (p*Abs(q)).matches(t): list(p, q) or NULL
.sset_match_pabs <- function(t) {
  if (.sset_is_fn(t, "Abs")) return(list(S_One, t$args[[1L]]))
  if (t$type == "Mul") {
    fac <- .sympy_args(t)
    ai <- which(vapply(fac, function(m) .sset_is_fn(m, "Abs"), NA))
    if (length(ai)) {
      k <- ai[length(ai)]
      return(list(.mul(fac[-k]), fac[[k]]$args[[1L]]))
    }
  }
  if (t$type == "Pow" && .sset_is_fn(t$args[[1L]], "Abs") && t$args[[2L]]$type == "Integer" && t$args[[2L]]$value > 1)
    return(list(.pow(t$args[[1L]], .int(t$args[[2L]]$value - 1L)), t$args[[1L]]$args[[1L]]))
  NULL
}

# ---------------------------------------------------------------------------
# _solveset

.sset_solveset <- function(f, symbol, domain, check = FALSE) {
  if (.is_true(f)) return(domain)
  orig_f <- f
  if (f$type == "Mul") {
    cf <- .sp_indep(f, list(symbol), as_add = FALSE)
    coeff <- cf[[1L]]; f <- cf[[2L]]
    if (is_const(coeff, "zoo") || is_const(coeff, "-oo") || is_const(coeff, "oo")) f <- together(orig_f)
  } else if (f$type == "Add") {
    ah <- .sp_indep(f, list(symbol))
    a <- ah[[1L]]; h <- ah[[2L]]
    mh <- .sp_indep(h, list(symbol), as_add = FALSE)
    m <- mh[[1L]]; h <- mh[[2L]]
    if (!(is_const(m, "zoo") || is_zero_num(m) || is_const(m, "oo") || is_const(m, "-oo")))
      f <- .add(list(.div(a, m), h))
  }
  solver <- function(f, x, dom = domain) .sset_solveset(f, x, dom)
  inverter <- function(f, rhs, symbol) .sset_invert(f, rhs, symbol, domain)
  result <- S_EmptySet
  expanded <- expand(f)
  if (is_relational(f) && !.has1(expanded, symbol)) {
    expanded <- simplify(expanded)
    if (.is_true(expanded)) return(domain)
    if (.is_false(expanded)) return(S_EmptySet)
    return(ConditionSet(symbol, expanded, domain))
  }
  if (!is_relational(f) && isTRUE(.sset_fact(expanded, "zero"))) return(domain)
  if (!.has1(expanded, symbol)) return(S_EmptySet)
  if (f$type == "Mul" && all(vapply(.sympy_args(f), function(m) isTRUE(.sset_is_finite_with_finite_vars(m, domain)), NA))) {
    result <- .sset_union(lapply(.sympy_args(f), function(m) solver(m, symbol)))
  } else if (!is_relational(f) && (.sset_is_function_class_equation("TrigonometricFunction", f, symbol) ||
             .sset_is_function_class_equation("HyperbolicFunction", f, symbol))) {
    result <- .sset_solve_trig(f, symbol, domain)
  } else if (.sset_is_fn(f, "arg")) {
    a <- f$args[[1L]]
    result <- Intersection(.sset_solveset(Gt(re(a), 0), symbol, domain), .sset_solveset(im(a), symbol, domain))
  } else if (f$type == "Piecewise") {
    for (pr in .sset_expr_set_pairs(f, domain, symbol)) {
      in_set <- pr[[2L]]
      if (is_relational(in_set)) in_set <- .sset_as_set(in_set)
      result <- Union(result, solver(pr[[1L]], symbol, in_set))
    }
  } else if (f$type == "Equality") {
    result <- solver(.sub(f$args[[1L]], f$args[[2L]]), symbol, domain)
  } else if (is_relational(f)) {
    result <- tryCatch(solve_univariate_inequality(f, symbol, domain = domain, relational = FALSE),
                       NotImplementedError = function(e) ConditionSet(symbol, f, domain))
    return(result)
  } else if (.sset_is_modular(f, symbol)) {
    result <- .sset_solve_modular(f, symbol, domain)
  } else {
    inv <- inverter(f, S_Zero, symbol)
    lhs <- inv[[1L]]; rhs_s <- inv[[2L]]
    if (.sset_eq(lhs, symbol)) {
      if (.sset_is_fs(rhs_s))
        rhs_s <- FiniteSet(lapply(rhs_s$args, function(i) {
          cp <- as_content_primitive(signsimp(i))
          .mul(list(cp[[1L]], cp[[2L]]))
        }))
      result <- rhs_s
    } else if (.sset_is_fs(rhs_s)) {
      for (rhs in .sset_elems(rhs_s)) {
        equation <- .sub(lhs, rhs)
        if (.sset_eq(equation, f)) {
          u <- .sp_unrad(f, list(symbol))
          if (!is.null(u)) {
            result <- Union(result, .sset_solve_radical(equation, u, symbol, solver))
          } else if (.sset_has_type(equation, "Abs")) {
            result <- Union(result, .sset_solve_abs(f, symbol, domain))
          } else {
            result_rational <- .sset_solve_as_rational(equation, symbol, domain)
            if (result_rational$type != "ConditionSet") {
              result <- Union(result, result_rational)
            } else {
              t_result <- .sset_transolve(equation, symbol, domain)
              if (t_result$type == "ConditionSet") {
                factored <- factor(equation)
                if (factored$type == "Mul" && !.sset_eq(equation, factored)) {
                  dep <- .sp_indep(factored, list(symbol))[[2L]]
                  if (dep$type != "Add") {
                    t_results <- list()
                    for (fac in .sset_mul_args(factored))
                      if (.has1(fac, symbol)) t_results[[length(t_results) + 1L]] <- solver(fac, symbol)
                    t_result <- .sset_union(t_results)
                  }
                }
              }
              result <- Union(result, t_result)
            }
          }
        } else {
          result <- Union(result, solver(equation, symbol))
        }
      }
    } else if (!.sset_is_empty(rhs_s)) {
      result <- ConditionSet(symbol, Eq(f, 0), domain)
    }
  }
  if (result$type == "ConditionSet" && .sset_is_expr(f)) {
    nd <- .as_numer_denom(f)
    if (.has1(nd[[2L]], symbol)) {
      r2 <- .sset_solveset(nd[[1L]], symbol, domain)
      if (r2$type != "ConditionSet") {
        singularities <- .sset_solveset(nd[[2L]], symbol, domain)
        result <- Complement(r2, singularities)
      }
    }
  }
  if (check) {
    if (result$type == "ConditionSet") return(result)
    if (.sset_is_expr(orig_f)) {
      fx <- .sp_indep(orig_f, list(symbol), as_add = TRUE)[[2L]]
      fx <- .sp_indep(fx, list(symbol), as_add = FALSE)[[2L]]
    } else fx <- orig_f
    if (.sset_is_fs(result))
      result <- FiniteSet(Filter(function(s) s$type == "CRootOf" || domain_check(fx, symbol, s), result$args))
  }
  result
}

# Piecewise.as_expr_set_pairs(domain)
.sset_expr_set_pairs <- function(f, domain, symbol) {
  exp_sets <- list()
  U <- domain
  cplx <- !.sset_sub_reals(domain)
  cond_free <- character()
  for (pc in .sps_pw_pairs(f)) {
    expr <- pc[[1L]]; cond <- pc[[2L]]
    cond_free <- union(cond_free, vapply(free_symbols(cond), .k, ""))
    if (length(cond_free) > 1L)
      .sset_err("NotImplementedError", .filldedent("multivariate conditions are not handled."))
    if (cplx) {
      rels <- Filter(is_relational, .sset_subexprs(cond))
      if (any(!vapply(rels, function(i) i$type %in% c("Equality", "Unequality"), NA)))
        .sset_err("ValueError", .filldedent("Inequalities in the complex domain are not supported. Try the real domain by setting domain=S.Reals"))
    }
    cs <- .sset_as_set(cond)
    cond_int <- Intersection(U, cs)
    U <- Complement(U, cond_int)
    if (!.sset_is_empty(cond_int)) exp_sets[[length(exp_sets) + 1L]] <- list(expr, cond_int)
  }
  exp_sets
}
.sset_subexprs <- function(e) {
  out <- list(e)
  for (a in e$args) out <- c(out, .sset_subexprs(a))
  out
}

# ---------------------------------------------------------------------------
# Modular equations

.sset_term_factors <- function(f) {
  out <- list()
  for (a in .sset_add_args(f)) out <- c(out, .sset_mul_args(a))
  out
}

.sset_is_modular <- function(f, symbol) {
  f <- as_sym(f); symbol <- as_sym(symbol)
  if (!.sset_has_type(f, "Mod")) return(FALSE)
  modterms <- Filter(function(a) .sset_is_fn(a, "Mod"), atoms(f, "Function"))
  length(modterms) == 1L && .has1(modterms[[1L]]$args[[1L]], symbol) &&
    isTRUE(.sset_fact(modterms[[1L]]$args[[2L]], "integer")) &&
    any(vapply(.sset_term_factors(f), function(t) .sset_is_fn(t, "Mod"), NA))
}

.sset_invert_modular <- function(modterm, rhs, n, symbol) {
  a <- modterm$args[[1L]]; m <- modterm$args[[2L]]
  if (identical(.sset_fact(rhs, "integer"), FALSE)) return(.sset_pair(symbol, S_EmptySet))
  if (identical(.sset_fact(rhs, "real"), FALSE) ||
      any(vapply(.sset_term_factors(a), function(t) identical(.sset_fact(t, "real"), FALSE), NA)))
    return(.sset_pair(modterm, rhs))
  ge <- Ge(Abs(rhs), Abs(m))
  if (.is_true(ge)) return(.sset_pair(symbol, S_EmptySet))
  if (!.is_false(ge)) .sset_err("TypeError", paste0("cannot determine truth value of Relational: ", .str(ge)))
  if (.sset_eq(a, symbol)) return(.sset_pair(symbol, ImageSet(Lambda(n, .add(list(.mul(list(m, n)), rhs))), Integers)))
  if (a$type == "Add") {
    gh <- .sp_indep(a, list(symbol))
    if (!is_zero_num(gh[[1L]])) {
      x_indep_term <- .sub(rhs, Mod(gh[[1L]], m))
      return(.sset_invert_modular(Mod(gh[[2L]], m), Mod(x_indep_term, m), n, symbol))
    }
  }
  if (a$type == "Mul") {
    gh <- .sp_indep(a, list(symbol))
    if (!is_one_num(gh[[1L]])) {
      x_indep_term <- .mul(list(rhs, invert(gh[[1L]], m)))
      return(.sset_invert_modular(Mod(gh[[2L]], m), Mod(x_indep_term, m), n, symbol))
    }
  }
  if (a$type == "Pow") {
    base <- a$args[[1L]]; expo <- a$args[[2L]]
    if (.has1(expo, symbol) && !.has1(base, symbol)) {
      if (m$type != "Integer" && rhs$type == "Integer" && base$type == "Integer") return(.sset_pair(modterm, rhs))
      mp <- gmp::as.bigz(m$value); rp <- gmp::as.bigz(rhs$value); bp <- gmp::as.bigz(base$value)
      mdiv <- mp %/% gmp::gcd(mp, rp)
      remainder <- tryCatch(discrete_log(.int(mdiv), .int(rp), .int(bp)), error = function(e) NULL)
      if (is.null(remainder)) return(.sset_pair(modterm, rhs))
      remainder <- as_sym(remainder)
      period <- as_sym(totient(m))
      for (p in divisors(period)) {
        p <- as_sym(p)
        md <- mp %/% gmp::gcd(mp, bp)
        if (gmp::powm(bp, gmp::as.bigz(p$value), md) == 1L) { period <- p; break }
      }
      return(.sset_pair(expo, ImageSet(Lambda(n, .add(list(.mul(list(period, n)), remainder))), Naturals0)))
    } else if (.has1(base, symbol) && !.has1(expo, symbol)) {
      rl <- tryCatch(nthroot_mod(rhs, expo, m, all_roots = TRUE), error = function(e) "fail")
      if (identical(rl, "fail")) return(.sset_pair(modterm, rhs))
      if (is.null(rl) || !length(rl)) return(.sset_pair(symbol, S_EmptySet))
      g_n <- S_EmptySet
      for (rem in rl) g_n <- Union(g_n, ImageSet(Lambda(n, .add(list(.mul(list(m, n)), as_sym(rem)))), Integers))
      return(.sset_pair(base, g_n))
    }
  }
  .sset_pair(modterm, rhs)
}

.sset_solve_modular <- function(f, symbol, domain) {
  unsolved_result <- ConditionSet(symbol, Eq(f, 0), domain)
  modterm <- Filter(function(a) .sset_is_fn(a, "Mod"), atoms(f, "Function"))[[1L]]
  rhs <- .neg(subs(f, modterm, S_Zero))
  cd <- .sset_coeff_of_term(f, modterm)
  if (isTRUE(.sset_fact(cd, "negative"))) rhs <- .neg(rhs)
  if (!isTRUE(.sset_subset(domain, Integers))) return(unsolved_result)
  if (.has1(rhs, symbol)) return(unsolved_result)
  n <- Dummy("n", integer = TRUE)
  r <- .sset_invert_modular(modterm, rhs, n, symbol)
  f_x <- r[[1L]]; g_n <- r[[2L]]
  if (.sset_eq(f_x, modterm) && .sset_eq(g_n, rhs)) return(unsolved_result)
  if (.sset_eq(f_x, symbol)) {
    if (!.sset_set_is(domain, "Integers")) return(Intersection(domain, g_n))
    return(g_n)
  }
  if (g_n$type == "ImageSet") {
    lam <- g_n$args[[1L]]
    lamda_expr <- lam$args[[2L]]
    lamda_vars <- lam$args[[1L]]$args
    base_sets <- g_n$args[-1L]
    sol_set <- .sset_solveset(.sub(f_x, lamda_expr), symbol, Integers)
    if (.sset_is_fs(sol_set)) {
      tmp <- S_EmptySet
      for (sol in .sset_elems(sol_set)) tmp <- Union(tmp, do.call(ImageSet, c(list(Lambda(lamda_vars, sol)), base_sets)))
      sol_set <- tmp
    } else {
      sol_set <- do.call(ImageSet, c(list(Lambda(lamda_vars, sol_set)), base_sets))
    }
    return(Intersection(domain, sol_set))
  }
  unsolved_result
}

# f.as_coefficients_dict()[term]
.sset_coeff_of_term <- function(f, term) {
  for (t in .sset_add_args(f)) {
    cm <- .sset_coeff_mul_num(t)
    if (.sset_eq(cm[[2L]], term)) return(cm[[1L]])
  }
  S_Zero
}
.sset_coeff_mul_num <- function(t) {
  if (is_number_type(t)) return(list(t, S_One))
  if (t$type == "Mul") {
    a <- .sympy_args(t)
    if (is_number_type(a[[1L]])) return(list(a[[1L]], .mul(a[-1L])))
  }
  list(S_One, t)
}

# ---------------------------------------------------------------------------
# Transcendental equations

.sset_solve_exponential <- function(lhs, rhs, symbol, domain) {
  unsolved_result <- ConditionSet(symbol, Eq(.sub(lhs, rhs), 0), domain)
  newlhs <- powdenest(lhs)
  if (!.sset_eq(lhs, newlhs)) {
    neweq <- factor(.sub(newlhs, rhs))
    if (!.sset_eq(neweq, .sub(lhs, rhs))) return(.sset_solveset(neweq, symbol, domain))
  }
  if (!(lhs$type == "Add" && length(lhs$args) == 2L)) return(unsolved_result)
  if (!is_zero_num(rhs)) return(unsolved_result)
  ab <- .sp_ordered(lhs$args)
  a <- ab[[1L]]; b <- ab[[2L]]
  a_term <- .sp_indep(a, list(symbol))[[2L]]
  b_term <- .sp_indep(b, list(symbol))[[2L]]
  abe <- .sset_as_base_exp(a_term); bbe <- .sset_as_base_exp(b_term)
  if (.sset_sub_reals(domain)) {
    conditions <- And(Gt(abe[[1L]], 0), Gt(bbe[[1L]], 0), Eq(im(abe[[2L]]), 0), Eq(im(bbe[[2L]]), 0))
  } else {
    conditions <- And(Ne(abe[[1L]], 0), Ne(bbe[[1L]], 0))
  }
  L <- expand_log(log(a), force = TRUE)
  R <- expand_log(log(.neg(b)), force = TRUE)
  solutions <- .sset_solveset(.sub(L, R), symbol, domain)
  ConditionSet(symbol, conditions, solutions)
}

# Expr.as_base_exp
.sset_as_base_exp <- function(e) {
  if (e$type == "Pow") return(e$args)
  if (.sset_is_fn(e, "exp")) return(list(S_E, e$args[[1L]]))
  list(e, S_One)
}

.sset_is_exponential <- function(f, symbol) {
  f <- as_sym(f); symbol <- as_sym(symbol)
  rv <- FALSE
  for (expr_arg in .sset_term_factors(f)) {
    if (!.has1(expr_arg, symbol)) next
    if ((expr_arg$type == "Pow" && !.has1(expr_arg$args[[1L]], symbol)) || .sset_is_fn(expr_arg, "exp")) rv <- TRUE
    else return(FALSE)
  }
  rv
}

.sset_solve_logarithm <- function(lhs, rhs, symbol, domain) {
  new_lhs <- logcombine(lhs, force = TRUE)
  .sset_solveset(.sub(new_lhs, rhs), symbol, domain)
}

.sset_is_logarithmic <- function(f, symbol) {
  f <- as_sym(f); symbol <- as_sym(symbol)
  rv <- FALSE
  for (term in .sset_add_args(f)) {
    saw_log <- FALSE
    for (term_arg in .sset_mul_args(term)) {
      if (!.has1(term_arg, symbol)) next
      if (.sset_is_fn(term_arg, "log")) {
        if (saw_log) return(FALSE)
        saw_log <- TRUE
      } else return(FALSE)
    }
    if (saw_log) rv <- TRUE
  }
  rv
}

.sset_is_lambert <- function(f, symbol) {
  f <- as_sym(f); symbol <- as_sym(symbol)
  term_factors <- .sset_term_factors(expand(f))
  dep <- Filter(function(a) .has1(a, symbol), term_factors)
  no_of_symbols <- length(dep)
  no_of_trig <- length(Filter(function(a) .has_function(a, c(.SSET_TRIG, .SSET_HYP)), term_factors))
  if (f$type == "Add" && no_of_symbols >= 2L) {
    if (any(vapply(dep, function(a) .sset_is_fn(a, c("log", .SSET_TRIG, .SSET_HYP)), NA))) {
      if (no_of_trig < no_of_symbols) return(TRUE)
    } else if (any(vapply(dep, function(a) (a$type == "Pow" || .sset_is_fn(a, "exp")) &&
                          .has1(.sset_as_base_exp(a)[[2L]], symbol), NA))) {
      return(TRUE)
    }
  }
  FALSE
}

.sset_transolve <- function(f, symbol, domain) {
  f <- as_sym(f); symbol <- as_sym(symbol); domain <- as_sym(domain)
  add_type <- function(lhs, rhs) {
    result <- ConditionSet(symbol, Eq(.sub(lhs, rhs), 0), domain)
    if (.sset_is_exponential(lhs, symbol)) result <- .sset_solve_exponential(lhs, rhs, symbol, domain)
    else if (.sset_is_logarithmic(lhs, symbol)) result <- .sset_solve_logarithm(lhs, rhs, symbol, domain)
    result
  }
  result <- ConditionSet(symbol, Eq(f, 0), domain)
  inv <- .sset_invert(f, S_Zero, symbol, domain)
  lhs <- inv[[1L]]; rhs_s <- inv[[2L]]
  if (.sset_is_fs(rhs_s)) {
    rhs <- rhs_s$args[[1L]]
    if (lhs$type == "Add") result <- add_type(lhs, rhs)
  } else result <- rhs_s
  result
}

# ---------------------------------------------------------------------------
# solveset

#' Solve an equation or inequality as a set
#'
#' `solveset(f, symbol, domain)` returns the set of values of `symbol` in
#' `domain` (`Complexes` by default) for which `f = 0` holds, or for which a
#' relational `f` is true, as SymPy's `solveset()`. Periodic solutions are
#' `ImageSet`s over the integers, and equations it cannot solve give a
#' `ConditionSet`. `solveset_real()` and `solveset_complex()` fix the domain.
#'
#' @param f An expression, equation or inequality.
#' @param symbol The variable; may be omitted when `f` has one free symbol.
#' @param domain The domain, a set.
#' @return A set.
#' @examples
#' x <- Symbol("x")
#' solveset(x^2 - 1, x)
#' solveset(exp(x) - 1, x)
#' solveset(cos(x), x, Reals)
#' solveset(x^2 > 1, x, Reals)
#' @export
solveset <- function(f, symbol = NULL, domain = Complexes) {
  if (is.list(f) && !is_sym(f))
    .sset_err("ValueError", paste0("[", paste(vapply(f, function(z) .str(as_sym(z)), ""), collapse = ", "),
                                   "] is not a valid SymPy expression"))
  if (is.logical(f) && length(f) == 1L) f <- if (f) S_true else S_false
  f <- as_sym(f)
  if (!is.null(symbol)) symbol <- as_sym(symbol)
  if (.is_true(f)) return(domain)
  if (.is_false(f)) return(S_EmptySet)
  if (!(.sset_is_expr(f) || is_relational(f)))
    .sset_err("ValueError", paste0(.str(f), " is not a valid SymPy expression"))
  if (!is.null(symbol) && !(.sset_is_expr(symbol) || is_relational(symbol)))
    .sset_err("ValueError", paste0(.str(symbol), " is not a valid SymPy symbol"))
  domain <- as_sym(domain)
  if (!is_set(domain)) .sset_err("ValueError", paste0(.str(domain), " is not a valid domain"))
  fsyms <- free_symbols(f)
  if (.has_type(f, "Piecewise")) f <- piecewise_fold(f)
  if (is.null(symbol) && !length(fsyms)) {
    b <- Eq(f, 0)
    if (.is_true(b)) return(domain)
    if (.is_false(b)) return(S_EmptySet)
    .sset_err("NotImplementedError", .filldedent(paste0("relationship between value and 0 is unknown: ", .str(b))))
  }
  if (is.null(symbol)) {
    if (length(fsyms) == 1L) symbol <- fsyms[[1L]]
    else .sset_err("ValueError", .filldedent("The independent variable must be specified for a multivariate equation."))
  } else if (!(symbol$type %in% c("Symbol", "Dummy"))) {
    rc <- .sp_recast(list(f), list(symbol))
    res <- solveset(rc$f[[1L]], rc$symbols[[1L]], domain)
    return(if (length(rc$restore)) xreplace(res, rc$restore) else res)
  }
  newsym <- NULL
  if (.sset_sub_reals(domain)) {
    if (!.sset_orig_is(symbol, "real")) newsym <- Dummy("R", real = TRUE)
  } else if (isTRUE(.sset_subset(domain, Complexes))) {
    if (!.sset_orig_is(symbol, "complex")) newsym <- Dummy("C", complex = TRUE)
  }
  if (!is.null(newsym)) {
    rv <- solveset(xreplace(f, list(list(symbol, newsym))), newsym, domain)
    rv2 <- tryCatch(xreplace(rv, list(list(newsym, symbol))), error = function(e) NULL)
    if (!is.null(rv2) && .sset_dummy_eq_any(rv, rv2)) rv <- rv2
    return(rv)
  }
  mk <- .sset_masked(f, "Abs")
  f <- rewrite(mk[[1L]], "Piecewise")
  for (de in mk[[2L]]) {
    e <- de[[2L]]
    e <- fcall("Abs", list(rewrite(e$args[[1L]], "Piecewise")))
    f <- xreplace(f, list(list(de[[1L]], e)))
  }
  f <- piecewise_fold(f)
  .sset_solveset(f, symbol, domain, check = TRUE)
}

#' @rdname solveset
#' @export
solveset_real <- function(f, symbol) {
  if (missing(symbol)) .sset_err("TypeError", "solveset_real() missing 1 required positional argument: 'symbol'")
  solveset(f, symbol, Reals)
}

#' @rdname solveset
#' @export
solveset_complex <- function(f, symbol) {
  if (missing(symbol)) .sset_err("TypeError", "solveset_complex() missing 1 required positional argument: 'symbol'")
  solveset(f, symbol, Complexes)
}

# ---------------------------------------------------------------------------
# solvify

#' Solve an equation as a list, as solve() would
#'
#' `solvify(f, symbol, domain)` solves `f = 0` with [solveset()] and returns
#' the solutions as a list; for periodic equations it lists the solutions in
#' one period, as SymPy's `solvify()`. It returns `NULL` when the solutions
#' cannot be listed.
#'
#' @param f An expression.
#' @param symbol The variable.
#' @param domain The domain, a set.
#' @return A list of expressions, or `NULL`.
#' @examples
#' x <- Symbol("x")
#' solvify(x^2 - 9, x, Reals)
#' solvify(sin(x) + 1, x, Reals)
#' @export
solvify <- function(f, symbol, domain) {
  f <- as_sym(f); symbol <- as_sym(symbol); domain <- as_sym(domain)
  solution_set <- solveset(f, symbol, domain)
  result <- NULL
  if (.sset_is_empty(solution_set)) {
    result <- list()
  } else if (solution_set$type == "ConditionSet") {
    .sset_err("NotImplementedError", "solveset is unable to solve this equation.")
  } else if (.sset_is_fs(solution_set)) {
    result <- .sset_elems(solution_set)
  } else {
    period <- periodicity(f, symbol)
    if (!is.null(period)) {
      solutions <- S_EmptySet
      iter_solutions <- list()
      if (solution_set$type == "ImageSet") iter_solutions <- list(solution_set)
      else if (solution_set$type == "Union" && all(vapply(solution_set$args, function(i) i$type == "ImageSet", NA)))
        iter_solutions <- .set_print_order(solution_set$args)
      for (solution in iter_solutions) solutions <- Union(solutions, Intersection(solution, Interval(0, period, FALSE, TRUE)))
      if (.sset_is_fs(solutions)) result <- .sset_elems(solutions)
    } else {
      solution <- Intersection(solution_set, domain)
      if (solution$type == "Union") {
        if (any(vapply(solution$args, .sset_is_fs, NA))) {
          result <- list()
          for (soln in .set_print_order(solution$args)) if (.sset_is_fs(soln)) result <- c(result, .sset_elems(soln))
        } else return(NULL)
      } else if (.sset_is_fs(solution)) {
        .sset_err("TypeError", "unsupported operand type(s) for +=: 'NoneType' and 'FiniteSet'")
      }
    }
  }
  result
}

# ---------------------------------------------------------------------------
# solve_decomposition

#' Solve a composite equation by decomposition
#'
#' `solve_decomposition(f, symbol, domain)` writes `f` as a composition of
#' simpler functions with [decompogen()] and solves the equations from the
#' outside in, as SymPy's `solve_decomposition()`.
#'
#' @param f An expression.
#' @param symbol The variable.
#' @param domain The domain, a set.
#' @return A set.
#' @examples
#' x <- Symbol("x")
#' solve_decomposition(exp(2*x) - 3*exp(x) + 2, x, Reals)
#' @export
solve_decomposition <- function(f, symbol, domain) {
  f <- as_sym(f); symbol <- as_sym(symbol); domain <- as_sym(domain)
  g_s <- decompogen(f, symbol)
  y_s <- FiniteSet(S_Zero)
  for (g in g_s) {
    frange <- function_range(g, symbol, domain)
    y_s <- Intersection(frange, y_s)
    result <- S_EmptySet
    if (.sset_is_fs(y_s)) {
      for (y in .sset_elems(y_s)) {
        solutions <- solveset(Eq(g, y), symbol, domain)
        if (solutions$type != "ConditionSet") result <- Union(result, solutions)
      }
    } else {
      iter_iset <- NULL
      if (y_s$type == "ImageSet") iter_iset <- list(y_s)
      else if (y_s$type == "Union") iter_iset <- .set_print_order(y_s$args)
      else if (.sset_is_empty(y_s)) return(S_EmptySet)
      for (iset in iter_iset) {
        lam <- iset$args[[1L]]
        new_solutions <- solveset(Eq(lam$args[[2L]], g), symbol, domain)
        dummy_var <- free_symbols(lam$args[[2L]])[[1L]]
        base_set <- iset$args[[2L]]
        new_exprs <- list()
        if (.sset_is_fs(new_solutions)) new_exprs <- .sset_elems(new_solutions)
        else if (new_solutions$type == "Intersection") {
          a2 <- new_solutions$args[[2L]]
          if (.sset_is_fs(a2)) new_exprs <- .sset_elems(a2)
        }
        for (new_expr in new_exprs) result <- Union(result, ImageSet(Lambda(dummy_var, new_expr), base_set))
      }
    }
    if (.sset_is_empty(result)) return(ConditionSet(symbol, Eq(f, 0), domain))
    y_s <- result
  }
  y_s
}

# ---------------------------------------------------------------------------
# _solveset_multi

.sset_solveset_multi <- function(eqs, syms, domains) {
  eqs <- lapply(eqs, as_sym); syms <- lapply(syms, as_sym); domains <- lapply(domains, as_sym)
  rep <- list()
  for (i in seq_along(syms)) if (.sset_set_is(domains[[i]], "Reals"))
    rep[[length(rep) + 1L]] <- list(syms[[i]], Symbol(syms[[i]]$name, real = TRUE))
  if (length(rep)) {
    eqs <- lapply(eqs, function(e) subs(e, rep))
    syms <- lapply(syms, function(s) subs(s, rep))
  }
  if (!length(eqs)) return(do.call(ProductSet, domains))
  if (length(syms) == 1L) {
    sym <- syms[[1L]]; domain <- domains[[1L]]
    solsets <- lapply(eqs, function(eq) solveset(eq, sym, domain))
    solset <- do.call(Intersection, solsets)
    return(.imageset_doit(Lambda(sym, Tuple(sym)), solset))
  }
  nfree <- vapply(eqs, function(eq) length(Filter(function(s) .has1(eq, s), syms)), 0L)
  eqs <- eqs[order(nfree, method = "radix")]
  for (n in seq_along(eqs)) {
    eq <- eqs[[n]]
    sols <- list()
    all_handled <- TRUE
    for (j in seq_along(syms)) {
      sym <- syms[[j]]
      if (!.has1(eq, sym)) next
      sol <- solveset(eq, sym, domains[[j]])
      if (.sset_is_fs(sol)) {
        symsp <- syms[-j]; domainsp <- domains[-j]; eqsp <- eqs[-n]
        for (s in .sset_elems(sol)) {
          eqsp_sub <- lapply(eqsp, function(e) subs(e, sym, s))
          sol_others <- .sset_solveset_multi(eqsp_sub, symsp, domainsp)
          body <- append(symsp, list(s), after = j - 1L)
          fun <- new_sym("Lambda", list(Tuple(do.call(Tuple, symsp)), do.call(Tuple, body)))
          sols[[length(sols) + 1L]] <- .imageset_doit(fun, sol_others)
        }
      } else all_handled <- FALSE
    }
    if (all_handled) return(.sset_union(sols))
  }
  NULL
}
